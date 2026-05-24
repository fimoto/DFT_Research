#!/usr/bin/env python3
"""Direct Monte Carlo/QMC implementation of Benites--Rosado--Manousakis
Eqs. (28)--(32), Phys. Rev. B 110, 195151 (2024).

This module is deliberately independent of the momentum-resolved imaginary-axis
kernel used by ``mc/xc_mc.py``.  It evaluates the BRM post-frequency-integration
representation directly:

  Delta epsilon_2b = Delta epsilon_2b,1 + Delta epsilon_2b,2,

with the coupling constant lambda and the mapped frequency variable u sampled
explicitly.  A global rotation fixes q1 along z, reducing the nominal 11-D
integral to an exactly equivalent 9-D integral.  The q2 proposal is centered on
the Coulomb singularity and uses a radial-uniform density; its Jacobian cancels
the explicit 1/|...|^2 factor.  k is sampled uniformly in the particle-hole
exit cap, so both Heaviside factors in BRM Eq. (30) are enforced exactly.

The output is in Ry per electron, as in BRM.  Adding ONSAGER_RY gives the full
kite value tabulated in their Table IX.

This evaluator generates the locked records in production/raw through
production/compute.py. Its randomized Sobol error bar measures
replicate-to-replicate QMC scatter; numerical systematics require the separate
checks described in REPRODUCE_SPIN.md and are not bounded by that error bar.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from scipy.stats import qmc

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '..', 'src'))
from ueg_exchange.lindhard import lindhard_shape as _lindhard_G  # legacy BRM reference boundary

ALPHA = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)
ZETA3 = 1.2020569031595942854
ONSAGER_RY = np.log(2.0) / 3.0 - 3.0 * ZETA3 / (2.0 * np.pi**2)
VBALL = 4.0 * np.pi / 3.0

# GL nodes used only in cancellation-prone regions of the BRM Lindhard g.
_GL_X, _GL_W = np.polynomial.legendre.leggauss(64)
_GL_K = 0.5 * (_GL_X + 1.0)
_GL_W = 0.5 * _GL_W
# 96-point rule on [0, 1] for the integral form of dg/dnu at q >= DG_INTEGRAL_Q
_GL96_X, _GL96_W = np.polynomial.legendre.leggauss(96)
_GL96_K = 0.5 * (_GL96_X + 1.0)
_GL96_W = 0.5 * _GL96_W
DG_INTEGRAL_Q = 3.0


def _clip_unit(x: np.ndarray) -> np.ndarray:
    """Keep inverse maps away from exact endpoints without measurable bias."""
    eps = np.finfo(float).eps
    return np.clip(np.asarray(x, float), eps, 1.0 - eps)


def _antideriv_P(w, xi, a, sgn):
    """Antiderivative of (w - sgn*a) log(w^2 + xi^2); used for the closed-form g."""
    r2 = w * w + xi * xi
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.where(r2 > 0, np.log(np.where(r2 > 0, r2, 1.0)), 0.0)
        at = np.where(xi > 0, xi * np.arctan(np.where(xi > 0, w / np.where(xi > 0, xi, 1.0), 0.0)), 0.0)
    return 0.5 * (r2 * L - w * w) - sgn * a * (w * L - 2 * w + 2 * at)


def _g_closed(q, nu):
    """g(q, i nu) = 2 G(q, nu) from the exact antiderivative of the k integral."""
    q = np.maximum(np.asarray(q, float), 1e-14)
    nu = np.abs(np.asarray(nu, float))
    a = 0.5 * q * q
    val = ((_antideriv_P(q + a, nu, a, +1.0) - _antideriv_P(a, nu, a, +1.0))
           - (_antideriv_P(q - a, nu, a, -1.0) - _antideriv_P(-a, nu, a, -1.0))) / (2 * q ** 3)
    return 2.0 * val


def g_brm(q: np.ndarray, nu: np.ndarray) -> np.ndarray:
    r"""BRM imaginary-axis Lindhard function g(q,i nu), Eq. (14).

    This normalization is twice ``ueg_exchange.lindhard.lindhard_shape``, to which the
    evaluation is delegated.
    """
    q, nu = np.broadcast_arrays(np.asarray(q, float), np.abs(np.asarray(nu, float)))
    out = np.empty_like(q)
    qsafe = np.maximum(q, 1.0e-14)

    # g = 2G is delegated to the validated evaluation of ueg_exchange.lindhard,
    # which switches between the small-q series, the exact antiderivative and a
    # log1p quadrature.  The closed BRM form below loses all its digits at large
    # q (at q = 1e8 it returns a negative value for a positive quantity), and a
    # fixed quadrature cannot resolve the k = q/2 logarithmic peak at small nu.
    zero_q = (qsafe < 1.0e-13) & (nu < 1.0e-13)   # removable endpoint: g(0,0) = 2
    out[zero_q] = 2.0
    rest = ~zero_q
    if np.any(rest):
        out[rest] = 2.0 * _lindhard_G(qsafe[rest], nu[rest])
    return out


def _g_brm_legacy(q, nu):
    """The original branch structure, kept for reference; see g_brm."""
    q = np.atleast_1d(np.asarray(q, float))
    nu = np.atleast_1d(np.abs(np.asarray(nu, float)))
    q, nu = np.broadcast_arrays(q, nu)
    out = np.empty(q.shape)
    qsafe = np.maximum(q, 1.0e-300)
    kp = qsafe + 0.5 * qsafe**2
    km = qsafe - 0.5 * qsafe**2
    ap = kp * kp + nu * nu
    am = km * km + nu * nu
    use_closed = (qsafe < 2.0e-3) & (nu <= 2.0 * (qsafe + 0.5 * qsafe**2))
    use_quad = (nu > 2.0 * (qsafe + 0.5 * qsafe**2))
    zero_q = (qsafe < 1.0e-13) & (nu < 1.0e-13)   # removable endpoint: g(0,0) = 2
    branch = (np.abs(qsafe - 2.0) < 1.0e-14) & (nu < 1.0e-14)
    use_quad = use_quad & ~branch
    use_closed = use_closed & ~branch
    use_quad = use_quad & ~zero_q
    use_closed = use_closed & ~zero_q
    direct = ~(use_quad | branch | use_closed | zero_q)
    out[branch] = 1.0
    out[zero_q] = 2.0
    if np.any(use_closed):
        out[use_closed] = _g_closed(qsafe[use_closed], nu[use_closed])
    if np.any(direct):
        qq = qsafe[direct]
        nn = nu[direct]
        kpp = kp[direct]
        kmm = km[direct]
        # log(ap/am) via log1p when the ratio is close to one.
        amm = kmm * kmm + nn * nn
        delta = kpp * kpp - kmm * kmm
        logratio = np.log1p(delta / amm)
        coeff = (nn * nn + kpp * kmm) / (2.0 * qq**3)
        angles = np.arctan2(kpp, nn) + np.arctan2(kmm, nn)
        # atan2(kappa,0) gives the correct one-sided nu->0 limit;
        # the prefactor nu makes this term vanish exactly at nu=0.
        out[direct] = 1.0 + coeff * logratio - (nn / qq) * angles

    if np.any(use_quad):
        idx = np.flatnonzero(use_quad.ravel())
        qflat = qsafe.ravel()
        nflat = nu.ravel()
        oflat = out.ravel()
        # Chunk to bound temporary memory.
        for start in range(0, len(idx), 4096):
            ii = idx[start : start + 4096]
            qq = qflat[ii][None, :]
            nn = nflat[ii][None, :]
            kk = _GL_K[:, None]
            a = 0.5 * qq * qq
            den = (qq * kk - a) ** 2 + nn * nn
            # numerator/denominator = 1 + 4*a*q*k/den
            val = (_GL_W[:, None] * kk * np.log1p(4.0 * a * qq * kk / den)).sum(axis=0) / qq[0]
            oflat[ii] = val
    return out


def dg_dnu_brm(q: np.ndarray, nu: np.ndarray) -> np.ndarray:
    """Derivative d g(q,i nu)/d nu, with stable finite-difference fallback."""
    q, nu = np.broadcast_arrays(np.asarray(q, float), np.abs(np.asarray(nu, float)))
    qsafe = np.maximum(q, 1.0e-14)
    # For q >= 3 the closed form below subtracts nearly equal terms (at q ~ 1e3-1e4 and
    # small nu it can return the wrong sign). There the integral form
    #   dg/dnu = -8 a nu int_0^1 k^2 dk / {[(qk+a)^2+nu^2][(qk-a)^2+nu^2]},  a = q^2/2,
    # has a smooth positive integrand (qk - a < 0 on [0, 1] for q > 2) and is used instead.
    large = qsafe >= DG_INTEGRAL_Q
    kp = qsafe + 0.5 * qsafe**2
    km = qsafe - 0.5 * qsafe**2
    ap = kp * kp + nu * nu
    am = km * km + nu * nu
    use_fd = (
        (qsafe < 2.0e-3)
        | (nu > 2.0 * (qsafe + 0.5 * qsafe**2))
    ) & ~large
    out = np.empty_like(qsafe)
    # At the 2k_F edge the analytic expression has removable 0/0 terms.
    # Its right derivative at nu=0 is half the q<2 cusp limit: -pi/4.
    edge = (qsafe == 2.0) & (nu == 0.0)
    out[edge] = -np.pi / 4.0
    direct = ~use_fd & ~large & ~edge
    if np.any(large):
        qq = qsafe[large][:, None]
        nn = nu[large][:, None]
        aa = 0.5 * qq * qq
        kk = _GL96_K[None, :]
        den = ((qq * kk + aa) ** 2 + nn * nn) * ((qq * kk - aa) ** 2 + nn * nn)
        out[large] = -8.0 * aa[:, 0] * nn[:, 0] * (_GL96_W[None, :] * kk * kk / den).sum(axis=1)
    if np.any(direct):
        qq = qsafe[direct]
        nn = nu[direct]
        kpp = kp[direct]
        kmm = km[direct]
        app = ap[direct]
        amm = am[direct]
        logratio = np.log1p((kpp * kpp - kmm * kmm) / amm)
        c = nn * nn + kpp * kmm
        angles = np.arctan2(kpp, nn) + np.arctan2(kmm, nn)
        dlog = 2.0 * nn / app - 2.0 * nn / amm
        dfirst = (2.0 * nn * logratio + c * dlog) / (2.0 * qq**3)
        datan_part = -angles / qq + (nn / qq) * (kpp / app + kmm / amm)
        out[direct] = dfirst + datan_part
    if np.any(use_fd):
        qq = qsafe[use_fd]
        nn = nu[use_fd]
        # g is even in nu but has a cusp at nu = 0, so the one-sided derivative
        # is not zero there.  Below q = 0.15 the series of ueg_exchange.lindhard is
        # differentiated analytically, which is exact at the cusp and uniformly
        # accurate away from it; replacing the whole region nu < 1e-3 q by the
        # endpoint value -pi/q is wrong by 0.1 % already at nu/q = 9e-4.
        # Elsewhere a central difference is used, with a step tied to
        # max(q, nu) because g varies on the scale nu ~ q.
        small = qq < 0.15
        fd = np.empty_like(qq)
        if np.any(small):
            z = nn[small] / qq[small]
            z2 = z * z
            # h'(z) = -arctan(1/z) + z/(1+z^2). Written directly this subtracts
            # two O(1/z) terms to expose an O(1/z^3) result, so it is evaluated
            # from its series at large z and from its Taylor form at small z.
            hp = np.empty_like(z)
            big = z >= 30.0
            tiny = z <= 1e-4
            mid = ~(big | tiny)
            if np.any(big):
                zb = z[big]
                hp[big] = -(2.0 / (3.0 * zb ** 3) - 4.0 / (5.0 * zb ** 5)
                            + 6.0 / (7.0 * zb ** 7) - 8.0 / (9.0 * zb ** 9)
                            + 10.0 / (11.0 * zb ** 11) - 12.0 / (13.0 * zb ** 13))
            if np.any(tiny):
                zt = z[tiny]
                hp[tiny] = -0.5 * np.pi + 2.0 * zt - 4.0 * zt ** 3 / 3.0
            if np.any(mid):
                zm = z[mid]
                hp[mid] = -np.arctan(1.0 / zm) + zm / (1.0 + zm * zm)
            c2p = z / (3.0 * (1.0 + z2) ** 3)
            c4p = z * (3.0 - 5.0 * z2) / (40.0 * (1.0 + z2) ** 5)
            qs = qq[small]
            fd[small] = 2.0 * (hp + qs ** 2 * c2p + qs ** 4 * c4p) / qs
        if np.any(~small):
            qb, nb = qq[~small], nn[~small]
            h = 2.0e-5 * np.maximum(qb, nb)
            lo = np.maximum(0.0, nb - h)
            hi = nb + h
            fd[~small] = (g_brm(qb, hi) - g_brm(qb, lo)) / (hi - lo)
        out[use_fd] = fd
    return out


def spin_radii(zeta: float) -> tuple[float, float]:
    if not (-1.0 <= zeta <= 1.0):
        raise ValueError("zeta must lie in [-1,1]")
    return ((1.0 + zeta) ** (1.0 / 3.0), (1.0 - zeta) ** (1.0 / 3.0))


def f_lambda_sigma(
    q: np.ndarray,
    nu: np.ndarray,
    lam: np.ndarray,
    rs: float,
    zeta: float,
    sigma: int,
) -> np.ndarray:
    r"""BRM f^sigma_lambda(q,nu,zeta), Eq. (32)."""
    q, nu, lam = np.broadcast_arrays(
        np.asarray(q, float), np.abs(np.asarray(nu, float)), np.asarray(lam, float)
    )
    xs = spin_radii(zeta)
    x = xs[0 if sigma > 0 else 1]
    if x < 1.0e-14:
        return np.zeros_like(q)
    ssum = np.zeros_like(q)
    for xp in xs:
        if xp < 1.0e-14:
            continue
        theta = x / xp
        ssum += xp * g_brm(theta * q, theta * theta * nu)
    den = x * x + lam * ALPHA * rs * ssum / (np.pi * np.maximum(q, 1.0e-300) ** 2)
    return x**3 * ssum / den


def df_dnu_lambda_sigma(
    q: np.ndarray,
    nu: np.ndarray,
    lam: np.ndarray,
    rs: float,
    zeta: float,
    sigma: int,
) -> np.ndarray:
    """Analytic nu derivative of Eq. (32), used for its divided difference."""
    q, nu, lam = np.broadcast_arrays(
        np.asarray(q, float), np.abs(np.asarray(nu, float)), np.asarray(lam, float)
    )
    xs = spin_radii(zeta)
    x = xs[0 if sigma > 0 else 1]
    if x < 1.0e-14:
        return np.zeros_like(q)
    ssum = np.zeros_like(q)
    dsum = np.zeros_like(q)
    for xp in xs:
        if xp < 1.0e-14:
            continue
        theta = x / xp
        qarg = theta * q
        narg = theta * theta * nu
        ssum += xp * g_brm(qarg, narg)
        dsum += xp * theta * theta * dg_dnu_brm(qarg, narg)
    beta = lam * ALPHA * rs / (np.pi * np.maximum(q, 1.0e-300) ** 2)
    den = x * x + beta * ssum
    # d[x^3 S/(x^2+beta S)] = x^5 S'/(x^2+beta S)^2
    return x**5 * dsum / den**2


def f_spin_sum(q: np.ndarray, nu: np.ndarray, lam: np.ndarray, rs: float, zeta: float) -> np.ndarray:
    return f_lambda_sigma(q, nu, lam, rs, zeta, +1) + f_lambda_sigma(q, nu, lam, rs, zeta, -1)


def df_spin_sum(q: np.ndarray, nu: np.ndarray, lam: np.ndarray, rs: float, zeta: float) -> np.ndarray:
    return df_dnu_lambda_sigma(q, nu, lam, rs, zeta, +1) + df_dnu_lambda_sigma(q, nu, lam, rs, zeta, -1)


def exit_cap_volume(q: np.ndarray) -> np.ndarray:
    """Volume of {p: |p|<=1, |p+q zhat|>1}."""
    q = np.asarray(q, float)
    return np.where(q >= 2.0, VBALL, np.pi * q * (12.0 - q * q) / 12.0)


def sample_exit_cap(q: np.ndarray, u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Uniform sample in the BRM particle-hole exit cap for q along z.

    Parameters
    ----------
    q : (n,) positive array
    u : (n,3) unit-cube points
    """
    q = np.asarray(q, float)
    u = _clip_unit(u)
    n = len(q)
    p = np.empty((n, 3), float)
    vol = exit_cap_volume(q)

    full = q >= 2.0
    if np.any(full):
        uu = u[full]
        r = uu[:, 0] ** (1.0 / 3.0)
        ct = 2.0 * uu[:, 1] - 1.0
        ph = 2.0 * np.pi * uu[:, 2]
        st = np.sqrt(np.maximum(0.0, 1.0 - ct * ct))
        p[full, 0] = r * st * np.cos(ph)
        p[full, 1] = r * st * np.sin(ph)
        p[full, 2] = r * ct

    part = ~full
    if np.any(part):
        qq = q[part]
        uu = u[part]
        v1 = np.pi * qq * (1.0 - 0.5 * qq) ** 2
        vt = vol[part]
        prob1 = np.clip(v1 / vt, 0.0, 1.0)
        choose1 = uu[:, 0] < prob1
        z = np.empty_like(qq)

        if np.any(choose1):
            t = uu[choose1, 0] / np.maximum(prob1[choose1], 1.0e-300)
            z[choose1] = -0.5 * qq[choose1] + (1.0 - 0.5 * qq[choose1]) * np.sqrt(t)

        choose2 = ~choose1
        if np.any(choose2):
            q2 = qq[choose2]
            z0 = 1.0 - q2
            t = (uu[choose2, 0] - prob1[choose2]) / np.maximum(1.0 - prob1[choose2], 1.0e-300)
            f0 = z0 - z0**3 / 3.0
            f1 = 2.0 / 3.0
            target = f0 + t * (f1 - f0)
            lo = z0.copy()
            hi = np.ones_like(lo)
            for _ in range(45):
                mid = 0.5 * (lo + hi)
                fm = mid - mid**3 / 3.0
                lowmask = fm < target
                lo = np.where(lowmask, mid, lo)
                hi = np.where(lowmask, hi, mid)
            z[choose2] = 0.5 * (lo + hi)

        upper2 = np.maximum(0.0, 1.0 - z * z)
        lower2 = np.maximum(0.0, 1.0 - (z + qq) ** 2)
        rho2 = lower2 + uu[:, 1] * np.maximum(0.0, upper2 - lower2)
        rho = np.sqrt(rho2)
        ph = 2.0 * np.pi * uu[:, 2]
        pp = np.column_stack([rho * np.cos(ph), rho * np.sin(ph), z])
        p[part] = pp
    return p, vol


def _sample_direction(u1: np.ndarray, u2: np.ndarray) -> np.ndarray:
    ct = 2.0 * u1 - 1.0
    ph = 2.0 * np.pi * u2
    st = np.sqrt(np.maximum(0.0, 1.0 - ct * ct))
    return np.column_stack([st * np.cos(ph), st * np.sin(ph), ct])


def _valid_exit_cap(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    pin = np.einsum("ij,ij->i", p, p) <= 1.0
    shifted = p.copy()
    shifted[:, 2] += q
    pout = np.einsum("ij,ij->i", shifted, shifted) > 1.0
    return pin & pout


def _divided_difference(
    q: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    tan_u: np.ndarray,
    lam: np.ndarray,
    rs: float,
    zeta: float,
) -> np.ndarray:
    na = a * tan_u
    nb = b * tan_u
    fa = f_spin_sum(q, na, lam, rs, zeta)
    fb = f_spin_sum(q, nb, lam, rs, zeta)
    den = a - b
    scale = np.maximum(1.0, np.abs(a) + np.abs(b))
    close = np.abs(den) < 2.0e-7 * scale
    out = np.empty_like(den)
    out[~close] = (fa[~close] - fb[~close]) / den[~close]
    if np.any(close):
        nm = 0.5 * (na[close] + nb[close])
        out[close] = tan_u[close] * df_spin_sum(q[close], nm, lam[close], rs, zeta)
    return out


def _q_map(t: np.ndarray, scale: float) -> tuple[np.ndarray, np.ndarray]:
    """q=scale*t/(1-t), returning q and d^3q probability weight."""
    t = _clip_unit(t)
    qv = scale * t / (1.0 - t)
    dqdt = scale / (1.0 - t) ** 2
    return qv, 4.0 * np.pi * qv**2 * dqdt


def evaluate_term(points: np.ndarray, rs: float, zeta: float, term: Literal[1, 2], q_scale: float) -> np.ndarray:
    """Vectorized unbiased importance-sampled BRM Eq. (28) or Eq. (29)."""
    if points.shape[1] != 9:
        raise ValueError("axis-reduced BRM sampler requires 9 Sobol coordinates")
    u = _clip_unit(points)
    lam = u[:, 0]
    qv, qweight = _q_map(u[:, 1], q_scale)
    k, kvol = sample_exit_cap(qv, u[:, 2:5])
    direction = _sample_direction(u[:, 6], u[:, 7])

    qvec = np.zeros_like(k)
    qvec[:, 2] = qv
    if term == 1:
        center = k
        radius = np.full_like(qv, 2.0)
    elif term == 2:
        center = -k - qvec
        radius = np.linalg.norm(center, axis=1) + 1.0
    else:
        raise ValueError("term must be 1 or 2")

    rho = radius * u[:, 5]
    q2 = center + rho[:, None] * direction
    valid = _valid_exit_cap(q2, qv)

    # a(q,p)=q^2/2+q*p_z when q is along z.
    aq2 = 0.5 * qv**2 + qv * q2[:, 2]
    ak = 0.5 * qv**2 + qv * k[:, 2]
    tan_u = np.tan(0.5 * np.pi * u[:, 8])

    pref = 3.0 * ALPHA * rs / (16.0 * np.pi**7)
    L = pref * lam**2 / np.maximum(qv, 1.0e-300) ** 4
    common = qweight * kvol * (4.0 * np.pi * radius) * (0.5 * np.pi) * L

    out = np.zeros_like(qv)
    good = valid & (aq2 > 0.0) & (ak > 0.0) & np.isfinite(common)
    if not np.any(good):
        return out

    if term == 1:
        dd = _divided_difference(
            qv[good], aq2[good], ak[good], tan_u[good], lam[good], rs, zeta
        )
        out[good] = common[good] * dd
    else:
        fa = f_spin_sum(qv[good], aq2[good] * tan_u[good], lam[good], rs, zeta)
        fb = f_spin_sum(qv[good], ak[good] * tan_u[good], lam[good], rs, zeta)
        denom = aq2[good] + ak[good]
        out[good] = -common[good] * (fa + fb) / denom
    return out


@dataclass
class Estimate:
    rs: float
    zeta: float
    m: int
    replicates: int
    seed: int
    q_scale: float
    delta1_ry: float
    delta1_err_ry: float
    delta2_ry: float
    delta2_err_ry: float
    correction_ry: float
    correction_err_ry: float
    full_kite_ry: float
    full_kite_err_ry: float
    full_kite_mry: float
    full_kite_err_mry: float
    replicate_values_delta1: list[float]
    replicate_values_delta2: list[float]


def estimate(
    rs: float,
    zeta: float,
    m: int = 16,
    replicates: int = 8,
    seed: int = 20260713,
    q_scale: float = 2.0,
) -> Estimate:
    """Randomized-Sobol estimate of BRM Eqs. (28)--(29)."""
    if m < 8:
        raise ValueError("m must be >= 8")
    if replicates < 2:
        raise ValueError("replicates must be >= 2 to estimate an error bar")
    # Input validation (audit F3, 1 Oct 2026): nonfinite or nonphysical inputs
    # previously produced a finite energy with zero error after the nonfinite
    # weights were masked out.
    if not (np.isfinite(rs) and rs > 0):
        raise ValueError("rs must be finite and positive")
    if not (np.isfinite(zeta) and -1 <= zeta <= 1):
        raise ValueError("zeta must be finite and within [-1, 1]")
    if not (np.isfinite(q_scale) and q_scale > 0):
        raise ValueError("q_scale must be finite and positive")
    vals1: list[float] = []
    vals2: list[float] = []
    for rep in range(replicates):
        sob1 = qmc.Sobol(d=9, scramble=True, seed=seed + 2 * rep)
        sob2 = qmc.Sobol(d=9, scramble=True, seed=seed + 2 * rep + 1)
        p1 = sob1.random_base2(m)
        p2 = sob2.random_base2(m)
        vals1.append(float(np.mean(evaluate_term(p1, rs, zeta, 1, q_scale))))
        vals2.append(float(np.mean(evaluate_term(p2, rs, zeta, 2, q_scale))))
    a1 = np.asarray(vals1)
    a2 = np.asarray(vals2)
    mu1 = float(a1.mean())
    mu2 = float(a2.mean())
    e1 = float(a1.std(ddof=1) / np.sqrt(replicates))
    e2 = float(a2.std(ddof=1) / np.sqrt(replicates))
    corr = mu1 + mu2
    ecorr = float(np.hypot(e1, e2))
    full = ONSAGER_RY + corr
    return Estimate(
        rs=rs,
        zeta=zeta,
        m=m,
        replicates=replicates,
        seed=seed,
        q_scale=q_scale,
        delta1_ry=mu1,
        delta1_err_ry=e1,
        delta2_ry=mu2,
        delta2_err_ry=e2,
        correction_ry=corr,
        correction_err_ry=ecorr,
        full_kite_ry=full,
        full_kite_err_ry=ecorr,
        full_kite_mry=1000.0 * full,
        full_kite_err_mry=1000.0 * ecorr,
        replicate_values_delta1=vals1,
        replicate_values_delta2=vals2,
    )


def self_test() -> None:
    """Fast analytic and sampling checks; raises AssertionError on failure."""
    # Closed-form g against its defining integral in non-pathological points.
    for qv, nv in [(0.2, 0.0), (1.0, 0.3), (2.0, 0.0), (3.0, 4.0), (0.01, 30.0)]:
        g = float(g_brm(np.array([qv]), np.array([nv]))[0])
        assert np.isfinite(g) and g >= 0.0, (qv, nv, g)

    # zeta=0 reduction of Eq. (32): f=2g/(1+lambda*2 alpha rs g/(pi q^2)).
    qv = np.array([0.4, 1.2, 2.7])
    nv = np.array([0.0, 0.5, 3.0])
    lv = np.array([0.2, 0.7, 1.0])
    rs = 3.0
    gg = g_brm(qv, nv)
    expected = 2.0 * gg / (1.0 + lv * 2.0 * ALPHA * rs * gg / (np.pi * qv**2))
    got = f_lambda_sigma(qv, nv, lv, rs, 0.0, +1)
    assert np.allclose(got, expected, rtol=3e-12, atol=3e-12)

    # Cap geometry and volume.
    rng = np.random.default_rng(4)
    qarr = np.repeat(np.array([0.1, 1.0, 1.9, 2.5]), 512)
    pts, vol = sample_exit_cap(qarr, rng.random((len(qarr), 3)))
    assert np.all(np.einsum("ij,ij->i", pts, pts) <= 1.0 + 1e-12)
    shifted = pts.copy(); shifted[:, 2] += qarr
    assert np.all(np.einsum("ij,ij->i", shifted, shifted) >= 1.0 - 1e-11)
    assert np.all(vol > 0.0)

    # rs->0 correction must vanish pointwise through the explicit rs prefactor.
    sob = qmc.Sobol(d=9, scramble=True, seed=8).random_base2(10)
    tiny = evaluate_term(sob, 1.0e-12, 0.0, 1, 2.0) + evaluate_term(sob, 1.0e-12, 0.0, 2, 2.0)
    assert abs(float(np.mean(tiny))) < 1.0e-10


def _main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rs", type=float, default=5.0)
    ap.add_argument("--zeta", type=float, default=0.0)
    ap.add_argument("--m", type=int, default=16, help="2**m Sobol points per term and replicate")
    ap.add_argument("--replicates", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260713)
    ap.add_argument("--q-scale", type=float, default=2.0)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        print("[PASS] Benites direct implementation self-test")
        return
    est = estimate(args.rs, args.zeta, args.m, args.replicates, args.seed, args.q_scale)
    print(f"BRM direct Eqs. (28)-(32): rs={args.rs:g}, zeta={args.zeta:g}")
    print(f"  Delta1 = {est.delta1_ry:+.8e} +- {est.delta1_err_ry:.2e} Ry")
    print(f"  Delta2 = {est.delta2_ry:+.8e} +- {est.delta2_err_ry:.2e} Ry")
    print(f"  correction = {est.correction_ry:+.8e} +- {est.correction_err_ry:.2e} Ry")
    print(f"  full kite  = {est.full_kite_mry:.4f} +- {est.full_kite_err_mry:.4f} mRy")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(asdict(est), indent=2) + "\n")
        print(f"  wrote {args.out}")


if __name__ == "__main__":
    _main()
