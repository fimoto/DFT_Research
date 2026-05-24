"""Imaginary-axis Lindhard function and RPA screening.

lindhard_shape(q, nu) = (1/2q) int_0^1 k dk log[((qk+q^2/2)^2+nu^2)/((qk-q^2/2)^2+nu^2)]
with lindhard_shape(q,0) = f_L(q/2) and lindhard_shape(0,0) = 1.  BRM's g equals 2 lindhard_shape.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from .constants import Q_TF_SQ_COEFF, spin_radii

_XK, _WK = leggauss(96)
kk = 0.5 * (_XK + 1.0)
wk = 0.5 * _WK


def h_shell(z):
    """h(z) = 1 - z arctan(1/z), the q -> 0 limit of lindhard_shape(q, qz).

    h(0) = 1; for large z the direct form is the difference of two numbers close
    to 1 and loses all digits (it returns 0 at z ~ 1e8), so the asymptotic
    series h = 1/(3z^2) - 1/(5z^4) + 1/(7z^6) - ... is used above z = 1e3.
    """
    z = np.asarray(z, float)
    out = np.ones(np.shape(z), float) if np.shape(z) else np.array(1.0)
    zz = np.atleast_1d(np.asarray(out, float) * 0 + z)
    res = np.empty(zz.shape, float)
    big = zz > 1e3
    zb = np.where(big, zz, 1.0)
    res_big = 1.0 / (3 * zb ** 2) - 1.0 / (5 * zb ** 4) + 1.0 / (7 * zb ** 6) - 1.0 / (9 * zb ** 8)
    zs = np.where(big, 1.0, np.maximum(zz, 1e-300))
    res_small = 1.0 - zs * np.arctan(1.0 / zs)
    res = np.where(big, res_big, res_small)
    res = np.where(zz == 0, 1.0, res)
    return res.reshape(np.shape(z)) if np.shape(z) else float(res[0])


def _P(w, nu, a, sgn):
    """Antiderivative of (w - sgn*a) log(w^2 + nu^2)."""
    r2 = w * w + nu * nu
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.where(r2 > 0, np.log(np.where(r2 > 0, r2, 1.0)), 0.0)
        at = np.where(nu > 0, nu * np.arctan(np.where(nu > 0, w / np.where(nu > 0, nu, 1.0), 0.0)), 0.0)
    return 0.5 * (r2 * L - w * w) - sgn * a * (w * L - 2 * w + 2 * at)


def _lindhard_series(q, z):
    """Small-q expansion of lindhard_shape(q, qz): h(z) + q^2 c2(z) + q^4 c4(z) + O(q^6).

    c2 = -1/(12(1+z^2)^2), c4 = (5z^2-1)/(240(1+z^2)^4), both derived from the
    k integral.  Relative error 5e-10 at q = 0.1 and 3e-7 at q = 0.3.
    """
    z2 = z * z
    return (h_shell(z)
            - q ** 2 / (12.0 * (1.0 + z2) ** 2)
            + q ** 4 * (5.0 * z2 - 1.0) / (240.0 * (1.0 + z2) ** 4))


def _lindhard_log1p(q, nu):
    """lindhard_shape by Gauss-Legendre on the log1p form, for q >= 3.

    lindhard_shape = (1/2q) int_0^1 k log1p(4 a q k/((qk-a)^2+nu^2)) dk, a = q^2/2.  For
    q > 2 the logarithmic peak k = q/2 lies outside the interval, so a fixed
    rule is exact to machine precision; log1p also removes the cancellation
    that destroys the closed form at large q (where it can even return a
    negative value).
    """
    q = np.asarray(q, float)[..., None]
    nu = np.asarray(nu, float)[..., None]
    a = 0.5 * q * q
    den = (q * kk - a) ** 2 + nu ** 2
    return 0.5 * np.sum(wk * kk * np.log1p(4.0 * a * q * kk / den), axis=-1) / q[..., 0]


def lindhard_shape(q, nu):
    """Imaginary-axis Lindhard function lindhard_shape(q, nu); lindhard_shape(q,0) = f_L(q/2), lindhard_shape(0,0) = 1.

    Three branches, each validated against a 30-digit reference to better than
    3e-7 and typically to 1e-9 or better:

      q < 0.15                      the small-q series (the closed form
                                    cancels here: at q = 1e-3 it is wrong by a
                                    factor 4e3),
      q >= 3, or nu >> q + q^2/2    Gauss-Legendre on the log1p form (the
                                    closed form loses all its digits: at
                                    q = 1e6 it returns -1.07e-9 for a quantity
                                    whose value is +1.33e-12),
      otherwise                     the exact antiderivative of the k integral.

    At q = 0 the two limits do not commute: lindhard_shape(0,0) = 1 is returned for nu = 0
    and lindhard_shape(0,nu) = 0 for nu > 0, which is the particle-number-conserving limit.

    Scalar and array input both work; the result has the broadcast shape.
    """
    qa = np.asarray(q, float)
    xa = np.abs(np.asarray(nu, float))
    scalar = (qa.ndim == 0 and xa.ndim == 0)
    qb, xb = np.broadcast_arrays(np.atleast_1d(qa), np.atleast_1d(xa))
    qb = np.array(qb, dtype=float, copy=True)
    xb = np.array(xb, dtype=float, copy=True)
    out = np.empty(qb.shape, float)

    # the log1p rule is used both for q >= 3 and, at any q, for nu far above
    # the peak scale, where the closed form again cancels
    # q -> 0 and nu -> 0 do not commute: lindhard_shape(0,0) = 1 is the static convention,
    # while at fixed nu > 0 particle-number conservation gives lindhard_shape ~ q^2/(3 nu^2)
    # -> 0.  The series branch below reproduces both through z = nu/q, but q = 0
    # exactly must be resolved explicitly.
    zero_q = qb <= 0.0
    if np.any(zero_q):
        out[zero_q] = np.where(xb[zero_q] > 0.0, 0.0, 1.0)
    small = (~zero_q) & (qb < 0.15)
    big = (~small) & (~zero_q) & ((qb >= 3.0) | (xb > 20.0 * (qb + 0.5 * qb ** 2)))
    mid = ~(small | big | zero_q)

    if np.any(small):
        qs = qb[small]
        z = np.where(qs > 0, xb[small] / np.where(qs > 0, qs, 1.0), 0.0)
        out[small] = _lindhard_series(qs, z)
    if np.any(mid):
        qm, xm = qb[mid], xb[mid]
        a = 0.5 * qm * qm
        out[mid] = ((_P(qm + a, xm, a, +1.0) - _P(a, xm, a, +1.0))
                    - (_P(qm - a, xm, a, -1.0) - _P(-a, xm, a, -1.0))) / (2 * qm ** 3)
    if np.any(big):
        out[big] = _lindhard_log1p(qb[big], xb[big])

    if scalar:
        return float(out[0])
    return out.reshape(np.broadcast(qa, xa).shape)


def f_L(x):
    """Static Lindhard shape, also the exchange self-energy shape F(x).

    f_L(x) = 1/2 + (1-x^2)/(4x) log|(1+x)/(1-x)| is even, with
    f_L(0)=1 and f_L(1)=1/2. Series avoid cancellation at either end;
    log1p and (1-x)(1+x) preserve accuracy near the cusp at x=1.
    """
    x = np.abs(np.asarray(x, float))
    out = np.empty_like(x)
    small = x < 1e-3
    xs = x[small]
    out[small] = 1.0 - xs ** 2 / 3 - xs ** 4 / 15 - xs ** 6 / 35
    # At large x the printed form subtracts 1/2 from a quantity tending to 1/2
    # to expose an O(x^-2) answer. The same function without that cancellation is
    #   f_L(x) = sum_{k>=1} x^{-2k}/((2k-1)(2k+1)) = 1/(3x^2) + 1/(15x^4) + ...
    big = x >= 30.0
    cusp = x == 1.0
    m = ~(small | big | cusp)
    xx = x[m]
    out[m] = (0.5 + (1 - xx) * (1 + xx) / (4 * xx)
              * np.log1p(2 * np.minimum(xx, 1.0) / np.abs(1 - xx)))
    if np.any(big):
        t = (1.0 / x[big]) ** 2
        acc = np.zeros_like(t)
        for k in range(10, 0, -1):
            acc = acc * t + 1.0 / ((2 * k - 1) * (2 * k + 1))
        out[big] = acc * t
    out[cusp] = 0.5
    return out


def brm_lindhard_g(q, nu):
    """BRM Eq. (14) normalization, g = 2 lindhard_shape, g(0,0) = 2."""
    return 2.0 * lindhard_shape(q, nu)


def lindhard_shape_spin(q, nu, zeta=0.0):
    """Two-component response in total-density units, Eq. (polarized-response)."""
    out = 0.0
    for s, x in zip((+1, -1), spin_radii(zeta)):
        if x < 1e-14:
            continue
        out = out + 0.5 * x * lindhard_shape(np.asarray(q) / x, np.asarray(nu) / x ** 2)
    return out


def screening_strength(q, nu, rs, zeta=0.0):
    """c = -v Pi_0 = (q_tf_sq/q^2) lindhard_shape_spin >= 0."""
    q = np.asarray(q, float)
    return Q_TF_SQ_COEFF * rs / q ** 2 * lindhard_shape_spin(q, nu, zeta)


def coupling_screening_factor(c):
    """R(c) = 2 int_0^1 lam dlam/(1+lam c), c >= 0.

    This dimensionless screening factor is 2 Wbar/v, not Wbar itself.
    R(0)=1 and R(inf)=0.
    """
    c = np.asarray(c, float)
    if np.any(c < 0):
        raise ValueError("coupling_screening_factor: c must be >= 0 (c = -v Pi_0)")
    small = (c >= 0) & (c < 1e-6)
    out = np.empty_like(c)
    cs = c[small]
    out[small] = 1.0 - 2 * cs / 3 + cs ** 2 / 2 - 2 * cs ** 3 / 5
    cb = c[~small]
    with np.errstate(invalid="ignore", divide="ignore"):
        out[~small] = np.where(np.isinf(cb), 0.0, 2.0 * (cb - np.log1p(cb)) / cb ** 2)
    return out


def noninteracting_response(q, nu, rs, zeta=0.0):
    """Pi_0 = chi_0 = -q_TF^2 L in the common inverse-interaction unit.

    Interaction unit v0=4*pi*e^2/k_F^2; response unit is 1/v0.
    This is the irreducible free response, not the positive Lindhard shape.
    """
    if np.any(np.asarray(rs) < 0):
        raise ValueError("rs must be nonnegative")
    return -Q_TF_SQ_COEFF * np.asarray(rs) * lindhard_shape_spin(q, nu, zeta)


def _screening_arguments(q, rs, coupling):
    """The Coulomb kernels are evaluated for q>0 and 0<=lambda<=1."""
    q = np.asarray(q, float)
    coupling = np.asarray(coupling, float)
    if np.any(q <= 0) or np.any(np.asarray(rs) < 0):
        raise ValueError("screening kernels require q>0 and rs>=0")
    if np.any((coupling < 0) | (coupling > 1)):
        raise ValueError("coupling must lie in [0,1]")
    return q, coupling


def rpa_dielectric(q, nu, rs, coupling=1.0, zeta=0.0):
    """epsilon_lambda = 1 - lambda*v*Pi_0 on the imaginary axis."""
    q, coupling = _screening_arguments(q, rs, coupling)
    return 1.0 + coupling * screening_strength(q, nu, rs, zeta)


def screened_interaction(q, nu, rs, coupling=1.0, zeta=0.0):
    """W_lambda = lambda*v/epsilon_lambda, including its bare coupling.

    v(q)=q^-2 in the declared interaction unit. In particular W_0=0,
    W_lambda/lambda -> v as lambda->0, and W_lambda->lambda*v at
    large imaginary frequency. The lambda=0 subscript here is a coupling
    endpoint, not the common G0W0 label for a non-self-consistent W.
    """
    q, coupling = _screening_arguments(q, rs, coupling)
    return coupling / q ** 2 / rpa_dielectric(q, nu, rs, coupling, zeta)


def coupling_averaged_interaction(q, nu, rs, zeta=0.0):
    """Wbar = int_0^1 W_lambda d lambda = v*R(c)/2; Wbar->v/2 unscreened."""
    q, _ = _screening_arguments(q, rs, 1.0)
    return coupling_screening_factor(screening_strength(q, nu, rs, zeta)) / (2.0 * q ** 2)


def rpa_density_response(q, nu, rs, coupling=1.0, zeta=0.0):
    """Reducible response chi_lambda = Pi_0/epsilon_lambda, in 1/v0 units."""
    return noninteracting_response(q, nu, rs, zeta) / rpa_dielectric(q, nu, rs, coupling, zeta)


def h_small_q(z):
    """Small-q limit lindhard_shape(q, qz) -> h(z); stable at large z (see h_shell)."""
    return h_shell(z)


def exchange_self_energy_shape(k):
    """Exchange self-energy shape F(k); Sigma_x = -(q_tf_sq/2) F(k).  F(0)=1, F(1)=1/2.

    This is mathematically the same shape as f_L, so it shares that function's
    small- and large-momentum series. In particular F(k) ~ 1/(3k^2) remains
    positive at large k, where direct subtraction loses all significant digits.
    """
    return f_L(k)
