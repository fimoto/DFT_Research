"""The constant-width Lorentzian kernel b_Lambda and its radial kernel.

b_Lambda(q, i nu) = (1/q^2) Lambda^2 / (Lambda^2 + z^2),   z = nu/q.

Radial energy density e_Lambda(q), including the energy prefactor (JPA Sec. 5 and Appendices G, H):
  (i)   small q, at fixed Lambda:  e_Lambda(q) = kappa(Lambda) q log(1/q) + beta q + o(q)
  (ii)  large q, s_lambda = Lambda/q fixed: e_Lambda(q) = (3/pi^3) p_loc(s_lambda)/q^2 + O(q^-3)
  (iii) plasmon regime (eta = q Lambda fixed): q log q demoted to q^2 log q
"""
import numpy as np
from scipy.integrate import quad
from .constants import A0


def p_loc(s_lambda):
    """Large-q local frequency kernel, p_loc(s_lambda) = 2 pi s_lambda / (1 + 2 s_lambda)^2; p_loc(0)=p_loc(inf)=0."""
    s_lambda = np.asarray(s_lambda, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = 2.0 * np.pi * s_lambda / (1.0 + 2.0 * s_lambda) ** 2
    return np.where(np.isinf(s_lambda), 0.0, out)


def radial_energy_large_q(q, Lam):
    """JPA large-momentum expansion (app-eq:largeq-block): (3/pi^3) p_loc(Lambda/q)/q^2 = (6/pi^2) Lam/[q (q+2Lam)^2]."""
    q = np.asarray(q, float)
    return 6.0 / np.pi ** 2 * Lam / (q * (q + 2.0 * Lam) ** 2)


def kappa(Lam):
    """JPA fixed-width small-momentum expansion (app-eq:kappa): kappa(Lambda) = (3/2pi^2) Lam int_0^1 c^2 dc/(c+Lam)^2; kappa(0)=kappa(inf)=0.

    For Lam > 1e3 the bracket is evaluated by its series to avoid cancellation.
    """
    Lam = np.asarray(Lam, float)
    out = np.zeros_like(Lam)
    fin = np.isfinite(Lam) & (Lam > 0)
    L = Lam[fin]
    big = L > 1e3
    Ls = np.where(big, 1.0, L)
    brk = 1.0 - 2.0 * Ls * np.log1p(1.0 / Ls) + Ls / (1.0 + Ls)
    Lb = np.where(big, L, 1.0)
    brk = np.where(big, 1.0 / (3 * Lb ** 2) - 1.0 / (2 * Lb ** 3) + 3.0 / (5 * Lb ** 4), brk)
    out[fin] = 3.0 / (2 * np.pi ** 2) * L * brk
    return out


def kappa_quadrature(Lam):
    """Independent check of kappa by direct quadrature."""
    return 3.0 / (2 * np.pi ** 2) * Lam * quad(
        lambda c: c ** 2 / (c + Lam) ** 2, 0.0, 1.0)[0]


def A_shell(z):
    """Coincidence kernel A(z) = 8 pi^2 [h(z) + z h'(z)], Eq. (Az).

    A(z) = 8 pi^2 [1 - 2 z arctan(1/z) + z^2/(1+z^2)];  int_0^inf A dz = 0.
    For z > 1e3 the large-z series alpha = -1/(3z^2) + 3/(5z^4) - 5/(7z^6) + ...
    is used, which follows from 2z arctan(1/z) = 2 - 2/(3z^2) + 2/(5z^4) - ...
    and z^2/(1+z^2) = 1 - 1/z^2 + 1/z^4 - ...
    """
    z = np.asarray(z, float)
    zz = np.maximum(z, 1e-300)
    zs = np.where(zz > 1e3, 1.0, zz)
    direct = 1.0 - 2.0 * zs * np.arctan(1.0 / zs) + zs ** 2 / (1.0 + zs ** 2)
    zb = np.where(zz > 1e3, zz, 1.0)
    series = -1.0 / (3.0 * zb ** 2) + 3.0 / (5.0 * zb ** 4) - 5.0 / (7.0 * zb ** 6) + 7.0 / (9.0 * zb ** 8)
    return 8.0 * np.pi ** 2 * np.where(zz > 1e3, series, direct)

def Xi_closed(z):
    """Xi(z) = z[Ci(z) sin z + (pi/2 - Si(z)) cos z], the Lorentzian screening kernel.

    Defining integral: Xi(z) = z int_0^inf e^{-zt}/(1+t^2) dt.
    """
    from scipy.special import sici
    z = np.asarray(z, float)
    out = np.empty_like(z, dtype=float)
    big = z > 1.0e6
    sml = z <= 0.0
    mid = ~(big | sml)
    if np.any(mid):
        si, ci = sici(z[mid])
        out[mid] = z[mid] * (ci * np.sin(z[mid]) + (np.pi / 2 - si) * np.cos(z[mid]))
    # Xi(0) = 0 by continuity; at large z the Si difference cancels, so use
    # Xi(z) = 1 - 2/z^2 + 24/z^4 - ...
    out[sml] = 0.0
    if np.any(big):
        zz = z[big]
        out[big] = 1.0 - 2.0 / zz ** 2 + 24.0 / zz ** 4
    return out if out.ndim else float(out)


def K_mellin(s):
    """K(s) = int_0^inf z^{s-1} Xi(z) dz = -(pi/2) Gamma(s+1)/sin(pi s/2).

    The defining integral converges only on -1 < Re s < 0, since Xi -> 1 at
    infinity; elsewhere the closed form is the analytic continuation.
    """
    from scipy.special import gamma
    return -0.5 * np.pi * gamma(s + 1.0) / np.sin(np.pi * s / 2.0)


def b0_closed():
    """b0 = 2 int_0^inf j1(t/2)^2 (1 - log t)/t dt = gamma/2 - 3/8.

    From the Weber-Schafheitlin formula (DLMF 10.22.57) at mu = nu = 3/2,
    differentiated at s = 0.
    """
    return 0.5 * np.euler_gamma - 0.375


def S_asymptote(Lam):
    """Large-Lambda form of S: 1 + (pi/(4 D0 Lam))(log Lam - 1/4).

    Uses C_0^Lambda = -pi/(16 D0) = -C_1^Lambda/4, which follows from
    b0 = gamma/2 - 3/8 together with K'(1) = -(pi/2)(1 - gamma).
    """
    D0 = np.pi ** 3 / 18 * np.log(2) - np.pi / 4 * 1.2020569031595942854
    return 1.0 + np.pi / (4 * D0 * Lam) * (np.log(Lam) - 0.25)
