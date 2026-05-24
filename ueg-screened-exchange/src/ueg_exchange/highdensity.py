"""High-density window analysis: the rs log^2 rs term and its cancellation (JPA Sec. 5.2 and Appendix G).

Window sqrt(q_tf_sq) << q << 1: R(c)-1 ~ -2c/3, X(q,qz) = A(z) log(1/q) + B(z).
    calA = int (dz/pi) h(z) A(z) = 4 pi int h^2 = (4 pi^2/3)(1 - log 2)
    C_ln2 = -(pi A0 / 3) a calA = -alpha_L (1 - log 2) / (6 pi^3)
Cancellation (the exchange-cancellation identity): Pi_vtx^(1)^log = -(q_tf_sq^2/4) alpha(z) log(1/q) and
Pi_SE^(1)^log = +(q_tf_sq^2/4) alpha(z) log(1/q), equal and opposite for every z.
Spin dependence: C_ln2(zeta) = C_ln2(0) * rho(zeta), rho in closed form.
"""
import numpy as np
from scipy.integrate import quad
from .constants import A0, Q_TF_SQ_COEFF, ALPHA_L, C_LN2
from .lindhard import h_small_q, f_L, exchange_self_energy_shape
from .lorentzian import A_shell


def int_h2():
    """int_0^inf h(z)^2 dz = (pi/3)(1 - log 2)."""
    return quad(lambda z: h_small_q(z) ** 2, 0, np.inf, limit=300)[0]


def calA():
    """calA = int (dz/pi) h A, numerically and in closed form."""
    num = quad(lambda z: h_small_q(z) * A_shell(z), 0, np.inf, limit=400)[0] / np.pi
    closed = 4 * np.pi ** 2 / 3 * (1 - np.log(2))
    return num, closed


def C_ln2_from_calA():
    """C_ln2 = -(pi A0/3) a calA."""
    return -np.pi * A0 / 3.0 * Q_TF_SQ_COEFF * calA()[1]


def eta_velocity(q, rs, c_p=1.0, lam=1.0):
    """delta -> q c_p (1 + eta), eta = lam q_tf_sq/4 log(2/(q|c_p|))."""
    q_tf_sq = Q_TF_SQ_COEFF * rs
    return lam * q_tf_sq / 4.0 * np.log(2.0 / (q * abs(c_p)))


def alpha_shell(z):
    """alpha(z) = h + z h' = 1 - 2 z arctan(1/z) + z^2/(1+z^2) = A(z)/(8 pi^2).

    Shares the stable large-z branch with lorentzian.A_shell (the direct form
    cancels to zero above z ~ 1e8).
    """
    from .lorentzian import A_shell
    return A_shell(z) / (8.0 * np.pi ** 2)


def response_log_coefficients(z, rs):
    """Coefficients of log(1/q) in Pi_vtx^(1) and Pi_SE^(1) at small q, fixed z (exchange cancellation).

    NOTE: this encodes the identification made in the manuscript (equal and
    opposite window coefficients); returning (-c, +c) is a statement of that
    result, not an independent diagrammatic evaluation of the two diagrams.
    """
    q_tf_sq = Q_TF_SQ_COEFF * rs
    c = q_tf_sq ** 2 / 4.0 * alpha_shell(z)
    return -c, +c


def Phi_spin(a):
    """Phi_s(a) = 2a^2 - a + (1-2a^3) log(1+a) + 2a^3 log a,  Phi_s(1) = 1 - log 2, Phi_s(0) = 0.

    The printed form cancels at both ends, so three branches are used.

    * 0 < a <= 1/2:  Phi_s = (3/2) a^2 + a^3/3 + 2 a^3 log a
                             + sum_{k>=4} (-1)^(k+1) (1/k + 2/(k-3)) a^k
    * 1/2 < a < 2:   the printed form.
    * a >= 2, t = 1/a:
          Phi_s = log a - 2/3 + sum_{n>=1} (-1)^(n+1) 3(n+1)/(n(n+3)) t^n
                = log a - 2/3 + 3/(2a) - 9/(10a^2) + 2/(3a^3) - 15/(28a^4)
                  + 9/(20a^5) - 7/(18a^6) + O(a^-7).

    Both series are exact expansions of log(1+a) and log(1+1/a), summed to
    _PHI_NTERMS terms (truncation below 2^-64 at the switch points). Checked
    against 90-digit evaluation of the printed form on 1e-12 <= a <= 1e15,
    including both sides of each switch point: relative error below 1e-14
    (tests/test_all.py::test_phi_spin_against_high_precision).
    """
    a = np.asarray(a, float)
    out = np.zeros_like(a)
    small = (a > 0) & (a <= 0.5)
    mid = (a > 0.5) & (a < 2.0)
    big = a >= 2.0
    if np.any(small):
        x = a[small]
        acc = np.zeros_like(x)
        for k in range(_PHI_NTERMS + 3, 3, -1):
            acc = acc * x + (-1) ** (k + 1) * (1.0 / k + 2.0 / (k - 3))
        out[small] = 1.5 * x ** 2 + x ** 3 / 3.0 + 2 * x ** 3 * np.log(x) + acc * x ** 4
    if np.any(mid):
        x = a[mid]
        out[mid] = 2 * x ** 2 - x + (1 - 2 * x ** 3) * np.log1p(x) + 2 * x ** 3 * np.log(x)
    if np.any(big):
        x = a[big]
        t = 1.0 / x
        acc = np.zeros_like(x)
        for n in range(_PHI_NTERMS, 0, -1):
            acc = acc * t + (-1) ** (n + 1) * 3.0 * (n + 1) / (n * (n + 3))
        out[big] = np.log(x) - 2.0 / 3.0 + acc * t
    return out


_PHI_NTERMS = 64


def I_theta(theta):
    """int_0^inf h(theta z) alpha(z) dz = (pi/(6a)) Phi_s(a), a = 1/theta."""
    a = 1.0 / theta
    return np.pi / (6 * a) * Phi_spin(a)


def rho_zeta(zeta):
    """C_ln2(zeta)/C_ln2(0) = (1/(4(1-log2))) sum_{s,s'} x_s^2 Phi_s(x_s'/x_s)."""
    from .constants import spin_radii
    xs = spin_radii(zeta)
    tot = 0.0
    for xa in xs:
        for xb in xs:
            if xa > 1e-14:
                tot += xa ** 2 * Phi_spin(xb / xa)
    return tot / (4 * (1 - np.log(2)))


def C_ln2_zeta(zeta):
    return C_LN2 * rho_zeta(zeta)


def self_energy_response_log_slope(z, q_tf_sq=1.0, qs=(1e-3, 3e-4, 1e-4, 3e-5), n=800):
    """Coefficient of log(1/q) in Pi_SE^(1) at fixed z, from its defining integral.

    This does NOT use the identification of the exchange-cancellation identity.  It evaluates

        Pi_SE^(1) = pi^2 q_tf_sq * 2 int d^3p/(2pi)^3 Delta(p;q) [Sigma_x(p+q)-Sigma_x(p)]
                                                 / (i nu - delta)^2

    on the small-q shell -- an approximation to the finite-q integral, not the
    integral itself -- where d^3p Delta f -> 2 pi q int dc c f, delta = q c and
    Sigma_x(k) = -(q_tf_sq/2) F(k) with F(1) = 1/2, and fits the result to
    a + b log(1/q). Multiplication by pi^2 q_tf_sq converts the
    loop measure to the common inverse-interaction response unit. The
    returned slope is compared in verify_analytics.py with +q_tf_sq^2 alpha(z)/4;
    the two agree to about 6e-7, which is an independent check of the frequency
    dependence and of the sign that fixes the cancellation against Pi_vtx^(1).
    """
    from numpy.polynomial.legendre import leggauss
    c_, w = leggauss(n)
    qs = np.asarray(qs, float)
    vals = []
    for q in qs:
        nu = q * z
        d = q * c_
        dS = -(q_tf_sq / 2.0) * (sigma_x_shape(np.abs(1.0 + q * c_)) - 0.5)
        vals.append((1.0 / (2.0 * np.pi ** 2))
                    * np.sum(w * np.sign(c_) * q * np.abs(c_) * dS
                             * np.real(1.0 / (1j * nu - d) ** 2)))
    A = np.stack([np.log(1.0 / qs), np.ones_like(qs)], 1)
    return float(np.pi ** 2 * q_tf_sq * np.linalg.lstsq(A, np.asarray(vals), rcond=None)[0][0])


def sigma_x_shape(k):
    """F(k) with Sigma_x(k) = -(q_tf_sq/2) F(k); F(0)=1, F(1)=1/2."""
    return exchange_self_energy_shape(k)


def c1_surface(n=200, nphi=160):
    """RETIRED diagnostic; do NOT use to verify a physical coefficient.

    This evaluates the angular integrand I_1 of an earlier draft, whose
    distributional expansion carried a sign error in the delta'(r-1) term and
    whose F(1,1) is valid on only one of the two cross-ball sectors.  Its
    vanishing therefore says nothing about the q^2 coefficient of X_c.  The
    current manuscript establishes c_1 = 0 using a regular fixed-domain
    representation instead; see JPA Appendix D and reviewer Checks 16, 25,
    36 and 40. Kept only to reproduce the earlier diagnostic.
    """
    from numpy.polynomial.legendre import leggauss
    c, wc = leggauss(n)
    ph, wph = leggauss(nphi)
    ph = np.pi * (ph + 1.0)
    wph = np.pi * wph
    tot = 0.0
    for cp, wp in zip(c, wc):
        sp = np.sqrt(max(1 - cp * cp, 0.0))
        for ck, wk in zip(c, wc):
            if cp * ck >= 0:
                continue
            sk = np.sqrt(max(1 - ck * ck, 0.0))
            cg = cp * ck + sp * sk * np.cos(ph)
            F11 = -0.5 / ((cg - 1.0) * (ck - cp))
            ap, bp, gp = cp, (1 - cp * cp) / 2, cp * cp / 2
            ak, bk, gk = ck, (1 - ck * ck) / 2, ck * ck / 2
            R2 = 2 - 2 * cg
            dFp = 2 * F11 - F11 * (2 - 2 * cg) / R2 - F11 * abs(cp) / (abs(cp) + abs(ck))
            dFk = 2 * F11 - F11 * (2 - 2 * cg) / R2 - F11 * abs(ck) / (abs(cp) + abs(ck))
            I1 = (ap * bk + ak * bp) * F11 + ap * gk * dFk + ak * gp * dFp
            tot += wp * wk * np.sum(wph * I1) * 2
    return float(tot)


def yukawa_pole(rs):
    """Auxiliary Yukawa pole of the static screened line.

    Continues W_stat(q) = 1/(q^2 + qTF^2 f_L(q/2)) to its screened pole
    q^2 = -kappa^2 and returns (kappa, Z) with Z the residue there, so that
    W_stat = Z/(q^2+kappa^2) + Delta(q); the remainder Delta carries the 2kF
    nonanalyticity of f_L.  On the imaginary axis
    f_L(i y) = 1/2 + (1+y^2) arctan(y)/(2y).
    """
    from scipy.optimize import brentq
    t = screening_wavevector_squared(rs)
    fLi = lambda y: 0.5 + (1 + y * y) * np.arctan(y) / (2 * y)
    root = lambda kk: -kk * kk + t * fLi(kk / 2.0)
    # kappa grows like pi t / 8 at strong screening, so a fixed upper bracket of
    # 50 fails above r_s ~ 1e2. Expand the bracket until the sign changes.
    lo, hi = 1.0e-8, max(50.0, np.pi * t / 2.0)
    for _ in range(80):
        if root(lo) * root(hi) < 0.0:
            break
        hi *= 2.0
    else:
        raise RuntimeError(f"could not bracket the screened pole at rs={rs}")
    k = brentq(root, lo, hi, xtol=1e-14, rtol=8.9e-16)
    # residue from the analytic derivative rather than a fixed-step difference:
    #   d f_L(i k/2)/d(k^2) = f'(k/2)/(4k),
    #   f'(y) = arctan y + 1/(2y) - (1+y^2) arctan y/(2 y^2)
    y = 0.5 * k
    fp = np.arctan(y) + 1.0 / (2 * y) - (1 + y * y) * np.arctan(y) / (2 * y * y)
    return float(k), float(1.0 / (1.0 - t * fp / (4.0 * k)))


def screening_wavevector_squared(rs):
    """qTF^2 in units of kF^2: 4 alpha_L rs / pi."""
    return 4.0 * (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0) * rs / np.pi


def friedel_share(rs, qgrid=None, qmax=60.0, n=400001):
    """Signed fraction of K[W_stat] carried by the remainder of eq:yukawa-friedel.

    This is the quantity quoted in the manuscript.  It uses the sampled kernel
    X_c (figs/xc_grid.npy, passed as qgrid=(q, X) or loaded by the caller) with
    the analytic endpoints, and integrates the remainder WITH ITS SIGN.  For a
    quick diagnostic that needs no grid, see friedel_share_proxy, which uses an
    endpoint-matched interpolation and the absolute value; the two differ by a
    few parts in a hundred and the proxy does not integrate to D_0.
    """
    from numpy import trapezoid as tz
    k, Z = yukawa_pole(rs)
    t = screening_wavevector_squared(rs)
    c0 = np.pi / 8 * (2 * np.log(2) - 1)
    Cinf = 2 * np.pi / 9
    if qgrid is None:
        import os
        g = np.load(os.path.join(os.path.dirname(__file__), "..", "..", "figs",
                                 "xc_grid.npy"))
        qgrid = (g[0], g[1])
    qg, Xg = qgrid
    q = np.linspace(1e-4, qmax, n)
    Xc = np.interp(q, qg, Xg, left=0.0, right=0.0)
    Xc = np.where(q < qg[0], c0 * q, Xc)
    Xc = np.where(q > qg[-1], Cinf / np.maximum(q, 1e-9) ** 4, Xc)
    W = 1.0 / (q * q + t * f_L_static(q / 2))
    Y = Z / (q * q + k * k)
    w = Xc * q * q
    return float(tz((W - Y) * w, q) / tz(W * w, q))


def friedel_share_proxy(rs, qmax=60.0, n=400001):
    """Grid-free diagnostic: absolute-value ratio with an interpolated kernel.

    The interpolation c0 q/(1+(c0/Cinf) q^5) matches both endpoints but
    integrates to 0.1846, not D_0 = 0.2499, so this is a diagnostic only.
    """
    from numpy import trapezoid as tz
    k, Z = yukawa_pole(rs)
    t = screening_wavevector_squared(rs)
    q = np.linspace(1e-4, qmax, n)
    c0 = np.pi / 8 * (2 * np.log(2) - 1)
    Cinf = 2 * np.pi / 9
    Xc = c0 * q / (1 + (c0 / Cinf) * q ** 5)
    W = 1.0 / (q * q + t * f_L_static(q / 2))
    Y = Z / (q * q + k * k)
    w = Xc * q * q
    return float(tz(np.abs(W - Y) * w, q) / tz(W * w, q))


def f_L_static(x):
    """Static Lindhard shape f_L, with f_L(0)=1 and f_L(1)=1/2."""
    return f_L(x)
