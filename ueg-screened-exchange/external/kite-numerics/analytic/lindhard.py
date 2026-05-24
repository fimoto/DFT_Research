"""
Imaginary-axis Lindhard function and the coupling-constant-integrated
RPA screening ratio, in the suite's units (k_F = 1, energies with
eps_k = k^2/2, so the particle-hole energies delta = q p_z + q^2/2 of the
Monte Carlo kernel and the frequency argument here share the same scale).

  chi0(q, i xi) = -(1/pi^2) * G(q, xi),
  G(q, xi) = (1/2q) Int_0^1 dk k ln[ ((qk+q^2/2)^2+xi^2)
                                    /((qk-q^2/2)^2+xi^2) ],
  G(q, 0) = f_L(q/2)                (static Lindhard, verified),
  G(q, xi) -> q^2/(3 xi^2) * ... at large xi (f-sum behaviour).

Dimensionless screening strength (density enters only through q_TF):

  c(q, xi; rs) = v_q |chi0| = (q_TF^2/q^2) G(q, xi),
  q_TF^2 = 4 alpha_L rs/pi,   alpha_L = (4/9pi)^(1/3),

which reduces at xi = 0 to the static factor q_TF^2 f_L(q/2)/q^2 used in
mc/run_scan.py.  The adiabatic-connection (coupling-constant) integral of
the screened-to-bare line ratio is closed form:

  R(q, xi; rs) = 2 Int_0^1 lam dlam / (1 + lam c) = (2/c^2)[c - ln(1+c)],

normalized so that R -> 1 as rs -> 0 (bare limit reproduced exactly).
"""
import numpy as np

ALPHA_L = (4/(9*np.pi))**(1/3.)
_KN, _KW = np.polynomial.legendre.leggauss(64)
_KN = 0.5*(_KN + 1.0)          # nodes on (0,1)
_KW = 0.5*_KW


def h_shell(z):
    """h(z) = 1 - z arctan(1/z), stable at large z (series above z = 1e3)."""
    z = np.asarray(z, float)
    zz = np.atleast_1d(z)
    big = zz > 1e3
    zb = np.where(big, zz, 1.0)
    zs = np.where(big, 1.0, np.maximum(zz, 1e-300))
    res = np.where(big,
                   1.0 / (3 * zb ** 2) - 1.0 / (5 * zb ** 4) + 1.0 / (7 * zb ** 6),
                   1.0 - zs * np.arctan(1.0 / zs))
    res = np.where(zz == 0, 1.0, res)
    return res.reshape(np.shape(z)) if np.shape(z) else float(res[0])


def f_lindhard(x):
    x = np.asarray(x, float)
    out = np.full_like(x, 0.5)
    small = x < 1e-8
    out[small] = 1.0                      # f_L(0) = 1 (the direct formula cancels to 0.5 at x=0)
    m = (np.abs(1 - x) > 1e-12) & ~small
    out[m] += (1 - x[m]**2)/(4*x[m])*np.log(np.abs((1 + x[m])/(1 - x[m])))
    return out


def G(q, xi):
    """Reduced imaginary-axis Lindhard, broadcasting over xi (array).
    Closed form of the k-integral (elementary antiderivatives), regular
    at xi = 0 including the log point k = q/2."""
    xi = np.atleast_1d(np.asarray(xi, float))
    a = q*q/2.0

    def I1(s):
        s2 = s*s + xi*xi
        return 0.5*((s2)*np.log(np.maximum(s2, 1e-300)) - s*s)

    def I0(s):
        s2 = s*s + xi*xi
        at = np.where(xi > 0, xi*np.arctan(s/np.maximum(xi, 1e-300)), 0.0)
        return s*np.log(np.maximum(s2, 1e-300)) - 2.0*s + 2.0*at

    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        Fp = (I1(a + q) - a*I0(a + q)) - (I1(a) - a*I0(a))
        Fm = (I1(q - a) + a*I0(q - a)) - (I1(-a) + a*I0(-a))
        out = (Fp - Fm)/(2.0*q**3)

    # large xi: the closed form suffers catastrophic cancellation
    # (differences of ~xi^2 ln xi terms); switch to a stable log1p
    # quadrature, which is smooth (and hence GL-exact) precisely there.
    # the same quadrature is smooth for every xi once q > 2 (no log point in k<=1);
    # the closed form loses all digits for q >~ 1e3, so switch to it for q > 100 as well.
    # The threshold lies above the largest node of the production grids (q = 60), so
    # the released grid records are unchanged bit for bit.
    big = (xi > 2.0*(q + a)) | np.full_like(xi, q > 100.0, dtype=bool)
    if np.any(big):
        xb = xi[big][None, :]
        k = _KN[:, None]
        den = (q*k - a)**2 + xb**2
        val = ((_KW[:, None]*k*np.log1p(4.0*a*q*k/den)).sum(axis=0))/(2*q)
        out[big] = val
    return out


def c_screen(q, xi, rs, zeta=0.0):
    """Dimensionless screening strength v_q |chi0|; for zeta != 0 the
    two-component dielectric (spin-sigma Fermi radii x = (1+s*zeta)^(1/3),
    chi0 scaling chi_x(q,xi) = x chi_1(q/x, xi/x^2))."""
    pref = (4*ALPHA_L*rs/np.pi)/q**2
    if zeta == 0.0:
        return pref*G(q, xi)
    tot = 0.0
    for sgn in (+1.0, -1.0):
        x = (1.0 + sgn*zeta)**(1/3.)
        if x > 1e-12:
            tot = tot + 0.5*x*G(q/x, np.asarray(xi, float)/x**2)
    return pref*tot


def R_ac(q, xi, rs, zeta=0.0):
    """lambda-integrated screened/bare ratio, R -> 1 as rs -> 0."""
    c = c_screen(q, xi, rs, zeta)
    out = np.ones_like(c)
    m = c > 1e-8
    out[m] = 2.0/c[m]**2*(c[m] - np.log1p(c[m]))
    out[~m] = 1.0 - 2*c[~m]/3.0
    return out


def R_static_noac(q, rs):
    """Static screened/bare ratio WITHOUT the AC lambda integral
    (the 'static RPA-screened kite' weight of run_scan.py)."""
    return 1.0/(1.0 + c_screen(q, np.array([0.0]), rs)[0])


if __name__ == "__main__":
    ok = True
    for q in (0.3, 1.0, 1.9, 3.0):
        g0 = G(q, [0.0])[0]
        ok &= abs(g0 - f_lindhard(np.array([q/2]))[0]) < 1e-10
    print(f"[{'PASS' if ok else 'FAIL'}] G(q,0) = f_L(q/2)")
    # f-sum-type large-xi decay: G ~ q^2/(3 xi^2) for q << xi.
    q, xi = 0.5, 200.0
    g = G(q, [xi])[0]
    pred = q**2/(3.0*xi**2)
    rel = abs(g/pred - 1.0)
    print(f"[{'PASS' if rel < 1e-4 else 'FAIL'}] "
          f"G coefficient at large xi (relative error {rel:.2e})")
    g2 = G(q, [2*xi])[0]
    print(f"[{'PASS' if abs(g/g2 - 4) < 0.05 else 'FAIL'}] "
          f"G ~ 1/xi^2 at large xi (ratio {g/g2:.3f})")
    r = R_ac(1.0, np.array([0.0, 1.0]), 1e-9)
    print(f"[{'PASS' if np.allclose(r, 1, atol=1e-8) else 'FAIL'}] "
          f"R -> 1 in the bare limit")
