"""SI Sec. S1 and SI Table S2: d-dimensional constants and the long-wavelength
Lindhard shape h_d(z).

Reproduces: the entries of SI Table S2 (exact), the hypergeometric form of h_d
against the angular average for d = 2, 3, 3.5, 4, 5 (18 digits), the closed forms
h_2, h_3, h_4 of main-text Eq. (10), the integrals of Eq. (S3), and the q^2
coefficient g_2(z) of the two-dimensional Lindhard function, Eq. (S4), with the
vanishing coefficients of q and q^3, from a 40-digit evaluation of the exact free
response at z = 0.3, 1 and 3.  Also the three-dimensional g_3 of main-text Sec. 2.9
from the shell representation, SI Eq. (S79).
"""
import mpmath as mp, sympy as sp
from common import Recorder

R = Recorder("c01_lindhard_d", "d-dimensional constants and the Lindhard shape h_d (SI Sec. S1, Table S2)")

# ---------------------------------------------------------------- SI Table S2 (exact)
d, rs = sp.symbols('d r_s', positive=True)
A = lambda d: 2 * sp.pi ** sp.Rational(d, 2) / sp.gamma(sp.Rational(d, 2))
C = lambda d: (4 * sp.pi) ** sp.Rational(d - 1, 2) * sp.gamma(sp.Rational(d - 1, 2))
alpha = lambda d: (2 ** (d - 1) * sp.gamma(sp.Rational(d, 2) + 1) ** 2) ** sp.Rational(-1, d)
B = lambda d: 2 * A(d) * C(d) / (2 * sp.pi) ** d
for dd, tab in [(2, dict(A=2 * sp.pi, C=2 * sp.pi, alpha=1 / sp.sqrt(2), B=2)),
                (3, dict(A=4 * sp.pi, C=4 * sp.pi, alpha=(sp.Rational(4, 9) / sp.pi) ** sp.Rational(1, 3), B=4 / sp.pi))]:
    for key, f in (('A', A), ('C', C), ('alpha', alpha), ('B', B)):
        R.exact("SI Table S2", f"{key}_{dd}", sp.simplify(f(dd) - tab[key]) == 0, True)
    qTF = (B(dd) * alpha(dd) * rs) ** sp.Rational(1, dd - 1)
    E0 = 2 / (alpha(dd) * rs) ** 2
    pref = sp.simplify(sp.Rational(dd, 8) / sp.pi * E0 * qTF ** (dd + 1))
    want = 2 * sp.sqrt(2) * rs / sp.pi if dd == 2 else 12 / sp.pi ** 3
    R.exact("SI Table S2", f"(d/8pi) E_0 q_TF^(d+1) for d={dd}", sp.simplify(pref - want) == 0, True)
    if dd == 2:
        R.exact("SI Table S2", "q_TF = sqrt(2) r_s, E_0 = 4/r_s^2 (d=2)", (sp.simplify(qTF - sp.sqrt(2) * rs) == 0, sp.simplify(E0 - 4 / rs ** 2) == 0), (True, True))
    else:
        R.exact("SI Table S2", "E_0 = 32/(pi^2 q_TF^4) (d=3)", sp.simplify(E0 - 32 / (sp.pi ** 2 * qTF ** 4)) == 0, True)

# ---------------------------------------------------------------- h_d: hypergeometric vs angular average
mp.mp.dps = 40
def h_hyp(d, z):
    return 1 - mp.hyp2f1(1, mp.mpf(1) / 2, mp.mpf(d) / 2, -1 / z ** 2)
def h_avg(d, z):
    d = mp.mpf(d)
    w = lambda m: (1 - m * m) ** ((d - 3) / 2) / mp.beta(mp.mpf(1) / 2, (d - 1) / 2)
    f = lambda m: w(m) * m * m / (m * m + z * z)
    # tanh-sinh handles the endpoint singularities of the weight for d < 3
    return 2 * mp.quad(f, [0, min(z, mp.mpf(1) / 2), 1])
worst = mp.mpf(0)
for dd in [2, 3, mp.mpf('3.5'), 4, 5]:
    for z in [mp.mpf('0.05'), mp.mpf('0.3'), mp.mpf(1), mp.mpf(3), mp.mpf(10)]:
        worst = max(worst, abs(h_hyp(dd, z) - h_avg(dd, z)))
R.bound("SI Sec. S1", "hypergeometric h_d vs direct angular average, d = 2, 3, 3.5, 4, 5 (25 points)", worst, mp.mpf('1e-18'), "agree to 18 digits")
dev = max(abs(h_hyp(2, z) - (1 - z / mp.sqrt(1 + z * z))) + abs(h_hyp(3, z) - (1 - z * mp.atan(1 / z))) + abs(h_hyp(4, z) - (mp.sqrt(1 + z * z) - z) ** 2)
          for z in [mp.mpf('0.1'), mp.mpf('0.7'), mp.mpf(4)])
R.bound("main Eq. (10)", "h_2, h_3, h_4 closed forms equal the hypergeometric form", dev, mp.mpf('1e-35'), "identity")
for dd in [2, 3, 4, 5]:
    z = mp.mpf(10) ** 4
    R.close("main Sec. 2.2", f"h_d(z) ~ 1/(d z^2) at large z, d={dd} (z = 1e4, relative)", h_hyp(dd, z) * dd * z * z, 1, mp.mpf('1e-7'))

# ---------------------------------------------------------------- SI Eq. (S3) and the 3D integral of h^2
r = sp.symbols('r', positive=True)
R.exact("SI Eq. (S3)", "Int h_2^3 dz = 8 Int_0^1 r^4/(1+r^2)^2 dr = 10 - 3 pi", sp.simplify(sp.integrate(8 * r ** 4 / (1 + r ** 2) ** 2, (r, 0, 1)) - (10 - 3 * sp.pi)) == 0, True)
R.exact("SI Eq. (S3)", "Int h_4^(5/3) dz = 60/91", sp.integrate(r ** sp.Rational(4, 3) + r ** sp.Rational(10, 3), (r, 0, 1)), sp.Rational(60, 91))
R.close("SI Eq. (S3)", "quadrature of Int h_2^3 dz over the real line", 2 * mp.quad(lambda z: h_hyp(2, z) ** 3, [0, 1, 10, mp.inf]), 10 - 3 * mp.pi, mp.mpf('1e-30'))
R.close("SI Eq. (S3)", "quadrature of Int h_4^(5/3) dz", 2 * mp.quad(lambda z: ((mp.sqrt(1 + z * z) - z) ** 2) ** (mp.mpf(5) / 3), [0, 1, 10, mp.inf]), mp.mpf(60) / 91, mp.mpf('1e-30'))
R.close("SI Sec. S1", "Int h^2 dz over the real line = (2 pi/3)(1 - log 2)", 2 * mp.quad(lambda z: (1 - z * mp.atan(1 / z)) ** 2, [0, 1, 10, mp.inf]), 2 * mp.pi / 3 * (1 - mp.log(2)), mp.mpf('1e-30'))

# ---------------------------------------------------------------- g_2(z), SI Eq. (S4): 40-digit evaluation of the exact response
mp.mp.dps = 40
def ReS(q, z):
    w = (q / 2 + 1j * z) ** 2 - 1
    s = 1j * mp.sqrt(-w)         # branch analytic near q = 0 that tends to i sqrt(1+z^2)
    return s.real
for zs, printed in [('0.3', '-0.0302318708'), ('1', None), ('3', '-0.0011858541')]:
    z = mp.mpf(zs)
    c = mp.taylor(lambda q: ReS(q, z), 0, 5)        # L = 1 - 2 sum_k c_k q^(k-1)
    h2 = 1 - z / mp.sqrt(1 + z * z)
    g2 = -z / (8 * (1 + z * z) ** mp.mpf(2.5))
    R.close("SI Sec. S1", f"coefficient of q^0 of L(q,qz) is h_2(z), z={zs}", 1 - 2 * c[1], h2, mp.mpf('1e-30'))
    R.close("SI Sec. S1", f"coefficient of q of L(q,qz) vanishes, z={zs}", -2 * c[2], 0, mp.mpf('1e-30'))
    R.close("SI Eq. (S4)", f"coefficient of q^2 equals g_2(z) = -z/[8(1+z^2)^(5/2)], z={zs}", -2 * c[3], g2, mp.mpf('1e-30'))
    R.close("SI Sec. S1", f"coefficient of q^3 vanishes, z={zs}", -2 * c[4], 0, mp.mpf('1e-30'))
    if printed:
        R.digits("SI Sec. S1", f"g_2({zs})", g2, printed)
    else:
        R.close("SI Sec. S1", "g_2(1) = -1/(32 sqrt 2)", g2, -1 / (32 * mp.sqrt(2)), mp.mpf('1e-35'))

# ---------------------------------------------------------------- g_3(z), main Sec. 2.9, from SI Eq. (S79)
def L3(q, z):          # SI Eq. (S79), 0 < q < 1
    m1 = lambda a: a
    m2 = lambda a: (1 - (a - q / 2) ** 2) / (2 * q)
    f = lambda a, m: a * m(a) / (a * a + z * z)
    return mp.quad(lambda a: f(a, m1), [0, 1 - q / 2]) + mp.quad(lambda a: f(a, m2), [1 - q / 2, 1 + q / 2])
for zs in ['0.3', '1', '3']:
    z = mp.mpf(zs)
    h = 1 - z * mp.atan(1 / z)
    qq = [mp.mpf(10) ** -k for k in (3, 4)]
    est = [(L3(q, z) - h) / q ** 2 for q in qq]
    R.close("main Sec. 2.9", f"[L(q,qz) - h(z)]/q^2 -> g_3(z) = -1/[12(1+z^2)^2], z={zs} (q = 1e-4)", est[1], -1 / (12 * (1 + z * z) ** 2), mp.mpf('1e-6'), rel=True)

raise SystemExit(R.write())
