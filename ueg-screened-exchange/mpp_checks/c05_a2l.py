"""Main-text Eqs. (39), (40), (59) and SI Sec. S5: the coefficient a_2l.

Reproduces: the r-substitution of SI Eq. (S2); the integrals A, B, C of SI Eq. (S11)
in closed form (symbolic for A, B and the antiderivative Phi; multiprecision
quadrature for C); the closed form of a_2l; the multiprecision quadratures of
(S11) in r at 40 and 60 digits and of main Eq. (39) in z at 100 digits; the
integer-relation search (PSLQ) for pi^2 a_2l over {1, pi, G, pi log 2}; the printed
digits in Ry and Hartree; and the size of the term at r_s = 0.1 and 1.
"""
import mpmath as mp, sympy as sp
from common import Recorder

R = Recorder("c05_a2l", "the coefficient a_2l (main Eqs. (39)-(40), SI Sec. S5)")
r, u = sp.symbols('r u', positive=True)
Z = (1 - r ** 2) / (2 * r)
R.exact("SI Eq. (S2)", "1+z^2 = [(1+r^2)/(2r)]^2, h_2 = 2r^2/(1+r^2), a = (1-r^2)/(1+r^2), |dz| = (1+r^2)/(2r^2) dr",
        [sp.simplify(1 + Z ** 2 - ((1 + r ** 2) / (2 * r)) ** 2) == 0, sp.simplify(1 - Z / ((1 + r ** 2) / (2 * r)) - 2 * r ** 2 / (1 + r ** 2)) == 0,
         sp.simplify(Z / ((1 + r ** 2) / (2 * r)) - (1 - r ** 2) / (1 + r ** 2)) == 0, sp.simplify(-sp.diff(Z, r) - (1 + r ** 2) / (2 * r ** 2)) == 0], [True] * 4)
A = sp.integrate(r ** 6 / (1 + r ** 2) ** 3, (r, 0, 1)); B = sp.integrate(r ** 4 * (1 - r ** 2) ** 2 / (1 + r ** 2) ** 4, (r, 0, 1))
R.exact("SI Sec. S5", "A = 3/2 - 15 pi/32", sp.simplify(A - (sp.Rational(3, 2) - 15 * sp.pi / 32)), 0)
R.exact("SI Sec. S5", "B = 13/6 - 11 pi/16", sp.simplify(B - (sp.Rational(13, 6) - 11 * sp.pi / 16)), 0)
Phi = sp.log(u) / 2 + sp.Rational(5, 2) / u - 2 / u ** 2 + sp.Rational(2, 3) / u ** 3
R.exact("SI Sec. S5", "Phi(1+r^2) is an antiderivative of r^3 (1-r^2)^2/(1+r^2)^4", sp.simplify(sp.diff(Phi.subs(u, 1 + r ** 2), r) - r ** 3 * (1 - r ** 2) ** 2 / (1 + r ** 2) ** 4), 0)
R.exact("SI Sec. S5", "Phi at r = 1 equals log(2)/2 + 5/6", sp.simplify(Phi.subs(u, 2) - (sp.log(2) / 2 + sp.Rational(5, 6))), 0)
G = mp.catalan
for dps, tol in [(40, mp.mpf('1e-41')), (60, mp.mpf('3e-61'))]:
    with mp.workdps(dps + 10):
        Cq = mp.quad(lambda x: x ** 3 * (1 - x * x) ** 2 * mp.atan(x) / (1 + x * x) ** 4, [0, mp.mpf(1) / 2, 1])
        Aq = mp.quad(lambda x: x ** 6 / (1 + x * x) ** 3, [0, 1]); Bq = mp.quad(lambda x: x ** 4 * (1 - x * x) ** 2 / (1 + x * x) ** 4, [0, 1])
        a2l_q = 2 / mp.pi ** 2 * (-32 * Aq + 16 * Bq - 16 * Cq)
        a2l_c = 7 / mp.pi + 4 * mp.log(2) / mp.pi - 160 / (9 * mp.pi ** 2) - 16 * G / mp.pi ** 2
        R.close("SI Sec. S5", f"quadrature of (S11) in r at {dps} digits vs closed form", a2l_q, a2l_c, tol)
        R.close("SI Sec. S5", f"C = -5/18 - pi log2/8 + pi/32 + G/2 ({dps} digits)", Cq, -mp.mpf(5) / 18 - mp.pi * mp.log(2) / 8 + mp.pi / 32 + G / 2, tol)
with mp.workdps(110):
    h2 = lambda z: 1 - z / mp.sqrt(1 + z * z)
    def R2(z):
        s = mp.sqrt(1 + z * z); rr = s - z; a = z / s
        return 2 / mp.pi * (h2(z) - z / s ** 3 + 2 * a * a * rr * mp.atan(rr) / (1 + rr * rr))
    I = 2 * mp.quad(lambda z: h2(z) ** 2 * R2(z), [0, mp.mpf(1) / 4, 1, 4, 16, mp.inf])
    a2l_z = -2 / mp.pi * I
    a2l_c = 7 / mp.pi + 4 * mp.log(2) / mp.pi - 160 / (9 * mp.pi ** 2) - 16 * G / mp.pi ** 2
    R.close("main Eq. (40), SI Sec. S5", "quadrature of main Eq. (39) in z vs closed form, 100 digits", a2l_z, a2l_c, mp.mpf('1e-100'))
    rel = mp.pslq([a2l_c * mp.pi ** 2, 1, mp.pi, G, mp.pi * mp.log(2)], maxcoeff=10 ** 6, maxsteps=10 ** 6)
R.exact("SI Sec. S5", "PSLQ of pi^2 a_2l over {1, pi, G, pi log 2}: relation 9 pi^2 a_2l = -160 + 63 pi - 144 G + 36 pi log 2", rel, [9, 160, -63, 144, -36])
mp.mp.dps = 40
a2l = 7 / mp.pi + 4 * mp.log(2) / mp.pi - 160 / (9 * mp.pi ** 2) - 16 * G / mp.pi ** 2
R.digits("main Eq. (40)", "a_2l (Ry)", a2l, "-0.1754613739208")
R.digits("main Eq. (40)", "a_2l (Hartree)", a2l / 2, "-0.0877306869604")
R.digits("abstract, main Eq. (59)", "a_2l (Ry)", a2l, "-0.17546")
R.digits("main Sec. 3.1", "a_2l r_s^2 log r_s at r_s = 0.1, in units of 1e-3 Ry", a2l * mp.mpf('0.01') * mp.log(mp.mpf('0.1')) * 1000, "4.0")
aL = -2 * mp.sqrt(2) * (10 / (3 * mp.pi) - 1)
R.digits("main Sec. 3.1", "ratio to a_l r_s log r_s at r_s = 0.1 ('about a tenth')", a2l * mp.mpf('0.1') / aL, "0.1")
R.exact("main Sec. 3.1", "the term vanishes at r_s = 1", mp.log(1) == 0, True)
raise SystemExit(R.write())
