"""Main-text Sec. 3.1 (properties (i)-(iii)), main Eqs. (38) and (58), SI Sec. S4:
closed form of R_2(z), its checks and the Fermi-surface factor of the vertex part.

Reproduces: the two forms of main-text Eq. (38) and the identity of SI Eq. (S9); the
Fourier coefficients of the Coulomb factor and the generating function of O_m;
R_2(0) = 2/pi and the large-z series 5/(12 pi z^4) - 11/(20 pi z^6) + 1283/(2240 pi z^8);
the direct angular quadrature of main Eq. (34) in the variables (sigma, delta) at
z = 0.3, 1, 3 (printed R_2(1) and the 7e-18 / 1e-13 agreements); the slopes of the
cut-off vertex part per unit log(1/epsilon) between epsilon = 1e-5 and 1e-7 at
z = 1e-4, 0.3, 1; the integration-by-parts identity behind SI Eq. (S10) and the
coefficient 0.03885 of main Eq. (58); and the printed values R_2(0.3), R_2(1) of the
caption of SI Table S3.
"""
import mpmath as mp, sympy as sp
from common import Recorder

R = Recorder("c04_R2", "closed form of R_2 and the vertex Fermi-surface factor (main Sec. 3.1, SI Sec. S4)")

# ---------------------------------------------------------------- symbolic: two forms, limits, large-z series
z = sp.symbols('z', positive=True)
s = sp.sqrt(1 + z ** 2); r = s - z; a = z / s; h2 = 1 - z / s; h2p = sp.diff(h2, z)
R2a = 2 / sp.pi * (h2 ** 2 + 2 * a ** 2 * (r * sp.atan(r) - r ** 2) / (1 + r ** 2))
R2b = 2 / sp.pi * (h2 + z * h2p + 2 * a ** 2 * r * sp.atan(r) / (1 + r ** 2))
R.exact("SI Eq. (S9)", "h_2 + z h_2' = h_2^2 - 2 a^2 r^2/(1+r^2)", sp.simplify(sp.radsimp(h2 + z * h2p - (h2 ** 2 - 2 * a ** 2 * r ** 2 / (1 + r ** 2)))) == 0, True)
R.exact("main Eq. (38)", "the two forms of R_2 agree (difference at z = 0.3, 1, 3 below 1e-40)", all(abs(sp.N((R2a - R2b).subs(z, v), 50)) < 1e-40 for v in (sp.Rational(3, 10), 1, 3)), True)
R.exact("main Sec. 3.1 (i)", "R_2(0) = 2/pi", sp.limit(R2b, z, 0), 2 / sp.pi)
w = sp.symbols('w', positive=True)
ser = sp.series(R2b.subs(z, 1 / w), w, 0, 10).removeO()
R.exact("main Sec. 3.1 (ii), SI Sec. S4", "large-z series 5/(12 pi z^4) - 11/(20 pi z^6) + 1283/(2240 pi z^8), no z^-2 term",
        [sp.nsimplify(sp.simplify(ser.coeff(w, k) * sp.pi)) for k in (2, 4, 6, 8)], [0, sp.Rational(5, 12), sp.Rational(-11, 20), sp.Rational(1283, 2240)])

# ---------------------------------------------------------------- Fourier ingredients
mp.mp.dps = 30
for m in [1, 2, 3, 5]:
    v = mp.quad(lambda ph: (mp.cos(m * ph) - 1) / (2 * abs(mp.sin(ph / 2))), [0, mp.pi]) / mp.pi
    R.close("SI Sec. S4", f"Fourier coefficient difference of 1/|2 sin(phi/2)| for m={m}: -(2/pi) O_m", v, -2 / mp.pi * sum(mp.mpf(1) / (2 * j - 1) for j in range(1, m + 1)), mp.mpf('1e-25'))
x = sp.symbols('x')
O = lambda m: sum(sp.Rational(1, 2 * j - 1) for j in range(1, m + 1))
gen = sp.series(sp.sqrt(x) * sp.atanh(sp.sqrt(x)) / (1 - x), x, 0, 10).removeO()
R.exact("SI Sec. S4", "generating function sum O_m x^m = sqrt(x) artanh(sqrt x)/(1-x) (through x^9)", sp.expand(gen - sum(O(m) * x ** m for m in range(1, 10))), 0)
zz = mp.mpf('0.7'); sz = mp.sqrt(1 + zz ** 2); rv = sz - zz
for m in [1, 2, 3]:
    v = mp.quad(lambda th: mp.cos(th) / (1j * zz - mp.cos(th)) * mp.exp(-1j * m * th), [0, mp.pi, 2 * mp.pi]) / (2 * mp.pi)
    R.close("SI Sec. S4", f"Fourier coefficient F_(1,{m}) of cos(theta) b, z = 0.7 (|difference|)", abs(v - (-1j * rv) ** (m - 1) * (1 - rv ** 2) / (2j * sz)), 0, mp.mpf('1e-25'))

def R2(z):
    s = mp.sqrt(1 + z * z); r = s - z; a = z / s
    return 2 / mp.pi * ((1 - z / s) - z / s ** 3 + 2 * a * a * r * mp.atan(r) / (1 + r * r))

# ---------------------------------------------------------------- direct angular quadrature, main Eq. (34) for d = 2
def R2_direct(z):
    """R_2 = -(1/2) <<cos t cos t' Re[(b-b')^2] / |2 sin((t-t')/2)|>> in sigma = (t+t')/2, delta = (t-t')/2:
    (mu-mu')^2 = 4 sin^2 sigma sin^2 delta, so the integrand 2 cos t cos t' sin^2 sigma |sin delta| Re[b^2 b'^2] is smooth."""
    def f(sg, dl):
        t, tp = sg + dl, sg - dl
        b, bp = 1 / (1j * z - mp.cos(t)), 1 / (1j * z - mp.cos(tp))
        return 2 * mp.cos(t) * mp.cos(tp) * mp.sin(sg) ** 2 * mp.sin(dl) * mp.re(b * b * bp * bp)
    S = mp.linspace(0, 2 * mp.pi, 9); D = mp.linspace(0, mp.pi, 5)
    I = mp.quad(f, S, D, method='gauss-legendre')
    return -mp.mpf(1) / 2 * 2 * I / (2 * mp.pi) ** 2
mp.mp.dps = 25
for zs, tol, pr in [('1', mp.mpf('7e-18'), '0.04977088289807023'), ('0.3', mp.mpf('1e-13'), None), ('3', mp.mpf('1e-13'), None)]:
    zv = mp.mpf(zs); q = R2_direct(zv); c = R2(zv)
    R.close("main Sec. 3.1, SI Sec. S4", f"angular quadrature of main Eq. (34) at z={zs} vs closed form (38)", q, c, mp.mpf('1e-20'))
    if pr:
        R.bound("SI Sec. S4", f"printed quadrature value R_2(1) = {pr} vs closed form", abs(mp.mpf(pr) - c), tol, "within 7e-18")
R.digits("SI Table S3 caption", "R_2(0.3)", R2(mp.mpf('0.3')), "0.3180676")
R.digits("SI Table S3 caption", "R_2(1)", R2(mp.mpf(1)), "0.0497709")

# ---------------------------------------------------------------- Ward slopes of the cut-off vertex part
mp.mp.dps = 25
def vertex_shell(z, e1, e2):
    """vertex part -(1/2)<<mu mu' Re(-2 b b')/|2 sin(phi/2)|>> restricted to e2 < |phi| < e1, phi = t' - t; phi = exp(u)."""
    def inner(u):
        ph = mp.exp(u)
        g = lambda t: sum(mp.cos(t) * mp.cos(t + sgn * ph) * mp.re(-2 / ((1j * z - mp.cos(t)) * (1j * z - mp.cos(t + sgn * ph)))) for sgn in (1, -1))
        th = mp.quad(g, mp.linspace(0, 2 * mp.pi, 9))
        return th * ph / (2 * abs(mp.sin(ph / 2)))
    I = mp.quad(inner, [mp.log(e2), mp.log(e1)])
    return -mp.mpf(1) / 2 * I / (2 * mp.pi) ** 2
for zs, pr in [('1e-4', '0.318246'), ('0.3', '0.1429308'), ('1', '-0.01930873')]:
    zv = mp.mpf(zs)
    slope = vertex_shell(zv, mp.mpf('1e-5'), mp.mpf('1e-7')) / mp.log(100)
    sz = mp.sqrt(1 + zv * zv)
    pred = ((1 - zv / sz) - zv / sz ** 3) / mp.pi
    R.digits("SI Sec. S4", f"slope per unit log(1/eps) of the vertex part, 1e-7 < eps < 1e-5, z={zs}", slope, pr)
    R.digits("SI Sec. S4", f"pi^-1 (h_2 + z h_2') at z={zs}", pred, pr)

# ---------------------------------------------------------------- isolated vertex in the ring chain, SI Eq. (S10), main Eq. (58)
mp.mp.dps = 30
hh = {2: (lambda z: 1 - z / mp.sqrt(1 + z * z), lambda z: -(1 + z * z) ** mp.mpf(-1.5)),
      3: (lambda z: 1 - z * mp.atan(1 / z), lambda z: -mp.atan(1 / z) + z / (1 + z * z))}
for dd in (2, 3):
    f, fp = hh[dd]
    for N in (1, 2, 3):
        L = mp.quad(lambda z: f(z) ** N * (f(z) + z * fp(z)), [0, 1, mp.inf]); Rr = mp.quad(lambda z: f(z) ** (N + 1), [0, 1, mp.inf]) * N / (N + 1)
        R.close("main Sec. 3.4, SI Sec. S4", f"Int h^N (h + z h') dz = [N/(N+1)] Int h^(N+1) dz, d={dd}, N={N}", L, Rr, mp.mpf('1e-25'))
vtx2 = 2 * (10 - 3 * mp.pi) / (3 * mp.pi ** 2)
R.digits("main Eq. (58)", "(2/(3 pi^2))(10 - 3 pi)", vtx2, "0.03885")
R.digits("main Table 4", "Pi^(1)_vtx in rings, d = 2", vtx2, "0.0389")
aL = (4 / (9 * mp.pi)) ** (mp.mpf(1) / 3)
chd = -aL * (1 - mp.log(2)) / (6 * mp.pi ** 3)
R.close("SI Sec. S4", "d = 3: (12 alpha_L/pi^3)(-1/(8 pi))(pi/3)(1 - log 2) = 3 c_hd", 12 * aL / mp.pi ** 3 * (-1 / (8 * mp.pi)) * (mp.pi / 3) * (1 - mp.log(2)), 3 * chd, mp.mpf('1e-28'))
raise SystemExit(R.write())
