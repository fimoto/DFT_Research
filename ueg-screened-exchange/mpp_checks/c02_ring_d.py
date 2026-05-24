"""Main-text Sec. 2.8, Table 2 and SI Sec. S2: ring series in d dimensions.

Reproduces: the scale-free ring integral, SI Eq. (S5); the Gell-Mann--Brueckner and
Rajagopal--Kimball coefficients, main-text Eqs. (30) and (31); A_d of main-text
Eq. (32) and Table 2 (d = 4-9, with the comparison column from Table I of
Schlesier et al. and the ratios); A_4 in closed form; (d-3)A_d/2 at d = 3.1, 3.01,
3.001; A_d/d at d = 100, 1000, 3000 and the limit -1/(2 e pi); the comparison of
main-text Eq. (32) with Eq. (21) of Schlesier et al. for 3.001 <= d <= 3000; and the
closed form of their series Sigma_D to 40 digits for d = 4, 5, 6, 9.
"""
import mpmath as mp
from common import Recorder

R = Recorder("c02_ring_d", "ring series in d dimensions (main Sec. 2.8, Table 2, SI Sec. S2)")
mp.mp.dps = 40

def alpha(d): return (2 ** (d - 1) * mp.gamma(d / 2 + 1) ** 2) ** (-1 / d)
def calA(d): return 2 * mp.pi ** (d / 2) / mp.gamma(d / 2)
def calC(d): return (4 * mp.pi) ** ((d - 1) / 2) * mp.gamma((d - 1) / 2)
def calB(d): return 2 * calA(d) * calC(d) / (2 * mp.pi) ** d
def hd(d, z):
    # main Eq. (9), h_d = 1 - 2F1(1, 1/2; d/2; -1/z^2), evaluated in the Pfaff-transformed form
    # 1 - z (1+z^2)^(-1/2) 2F1(d/2 - 1, 1/2; d/2; 1/(1+z^2)), which is stable for small z
    if z == 0: return mp.mpf(1)
    return 1 - z / mp.sqrt(1 + z * z) * mp.hyp2f1(d / 2 - 1, mp.mpf(1) / 2, d / 2, 1 / (1 + z * z))
def int_hN(d, N):
    pts = [0, 1 / mp.sqrt(d) / 4, 1 / mp.sqrt(d), 4 / mp.sqrt(d), 1, 10, mp.inf]
    pts = sorted(set(pts))
    return 2 * mp.quad(lambda z: hd(d, z) ** N, pts)
def Ad(d):
    d = mp.mpf(d); N = (d + 1) / (d - 1)
    return d * alpha(d) ** (N - 2) * calB(d) ** N / (4 * (d + 1) * mp.sin(mp.pi * N)) * int_hN(d, N)

# SI Eq. (S5)
for d in [mp.mpf(4), mp.mpf('5.5'), mp.mpf(9)]:
    N = (d + 1) / (d - 1); hh = mp.mpf('0.37')
    with mp.workdps(60):
        lhs = mp.quad(lambda t: t ** d * (mp.log1p(hh * t ** (1 - d)) - hh * t ** (1 - d)), [0, mp.mpf(1) / 4, 1, 4, 16, 100, mp.inf])
    R.close("SI Eq. (S5)", f"scale-free ring integral at d={mp.nstr(d,3)}, h=0.37", lhs, hh ** N / (d - 1) * mp.pi / (N * mp.sin(mp.pi * N)), mp.mpf('1e-20'))
# main Eqs. (30), (31)
h = lambda z: 1 - z * mp.atan(1 / z)
GB = 3 / mp.pi ** 3 * 2 * mp.quad(lambda z: h(z) ** 2, [0, 1, mp.inf])
R.close("main Eq. (30)", "(3/pi^3) Int h^2 dz = 2(1 - log 2)/pi^2", GB, 2 * (1 - mp.log(2)) / mp.pi ** 2, mp.mpf('1e-35'))
R.digits("main Eq. (30)", "Gell-Mann--Brueckner coefficient A", GB, "0.06218")
RK = -2 * mp.sqrt(2) / (3 * mp.pi) * 2 * mp.quad(lambda z: hd(mp.mpf(2), z) ** 3, [0, 1, mp.inf])
R.close("main Eq. (31)", "Rajagopal--Kimball coefficient -2 sqrt2 (10/(3 pi) - 1)", RK, -2 * mp.sqrt(2) * (10 / (3 * mp.pi) - 1), mp.mpf('1e-35'))
R.digits("main Eq. (2)", "a_l", RK, "-0.17263")

# main Eq. (32), Table 2
A4c = -24 * 2 ** (mp.mpf(5) / 12) / (91 * mp.sqrt(3))
R.close("main Sec. 2.8 (i)", "A_4 from Eq. (32) equals -24 2^(5/12)/(91 sqrt3)", Ad(4), A4c, mp.mpf('1e-30'))
schlesier_Ha = {4: '-0.0196', 5: '-0.0285', 6: '-0.0391', 7: '-0.0509', 8: '-0.0638', 9: '-0.0773'}   # Table I, xi = 0
tab2 = {4: ('-0.20325', '-0.0392', '5.19'), 5: ('-0.17975', '-0.0570', '3.15'), 6: ('-0.19753', '-0.0782', '2.53'),
        7: ('-0.22709', '-0.1018', '2.23'), 8: ('-0.26238', '-0.1276', '2.06'), 9: ('-0.30123', '-0.1546', '1.95')}
table = {}
for d, (pA, pS, pr) in tab2.items():
    a = Ad(d); s = 2 * mp.mpf(schlesier_Ha[d])
    table[d] = dict(A_d=mp.nstr(a, 12), schlesier_Ry=mp.nstr(s, 4), ratio=mp.nstr(a / s, 6))
    R.digits("main Table 2", f"A_{d} (Ry)", a, pA)
    R.digits("main Table 2", f"Table I of Schlesier et al. for d={d}, converted to Ry", s, pS)
    R.digits("main Table 2", f"ratio for d={d}", a / s, pr)

# (ii) d -> 3+
for ds, pr in [('3.1', '-0.0639390'), ('3.01', '-0.0623257'), ('3.001', '-0.0621955')]:
    d = mp.mpf(ds)
    R.digits("SI Sec. S2", f"(d-3) A_d/2 at d={ds}", (d - 3) * Ad(d) / 2, pr)
R.digits("SI Sec. S2", "-A = -2(1 - log 2)/pi^2", -2 * (1 - mp.log(2)) / mp.pi ** 2, "-0.0621814")
# (iii) large d
for d, pr in [(100, '-0.052756'), (1000, '-0.057680'), (3000, '-0.058216')]:
    R.digits("SI Sec. S2", f"A_d/d at d={d}", Ad(d) / d, pr)
R.digits("main Sec. 2.8 (iii)", "-1/(2 e pi)", -1 / (2 * mp.e * mp.pi), "-0.05855")
R.digits("SI Sec. S2", "-1/(2 e pi)", -1 / (2 * mp.e * mp.pi), "-0.058550")
# Gaussian limit of the z integral
for d in [mp.mpf(1000), mp.mpf(10000)]:
    I = 2 * mp.quad(lambda y: hd(d, y / mp.sqrt(d)), [0, 1, 10, mp.inf])
    R.close("SI Sec. S2", f"Int h_d(y/sqrt d) dy -> sqrt(2 pi), d={int(d)} (relative)", I, mp.sqrt(2 * mp.pi), mp.mpf(5) / mp.sqrt(d), rel=True)

# Schlesier et al., Eq. (21) and Sigma_D (their Appendix B), xi = 0, Hartree
def Sigma_series(D):
    a = 2 / (D - 1)
    return mp.nsum(lambda n: (-1) ** n / (n * (n + a + 1)) - (-1) ** n / (n * (a + 1 - n)), [1, mp.inf]) - 1 / (a + 1) ** 2
def R_sch(D, u):
    return mp.gamma(mp.mpf(3) / 2) * mp.gamma((D - 1) / 2) / (2 * u * u * mp.gamma(D / 2 + 1)) * mp.hyp2f1(1, mp.mpf(3) / 2, D / 2 + 1, -1 / u ** 2)
def c_sch(D, Sigma):
    N = (D + 1) / (D - 1); aD = 1 / alpha(D)   # their alpha_D = 1/alpha_d
    pts = sorted(set([0, 1 / mp.sqrt(D) / 4, 1 / mp.sqrt(D), 4 / mp.sqrt(D), 1, 10, mp.inf]))
    I = 2 * mp.quad(lambda u: R_sch(D, u) ** N, pts)
    return -2 * D * Sigma / (mp.pi ** 3 * (D - 1)) * (aD * mp.pi / 4) ** ((D - 3) / (D - 1)) * I
for d in [4, 5, 6, 9]:
    D = mp.mpf(d); N = (D + 1) / (D - 1)
    R.close("SI Sec. S2", f"Sigma_D series equals -pi/[N_r sin(pi N_r)], d={d} (40 digits)", Sigma_series(D), -mp.pi / (N * mp.sin(mp.pi * N)), mp.mpf('1e-40'))
mp.mp.dps = 30
worst = mp.mpf(0); grid = ['3.001', '3.01', '3.1', '3.5', '4', '5', '6', '7', '8', '9', '12', '20', '50', '100', '300', '1000', '3000']
for ds in grid:
    D = mp.mpf(ds); N = (D + 1) / (D - 1)
    ours = Ad(D); theirs = 2 * c_sch(D, Sigma_series(D))      # Hartree -> Ry
    worst = max(worst, abs(ours / theirs - 1))
R.bound("main Sec. 2.8, SI Sec. S2", f"main Eq. (32) vs Eq. (21) of Schlesier et al. (Hartree -> Ry), {len(grid)} values of d in [3.001, 3000], relative", worst, mp.mpf('1e-15'), "machine precision")
mp.mp.dps = 40
D = mp.mpf(3000)
R.close("SI Sec. S2", "their large-d limit -D/(2 e pi^2) Ha is smaller in magnitude than -d/(4 e pi) Ha by pi/2",
        (D / (4 * mp.e * mp.pi)) / (D / (2 * mp.e * mp.pi ** 2)), mp.pi / 2, mp.mpf('1e-35'))
raise SystemExit(R.write(dict(table2=table)))
