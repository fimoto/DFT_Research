"""SI Sec. S8 and SI Table S4: spin dependence of a_l and a_2l in two dimensions.

Reproduces SI Eq. (S12) by quadrature: the ratios of SI Table S4, the fully polarized
values 1/(4 sqrt 2) and 1/8, the small-zeta coefficients -0.289 and -0.0355, and
the limit zeta -> 1 of the a_2l ratio, 0.81649, with the printed approach values.
"""
import mpmath as mp
from common import Recorder

R = Recorder("c08_spin_2d", "spin dependence in two dimensions (SI Sec. S8, Table S4)")
mp.mp.dps = 30
h2 = lambda z: 1 - z / mp.sqrt(1 + z * z)
def R2(z):
    s = mp.sqrt(1 + z * z); rr = s - z; a = z / s
    return 2 / mp.pi * (h2(z) - z / s ** 3 + 2 * a * a * rr * mp.atan(rr) / (1 + rr * rr))
DEN1 = 8 * mp.quad(lambda z: h2(z) ** 3, [0, 1, 4, mp.inf])
DEN2 = 8 * mp.quad(lambda z: R2(z) * h2(z) ** 2, [0, 1, 4, mp.inf])
def ratios(zeta):
    xp = mp.sqrt(1 + zeta); xm = mp.sqrt(1 - zeta)
    pts = sorted(set([0, xm / 4, xm, 4 * xm, 1, 4, mp.inf] if xm > 0 else [0, 1, 4, mp.inf]))
    S = lambda z: h2(z / xp) + (h2(z / xm) if xm > 0 else 0)
    SR = lambda z: R2(z / xp) / xp + (R2(z / xm) / xm if xm > 0 else 0)
    return mp.quad(lambda z: S(z) ** 3, pts) / DEN1, mp.quad(lambda z: SR(z) * S(z) ** 2, pts) / DEN2
tab = {'0.2': ('0.98826', '0.99855'), '0.4': ('0.95076', '0.99377'), '0.6': ('0.87865', '0.98391'),
       '0.8': ('0.74375', '0.96281'), '0.9': ('0.61946', '0.93979')}
a, b = ratios(mp.mpf(0))
R.close("SI Table S4", "zeta = 0: both ratios equal 1", a + b, 2, mp.mpf('1e-25'))
for zs, (pa, pb) in tab.items():
    a, b = ratios(mp.mpf(zs))
    R.digits("SI Table S4", f"a_l(zeta)/a_l(0) at zeta={zs}", a, pa)
    R.digits("SI Table S4", f"a_2l(zeta)/a_2l(0) at zeta={zs}", b, pb)
a, b = ratios(mp.mpf(1))
R.close("SI Sec. S8", "zeta = 1: a_l ratio = 1/(4 sqrt 2)", a, 1 / (4 * mp.sqrt(2)), mp.mpf('1e-25'))
R.close("SI Sec. S8", "zeta = 1: a_2l ratio = 1/8", b, mp.mpf(1) / 8, mp.mpf('1e-25'))
# small zeta: (ratio - 1)/zeta^2 by Richardson extrapolation from zeta = 0.01, 0.005
ests = []
for zs in ['0.01', '0.005']:
    a, b = ratios(mp.mpf(zs)); z2 = mp.mpf(zs) ** 2
    ests.append(((a - 1) / z2, (b - 1) / z2))
ca = (4 * ests[1][0] - ests[0][0]) / 3; cb = (4 * ests[1][1] - ests[0][1]) / 3
R.digits("SI Sec. S8", "small-zeta coefficient of a_l(zeta)/a_l(0)", ca, "-0.289")
R.digits("SI Sec. S8", "small-zeta coefficient of a_2l(zeta)/a_2l(0)", cb, "-0.0355")
lim = mp.mpf(1) / 8 + mp.quad(lambda u: R2(u) * (1 + h2(u)) ** 2, [0, 1, 4, mp.inf]) / DEN2
R.digits("SI Sec. S8", "limit zeta -> 1 of a_2l(zeta)/a_2l(0)", lim, "0.81649")
for e, pr in [(2, '0.87373'), (4, '0.82392'), (6, '0.81727'), (8, '0.81657')]:
    a, b = ratios(1 - mp.mpf(10) ** -e)
    R.digits("SI Sec. S8", f"a_2l ratio at 1 - zeta = 1e-{e}", b, pr)
raise SystemExit(R.write())
