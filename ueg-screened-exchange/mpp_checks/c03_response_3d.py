"""Main-text Sec. 2.9 and SI Sec. S3: first-order exchange polarization in three
dimensions and the Carr--Maradudin coefficient.

Reproduces: R_3 = R^(1)/(4 pi), SI Eq. (S8), by quadrature of the reduced double
integral at z = 0.1, 0.3, 1, 3 (printed values and the 2.5e-17 agreement); the
azimuthal average 1/(2|mu - mu'|); the integrals Int h R_3 dz and Int h^2 g_3 dz;
the two parts of the Carr--Maradudin coefficient and their sum against the closed
form of C to 60 digits; and the quadratures of main-text Eqs. (30) and (37) at 50
digits against their closed forms (40 digits).
"""
import mpmath as mp
from common import Recorder

R = Recorder("c03_response_3d", "first-order exchange polarization in 3D and the Carr--Maradudin coefficient (main Sec. 2.9, SI Sec. S3)")
mp.mp.dps = 30
h = lambda z: 1 - z * mp.atan(1 / z)
R1 = lambda z: 4 / (1 + z * z) * ((1 + 3 * z * z) - (2 + 3 * z * z) * z * mp.atan(1 / z))
g3 = lambda z: -1 / (12 * (1 + z * z) ** 2)

def R3_quad(z):
    """SI Eq. (S8): -(1/8pi) Int Int mu mu' |mu-mu'| Re[(iz-mu)^-2 (iz-mu')^-2], on the triangle mu' < mu, doubled."""
    f = lambda m, n: m * n * (m - n) * mp.re(1 / ((1j * z - m) ** 2 * (1j * z - n) ** 2))
    br = sorted(set([-1, -z, 0, z, 1]) if z < 1 else [-1, 0, 1])
    inner = lambda m: mp.quad(lambda n: f(m, n), [x for x in br if x < m] + [m])
    return -2 * mp.quad(inner, br) / (8 * mp.pi)

printed = {'0.1': '0.23049451763139286', '0.3': '0.11645170899574007', '1': '0.011619772367581344', '3': '0.0002451548733833656'}
worst_pr = mp.mpf(0)
for zs, pv in printed.items():
    z = mp.mpf(zs); q = R3_quad(z); c = R1(z) / (4 * mp.pi)
    R.close("SI Sec. S3", f"quadrature of SI Eq. (S8) equals R^(1)/(4 pi), z={zs}", q, c, mp.mpf('1e-22'))
    worst_pr = max(worst_pr, abs(mp.mpf(pv) - c))
R.bound("main Sec. 2.9, SI Sec. S3", "printed quadrature values R_3(z) against R^(1)(z)/(4 pi)", worst_pr, mp.mpf('2.5e-17'), "within 2.5e-17")
mu, mu2 = mp.mpf('0.3'), mp.mpf('-0.5')
av = mp.quad(lambda ph: 1 / (2 * (1 - (mu * mu2 + mp.sqrt((1 - mu ** 2) * (1 - mu2 ** 2)) * mp.cos(ph)))), [0, mp.pi]) / mp.pi
R.close("SI Sec. S3", "azimuthal average of |p - p'|^-2 at fixed mu, mu' equals 1/(2|mu - mu'|)", av, 1 / (2 * abs(mu - mu2)), mp.mpf('1e-25'))

mp.mp.dps = 50
IhR3 = 2 * mp.quad(lambda z: h(z) * R1(z) / (4 * mp.pi), [0, 1, mp.inf])
Ih2g3 = 2 * mp.quad(lambda z: h(z) ** 2 * g3(z), [0, 1, mp.inf])
R.close("SI Sec. S3", "Int h R_3 dz = (pi^2 - 12 log 2)/12", IhR3, (mp.pi ** 2 - 12 * mp.log(2)) / 12, mp.mpf('1e-45'))
R.close("SI Sec. S3", "Int h^2 g_3 dz = -pi (pi^2 - 6)/288", Ih2g3, -mp.pi * (mp.pi ** 2 - 6) / 288, mp.mpf('1e-45'))
aL = (4 / (9 * mp.pi)) ** (mp.mpf(1) / 3)
# residue integral of main Eq. (37) and the ring part, both from quadrature
dub = 3 * aL / (2 * mp.pi ** 4) * 2 * mp.quad(lambda z: h(z) * R1(z), [0, 1, mp.inf])
ring = -12 / mp.pi ** 3 * Ih2g3 * 2 * aL / mp.pi
R.close("main Eq. (37), SI Sec. S19", "residue integral (50-digit quadrature) vs alpha_L (pi^2 - 12 log 2)/(2 pi^3), 40 digits", dub, aL * (mp.pi ** 2 - 12 * mp.log(2)) / (2 * mp.pi ** 3), mp.mpf('1e-40'))
R.digits("SI Sec. S19", "residue integral", dub, "0.0130393528")
R.digits("main Sec. 2.9", "exchange part of the Carr--Maradudin coefficient", dub, "0.01304")
R.digits("main Sec. 2.9", "ring part alpha_L(pi^2 - 6)/(12 pi^3)", ring, "0.00542")
Ctot = dub + ring
Cclosed = aL * (7 * mp.pi ** 2 / 6 - 12 * mp.log(2) - 1) / (2 * mp.pi ** 3)
with mp.workdps(80):
    dub80 = 3 * aL / (2 * mp.pi ** 4) * 2 * mp.quad(lambda z: h(z) * R1(z), [0, 1, mp.inf])
    ring80 = -12 / mp.pi ** 3 * (2 * mp.quad(lambda z: h(z) ** 2 * g3(z), [0, 1, mp.inf])) * 2 * aL / mp.pi
    C80 = aL * (7 * mp.pi ** 2 / 6 - 12 * mp.log(2) - 1) / (2 * mp.pi ** 3)
    R.close("SI Sec. S3", "sum of the two parts (quadrature) vs closed form of C, 60 digits", dub80 + ring80, C80, mp.mpf('1e-60'))
R.digits("SI Sec. S3", "sum of the two parts of C", Ctot, "0.0184584246606864")
R.digits("main Eq. (1)", "C", Cclosed, "0.0185")
GB = 3 / mp.pi ** 3 * 2 * mp.quad(lambda z: h(z) ** 2, [0, 1, mp.inf])
R.close("main Eq. (30), SI Sec. S19", "(3/pi^3) Int h^2 dz (50-digit quadrature) vs 2(1 - log 2)/pi^2, 40 digits", GB, 2 * (1 - mp.log(2)) / mp.pi ** 2, mp.mpf('1e-40'))
R.digits("SI Sec. S19", "(3/pi^3) Int h^2 dz", GB, "0.0621813817")
R.digits("SI Sec. S19", "Int_0^inf h^2 dz", GB * mp.pi ** 3 / 6, "0.32133552109539704")
raise SystemExit(R.write())
