"""SI Sec. S7 (screened internal exchange line) and main Sec. 3.1: replacing the
bare internal exchange line 1/|p-k| of SI Eq. (S6) by 1/(|p-k| + kappa_s) in both
the vertex and the self-energy parts changes R_2(1) linearly in kappa_s.

Reproduces [R_2(1; kappa_s) - R_2(1)]/kappa_s for kappa_s = 1e-1 ... 1e-5 and the
kappa_s -> 0 value -0.0321699 obtained from the integrand.  In the long-wavelength
limit at fixed z the line becomes 1/(|p^ - p^'| + kappa_s) with |p^ - p^'| = 2|sin delta|.
"""
import mpmath as mp
from common import Recorder

R = Recorder("c07_screened_line", "screened internal exchange line (SI Sec. S7)")
mp.mp.dps = 20
z = mp.mpf(1)
def kernel(sg, dl):
    """mu mu' Re[(b-b')^2] / D^2 times D^2, written as cos t cos t' 4 sin^2(sigma) sin^2(delta) Re[b^2 b'^2]."""
    t, tp = sg + dl, sg - dl
    b, bp = 1 / (1j * z - mp.cos(t)), 1 / (1j * z - mp.cos(tp))
    return mp.cos(t) * mp.cos(tp) * 4 * mp.sin(sg) ** 2 * mp.sin(dl) ** 2 * mp.re(b * b * bp * bp)
pref = mp.mpf(1) / 2 * 2 / (2 * mp.pi) ** 2      # -(1/2) << >> with the (sigma, delta) Jacobian, sign absorbed below
S = mp.linspace(0, 2 * mp.pi, 9)
def slope(kap):
    # [R(kappa) - R(0)]/kappa = +(1/2) << N / (D (D + kappa)) >>
    D = lambda dl: 2 * mp.sin(dl)
    br = sorted(set([0, kap / 10, kap, 10 * kap, mp.mpf(1) / 10, mp.pi / 2, mp.pi - mp.mpf(1) / 10, mp.pi - 10 * kap, mp.pi - kap, mp.pi - kap / 10, mp.pi]))
    br = [x for x in br if 0 <= x <= mp.pi]
    if kap == 0:
        f = lambda sg, dl: kernel(sg, dl) / D(dl) ** 2
    else:
        f = lambda sg, dl: kernel(sg, dl) / (D(dl) * (D(dl) + kap))
    return pref * mp.quad(f, S, br)
for ks, pr in [('1e-1', '-0.02951'), ('1e-2', '-0.03180'), ('1e-3', '-0.03212'), ('1e-4', '-0.03216'), ('1e-5', '-0.03217')]:
    R.digits("SI Sec. S7", f"[R_2(1; kappa_s) - R_2(1)]/kappa_s at kappa_s = {ks}", slope(mp.mpf(ks)), pr)
R.digits("SI Sec. S7", "kappa_s -> 0 value from the integrand", slope(mp.mpf(0)), "-0.0321699")
raise SystemExit(R.write())
