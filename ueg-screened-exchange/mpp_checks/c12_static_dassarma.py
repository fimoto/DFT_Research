"""SI Sec. S6, static comparison with Table I of Das Sarma (arXiv:2610.08309):
Delta R/R(0) = [R(q,0) - 2/pi]/(2/pi) at x = q/k_F = 0.1, 0.2, 0.3, 0.4.

finiteq2d_static_tips.py evaluates the crescent-coordinate product quadrature with the phi panels
ending at the crescent tips cos(phi0) = -q/2, where E_p = p.q + q^2/2 vanishes and the static
integrand is singular, with nodes graded towards the tips (phi = phi0 -/+ L u^4) and E_p and 1 - rho*
evaluated without cancellation.  The limit is estimated from
  * nphi = 800, 1200, 1600 at nt = 8, extrapolated with e(N) = L + A N^-p (p from the three values), and
  * the radial correction v(800, 24) - v(800, 8) plus the geometric tail of the nt = 8, 16, 24 sequence.
The four-digit values are compared with SI Sec. S6 (manuscript v4.21), which prints them, and with
Table I of Das Sarma (agreement within one unit of the fourth digit).

finiteq2d_static.py is the earlier crescent code, with the phi panels split at +-pi/2.  It is not
converged at z = 0.  It produced the row printed in manuscript v4.20 (4.454e-4, 2.020e-3, 5.096e-3,
9.625e-3, the values at nphi = 400, nt = 10 truncated to four digits), which is reproduced here as a
record.  About 55 min on one core.
"""
import numpy as np, mpmath as mp
from common import Recorder
import finiteq2d_static as orig
import finiteq2d_static_tips as tips

R = Recorder("c12_static_dassarma", "static finite-q comparison with Das Sarma (SI Sec. S6)")
R0 = 2 / np.pi
printed = {0.1: ('4.523', 1e-4), 0.2: ('2.104', 1e-3), 0.3: ('5.141', 1e-3), 0.4: ('9.699', 1e-3)}   # SI v4.21
ds = {0.1: '4.524', 0.2: '2.104', 0.3: '5.141', 0.4: '9.699'}                                       # Das Sarma, Table I
old_row = {0.1: '4.454', 0.2: '2.020', 0.3: '5.096', 0.4: '9.625'}                                  # SI v4.20
out = {}
for x in (0.1, 0.2, 0.3, 0.4):
    pr, sc = printed[x]
    f = lambda n, t: (tips.R(x, n, t, 4.0) - R0) / R0
    e8 = [f(n, 8) for n in (800, 1200, 1600)]
    d1, d2 = e8[0] - e8[1], e8[1] - e8[2]
    # e(N) = L + A N^-p: solve d1/d2 = (800^-p - 1200^-p)/(1200^-p - 1600^-p) for p
    g = lambda p: (800.0 ** -p - 1200.0 ** -p) / (1200.0 ** -p - 1600.0 ** -p) - d1 / d2
    p = float(mp.findroot(g, 2.5))
    tail_phi = d2 * 1600.0 ** -p / (1200.0 ** -p - 1600.0 ** -p)
    L8 = e8[2] - tail_phi
    t16, t24 = f(800, 16), f(800, 24)
    a, b = t16 - f(800, 8), t24 - t16
    r = b / a if a != 0 else 0.0
    tail_t = b * r / (1 - r) if 0 < r < 1 else 0.0
    est = L8 + (t24 - f(800, 8)) + tail_t
    unc = abs(tail_phi) + abs(tail_t)
    out[x] = dict(nt8=dict(zip(('800', '1200', '1600'), e8)), p=p, phi_tail=tail_phi, v800_16=t16, v800_24=t24,
                  nt_tail=tail_t, estimate=est, uncertainty=unc)
    R.note("SI Sec. S6", f"x={x}: nt=8 sequence {', '.join(f'{v:.6e}' for v in e8)}, p = {p:.2f}; (800,16) {t16:.6e}, (800,24) {t24:.6e}",
           f"estimate {est:.6e} +- {unc:.1e}")
    R.digits("SI Sec. S6", f"Delta R/R(0) at x={x} (units of {sc:g}), tip-adapted quadrature, extrapolated", mp.mpf(est) / sc, pr)
    R.close("SI Sec. S6", f"Delta R/R(0) at x={x} vs Table I of Das Sarma, within one unit of the fourth digit", est / sc, float(ds[x]), 1.0e-3 + 0.5e-3, printed=ds[x])
    v_old = (orig.R(x, 0.0, nphi=400, nt=10) - R0) / R0
    R.digits("record: SI Sec. S6 of v4.20", f"earlier crescent code (panels split at +-pi/2), nphi = 400, nt = 10, x={x} (units of {sc:g})", mp.mpf(v_old) / sc, old_row[x])
raise SystemExit(R.write(dict(values=out)))
