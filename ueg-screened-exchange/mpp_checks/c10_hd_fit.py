"""Main-text Sec. 3.3 and SI Sec. S19, SI Table S5: weighted least-squares fits of
the direct high-density values of e_kite.

Reads the stored direct values (validation/high_density/results/direct_high_density.json,
64 replicates per density, disjoint seeds) and the analytic coefficients
(validation/high_density/results/hd_coefficients.json).  Reproduces the entries of
SI Table S5 (rounded from the unrounded means), e_2x = 48.35832 mRy, the free
three-parameter fits c_fit r_s log^2 r_s + c_1,fit r_s log r_s + c_2 r_s for r_s <= 0.1
and r_s <= 0.03, the weighted sums of squared residuals with and without the
r_s log^2 r_s term, the fits with c_fit = c_hd fixed for r_s <= 0.1, 0.03, 0.01, and
the distance of the r_s <= 0.03 value from c_1 in standard errors.
"""
import json, os
import numpy as np, mpmath as mp
from common import Recorder, HERE

R = Recorder("c10_hd_fit", "fits of the direct high-density kite energies (main Sec. 3.3, SI Sec. S19, Table S5)")
root = os.path.dirname(HERE)
J = json.load(open(os.path.join(root, 'validation/high_density/results/direct_high_density.json')))
C = json.load(open(os.path.join(root, 'validation/high_density/results/hd_coefficients.json')))['results']['0.0']
rows = {float(k): v for k, v in J['rows'].items()}
rs = np.array([0.001, 0.003, 0.01, 0.03, 0.1])
e = np.array([rows[r]['mean_mry'] for r in rs]); se = np.array([rows[r]['se_mry'] for r in rs])
reps = [rows[r]['replicates'] for r in rs]
e2x = float(mp.log(2) / 3 - 3 * mp.zeta(3) / (2 * mp.pi ** 2)) * 1000
R.digits("SI Sec. S19, Table S5", "e_2x in mRy", e2x, "48.35832")
R.exact("SI Table S5", "64 replicates per density", reps, [64] * 5)
printed = {0.001: ('48.22031', '-0.13801', '0.00073'), 0.003: ('48.01855', '-0.33977', '0.00100'), 0.01: ('47.47046', '-0.88786', '0.00196'),
           0.03: ('46.29632', '-2.06200', '0.00326'), 0.1: ('43.43311', '-4.92520', '0.00614')}
for r, ei, si in zip(rs, e, se):
    pe, pd, ps = printed[r]
    R.digits("SI Table S5", f"e_kite at r_s={r} (mRy)", ei, pe)
    R.digits("SI Table S5", f"e_kite - e_2x at r_s={r} (mRy), from the unrounded values", ei - e2x, pd)
    R.digits("SI Table S5", f"standard error at r_s={r} (mRy)", si, ps)
R.digits("main Sec. 3.3", "direct value at r_s = 0.1 (mRy, 43.433(6))", e[-1], "43.433")
y = (e - e2x) / 1000; s = se / 1000; L = np.log(rs)
def fit(cols, mask, yy):
    A = np.vstack([c[mask] for c in cols]).T; W = 1 / s[mask] ** 2
    M = A.T @ (A * W[:, None]); c = np.linalg.solve(M, A.T @ (W * yy[mask])); cov = np.linalg.inv(M)
    chi = float(np.sum(W * (yy[mask] - A @ c) ** 2)); return c, np.sqrt(np.diag(cov)), chi
allm = np.ones(5, bool)
c, ce, chi = fit([rs * L ** 2, rs * L, rs], allm, y)
R.digits("main Sec. 3.3, SI Sec. S19", "c_fit, all five values (Ry)", c[0], "-0.00091")
R.digits("main Sec. 3.3, SI Sec. S19", "its standard error (Ry)", ce[0], "0.00005")
R.digits("main Sec. 3.3, SI Sec. S19", "weighted sum of squared residuals, two degrees of freedom", chi, "0.02")
c4, ce4, _ = fit([rs * L ** 2, rs * L, rs], rs <= 0.03 + 1e-12, y)
R.digits("main Sec. 3.3, SI Sec. S19", "c_fit, four values with r_s <= 0.03 (Ry)", c4[0], "-0.00092")
R.digits("main Sec. 3.3, SI Sec. S19", "its standard error (Ry)", ce4[0], "0.00012")
_, _, chi0 = fit([rs * L, rs], allm, y)
R.digits("main Sec. 3.3, SI Sec. S19", "weighted sum of squared residuals without the r_s log^2 r_s term, three degrees of freedom", chi0, "380")
chd = C['c_hd']; c1 = C['c_1']
R.digits("main Eq. (52)", "c_hd (Ry)", chd, "-0.0008594457")
R.digits("main Sec. 3.3", "c_1 (Ry), companion paper", c1, "0.01150")
for cut, pv, pe in [(0.1, '0.01125', '0.00006'), (0.03, '0.01134', '0.00012'), (0.01, '0.01143', '0.00025')]:
    cc, cce, _ = fit([rs * L, rs], rs <= cut + 1e-12, y - chd * rs * L ** 2)
    R.digits("main Sec. 3.3, SI Sec. S19", f"c_1,fit with c_fit = c_hd, r_s <= {cut} (Ry)", cc[0], pv)
    R.digits("main Sec. 3.3, SI Sec. S19", f"its standard error, r_s <= {cut} (Ry)", cce[0], pe)
    if cut == 0.03:
        R.digits("main Sec. 3.3, SI Sec. S19", "distance of the r_s <= 0.03 value from c_1, in standard errors ('within 1.4 sigma')", abs(cc[0] - c1) / cce[0], "1.4")
raise SystemExit(R.write())
