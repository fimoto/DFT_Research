"""Compare the shipped Phi grid with independent quadrature at selected abscissae.

Writes figs/phi_grid_point_checks.csv. The five abscissae are grid points; the
reference is scripts/phi_point_crosscheck.py at U = 8000 (tail included), and the
change from U = 4000 is recorded as a convergence indicator. These are selected
points, not a uniform accuracy bound on the 900-point grid.

Run from the archive root:  python scripts/phi_grid_point_checks.py
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import phi_point_crosscheck as X

TARGETS = [0.49759510651824584, 0.990693122989711, 2.021543164741066,
           5.021962607149256, 9.998538477652408]
ys, ph = np.load(os.path.join(HERE, '..', 'figs', 'phi_grid.npy'))
rows = []
for t in TARGETS:
    i = int(np.argmin(np.abs(ys - t)))
    v8, v4 = X.phi(ys[i], U=8000.0), X.phi(ys[i], U=4000.0)
    rows.append((ys[i], ph[i], v8, abs(ph[i] - v8), abs(ph[i] - v8) / abs(v8), abs(v8 - v4)))
    print('y=%.17g shipped %+.14e recomputed %+.14e abs %.3e rel %.3e U-change %.2e' % rows[-1])
out = os.path.join(HERE, '..', 'figs', 'phi_grid_point_checks.csv')
with open(out, 'w') as f:
    f.write('y,shipped,recomputed_U8000,abs_diff,rel_diff,abs_change_U4000_to_U8000\n')
    for r in rows:
        f.write('%.17g,%.17g,%.17g,%.4e,%.4e,%.3e\n' % r)
print('wrote figs/phi_grid_point_checks.csv; max relative difference %.2e' % max(r[4] for r in rows))
