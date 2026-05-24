"""Panel-resolution check for scripts/sample_phi_grid.py.

1. Every u panel produced by _u_edges(y) is no wider than min(pi/y, 1) for
   y > pi, and no wider than min(pi/y, max(1, GROWTH*u)) for y <= pi; the panel
   count stays small at the smallest grid abscissa, y = 1e-7.
2. The contribution of u in [5, 10] to Phi (prefactor y included) at y = 250 and
   400, with the v rule of sample_phi_grid.py: the former single 48-node panel
   against the current split panels and a further 2x refinement.
3. Regression at two grid abscissae against scripts/phi_point_crosscheck.py.

Run from the archive root:  python scripts/check_phi_panels.py
"""
import os, sys
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sample_phi_grid as S
import phi_point_crosscheck as X


def piece(y, edges, nu):
    vv, wv = S._v_nodes(y)
    x, w = leggauss(nu)
    tot = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        u = 0.5 * (b - a) * x + 0.5 * (a + b); wu = 0.5 * (b - a) * w
        f = np.log((u + 1) / (u - 1))[:, None] * np.log(1 / (u[:, None] * vv[None, :])) \
            * S._j1(0.5 * y * (u[:, None] + vv[None, :])) * S._j1(0.5 * y * (u[:, None] - vv[None, :]))
        tot += float(wu @ f @ wv)
    return y * tot


ok = True
for y in (1e-7, 0.3, 3.0, 9.9985, 250.0, 400.0):
    e = S._u_edges(y); w = np.diff(e); a = e[:-1]
    if y > np.pi:
        lim = np.full_like(w, min(np.pi / y, 1.0))
    else:
        lim = np.where(a < 10.0 - 1e-12, min(np.pi / y, 1.0),
                       np.minimum(np.pi / y, np.maximum(1.0, S.GROWTH * a)))
    good = bool(np.all(w <= lim * (1 + 1e-12))) and (y > np.pi or len(w) < 1e4)
    ok &= good
    print(f"y={y:<8g} {len(w):7d} u panels, max width/limit {np.max(w / lim):.6f}  {'OK' if good else 'FAIL'}")

for y in (250.0, 400.0):
    e = S._u_edges(y); cur = e[(e >= 5.0 - 1e-12) & (e <= 10.0 + 1e-12)]
    fine = np.linspace(5.0, 10.0, 2 * (len(cur) - 1) + 1)
    old, new, ref = piece(y, np.array([5.0, 10.0]), 48), piece(y, cur, S.NU), piece(y, fine, S.NU)
    good = abs(new - ref) < 1e-14
    ok &= good
    print(f"y={y:g} [5,10]: former single panel {old:+.10e}  current {new:+.12e}  2x refined {ref:+.12e}  {'OK' if good else 'FAIL'}")

for y in (0.49759510651824584, 9.998538477652408):
    a = S.phi(y, add_tail=True); b = X.phi(y)
    good = abs(a - b) < 1e-9
    ok &= good
    print(f"y={y:.6f}: sample_phi_grid+tail {a:+.12e}  crosscheck {b:+.12e}  diff {abs(a-b):.1e}  {'OK' if good else 'FAIL'}")
print("ALL OK" if ok else "FAILED")
sys.exit(0 if ok else 1)
