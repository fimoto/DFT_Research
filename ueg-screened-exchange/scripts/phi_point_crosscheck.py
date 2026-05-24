"""Second independent evaluation of Phi(y) at one abscissa (default: the grid
point y = 0.49759510651824584 used in LIMITATIONS.md).

    Phi(y) = y int_1^inf du L(u) int_0^1 dv log(1/(uv)) j1((u+v)y/2) j1((u-v)y/2)

Independent of scripts/sample_phi_grid.py in its variables and panels:
  * v = exp(-t), t in [0, 45]: removes the log singularity at v = 0;
  * u = 1 + exp(s), s in [-40, log 9]: removes the log singularity of L at u = 1,
    with L = log(2 + e^s) - s evaluated without cancellation;
  * u in [10, U]: unit panels (period of cos(u y) is 2 pi / y);
  * u > U: leading non-oscillatory tail of Phi itself (prefactor y included),
        T(U) = [2 Si(y) - (2 log U + 1) sin y] / (y^2 U^2),
    from L(u) ~ 2/u, j1(x) ~ -cos x / x and the product-to-sum identity.

VALID RANGE: y <= 10. The u panels of width 1 carry nu = 24 nodes, and the s and
t panels are fixed, so the oscillation of period 2 pi / y is not resolved at large
y. At y = 10.0 this script agrees with sample_phi_grid.py (fully split panels,
tail added) to 2.4e-11; at y = 30.24 it is off by 4.0e-4 and at y = 57.3 by
2.9e-4, while doubling U changes it by only 1e-14 -- a stable cut is not
evidence of resolved panels. phi() refuses y > 10; use sample_phi_grid.py there.

Recorded result (U = 32000 and U + pi/y, which agree to 1e-14):
    Phi(0.49759510651824584) = -0.21328963265(1)
    shipped figs/phi_grid.npy : -0.21329183661294945  (difference 2.20e-6, 1.03e-5 relative)
This agrees with sample_phi_grid.py plus the same tail to about 1e-11.

    python scripts/phi_point_crosscheck.py [y] [U]      # about 3 minutes at U = 32000
"""
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sici


def _j1(x):
    return np.where(np.abs(x) < 1e-6, x / 3 - x ** 3 / 30, (np.sin(x) / x - np.cos(x)) / x)


def _composite(a, b, npan, n):
    x, w = leggauss(n)
    e = np.linspace(a, b, npan + 1)
    lo, hi = e[:-1, None], e[1:, None]
    return (0.5 * (hi - lo) * x + 0.5 * (hi + lo)).ravel(), (0.5 * (hi - lo) * w).ravel()


def _inner_v(u, y, t, wt):
    v = np.exp(-t)
    return ((t[None, :] - np.log(u)[:, None]) * _j1(0.5 * y * (u[:, None] + v))
            * _j1(0.5 * y * (u[:, None] - v)) * v) @ wt


Y_VALID = 10.0


def phi(y, U=32000.0, ns=32, nt=24, nu=24, block=500, allow_large_y=False):
    if y > Y_VALID and not allow_large_y:
        raise ValueError(f'phi_point_crosscheck is validated only for y <= {Y_VALID}; '
                         'the fixed panels do not resolve the oscillation at larger y '
                         '(use scripts/sample_phi_grid.py)')
    t, wt = _composite(0.0, 45.0, 45, nt)
    s, ws = _composite(-40.0, np.log(9.0), 60, ns)
    u1 = 1 + np.exp(s)
    total = np.sum(ws * np.exp(s) * (np.log(2 + np.exp(s)) - s) * _inner_v(u1, y, t, wt))
    edges = np.append(np.arange(10.0, U, 1.0), U)
    x, w = leggauss(nu)
    for k in range(0, len(edges) - 1, block):
        hi = edges[k + 1:k + block + 1]
        lo = edges[k:k + len(hi)]
        uu = (0.5 * (hi - lo)[:, None] * x + 0.5 * (hi + lo)[:, None]).ravel()
        wu = (0.5 * (hi - lo)[:, None] * w).ravel()
        total += np.sum(wu * np.log((uu + 1) / (uu - 1)) * _inner_v(uu, y, t, wt))
    tail = (2 * sici(y)[0] - (2 * np.log(U) + 1) * np.sin(y)) / (y ** 2 * U ** 2)
    return y * total + tail


if __name__ == "__main__":
    y = float(sys.argv[1]) if len(sys.argv) > 1 else 0.49759510651824584
    U = float(sys.argv[2]) if len(sys.argv) > 2 else 32000.0
    print(f"Phi({y!r}) = {phi(y, U):.13f}  (U = {U:g})")
