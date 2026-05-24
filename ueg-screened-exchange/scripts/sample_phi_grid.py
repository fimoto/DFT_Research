"""Independent recomputation of the geometric profile Phi(y).

    Phi(y) = y int_1^inf du L(u) int_0^1 dv l(u,v)
                 j1((u+v)y/2) j1((u-v)y/2),
    L(u) = log((u+1)/(u-1)),   l(u,v) = log(1/(uv)),

which is Eq. (Phi) of the JPA paper.

WHAT THIS SCRIPT IS
    An independent evaluation of the same integral, so that the shipped grid
    figs/phi_grid.npy can be checked rather than taken on trust.  It is NOT a
    bit-for-bit reproduction: the panel layout used to produce that grid was not
    recorded, and the two evaluations differ at the level quoted below.  Output
    goes to figs/phi_grid_recomputed.npy; the shipped grid is not overwritten,
    and the figures continue to use the shipped grid.

QUADRATURE
    Both integrands oscillate on scale 1/y, so the panels are chosen per y:
      * u: graded edges 1, 1+1e-8, ..., 5, 10 to resolve the logarithmic
        singularity of L(u) at u = 1, out to u_max = UMAX_SCALE * max(200, 400/y);
        below u = 10 every interval is split into panels no wider than
        min(pi/y, 1). Above u = 10: for y > pi, panels of equal width at most
        pi/y (half a period of cos(u y)); for y <= pi, where u_max reaches 1e10,
        widths grow geometrically, min(pi/y, max(1, GROWTH*u)), which follows
        the algebraic decay while u*y is small and is capped at half a period
        once the oscillation sets in (505 panels at y = 1e-7; the equal-width
        rule would need 1e10). NU Gauss-Legendre nodes per panel; layouts with
        more than MAX_U_PANELS panels are refused;
      * v: graded panels 0, 1e-10, ..., 1e-1 for the log singularity of
        log(1/(uv)) at v = 0, then panels of width (2 pi/y)/8, NV nodes each.
    The u integral is accumulated in blocks of panels whose size is set from a
    memory budget (MEM_BUDGET_MIB, default 256 MiB), so memory does not grow
    with the number of u nodes.  phi() returns the integral up to u_max; the
    tail T(U) below is added only with add_tail=True (or --tail).

MEASURED AGREEMENT WITH THE SHIPPED GRID
    At the grid point y = 0.49759510651824584 the shipped value is
    -0.21329183661294945.  This script gives -0.2132527211, -0.2132891275 and
    -0.2132895959 at UMAX_SCALE = 1, 10 and 40.  Adding the leading
    non-oscillatory tail beyond u_max,
        T(U) = [2 Si(y) - (2 log U + 1) sin y] / (y^2 U^2),
    gives -0.2132896325 at UMAX_SCALE = 10 and -0.21328963265 at 40 and 80,
    where the two agree to about 1e-12.
    scripts/phi_point_crosscheck.py, which uses different variables, panels and
    cut, gives -0.21328963265 as well.  The two independent evaluations agree to
    about 1e-11 and both differ from the shipped value by 2.20e-6
    (1.03e-5 relative).  At this point the residual is therefore an error of the
    shipped grid, not of the recomputation.  Only this abscissa has been resolved
    to that level; the grid has not been recomputed at full precision everywhere.
    Differences of this size are invisible in the Phi and S(Lambda) figures (JPA Figs. 5 and 6).

    Large y: earlier layouts did not resolve cos(u y) at large y. Width-1 panels
    above u = 10 were replaced by half-period panels in the 2026-09-16 revision,
    but the graded panels below u = 10 were left unsplit, so at y = 250 the
    single panel [5, 10] gave -3.057e-5 against the converged 3.2707e-6. All
    intervals are now split; see scripts/check_phi_panels.py.

    On the integral over the shipped grid: the grid itself was not changed, and
    the value depends on the rule used on the 900 logarithmically spaced points:
        sum Phi * gradient(log y)       0.249891331
        trapezoid of Phi/y in y         0.249914892
        trapezoid of Phi in log y       0.249889720
        exact D_0                       0.249901897
    The spread (about 1e-4 relative) comes from the integration rule and the
    spacing in the oscillatory tail, where the grid has about 6 points per period
    of sin y at y = 40, 2.5 at y = 100 and 1 at y = 250.  None of these is a
    high-precision evaluation of D_0, and none is used as one; D_0 is known in
    closed form and checked in the reviewer notebook.

    The panel layout that produced the shipped grid was not recorded, so this is
    a validation, not a bit-for-bit reproduction.

USAGE
    python scripts/sample_phi_grid.py --check [--tail]   # 45 points, prints the table
    python scripts/sample_phi_grid.py            # full 900-point grid (slow)
"""
import os
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

NY, YLO, YHI = 900, 1.0e-7, 400.0
NU, NV = 48, 24
UMAX_SCALE = 10.0   # multiplies the u cut; see MEASURED AGREEMENT above
GROWTH = 0.1        # relative u-panel width in the non-oscillatory region (y <= pi only)
MAX_U_PANELS = 2_000_000  # refuse layouts beyond this instead of exhausting memory
V_PANELS_PER_PERIOD = 8   # v panels per period 2 pi/y beyond v = 0.1


def _panel(a, b, n):
    x, w = leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (a + b), 0.5 * (b - a) * w


def _j1(x):
    return np.where(np.abs(x) < 1e-8, x / 3.0, (np.sin(x) / x - np.cos(x)) / x)


def _u_edges(y):
    umax = min(max(200.0 * UMAX_SCALE, 400.0 * UMAX_SCALE / max(y, 1e-12)), 1e10)
    graded = [1.0, 1 + 1e-8, 1 + 1e-6, 1 + 1e-4, 1.01, 1.1, 1.5, 2.0, 3.0, 5.0,
              10.0, umax]
    # cos(u y) from the product-to-sum form of j1 j1 oscillates in u with period
    # 2 pi / y. Every interval between the graded edges, including those below
    # u = 10, is subdivided so that no panel is wider than half a period (and
    # never wider than 1). Earlier versions subdivided only u > 10, which left
    # e.g. the single panel [5, 10] spanning ~200 periods at y = 250.
    step = min(np.pi / max(y, 1e-12), 1.0)
    edges = [graded[0]]
    if y <= np.pi:
        # Small y: u_max ~ 4000/y, so width-1 panels out to u_max are unaffordable
        # (1e10 panels at y = 1e-7). Below u = 10 the width-1 split is kept; above
        # it the panel width grows geometrically, w = max(1, GROWTH*u), which
        # follows the non-oscillatory 1/u-type decay while u*y << 1, and is capped
        # at half a period, pi/y, once the oscillation sets in. For y > pi the
        # layout below is unchanged.
        for a, b in zip(graded[:-2], graded[1:-1]):
            n = max(1, int(np.ceil((b - a) / step - 1e-12)))
            edges.extend(a + (b - a) * np.arange(1, n + 1) / n)
        u = graded[-2]
        while u < umax:
            u = min(u + min(np.pi / y, max(1.0, GROWTH * u)), umax)
            edges.append(u)
            if len(edges) > MAX_U_PANELS:
                raise ValueError(f'more than {MAX_U_PANELS} u panels at y={y}')
        return np.array(edges)
    n_est = (umax - graded[0]) / step
    if n_est > MAX_U_PANELS:
        raise ValueError(f'about {n_est:.2e} u panels at y={y} exceed {MAX_U_PANELS}')
    for a, b in zip(graded[:-1], graded[1:]):
        n = max(1, int(np.ceil((b - a) / step - 1e-12)))
        edges.extend(a + (b - a) * np.arange(1, n + 1) / n)
    return np.array(edges)


def _v_nodes(y, nv=NV):
    # l(u,v) = log(1/(uv)) is logarithmically singular at v = 0, so the first
    # decades are resolved by graded panels; Gauss-Legendre on a single panel
    # reaching to v = 0 converges far too slowly and was the dominant error of
    # the earlier layout.
    edges = [0.0, 1e-10, 1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1]
    step = min(2 * np.pi / max(y, 1e-12) / V_PANELS_PER_PERIOD, 1.0)
    t = edges[-1]
    while t < 1.0:
        t = min(t + max(step, 1e-3), 1.0)
        edges.append(t)
    xs, ws = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        x, w = _panel(a, b, nv)
        xs.append(x)
        ws.append(w)
    return np.concatenate(xs), np.concatenate(ws)


MEM_BUDGET_MIB = 256   # working-array budget for one block of u panels


def tail(y, umax):
    """Leading non-oscillatory tail of Phi beyond u = umax (prefactor y included)."""
    from scipy.special import sici
    return (2 * sici(y)[0] - (2 * np.log(umax) + 1) * np.sin(y)) / (y ** 2 * umax ** 2)


def phi(y, nu=NU, nv=NV, block=None, add_tail=False, mem_mib=MEM_BUDGET_MIB):
    """Phi(y) by the rule described in the module docstring.

    The u integral is accumulated over blocks of panels. With block=None the
    block size is chosen so that the float64 working arrays of one block
    (about 8 arrays of nu * len(v nodes) entries per panel) stay within
    `mem_mib` MiB; the nodes and weights are the same as an unblocked
    evaluation. add_tail=True adds tail(y, u_max); without it the value is the
    integral up to u_max only.
    """
    edges = _u_edges(y)
    vv, wv = _v_nodes(y, nv)
    x, w = leggauss(nu)
    if block is None:
        per_panel = 8 * 8 * nu * len(vv)
        block = max(1, int(mem_mib * 2 ** 20 // per_panel))
    total = 0.0
    for k in range(0, len(edges) - 1, block):
        hi = edges[k + 1:k + block + 1]
        lo = edges[k:k + len(hi)]
        uu = (0.5 * (hi - lo)[:, None] * x[None, :] + 0.5 * (hi + lo)[:, None]).ravel()
        wu = (0.5 * (hi - lo)[:, None] * w[None, :]).ravel()
        Lu = np.log((uu + 1) / (uu - 1))
        ell = np.log(1.0 / (uu[:, None] * vv[None, :]))
        up = 0.5 * y * (uu[:, None] + vv[None, :])
        um = 0.5 * y * (uu[:, None] - vv[None, :])
        total += float(wu @ (Lu[:, None] * ell * _j1(up) * _j1(um)) @ wv)
    out = y * total
    return out + tail(y, edges[-1]) if add_tail else out


def _shipped():
    here = os.path.dirname(os.path.abspath(__file__))
    return np.load(os.path.join(here, "..", "figs", "phi_grid.npy"))


def check(stride=20):
    """Compare with the shipped grid on every `stride`-th abscissa."""
    g = _shipped()
    ys, ph = g[0], g[1]
    idx = list(range(0, len(ys), stride))
    new = np.array([phi(ys[i], add_tail="--tail" in sys.argv) for i in idx])
    old = ph[idx]
    env = np.array([np.max(np.abs(ph[max(0, i - 25):i + 25])) for i in idx])
    scaled = np.abs(new - old) / np.maximum(env, 1e-30)
    print("  y            shipped          recomputed       scaled diff")
    for i, a, b, s in zip(idx, old, new, scaled):
        print("%12.4e %+16.8e %+16.8e %12.2e" % (ys[i], a, b, s))
    print("\nscaled by the local envelope: median %.2e, max %.2e at y = %.3e"
          % (np.median(scaled), scaled.max(),
             ys[idx[int(np.argmax(scaled))]]))


if __name__ == "__main__":
    if "--check" in sys.argv:
        check()
    else:
        ys = np.exp(np.linspace(np.log(YLO), np.log(YHI), NY))
        out = np.array([phi(y) for y in ys])
        here = os.path.dirname(os.path.abspath(__file__))
        np.save(os.path.join(here, "..", "figs", "phi_grid_recomputed.npy"),
                np.stack([ys, out]))
        w = np.gradient(np.log(ys))
        exact = np.pi ** 3 / 18 * np.log(2) - np.pi / 4 * 1.2020569031595943
        print("wrote figs/phi_grid_recomputed.npy  (%d points)" % NY)
        print("  int Phi dy/y = %.9f (sum Phi*gradient(log y)), exact D_0 = %.9f, relative %.2e"
              % (np.sum(out * w), exact, abs(np.sum(out * w) / exact - 1)))
