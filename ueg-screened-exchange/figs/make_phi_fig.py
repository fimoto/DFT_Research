"""fig_Phi.pdf : the geometric profile Phi(y) and its two endpoints.

Phi is the y-space profile of the master formula; its small-y endpoint fixes the
constants a_0 = 1/2 and b_0 = gamma/2 - 3/8 that produce the double pole of
F(-s) at s=1, and its large-y behaviour carries the oscillatory tail whose
non-oscillatory part gives the simple pole at s=2.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dataout import write as write_dat

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm",
    "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 0.6, "lines.linewidth": 1.2,
})

# phi_grid_corrected.npy is the supplied phi_grid.npy with every abscissa below
# y = 60 recomputed (scripts/build_phi_grid_corrected.py merges the y < 4 and
# 4 < y < 60 runs); the supplied grid was off by up to 3.9e-5 in 4 < y < 60, which
# y^2 turns into 0.12 in the right panel. The 78 points above y = 60 were recomputed too.
g = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), "phi_grid_corrected.npy"))
ys, Phi = g[0], g[1]
a0 = 0.5
b0 = 0.5 * np.euler_gamma - 0.375

fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.9))

m = (ys > 3e-3) & (ys < 1.5)
ax[0].plot(ys[m], Phi[m], "o", ms=2.2, color="C0", markevery=6,
           label=r"$\Phi(y)$")
ax[0].plot(ys[m], a0 * ys[m] * np.log(ys[m]) + b0 * ys[m], "-", color="C3",
           lw=1.0, label=r"$a_0\,y\log y+b_0\,y$")
ax[0].set_xscale("log")
ax[0].set_xlabel(r"$y$")
ax[0].set_ylabel(r"$\Phi(y)$")
ax[0].legend(frameon=False, loc="upper left")

# The grid has only 4-6 points per period of the sin y oscillation here, and the
# oscillation of y^2 Phi has a half-range of 11-14, so a plain average over the
# grid points inside one period leaves a visible residue (up to 1.7). Instead the
# grid values are joined by a cubic spline s and averaged over one period twice
# (a triangular window of half-width 2 pi), which also removes the residue of an
# oscillation whose amplitude grows linearly in y. With H the second
# antiderivative of s, the double mean is exactly
#     [H(t + 2 pi) - 2 H(t) + H(t - 2 pi)] / (2 pi)^2 .
# The curve ends at the last grid point below y = 60 and the average is drawn
# only where its window lies inside the data, so nothing is extrapolated.
# Interpolation error measured at grid midpoints: 0.06 at y = 50.1 and 0.19 at
# y = 56.6 in y^2 Phi.
from scipy.interpolate import CubicSpline
ms = (ys > 0.0) & (ys < 60.0)
spl = CubicSpline(ys[ms], ys[ms] ** 2 * Phi[ms])
H = spl.antiderivative(2)
per = 2 * np.pi   # period of the sin y oscillation of y^2 Phi
y_last = ys[ms][-1]

y2 = np.linspace(4.0, y_last, 2000)
P2 = spl(y2)
ax[1].plot(y2, P2, "-", color="0.75", lw=0.7, label=r"$y^{2}\Phi(y)$")
ta = np.linspace(per, y_last - per, 600)
run_a = (H(ta + per) - 2 * H(ta) + H(ta - per)) / per ** 2
ax[1].plot(ta, run_a, "-", color="C0", lw=1.2,
           label=r"twice-averaged")
run = np.interp(y2, ta, run_a, left=np.nan, right=np.nan)
ax[1].axhline(2 * np.pi * np.log(2), color="C3", ls="--", lw=1.0,
              label=r"$2\pi\log 2$")
ax[1].set_ylim(-14, 22)
ax[1].set_xlabel(r"$y$")
ax[1].set_ylabel(r"$y^{2}\Phi(y)$")
ax[1].legend(frameon=False, loc="upper center", ncol=3,
             bbox_to_anchor=(0.5, 1.22), columnspacing=1.0,
             handletextpad=0.4, borderaxespad=0.0, fontsize=7.5)

ax[0].set_title("(a) Small-$y$ expansion", fontsize=9, pad=28)
ax[1].set_title("(b) Oscillatory large-$y$ tail", fontsize=9, pad=28)
fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig_Phi.pdf"))
print("wrote figs/fig_Phi.pdf")

# ------------------------------------------------------------ ASCII export
write_dat("fig2_left_profile.dat",
          [("y", ys[m]), ("Phi", Phi[m]),
           ("a0*y*log(y)+b0*y", a0 * ys[m] * np.log(ys[m]) + b0 * ys[m])],
          "JPA Fig. 5 left: the geometric profile Phi(y) and its small-y endpoint",
          ["a0 = 1/2, b0 = gamma/2 - 3/8 = %.12f" % b0,
           "grid: figs/phi_grid_corrected.npy, restricted to the plotted window"])
write_dat("fig2_right_tail.dat",
          [("y", y2), ("y^2*Phi (cubic spline through the grid)", P2),
           ("period_average (nan where the window leaves the data)", run)],
          "JPA Fig. 5 right: the large-y tail of y^2 Phi and its twice-applied period average",
          ["the average is the one-period mean applied twice (triangular window, half-width 2 pi)",
           "constant line: 2 pi log 2 = %.12f" % (2 * np.pi * np.log(2)),
           "grid: figs/phi_grid_corrected.npy (all points below y = 60 recomputed)"])
