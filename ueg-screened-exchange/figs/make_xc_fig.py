"""fig_xckernel.pdf : the static geometric kernel X_c(q).

Plots the saved kernel figs/xc_grid.npy against the two closed-form endpoints.
The grid was produced by scripts/sample_xc_grid.py in the separate code archive,
which records the sampler settings; this script only draws it. Endpoints: the
small-q slope c_0 = pi(2 log 2 - 1)/8 and the large-q tail C_inf/q^4 with
C_inf = 2 pi/9.  The linearity at small q is the visible content of c_1 = 0.
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

c0 = np.pi / 8 * (2 * np.log(2) - 1)
Cinf = 2 * np.pi / 9

g = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), "xc_grid.npy"))
q, X, E = g[0], g[1], g[2]

fig, (ax, bx) = plt.subplots(2, 1, figsize=(3.4, 4.7))

# Upper: both closed-form endpoints on logarithmic axes, where they are visible
# over the whole range.
ax.errorbar(q, X, yerr=E, fmt="o", ms=3.0, lw=0.8, capsize=1.5, color="C0",
            mfc="w", label=r"$X_c(q)$, Monte Carlo")
qq = np.linspace(0.03, 1.0, 60)
ax.plot(qq, c0 * qq, "--", color="C3", lw=1.0, label=r"$c_0q$")
qt = np.linspace(2.2, 5.5, 60)
ax.plot(qt, Cinf / qt ** 4, "-", color="0.35", lw=1.0,
        label=r"$C_\infty/q^{4}$")
ax.axvline(2.0, color="0.4", ls=":", lw=0.8)
ax.set_xscale("log")
ax.set_yscale("log")
from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
ax.xaxis.set_major_locator(FixedLocator([0.05, 0.1, 0.5, 1, 2, 5]))
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel(r"$q$ (in units of $k_F$)")
ax.set_ylabel(r"$X_c(q)$")
ax.legend(frameon=False, loc="lower center", ncol=3,
          bbox_to_anchor=(0.5, 1.02), columnspacing=1.0,
          handletextpad=0.4, borderaxespad=0.0, fontsize=7.5)

# Lower: compensated view.  Above, the error bars are 1-5 % of the values;
# dividing by q puts the intercept at c_0 and turns a quadratic term into a
# slope, so both are read directly with visible errors.
m = q <= 1.0
bx.errorbar(q[m], X[m] / q[m], yerr=E[m] / q[m], fmt="o", ms=3.0, lw=0.9,
            capsize=1.5, color="C0", mfc="w")
bx.axhline(c0, color="C3", ls="--", lw=1.0)
bx.text(0.72, c0 + 0.003, r"$c_0=\frac{\pi}{8}(2\log 2-1)$", fontsize=7.5,
        color="C3")
bx.set_xlim(0, 1.05)
bx.set_xlabel(r"$q$ (in units of $k_F$)")
bx.set_ylabel(r"$X_c(q)/q$")

ax.set_title("(a) Static kernel and endpoint forms", fontsize=9, pad=30)
bx.set_title("(b) Small-momentum slope", fontsize=9)
fig.tight_layout(h_pad=0.8)
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig_xckernel.pdf"))
print("wrote figs/fig_xckernel.pdf")

# ------------------------------------------------------------ ASCII export
write_dat("fig6_xckernel.dat",
          [("q (dimensionless)", q), ("X_c", X), ("error", E),
           ("X_c/q", X / q), ("error/q", E / q)],
          "JPA Fig. 8: the static geometric kernel X_c(q), Monte Carlo",
          ["grid: figs/xc_grid.npy; original code archive: scripts/sample_xc_grid.py (seed 3)",
           "errors are standard errors of the mean over 4e6 points per q",
           "upper panel plots columns 1:2:3, lower panel 1:4:5 for q <= 1",
           "endpoints, drawn as lines in the figure:",
           "  c_0 q with c_0 = pi(2 log 2 - 1)/8 = %.12f  (q < 1)" % c0,
           "  C_inf/q^4 with C_inf = 2 pi/9 = %.12f      (q > 2)" % Cinf])
