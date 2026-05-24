#!/usr/bin/env python3
"""Figures for the spectral-kernel paper.

fig_dispersion.pdf : the dispersion of the approximate screening pole and the
                     frequency kernel p_loc(s_lambda) it is integrated against.
fig_Slambda.pdf    : the one-dimensional reduction S(Lambda) and its large-Lambda
                     Mellin asymptote.

Run from any directory; writes next to this script in figs/.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import sici

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dataout import write as write_dat

plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 0.7, "lines.linewidth": 1.2,
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

# --------------------------------------------------------------- dispersion
s_lambda = lambda x, lam: 0.5*np.sqrt(1.0 + lam*x**-4)
p_loc  = lambda m: 2*np.pi*m/(1 + 2*m)**2

fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.5))

x = np.logspace(-1.0, 1.3, 400)
# the pole frequency itself, omega = Lambda q, in units of q_*^2: a dispersion
# running from the plasma scale to the free-electron limit q^2/2
for lam, ls in [(1.0, "-"), (0.5, "--"), (0.25, ":")]:
    ax[0].loglog(x, 0.5*x**2*np.sqrt(1.0 + lam*x**-4), ls, color="C0",
                 label=r"$\lambda=%g$" % lam)
ax[0].loglog(x, 0.5*x**2, color="0.35", lw=0.9)
ax[0].text(0.98, 0.40, "free-electron limit\n" + r"$\omega=q^{2}/2$",
           transform=ax[0].transAxes, ha="right", fontsize=7.5, color="0.25")
ax[0].axhline(0.5, color="C3", lw=0.8, ls="-.")
ax[0].text(0.105, 0.135, r"plasma limit ($\lambda=1$): $\omega/q_*^2=1/2$",
           fontsize=7.5, color="C3")
ax[0].set_xlabel(r"$q/q_*$")
ax[0].set_ylabel(r"approximate pole $\omega_\lambda/q_*^{2}$")
ax[0].set_ylim(0.1, 40)
ax[0].legend(frameon=False, loc="upper left")

m = np.logspace(-1.3, 1.6, 400)
ax[1].semilogx(m, p_loc(m), color="C0")
ax[1].axvline(0.5, color="0.35", lw=0.8)
ax[1].plot([0.5], [p_loc(0.5)], "o", ms=3.5, color="C3")
ax[1].text(0.55, p_loc(0.5)*0.93, r"$s_\Lambda=\frac{1}{2}$", fontsize=7.5, color="0.25")
ax[1].set_xlabel(r"$s_\Lambda=\Lambda/q$")
ax[1].set_ylabel(r"$p_{\rm loc}(s_\Lambda)=2\pi s_\Lambda/(1+2s_\Lambda)^{2}$")
ax[1].set_ylim(0, 0.95)
ax[0].set_title("(a) Approximate screening pole", fontsize=9)
ax[1].set_title("(b) Local frequency kernel", fontsize=9)
fig.subplots_adjust(wspace=0.38)
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig_dispersion.pdf"))
plt.close(fig)

# ------------------------------------------------------------ ASCII export
write_dat("fig5_left_dispersion.dat",
          [("q/q_*", x)]
          + [("omega/q_*^2 (lambda=%g)" % lam,
              0.5 * x ** 2 * np.sqrt(1.0 + lam * x ** -4))
             for lam in (1.0, 0.5, 0.25)]
          + [("free electron q^2/2", 0.5 * x ** 2)],
          "JPA Fig. 7 left: pole frequency omega/q_*^2 against q/q_*",
          ["omega_lambda(q)/q_*^2 = (q/q_*)^2 sqrt(1 + lambda (q_*/q)^4)/2",
           "plasma-frequency end: omega -> sqrt(lambda) q_*^2/2"])
write_dat("fig5_right_weight.dat",
          [("s_lambda=Lambda/q", m), ("p_loc(s_lambda)", p_loc(m))],
          "JPA Fig. 7 right: the large-q kernel weight p_loc(s_lambda) = 2 pi s_lambda/(1+2 s_lambda)^2",
          ["maximum at s_lambda = 1/2, p_loc(1/2) = pi/4 = %.12f" % (np.pi / 4)])

# ----------------------------------------------------------------- S(Lambda)
g = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), "phi_grid_corrected.npy"))
ys, Phi_y = g[0], g[1]
w = np.gradient(np.log(ys))
D0 = np.pi**3/18*np.log(2) - np.pi/4*1.2020569031595943

def Xi(z):
    si, ci = sici(z)
    return z*(ci*np.sin(z) + (np.pi/2 - si)*np.cos(z))

def S(L):
    return np.sum(Phi_y*Xi(L*ys)*w)/D0

L = np.logspace(-2, 2.6, 260)
Sv = np.array([S(l) for l in L])
C1 = np.pi/(4*D0); C0 = -C1/4        # JPA eq:largeMu (Appendix F)

fig, ax = plt.subplots(figsize=(3.4, 2.5))
ax.semilogx(L, Sv, color="C0", label=r"$S(\Lambda)$")
Lb = L[L > 3]
ax.semilogx(Lb, 1 + (C1*np.log(Lb) + C0)/Lb, "--", color="C3",
            label=r"$1+(C_1^\Lambda\log\Lambda+C_0^\Lambda)/\Lambda$")
ax.axhline(1.0, color="0.35", lw=0.8)
ax.text(0.012, 1.06, r"$S(\infty)=1$: Onsager", fontsize=7.5, color="0.25")
ax.set_xlabel(r"$\Lambda$")
ax.set_ylabel(r"$S(\Lambda)=\mathcal{K}[b_\Lambda]\,/\,e_{2x}$")
ax.set_ylim(0, 2.6)
ax.legend(frameon=False, loc="upper left", fontsize=7.5)
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig_Slambda.pdf"))
plt.close(fig)
print("wrote figs/fig_dispersion.pdf and figs/fig_Slambda.pdf")

write_dat("fig4_slambda.dat",
          [("Lambda", L), ("S(Lambda)", Sv),
           ("1+(C1 log L + C0)/L", 1 + (C1 * np.log(L) + C0) / L)],
          "JPA Fig. 6: the reduced kernel energy S(Lambda) and its large-Lambda form",
          ["S = K[b_Lambda]/e_2x, evaluated from figs/phi_grid_corrected.npy",
           "C1 = pi/(4 D_0) = %.12f, C0 = -C1/4 = %.12f" % (C1, C0),
           "the paper draws the asymptote only for Lambda > 3",
           "S(infinity) = 1 is the Onsager value"])
