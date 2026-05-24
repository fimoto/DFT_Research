#!/usr/bin/env python3
"""
Regenerate the Gaussian-process figure (manuscript Fig. 6) from the outputs of
gp_fit.py:

  left  : fixed-active-set MAP estimate of delta_g(y), with the local Gaussian
          +/-2-sigma parameter band                         (eg_GP_band.dat)
  right : fixed-density A=140 NLF prediction at the KS density, with the same
          covariance propagated to an approximate +/-2-sigma NLF band
                                                              (gp_holdout_A140_nlf.dat)

Run gp_fit.py first (in scripts/, it writes these files to the working
directory); then point --gp at that directory.
"""
import numpy as np, os, argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gp", default=".", help="directory with eg_GP_band.dat, gp_holdout_A140_nlf.dat, gp_shift.txt")
    ap.add_argument("--out", default="fig6.png")
    args = ap.parse_args()

    band = np.loadtxt(os.path.join(args.gp, "eg_GP_band.dat"))   # y, dg, sigma
    y, dg, sd = band[:, 0], band[:, 1], band[:, 2]
    shift = None
    sp = os.path.join(args.gp, "gp_shift.txt")
    if os.path.exists(sp):
        shift = float(open(sp).read().strip())

    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.0))

    a = ax[0]
    a.fill_between(
        y, dg - 2 * sd, dg + 2 * sd,
        color="tab:blue", alpha=0.25,
        label=r"fixed-active-set $\pm2\sigma$",
    )
    a.plot(y, dg, "b-", lw=1.8, label="fixed-active-set MAP")
    a.axhline(0, color="k", lw=0.6)
    a.set_xlim(0, 12); a.set_xlabel("y"); a.set_ylabel(r"$\delta g(y)$")
    ttl = "GP kernel correction"
    if shift is not None:
        ttl += r"  ($s=%.3f$)" % shift
    a.set_title(ttl); a.legend(loc="upper right", fontsize=9, frameon=False)

    b = ax[1]
    ho = os.path.join(args.gp, "gp_holdout_A140_nlf.dat")
    if os.path.exists(ho):
        d = np.loadtxt(ho)          # r, NLF_model, NLF_KS, 2sigma
        r, nm, nk, band2 = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
        b.fill_between(
            r, nm - band2, nm + band2,
            color="tab:red", alpha=0.2,
            label=r"fixed-active-set $\pm2\sigma$",
        )
        b.plot(r, nk, "k-", lw=1.4, label="KS reference")
        b.plot(r, nm, "r-", lw=1.8, label="fixed-density GP prediction")
        b.set_xlim(0, 9); b.set_ylim(0, 1.05)
        b.set_xlabel("r [fm]"); b.set_ylabel("species NLF")
        b.set_title(r"$A=140$: fixed-density prediction")
        b.legend(loc="lower left", fontsize=9, frameon=False)
    else:
        b.text(0.5, 0.5, "run gp_fit.py first", ha="center", va="center", transform=b.transAxes)

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
