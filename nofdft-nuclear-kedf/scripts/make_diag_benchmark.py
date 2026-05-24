#!/usr/bin/env python3
"""
Regenerate the constant-kF ^16O diagonalization-benchmark numbers (manuscript
Table `tab:diag` and the ITE-vs-diagonalization comparison) from the raw data in
results/raw/diag_benchmark/.

Reports, per (lambda, alpha) diagonalization run: outer iterations (history
rows), converged energy, and the auxiliary gap; and for the fixed-step ITE runs:
the energy offset and relative L2 density difference from the (3,1.0)
diagonalization endpoint.
"""
import numpy as np, os, re, argparse

RUNS = [("l1_a1", 1.0, 1.0), ("l1_a03", 1.0, 0.3), ("l3_a1", 3.0, 1.0), ("l10_a1", 10.0, 1.0)]


def energy(path):
    m = re.search(r"E_total_NeqZ_MeV\s+(-?\d+\.\d+E[+-]\d+)", open(path).read())
    return float(m.group(1))


def n_iter(hist):
    return sum(1 for ln in open(hist) if ln.strip() and not ln.startswith("#"))


def final_gap(hist):
    import numpy as np
    h = np.loadtxt(hist)
    return float(h[-1, 5])            # col 6 (0-indexed 5) = gap(eps1-eps0)


def converged(summ):
    for ln in open(summ):
        if ln.strip().startswith("converged"):
            return ln.split()[-1] == "T"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="results/raw/diag_benchmark")
    args = ap.parse_args()
    R = args.raw

    print("Rearranged diagonalization, constant-kF ^16O (nr=250, kF=0.80, nmu=24):")
    print("  %-8s %5s %5s %8s %14s %10s %6s" % ("tag", "lam", "alp", "iters", "E [MeV]", "gap[MeV]", "conv"))
    for tag, lam, alp in RUNS:
        h = os.path.join(R, "diag_%s_h.dat" % tag)
        summ = os.path.join(R, "diag_%s_s.dat" % tag)
        conv = converged(summ)
        if conv:
            print("  %-8s %5.0f %5.1f %8d %14.7f %10.3f %6s"
                  % (tag, lam, alp, n_iter(h), energy(summ), final_gap(h), "yes"))
        else:
            print("  %-8s %5.0f %5.1f %8s %14s %10s %6s"
                  % (tag, lam, alp, ">=%d" % n_iter(h), "--", "--", "no"))

    # ITE vs (3,1.0) diagonalization
    Ed = energy(os.path.join(R, "diag_l3_a1_s.dat"))
    dd = np.loadtxt(os.path.join(R, "diag_l3_a1_d.dat"))
    r = dd[:, 0]
    print("\nFixed-step ITE vs (3,1.0) diagonalization (E_diag = %.7f MeV):" % Ed)
    for nit in (8000, 20000):
        Ei = energy(os.path.join(R, "ite_%d_s.dat" % nit))
        di = np.loadtxt(os.path.join(R, "ite_%d_d.dat" % nit))
        dii = np.interp(r, di[:, 0], di[:, 1])
        rel = np.sqrt(np.sum(r**2 * (dd[:, 1] - dii) ** 2) / np.sum(r**2 * dd[:, 1] ** 2))
        print("  %5d steps:  dE = %8.3f keV    rel-L2 density = %.2e" % (nit, (Ei - Ed) * 1e3, rel))

    # independently monitored EL residual at the (3,1.0) diagonalization density
    er = os.path.join(R, "elres_diag_e.dat")
    if os.path.exists(er):
        m = re.search(r"residual_phi_L2\s+(-?\d+\.\d+E[+-]\d+)", open(er).read())
        if m:
            print("\nIndependently monitored EL residual at the (3,1.0) density: %.3f MeV"
                  % float(m.group(1)))


if __name__ == "__main__":
    main()
