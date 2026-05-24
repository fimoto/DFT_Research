#!/usr/bin/env python3
"""
Generate the selected-path trial densities for the shared-fit energy-rise
penalty (manuscript Sec. "Learning the kernel correction", item 6).

The paths follow the density change that produced the central-density artifact
in preliminary self-consistent calculations: the A=80 hollow-bubble endpoint of
the preliminary no-path shared fit eg_SHNP14 (scripts/path_data/a80_artifact_endpoint.dat, imaginary-time
relaxation from the KS density at the production ledger). With

    delta rho(r) = rho_bubble(r) - rho_KS(r)          (particle-number neutral)
    rho_j(r)     = rho_KS(r) + xi_j * delta rho(r),   xi_j in {0.35, 0.60, 1.00}

each trial density is a convex combination of two nonnegative densities, hence
nonnegative, and conserves the particle number exactly.

The trial densities are written to scripts/path_data/path_A80_xi<...>.dat as
two-column (r, rho_species) files. kernel_fit.py recomputes the affine
energy-rise rows (q, dE_base) from these files at run time, so the rows are
reproducible for any kernel basis.
"""
import os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "scripts", "path_data")
XIS = (0.35, 0.60, 1.00)


def main():
    ks = np.loadtxt(os.path.join(ROOT, "inputs", "ks", "ks_A80_nr250.dat"))
    bb = np.loadtxt(os.path.join(ROOT, "scripts", "path_data", "a80_artifact_endpoint.dat"))
    r = ks[:, 0]
    if not np.allclose(r, bb[:, 0]):
        raise SystemExit("grid mismatch between KS and bubble densities")
    drho = bb[:, 1] - ks[:, 1]
    os.makedirs(OUT, exist_ok=True)
    for xi in XIS:
        rho = ks[:, 1] + xi * drho
        if rho.min() < -1e-14:
            raise SystemExit("negative trial density at xi=%.2f" % xi)
        fn = os.path.join(OUT, "path_A80_xi%.2f.dat" % xi)
        np.savetxt(fn, np.column_stack([r, np.maximum(rho, 0.0)]),
                   header="r[fm] rho_species[fm^-3]  (A80 trial path, xi=%.2f; "
                          "see make_paths.py)" % xi)
        dr = r[1] - r[0]
        N = np.sum(4 * np.pi * r ** 2 * dr * rho)
        print("wrote %s  (N_species = %.6f)" % (os.path.basename(fn), N))


if __name__ == "__main__":
    main()
