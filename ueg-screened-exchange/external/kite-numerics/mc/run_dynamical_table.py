"""Pipeline for the dynamical kite table:
per-q raw CSV + seeds + metadata for each (rs, zeta) point.

Usage:
  python mc/run_dynamical_table.py --rs 1 5 --zeta 0 --samples 150000 \\
         --seed 41 --out outputs/dynamical_table
The energy assembly uses trapezoidal integration plus the analytic lower
endpoint term. Grid, tail, and frequency systematics must be assessed with
separate convergence runs.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"mc"))
sys.path.insert(0, str(ROOT/"analytic"))
from xc_mc import estimate, D0_EXACT  # noqa: E402

C_MRY = 6e3/np.pi**3


def grid(dense=True):
    if dense:
        return np.unique(np.concatenate([np.geomspace(0.006, 6, 38),
                                         np.linspace(1.85, 2.3, 8),
                                         np.geomspace(6.8, 10, 3)]))
    return np.geomspace(0.02, 6, 20)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rs", nargs="+", type=float, required=True)
    ap.add_argument("--zeta", type=float, default=0.0)
    ap.add_argument("--samples", type=int, default=150_000)
    ap.add_argument("--n-batch", type=int, default=32,
                    help="number of Monte Carlo batches (default: 32)")
    ap.add_argument("--seed", type=int, default=41)
    ap.add_argument("--out", default="data/dynamical_table")
    ap.add_argument("--coarse", action="store_true")
    ap.add_argument(
        "--save-batch-means", action="store_true",
        help=("write batch_means_rs*_z*.csv for covariance and heavy-tail "
              "diagnostics; this does not alter the estimator"))
    a = ap.parse_args()
    if a.n_batch < 2:
        ap.error("--n-batch must be at least 2")
    os.makedirs(a.out, exist_ok=True)
    g = grid(dense=not a.coarse)
    samples_used = max(1, a.samples//a.n_batch)*a.n_batch
    if a.samples//a.n_batch < 100_000:
        print(f"WARNING: {a.samples} samples -> <1e5 per batch; "
              "demonstration statistics only (production rule: n/n_b>=1e5)")
    # Spin loop: identical treatment to run_scan.py dyn --
    # spin-sigma Fermi radius x=(1+sigma*zeta)^(1/3), weight x^3/2, and
    # xspin=x passed to the estimator; independent seed stream per spin.
    spins = ([(0, 1.0)] if a.zeta == 0.0 else
             [(i, (1 + s*a.zeta)**(1/3.))
              for i, s in enumerate((+1.0, -1.0)) if (1 + s*a.zeta) > 1e-12])
    for rs in a.rs:
        acc = np.zeros((len(g), 2))
        batch_acc = (np.zeros((len(g), a.n_batch))
                     if a.save_batch_means else None)
        for isp, x in spins:
            wsp = 0.5*x**3 if a.zeta != 0.0 else 1.0
            rng = np.random.default_rng(a.seed + isp*10_000_019)
            for iq, q in enumerate(g):
                result = estimate(
                    q, a.samples, rng, n_batch=a.n_batch,
                    return_batch_means=a.save_batch_means, kernel="diff",
                    rs=rs, zeta=a.zeta, xspin=x)
                if a.save_batch_means:
                    mu, e, batch_means = result
                    batch_acc[iq] += wsp*batch_means
                else:
                    mu, e = result
                acc[iq, 0] += wsp*mu
                acc[iq, 1] += (wsp*e)**2
        rows = np.column_stack([g, acc[:, 0], np.sqrt(acc[:, 1])])
        tag = f"rs{rs:g}_z{a.zeta:g}"
        np.savetxt(f"{a.out}/per_q_{tag}.csv", rows, delimiter=",",
                   header="q,dX,err", comments="")
        batch_name = None
        if a.save_batch_means:
            batch_name = f"batch_means_{tag}.csv"
            batch_rows = np.column_stack([
                np.repeat(g, a.n_batch),
                np.tile(np.arange(1, a.n_batch + 1), len(g)),
                batch_acc.reshape(-1),
            ])
            np.savetxt(
                f"{a.out}/{batch_name}", batch_rows, delimiter=",",
                header="q,batch,dX_batch", comments="",
                fmt=["%.17g", "%d", "%.17g"])
        ys, es = rows[:, 1], rows[:, 2]
        eK = C_MRY*(D0_EXACT + np.trapezoid(ys, g) + 0.5*ys[0]*g[0])
        tw = np.zeros_like(g)
        tw[1:] += 0.5*np.diff(g); tw[:-1] += 0.5*np.diff(g)
        tw[0] += 0.5*g[0]
        sK = C_MRY*np.sqrt(np.sum((tw*es)**2))
        print(f"rs={rs} zeta={a.zeta}: e_kite = {eK:.2f}({sK:.2f}) mRy "
              f"[MC component only]", flush=True)
        metadata = {"command": " ".join(sys.argv),
                    "samples_requested_per_q_per_spin": a.samples,
                    "samples_used_per_q_per_spin": samples_used,
                    "samples_per_batch": samples_used//a.n_batch,
                    "n_batch": a.n_batch,
                    "seed": a.seed,
                    "seed_stream": "default_rng(seed + spin_index*10000019)",
                    "cross_rs_correlation": (
                        "The RNG is reset for every rs, so different-density "
                        "outputs use common random numbers and are correlated. "
                        "Do not fit them as independent observations unless "
                        "batch-level covariance is retained or distinct seeds are used."),
                    "zeta": a.zeta, "rs": rs,
                    "spin_loop": [list(sp) for sp in spins],
                    "grid_points": len(g),
                    "q_grid": [float(x) for x in g],
                    "batch_means_file": batch_name,
                    "batch_means_columns": (
                        ["q", "batch", "dX_batch"] if batch_name else None),
                    "batch_means_note": (
                        "Spin-weighted raw batch means before any "
                        "median-of-means reporting guard. Batch IDs are "
                        "aligned across density points generated by this "
                        "command because the RNG is reset for every rs."
                        if batch_name else
                        "Not saved; rerun with --save-batch-means for "
                        "covariance and heavy-tail diagnostics."),
                    "xi_grid": "corner-refined log panels [1e-7,1]*s0, "
                               "nseg=10(+8 corner)*24 GL nodes (Python)",
                    "note": "finite-sample; errors are MC batch-scatter only"}
        with open(f"{a.out}/metadata_{tag}.json", "w") as fh:
            json.dump(metadata, fh, indent=1)
    return a.out


if __name__ == "__main__":
    main()
