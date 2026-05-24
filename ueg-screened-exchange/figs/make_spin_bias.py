"""Reanalyse the stored spin map against endpoint exchange interpolation.

No integral is evaluated and no density asymptote is fitted. For every term
and independent scramble, the residual is formed from the three paired spin
nodes before its variance is computed. The two terms use independent scrambles.
Run from any directory: python figs/make_spin_bias.py
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'production'))
from common import check_source_compatibility


def exchange_interpolation(zeta):
    z = np.asarray(zeta, dtype=float)
    return ((1 + z)**(4 / 3) + (1 - z)**(4 / 3) - 2) / (2**(4 / 3) - 2)


def raw_records():
    config = json.loads((ROOT / "production/config.json").read_text())
    lock = json.loads((ROOT / "production/reference_lock.json").read_text())["sha256"]
    rows = {}
    for rs in config["map"]["rs"]:
        for zeta in config["map"]["zeta"]:
            key = f"map/rs_{rs:g}_z_{zeta:g}.json"
            path = ROOT / "production/raw" / key
            if hashlib.sha256(path.read_bytes()).hexdigest() != lock[key]:
                raise ValueError(f"Changed raw reference: {key}")
            row = json.loads(path.read_text())
            check_source_compatibility(row["metadata"]["sources"])
            if (row["rs"], row["zeta"]) != (rs, zeta):
                raise ValueError(f"Wrong raw-node identity: {key}")
            rows[rs, zeta] = row
    settings = [row["settings"] for row in rows.values()]
    if any(setting != settings[0] for setting in settings):
        raise ValueError("Spin-map nodes do not share the same sampling settings")
    return config, rows


def summarize():
    config, rows = raw_records()
    output = []
    for rs in config["map"]["rs"]:
        unpolarized = np.asarray(rows[rs, 0]["term_replicates_ry"])
        polarized = np.asarray(rows[rs, 1]["term_replicates_ry"])
        for zeta in config["map"]["zeta"]:
            terms = np.asarray(rows[rs, zeta]["term_replicates_ry"])
            fx = float(exchange_interpolation(zeta))
            # Reuse the recorded common random numbers across spin nodes.
            residual = 1000 * (terms - unpolarized - fx * (polarized - unpolarized))
            term_means = residual.mean(axis=0)
            mean = term_means.sum()
            se = np.sqrt(residual.var(axis=0, ddof=1).sum() / len(residual))
            output.append([rs, zeta, fx, mean, se, *term_means])
    return np.asarray(output)


def main():
    data = summarize()
    data_dir = ROOT / "figs/data"
    data_dir.mkdir(exist_ok=True)
    np.savetxt(data_dir / "fig_spin_bias.csv", data, delimiter=",", fmt="%.17g",
               header="rs,zeta,fx,residual_mry,se_mry,term1_residual_mry,term2_residual_mry",
               comments="")

    plt.rcdefaults()
    plt.rcParams.update({"font.family": "serif", "font.size": 8.5,
                         "mathtext.fontset": "dejavuserif", "pdf.fonttype": 42,
                         "axes.linewidth": 0.7, "xtick.direction": "in",
                         "ytick.direction": "in"})
    fig, ax = plt.subplots(figsize=(3.4, 2.8), layout="constrained")
    for zeta, color, marker in [(0.4, "#0072B2", "s"), (0.8, "#D55E00", "o")]:
        series = data[data[:, 1] == zeta]
        ax.errorbar(series[:, 0], series[:, 3], yerr=series[:, 4], color=color,
                    marker=marker, markersize=3.2, lw=1.1, elinewidth=0.8,
                    capsize=2, label=rf"$\zeta={zeta:g}$")
    ax.axhline(0, color="0.55", lw=0.65, zorder=0)
    ax.set(xscale="log", xlim=(0.8, 1650), xlabel=r"Density parameter $r_s$",
           ylabel=r"$\delta e_{\rm kite}$ (mRy)")
    ax.legend(loc="lower right", frameon=False, fontsize=8)
    ax.tick_params(top=True, right=True)
    stamp = datetime(2026, 9, 27, tzinfo=timezone.utc)
    fig.savefig(ROOT / "figs/fig_spin_bias.pdf",
                metadata={"Creator": "Paired spin-map interpolation residual",
                          "CreationDate": stamp, "ModDate": stamp})
    plt.close(fig)

    lines = ["% Generated from stored paired spin-map replicates; uncertainties are replicate standard errors."]
    for rs, name in [(1, "One"), (10, "Ten"), (60, "Sixty")]:
        row, = data[(data[:, 0] == rs) & (data[:, 1] == 0.8)]
        lines.append(rf"\newcommand{{\BiasRs{name}}}{{{row[3]:+.3f}}}")
        lines.append(rf"\newcommand{{\BiasRs{name}SE}}{{{row[4]:.3f}}}")
    target = ROOT / "spin/generated/spin_bias_numbers.tex"
    target.parent.mkdir(exist_ok=True)
    target.write_text("\n".join(lines) + "\n")
    print(f"Wrote {len(data)} spin-map residual rows, one figure and six numerical macros.")
    print("Reanalysis of saved replicate means; no new integral or HPC evaluation.")


if __name__ == "__main__":
    main()
