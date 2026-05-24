#!/usr/bin/env python3
"""Regenerate all six manuscript figures from the audited release data."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results" / "raw"
KS = ROOT / "results" / "reference"
KERNELS = ROOT / "kernels"
DIAGNOSTICS = ROOT / "results" / "reference" / "diagnostics_v2"


def save_png(fig: plt.Figure, output: Path) -> None:
    """Write a complete PNG before atomically replacing the public output."""
    temporary = output.with_name(f".{output.name}.tmp")
    fig.savefig(temporary, format="png", dpi=220)
    if not temporary.read_bytes().endswith(b"\x00\x00\x00\x00IEND\xaeB`\x82"):
        raise OSError(f"incomplete PNG generated at {temporary}")
    temporary.replace(output)


def _lindhard_l(x: np.ndarray) -> np.ndarray:
    """Static Lindhard function with the same limiting branches as Fortran."""
    ax = np.abs(x)
    out = np.empty_like(ax)
    small = ax < 1.0e-8
    near_one = np.abs(ax - 1.0) < 1.0e-7
    regular = ~(small | near_one)
    out[small] = 1.0 - x[small] ** 2 / 3.0 - x[small] ** 4 / 15.0
    out[near_one] = 0.5
    ar = ax[regular]
    out[regular] = 0.5 * (
        1.0 + (1.0 - ar * ar) / (2.0 * ar)
        * np.log(np.abs((1.0 + ar) / (1.0 - ar)))
    )
    return out


def base_dimensionless_kernel(y: np.ndarray, xqmax: float = 80.0, nq: int = 6000) -> np.ndarray:
    """Finite-table base g(y), using the production Simpson quadrature."""
    q = np.linspace(0.0, xqmax, nq + 1)
    weights = np.ones(nq + 1)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    x = q / 2.0
    l1 = (5.0 / 8.0) * (1.0 / _lindhard_l(x) - 3.0 * x * x + 3.0 / 5.0)
    prefactor = (xqmax / nq) / 3.0 / (2.0 * np.pi * np.pi)
    result = np.empty_like(y)
    zero = np.isclose(y, 0.0)
    result[zero] = prefactor * np.sum(weights * q * q * l1)
    indices = np.flatnonzero(~zero)
    for start in range(0, len(indices), 128):
        take = indices[start:start + 128]
        yy = y[take, None]
        result[take] = prefactor * np.sum(
            weights[None, :] * q[None, :] * l1[None, :] * np.sin(yy * q[None, :]),
            axis=1,
        ) / y[take]
    return result


def nlf_distance(reference: np.ndarray, model: np.ndarray, cutoff: float) -> float:
    """Paper definition: r^2 rho_KS-weighted NLF distance below cutoff."""
    r = model[:, 0]
    rho_ks = np.interp(r, reference[:, 0], reference[:, 1])
    nlf_ks = np.interp(r, reference[:, 0], reference[:, 5])
    mask = r < cutoff
    weights = r[mask] ** 2 * rho_ks[mask]
    return float(np.sqrt(np.sum(weights * (model[mask, 6] - nlf_ks[mask]) ** 2) / np.sum(weights)))


def contrast_figure(output: Path) -> None:
    ks = np.loadtxt(KS / "ks_nlf_O16.dat")
    value_density = np.loadtxt(DIAGNOSTICS / "O16_value_only_density.dat")
    value_nlf = np.loadtxt(DIAGNOSTICS / "O16_value_only_nlf.dat")
    final_density = np.loadtxt(RAW / "dens_pernuc_O16.dat")
    final_nlf = np.loadtxt(RAW / "nlf_pernuc_O16.dat")

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.0), constrained_layout=True)
    axes[0].plot(ks[:, 0], 2.0 * ks[:, 1], "k-", lw=2, label="KS reference")
    axes[0].plot(value_density[:, 0], 2.0 * value_density[:, 1], "b--", lw=1.8,
                 label="value-only fit endpoint")
    axes[0].plot(final_density[:, 0], 2.0 * final_density[:, 1], "r-", lw=1.8,
                 label="functional-derivative + energy fit endpoint")
    axes[0].set_xlim(0, 7)
    axes[0].set_xlabel(r"$r$ [fm]")
    axes[0].set_ylabel(r"$\rho_{\rm tot}$ [fm$^{-3}$]")
    axes[0].set_title(r"(a) $^{16}$O density after 8000 ITE steps")
    axes[0].legend(fontsize=8)

    axes[1].plot(ks[:, 0], ks[:, 5], "k-", lw=2, label="KS reference")
    axes[1].plot(value_nlf[:, 0], value_nlf[:, 6], "b--", lw=1.8,
                 label="value-only fit endpoint")
    axes[1].plot(final_nlf[:, 0], final_nlf[:, 6], "r-", lw=1.8,
                 label="functional-derivative + energy fit endpoint")
    axes[1].set_xlim(0, 7)
    axes[1].set_ylim(0, 1.05)
    axes[1].set_xlabel(r"$r$ [fm]")
    axes[1].set_ylabel("species NLF")
    axes[1].set_title("(b) Species NLF at the 8000-step endpoints")
    axes[1].legend(fontsize=8)
    save_png(fig, output)
    plt.close(fig)


def light_nuclei_figure(output: Path) -> None:
    nuclei = ("O16", "Ca40")
    labels = (r"$^{16}$O", r"$^{40}$Ca")
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 6.8), constrained_layout=True)
    for j, (nucleus, label) in enumerate(zip(nuclei, labels)):
        ks = np.loadtxt(KS / f"ks_nlf_{nucleus}.dat")
        density = np.loadtxt(RAW / f"dens_pernuc_{nucleus}.dat")
        final_nlf = np.loadtxt(RAW / f"nlf_pernuc_{nucleus}.dat")
        before = np.loadtxt(DIAGNOSTICS / f"{nucleus}_eg_NT14_{nucleus}_nlf_fixed.dat")
        distance = nlf_distance(ks, final_nlf, 7.0)

        axes[0, j].plot(ks[:, 0], 2.0 * ks[:, 1], "k-", lw=2, label="KS reference")
        axes[0, j].plot(density[:, 0], 2.0 * density[:, 1], "r--", lw=1.8,
                        label="final OFDFT density")
        axes[0, j].set_xlim(0, 10)
        axes[0, j].set_ylabel(r"$\rho_{\rm tot}$ [fm$^{-3}$]")
        axes[0, j].set_title(label)
        axes[0, j].legend(fontsize=8)

        axes[1, j].plot(ks[:, 0], ks[:, 5], "k-", lw=2, label="KS reference")
        axes[1, j].plot(before[:, 0], before[:, 1], "b:", lw=1.6,
                        label="without tail penalties (fixed KS density)")
        axes[1, j].plot(final_nlf[:, 0], final_nlf[:, 6], "r--", lw=1.8,
                        label=rf"final OFDFT NLF ($d_{{\rm NLF}}={distance:.3f}$)")
        axes[1, j].set_xlim(0, 10)
        axes[1, j].set_ylim(0, 1.05)
        axes[1, j].set_xlabel(r"$r$ [fm]")
        axes[1, j].set_ylabel("species NLF")
        axes[1, j].legend(fontsize=8)
    fig.suptitle("Per-nucleus fits with the quadratic tail-margin penalty")
    save_png(fig, output)
    plt.close(fig)


def heavy_nuclei_figure(output: Path) -> None:
    nuclei = ("A80", "A140")
    labels = (r"$A=80$ ($Z=N=40$)", r"$A=140$ ($Z=N=70$)")
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 6.8), constrained_layout=True)
    for j, (nucleus, label) in enumerate(zip(nuclei, labels)):
        ks = np.loadtxt(KS / f"ks_nlf_{nucleus}.dat")
        density = np.loadtxt(RAW / f"dens_pernuc_{nucleus}.dat")
        nlf = np.loadtxt(RAW / f"nlf_pernuc_{nucleus}.dat")
        distance = nlf_distance(ks, nlf, 9.0)
        axes[0, j].plot(ks[:, 0], 2.0 * ks[:, 1], "k-", lw=2, label="KS reference")
        axes[0, j].plot(density[:, 0], 2.0 * density[:, 1], "r--", lw=1.8,
                        label="final OFDFT density")
        axes[0, j].set_xlim(0, 11)
        axes[0, j].set_ylabel(r"$\rho_{\rm tot}$ [fm$^{-3}$]")
        axes[0, j].set_title(label)
        axes[0, j].legend(fontsize=8)
        axes[1, j].plot(ks[:, 0], ks[:, 5], "k-", lw=2, label="KS reference")
        axes[1, j].plot(nlf[:, 0], nlf[:, 6], "r--", lw=1.8,
                        label=rf"final OFDFT NLF ($d_{{\rm NLF}}={distance:.3f}$)")
        axes[1, j].set_xlim(0, 11)
        axes[1, j].set_ylim(0, 1.05)
        axes[1, j].set_xlabel(r"$r$ [fm]")
        axes[1, j].set_ylabel("species NLF")
        axes[1, j].legend(fontsize=8)
    fig.suptitle("Per-nucleus fits for heavier closed-shell systems")
    save_png(fig, output)
    plt.close(fig)


def shared_multi_nucleus_figure(output: Path) -> None:
    nuclei = ("O16", "Ca40", "A80", "A140")
    labels = (r"$^{16}$O", r"$^{40}$Ca", r"$A=80$", r"$A=140$ (not used in fit)")
    radial_cutoffs = (7.0, 7.0, 9.0, 9.0)
    plot_limits = (8.0, 8.0, 11.0, 11.0)
    fig, axes = plt.subplots(2, 4, figsize=(13.0, 5.8), constrained_layout=True)
    for j, (nuc, label, cutoff, xmax) in enumerate(zip(nuclei, labels, radial_cutoffs, plot_limits)):
        ks = np.loadtxt(KS / f"ks_nlf_{nuc}.dat")
        den = np.loadtxt(RAW / f"v5_{nuc}_d.dat")
        nlf = np.loadtxt(RAW / f"v5nlf_{nuc}.dat")
        ax = axes[0, j]
        ax.plot(ks[:, 0], 2 * ks[:, 1], "k-", lw=1.8, label="KS reference")
        ax.plot(den[:, 0], 2 * den[:, 1], "r--", lw=1.6,
                label="shared multi-nucleus fit")
        ax.set_title(label, fontweight="bold" if nuc == "A140" else "normal")
        ax.set_ylabel(r"$\rho_{\rm tot}$ [fm$^{-3}$]" if j == 0 else "")
        ax.set_xlim(0, xmax)
        ax.legend(fontsize=8)
        ax = axes[1, j]
        ax.plot(ks[:, 0], ks[:, 5], "k-", lw=1.8, label="KS reference")
        ax.plot(nlf[:, 0], nlf[:, 6], "r--", lw=1.6,
                label="shared multi-nucleus fit")
        rho_ks = np.interp(nlf[:, 0], ks[:, 0], ks[:, 1])
        nlf_ks = np.interp(nlf[:, 0], ks[:, 0], ks[:, 5])
        mask = nlf[:, 0] < cutoff
        weights = nlf[mask, 0] ** 2 * rho_ks[mask]
        distance = np.sqrt(np.sum(weights * (nlf[mask, 6] - nlf_ks[mask]) ** 2) / np.sum(weights))
        ax.text(0.04, 0.94, rf"$d_{{\rm NLF}}={distance:.3f}$", transform=ax.transAxes,
                ha="left", va="top", fontsize=9,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.5})
        ax.set_ylim(0, 1.03)
        ax.set_xlim(0, xmax)
        ax.set_xlabel(r"$r$ [fm]")
        ax.set_ylabel("species NLF" if j == 0 else "")
        ax.legend(fontsize=8)
    fig.suptitle("Shared multi-nucleus fit; $A=140$ was not used in fitting")
    save_png(fig, output)
    plt.close(fig)


def gp_figure(gp_dir: Path, output: Path) -> None:
    kernel = np.loadtxt(gp_dir / "eg_GP_band.dat")
    a140_prediction = np.loadtxt(gp_dir / "gp_holdout_A140_nlf.dat")
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.0), constrained_layout=True)

    y, mean, sigma = kernel[:, 0], kernel[:, 1], kernel[:, 2]
    axes[0].plot(y, mean, color="tab:blue", lw=2, label=r"GP MAP $\delta g$")
    axes[0].fill_between(y, mean - 2 * sigma, mean + 2 * sigma,
                         color="tab:blue", alpha=0.20,
                         label=r"fixed-active-set $\pm2\sigma$")
    axes[0].axhline(0, color="0.4", lw=0.7)
    axes[0].set_xlim(0, 12)
    axes[0].set_xlabel(r"$y$")
    axes[0].set_ylabel(r"$\delta g(y)$")
    axes[0].set_title("Kernel correction")
    axes[0].legend(fontsize=8)

    r, nlf_gp, nlf_ks, band_2sigma = (
        a140_prediction[:, 0], a140_prediction[:, 1],
        a140_prediction[:, 2], a140_prediction[:, 3]
    )
    axes[1].plot(r, nlf_ks, "k-", lw=2, label="KS reference")
    axes[1].plot(r, nlf_gp, color="tab:blue", lw=2,
                 label="GP prediction at fixed KS density")
    axes[1].fill_between(r, np.clip(nlf_gp - band_2sigma, 0, 1),
                         np.clip(nlf_gp + band_2sigma, 0, 1),
                         color="tab:blue", alpha=0.20,
                         label=r"fixed-active-set $\pm2\sigma$")
    axes[1].set_xlim(0, 10)
    axes[1].set_ylim(0, 1.03)
    axes[1].set_xlabel(r"$r$ [fm]")
    axes[1].set_ylabel("species NLF")
    axes[1].set_title(r"Fixed-density prediction for $A=140$")
    axes[1].legend(fontsize=8)
    save_png(fig, output)
    plt.close(fig)


def kernel_shapes_figure(gp_dir: Path, output: Path) -> None:
    per_nucleus = {
        "O16": KERNELS / "per_nucleus" / "eg_PN14_O16.dat",
        "Ca40": KERNELS / "per_nucleus" / "eg_PN14_Ca40.dat",
        "A80": KERNELS / "per_nucleus" / "eg_PN14_A80.dat",
        "A140": KERNELS / "per_nucleus" / "eg_PN14_A140.dat",
    }
    shared = np.loadtxt(KERNELS / "universal" / "eg_SH14.dat")
    gp = np.loadtxt(gp_dir / "eg_GP_band.dat")
    mask = shared[:, 0] <= 12.0
    y = shared[mask, 0]
    base = base_dimensionless_kernel(y)

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.0), constrained_layout=True)
    axes[0].plot(y, base, "k-", lw=2.0, label="finite-table Lindhard-derived base")
    axes[0].axhline(0.0, color="0.5", lw=0.7)
    axes[0].set_xlim(0, 12)
    axes[0].set_xlabel(r"$y=k_F R$")
    axes[0].set_ylabel(r"$g(y)$")
    axes[0].set_title("(a) Base dimensionless kernel")
    axes[0].legend(fontsize=8)

    for nucleus, path in per_nucleus.items():
        values = np.loadtxt(path)
        axes[1].plot(values[:, 0], values[:, 1], lw=1.2, label=f"per-nucleus {nucleus}")
    axes[1].plot(shared[:, 0], shared[:, 1], "k--", lw=1.8,
                 label="shared multi-nucleus fit")
    axes[1].plot(gp[:, 0], 10.0 * gp[:, 1], color="magenta", lw=1.8,
                 label=r"GP fixed-active-set MAP ($\times10$)")
    axes[1].fill_between(gp[:, 0], 10.0 * (gp[:, 1] - 2.0 * gp[:, 2]),
                         10.0 * (gp[:, 1] + 2.0 * gp[:, 2]),
                         color="magenta", alpha=0.16,
                         label=r"GP fixed-active-set $2\sigma$ ($\times10$)")
    axes[1].axhline(0.0, color="0.5", lw=0.7)
    axes[1].set_xlim(0, 12)
    axes[1].set_xlabel(r"$y$")
    axes[1].set_ylabel(r"$\delta g(y)$")
    axes[1].set_title("(b) Fitted kernel corrections")
    axes[1].legend(fontsize=7.5, ncol=2)
    save_png(fig, output)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gp-dir", type=Path, default=ROOT / "results" / "gp")
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "results" / "figures")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    light_nuclei_figure(args.output_dir / "fig1.png")
    contrast_figure(args.output_dir / "fig2.png")
    heavy_nuclei_figure(args.output_dir / "fig3.png")
    shared_multi_nucleus_figure(args.output_dir / "fig4.png")
    kernel_shapes_figure(args.gp_dir, args.output_dir / "fig5.png")
    gp_figure(args.gp_dir, args.output_dir / "fig6.png")
    print(f"Wrote figures to {args.output_dir}")


if __name__ == "__main__":
    main()
