#!/usr/bin/env python3
"""Compare one fresh solver run with the released production result.

The solver-reproduction driver writes the following files in each case
directory::

    ite_density.dat  final one-species density from imaginary-time evolution
    energy.dat       nmu=64 energy/functional-derivative summary
    nlf.dat          nmu=48 NLF table

This script compares those files with the corresponding immutable data under
``results/raw``.  It deliberately uses numerical tolerances rather than file
hash equality: floating-point roundoff can change across supported compilers
and BLAS/LAPACK implementations even when the calculation is reproduced.

Exit status is 0 for a numerical pass, 1 for a completed comparison with one
or more tolerance failures, and 2 for invalid input or an unreadable result.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


NUCLEI = {"o16": "O16", "ca40": "Ca40", "a80": "A80", "a140": "A140"}
RCUT_FM = {"O16": 7.0, "Ca40": 7.0, "A80": 9.0, "A140": 9.0}

DEFAULT_TOLERANCES = {
    # One keV on a total nuclear energy.  This is tight enough to catch a
    # changed trajectory while allowing harmless compiler-level roundoff.
    "energy_absolute_difference_mev": 1.0e-3,
    "density_relative_radial_l2": 1.0e-5,
    "central_total_density_difference_fm3": 1.0e-5,
    "rms_radius_difference_fm": 1.0e-4,
    "nlf_distance_difference": 1.0e-5,
    "nlf_curve_weighted_rms_difference": 1.0e-4,
}


class ComparisonInputError(RuntimeError):
    """Raised when a requested comparison cannot be evaluated safely."""


def _parse_case(value: str) -> tuple[str, str]:
    """Normalize a public case name to (family, nucleus)."""
    normalized = value.strip().replace(":", "/")
    parts = normalized.split("/")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            "case must be FAMILY/NUCLEUS, for example per-nucleus/O16"
        )

    family_aliases = {
        "per-nucleus": "per-nucleus",
        "per_nucleus": "per-nucleus",
        "pernucleus": "per-nucleus",
        "per": "per-nucleus",
        "shared": "shared",
        "universal": "shared",
    }
    family = family_aliases.get(parts[0].strip().lower())
    nucleus = NUCLEI.get(parts[1].strip().lower())
    if family is None:
        raise argparse.ArgumentTypeError(
            "family must be 'per-nucleus' (or 'per_nucleus') or 'shared'"
        )
    if nucleus is None:
        raise argparse.ArgumentTypeError(
            "nucleus must be one of O16, Ca40, A80, or A140"
        )
    return family, nucleus


def _load_table(path: Path, minimum_columns: int, description: str) -> np.ndarray:
    """Load and validate a whitespace-delimited numerical table."""
    if not path.is_file():
        raise ComparisonInputError(f"missing {description}: {path}")
    try:
        table = np.loadtxt(path, comments="#", ndmin=2)
    except (OSError, ValueError) as exc:
        raise ComparisonInputError(f"cannot read {description} {path}: {exc}") from exc
    if table.shape[0] < 2 or table.shape[1] < minimum_columns:
        raise ComparisonInputError(
            f"{description} {path} has shape {table.shape}; expected at least "
            f"2 rows and {minimum_columns} columns"
        )
    if not np.all(np.isfinite(table)):
        raise ComparisonInputError(f"{description} contains NaN or infinity: {path}")
    radius = table[:, 0]
    if np.any(np.diff(radius) <= 0.0):
        raise ComparisonInputError(f"radius grid is not strictly increasing in {path}")
    return table


def _read_summary_value(path: Path, key: str) -> float:
    """Read one named scalar from a Fortran summary file."""
    if not path.is_file():
        raise ComparisonInputError(f"missing energy summary: {path}")
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ComparisonInputError(f"cannot read energy summary {path}: {exc}") from exc
    for line in lines:
        fields = line.split()
        if len(fields) >= 2 and fields[0] == key:
            try:
                value = float(fields[1].replace("D", "E").replace("d", "e"))
            except ValueError as exc:
                raise ComparisonInputError(
                    f"invalid {key} value in {path}: {fields[1]!r}"
                ) from exc
            if not math.isfinite(value):
                raise ComparisonInputError(f"non-finite {key} value in {path}")
            return value
    raise ComparisonInputError(f"key {key!r} was not found in {path}")


def _require_grid_coverage(source_r: np.ndarray, target_r: np.ndarray, label: str) -> None:
    """Reject interpolation that would silently extrapolate at either end."""
    scale = max(1.0, float(abs(target_r[0])), float(abs(target_r[-1])))
    allowance = 1.0e-10 * scale
    if source_r[0] > target_r[0] + allowance or source_r[-1] < target_r[-1] - allowance:
        raise ComparisonInputError(
            f"{label} grid [{source_r[0]:.8g}, {source_r[-1]:.8g}] fm does not "
            f"cover target grid [{target_r[0]:.8g}, {target_r[-1]:.8g}] fm"
        )


def _relative_radial_l2(reference: np.ndarray, trial: np.ndarray) -> float:
    """Return sqrt(sum r^2(drho)^2 / sum r^2 rho_ref^2)."""
    r_ref, rho_ref = reference[:, 0], reference[:, 1]
    _require_grid_coverage(trial[:, 0], r_ref, "trial density")
    rho_trial = np.interp(r_ref, trial[:, 0], trial[:, 1])
    denominator = float(np.sum(r_ref**2 * rho_ref**2))
    if denominator <= 0.0:
        raise ComparisonInputError("reference density has zero radial L2 norm")
    numerator = float(np.sum(r_ref**2 * (rho_trial - rho_ref) ** 2))
    return math.sqrt(max(0.0, numerator / denominator))


def _rms_radius(table: np.ndarray) -> float:
    """Compute the rms radius used by the public table-generation script."""
    r, rho = table[:, 0], table[:, 1]
    denominator = float(np.sum(r**2 * rho))
    numerator = float(np.sum(r**4 * rho))
    if denominator <= 0.0 or numerator < 0.0:
        raise ComparisonInputError("density cannot define a finite rms radius")
    return math.sqrt(numerator / denominator)


def _nlf_distance(model: np.ndarray, ks: np.ndarray, cutoff_fm: float) -> float:
    """Paper d_NLF: KS-density-weighted NLF distance below the radial cutoff."""
    r = model[:, 0]
    _require_grid_coverage(ks[:, 0], r, "KS reference")
    mask = r < cutoff_fm
    if not np.any(mask):
        raise ComparisonInputError(f"NLF grid has no points below {cutoff_fm:g} fm")
    rho_ks = np.interp(r, ks[:, 0], ks[:, 1])
    nlf_ks = np.interp(r, ks[:, 0], ks[:, 5])
    weights = r[mask] ** 2 * rho_ks[mask]
    denominator = float(np.sum(weights))
    if denominator <= 0.0:
        raise ComparisonInputError("KS density gives zero NLF weight")
    squared = float(np.sum(weights * (model[mask, 6] - nlf_ks[mask]) ** 2))
    return math.sqrt(max(0.0, squared / denominator))


def _nlf_curve_difference(
    reference: np.ndarray, trial: np.ndarray, ks: np.ndarray, cutoff_fm: float
) -> tuple[float, float]:
    """Compare two NLF curves using the paper's KS-density radial weight."""
    r = reference[:, 0]
    _require_grid_coverage(trial[:, 0], r, "trial NLF")
    _require_grid_coverage(ks[:, 0], r, "KS reference")
    mask = r < cutoff_fm
    model_trial = np.interp(r, trial[:, 0], trial[:, 6])
    rho_ks = np.interp(r, ks[:, 0], ks[:, 1])
    weights = r[mask] ** 2 * rho_ks[mask]
    difference = model_trial[mask] - reference[mask, 6]
    denominator = float(np.sum(weights))
    if denominator <= 0.0:
        raise ComparisonInputError("KS density gives zero NLF curve weight")
    weighted_rms = math.sqrt(
        max(0.0, float(np.sum(weights * difference**2)) / denominator)
    )
    maximum_absolute = float(np.max(np.abs(difference)))
    return weighted_rms, maximum_absolute


def _canonical_paths(root: Path, family: str, nucleus: str) -> dict[str, Path]:
    raw = root / "results" / "raw"
    if family == "per-nucleus":
        density = raw / f"dens_pernuc_{nucleus}.dat"
        energy = raw / f"pnE_{nucleus}.dat"
        nlf = raw / f"nlf_pernuc_{nucleus}.dat"
    else:
        density = raw / f"v5_{nucleus}_d.dat"
        energy = raw / f"v5e64_{nucleus}.dat"
        nlf = raw / f"v5nlf_{nucleus}.dat"
    return {
        "density": density,
        "energy": energy,
        "nlf": nlf,
        "ks_nlf": root / "results" / "reference" / f"ks_nlf_{nucleus}.dat",
    }


def _display_path(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def _check(name: str, value: float, tolerance: float) -> dict[str, Any]:
    return {
        "name": name,
        "value": value,
        "tolerance": tolerance,
        "passed": bool(value <= tolerance),
    }


def compare(
    root: Path,
    run_directory: Path,
    family: str,
    nucleus: str,
    tolerances: dict[str, float],
) -> dict[str, Any]:
    """Evaluate all numerical metrics and return a JSON-serializable report."""
    canonical = _canonical_paths(root, family, nucleus)
    run_paths = {
        "density": run_directory / "ite_density.dat",
        "energy": run_directory / "energy.dat",
        "nlf": run_directory / "nlf.dat",
    }

    density_ref = _load_table(canonical["density"], 2, "reference density")
    density_run = _load_table(run_paths["density"], 2, "trial density")
    nlf_ref = _load_table(canonical["nlf"], 7, "reference NLF")
    nlf_run = _load_table(run_paths["nlf"], 7, "trial NLF")
    ks = _load_table(canonical["ks_nlf"], 6, "KS NLF reference")

    energy_ref = _read_summary_value(canonical["energy"], "E_total_NeqZ_MeV")
    energy_run = _read_summary_value(run_paths["energy"], "E_total_NeqZ_MeV")
    energy_difference = abs(energy_run - energy_ref)

    density_l2 = _relative_radial_l2(density_ref, density_run)
    # rho_q is a one-species density.  The systems are N=Z, hence rho_tot=2 rho_q.
    central_radius = float(density_ref[0, 0])
    central_ref = 2.0 * float(density_ref[0, 1])
    central_run = 2.0 * float(
        np.interp(central_radius, density_run[:, 0], density_run[:, 1])
    )
    central_difference = abs(central_run - central_ref)
    rms_ref = _rms_radius(density_ref)
    rms_run = _rms_radius(density_run)
    rms_difference = abs(rms_run - rms_ref)

    cutoff = RCUT_FM[nucleus]
    dnlf_ref = _nlf_distance(nlf_ref, ks, cutoff)
    dnlf_run = _nlf_distance(nlf_run, ks, cutoff)
    dnlf_difference = abs(dnlf_run - dnlf_ref)
    nlf_curve_rms, nlf_curve_max = _nlf_curve_difference(nlf_ref, nlf_run, ks, cutoff)

    checks = [
        _check(
            "energy_absolute_difference_mev",
            energy_difference,
            tolerances["energy_absolute_difference_mev"],
        ),
        _check(
            "density_relative_radial_l2",
            density_l2,
            tolerances["density_relative_radial_l2"],
        ),
        _check(
            "central_total_density_difference_fm3",
            central_difference,
            tolerances["central_total_density_difference_fm3"],
        ),
        _check(
            "rms_radius_difference_fm",
            rms_difference,
            tolerances["rms_radius_difference_fm"],
        ),
        _check(
            "nlf_distance_difference",
            dnlf_difference,
            tolerances["nlf_distance_difference"],
        ),
        _check(
            "nlf_curve_weighted_rms_difference",
            nlf_curve_rms,
            tolerances["nlf_curve_weighted_rms_difference"],
        ),
    ]
    passed = all(item["passed"] for item in checks)

    return {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "case": {"family": family, "nucleus": nucleus},
        "status": "PASS" if passed else "FAIL",
        "paths": {
            "run_directory": _display_path(run_directory, root),
            "run": {key: _display_path(path, root) for key, path in run_paths.items()},
            "reference": {
                key: _display_path(path, root) for key, path in canonical.items()
            },
        },
        "definitions": {
            "density": "one-species density rho_q; central total density is 2*rho_q(r_min)",
            "density_relative_radial_l2": (
                "sqrt[sum_i r_i^2 (rho_run-rho_ref)^2 / "
                "sum_i r_i^2 rho_ref^2]"
            ),
            "rms_radius": "sqrt[sum_i r_i^4 rho_q / sum_i r_i^2 rho_q]",
            "nlf_cutoff_fm": cutoff,
            "nlf_weight": "r^2 rho_KS below the nucleus-dependent cutoff",
        },
        "tolerances": tolerances,
        "metrics": {
            "energy": {
                "run_mev": energy_run,
                "reference_mev": energy_ref,
                "absolute_difference_mev": energy_difference,
            },
            "density": {
                "relative_radial_l2": density_l2,
                "central_grid_radius_fm": central_radius,
                "central_total_run_fm3": central_run,
                "central_total_reference_fm3": central_ref,
                "central_total_absolute_difference_fm3": central_difference,
                "rms_radius_run_fm": rms_run,
                "rms_radius_reference_fm": rms_ref,
                "rms_radius_absolute_difference_fm": rms_difference,
            },
            "nlf": {
                "distance_to_ks_run": dnlf_run,
                "distance_to_ks_reference": dnlf_ref,
                "distance_to_ks_absolute_difference": dnlf_difference,
                "curve_weighted_rms_difference": nlf_curve_rms,
                "curve_maximum_absolute_difference_below_cutoff": nlf_curve_max,
            },
        },
        "checks": checks,
    }


def _print_human_report(report: dict[str, Any], json_path: Path) -> None:
    case = report["case"]
    metrics = report["metrics"]
    checks = {item["name"]: item for item in report["checks"]}

    print(f"Reproduction comparison: {case['family']}/{case['nucleus']}")
    print(
        "  energy       run={:.9f}  reference={:.9f}  |delta|={:.3e} MeV".format(
            metrics["energy"]["run_mev"],
            metrics["energy"]["reference_mev"],
            metrics["energy"]["absolute_difference_mev"],
        )
    )
    print(
        "  density      relative radial L2={:.3e}; |delta rho_tot,c|={:.3e} fm^-3; "
        "|delta r_rms|={:.3e} fm".format(
            metrics["density"]["relative_radial_l2"],
            metrics["density"]["central_total_absolute_difference_fm3"],
            metrics["density"]["rms_radius_absolute_difference_fm"],
        )
    )
    print(
        "  NLF          d_NLF(run)={:.7f}; d_NLF(ref)={:.7f}; |delta|={:.3e}; "
        "curve weighted RMS={:.3e}".format(
            metrics["nlf"]["distance_to_ks_run"],
            metrics["nlf"]["distance_to_ks_reference"],
            metrics["nlf"]["distance_to_ks_absolute_difference"],
            metrics["nlf"]["curve_weighted_rms_difference"],
        )
    )
    for item in report["checks"]:
        marker = "PASS" if item["passed"] else "FAIL"
        print(
            f"  [{marker}] {item['name']}: {item['value']:.3e} "
            f"<= {item['tolerance']:.3e}"
        )
    print(f"Overall: {report['status']}")
    print(f"JSON: {json_path}")


def _positive_tolerance(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("tolerance must be a number") from exc
    if not math.isfinite(parsed) or parsed < 0.0:
        raise argparse.ArgumentTypeError("tolerance must be finite and nonnegative")
    return parsed


def _build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Compare one fresh ks.x reproduction with released raw results."
    )
    parser.add_argument(
        "--case",
        required=True,
        type=_parse_case,
        metavar="FAMILY/NUCLEUS",
        help="per-nucleus/O16 or shared/O16 (four supported nuclei)",
    )
    parser.add_argument(
        "--run",
        required=True,
        type=Path,
        metavar="DIR",
        help="case directory containing ite_density.dat, energy.dat, and nlf.dat",
    )
    parser.add_argument(
        "--repository-root",
        type=Path,
        default=root,
        help="repository root containing results/raw (default: inferred from script)",
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=None,
        metavar="FILE",
        help="JSON report path (default: RUN/comparison.json)",
    )
    parser.add_argument(
        "--energy-tol-mev",
        type=_positive_tolerance,
        default=DEFAULT_TOLERANCES["energy_absolute_difference_mev"],
    )
    parser.add_argument(
        "--density-rel-l2-tol",
        type=_positive_tolerance,
        default=DEFAULT_TOLERANCES["density_relative_radial_l2"],
    )
    parser.add_argument(
        "--central-density-tol-fm3",
        type=_positive_tolerance,
        default=DEFAULT_TOLERANCES["central_total_density_difference_fm3"],
    )
    parser.add_argument(
        "--rms-radius-tol-fm",
        type=_positive_tolerance,
        default=DEFAULT_TOLERANCES["rms_radius_difference_fm"],
    )
    parser.add_argument(
        "--nlf-distance-tol",
        type=_positive_tolerance,
        default=DEFAULT_TOLERANCES["nlf_distance_difference"],
    )
    parser.add_argument(
        "--nlf-curve-rms-tol",
        type=_positive_tolerance,
        default=DEFAULT_TOLERANCES["nlf_curve_weighted_rms_difference"],
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    family, nucleus = args.case
    root = args.repository_root.resolve()
    run_directory = args.run.resolve()
    json_path = args.json.resolve() if args.json is not None else run_directory / "comparison.json"
    tolerances = {
        "energy_absolute_difference_mev": args.energy_tol_mev,
        "density_relative_radial_l2": args.density_rel_l2_tol,
        "central_total_density_difference_fm3": args.central_density_tol_fm3,
        "rms_radius_difference_fm": args.rms_radius_tol_fm,
        "nlf_distance_difference": args.nlf_distance_tol,
        "nlf_curve_weighted_rms_difference": args.nlf_curve_rms_tol,
    }

    try:
        report = compare(root, run_directory, family, nucleus, tolerances)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    except (ComparisonInputError, OSError) as exc:
        print(f"comparison error: {exc}", file=sys.stderr)
        return 2

    _print_human_report(report, json_path)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
