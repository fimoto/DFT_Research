#!/usr/bin/env python3
"""Validate the eight-case index and all 24 production namelists."""
from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STAGES = ("ite", "energy", "nlf")


def parse_value(raw: str) -> object:
    value = raw.strip()
    if len(value) >= 2 and value[0] == value[-1] == "'":
        return value[1:-1]
    lower = value.lower()
    if lower == ".true.":
        return True
    if lower == ".false.":
        return False
    number = value.replace("D", "e").replace("d", "e")
    try:
        if re.fullmatch(r"[+-]?\d+", number):
            return int(number)
        return float(number)
    except ValueError:
        return value


def parse_namelist(path: Path) -> tuple[dict[str, object], list[str]]:
    lines = []
    for line in path.read_text(encoding="ascii").splitlines():
        lines.append(line.split("!", 1)[0])
    text = "\n".join(lines)
    match = re.search(r"&input\b(.*?)(?:^|\n)\s*/\s*(?:\n|$)", text, re.I | re.S)
    if not match:
        raise ValueError("missing or unterminated &input group")
    assignments = re.findall(
        r"\b([A-Za-z][A-Za-z0-9_]*)\s*=\s*(" 
        r"'(?:[^']|'')*'|\.true\.|\.false\.|[^,\n/]+)",
        match.group(1),
        re.I,
    )
    names = [name.lower() for name, _ in assignments]
    duplicates = [name for name, count in Counter(names).items() if count > 1]
    values = {name.lower(): parse_value(value) for name, value in assignments}
    return values, sorted(duplicates)


def declared_input_names(source: Path) -> set[str]:
    text = source.read_text(encoding="utf-8")
    match = re.search(r"namelist\s*/input/(.*?)(?=\n\s*\n)", text, re.I | re.S)
    if not match:
        raise ValueError(f"cannot locate /input/ declaration in {source}")
    block = re.sub(r"!.*", "", match.group(1)).replace("&", " ")
    return {item.strip().lower() for item in block.split(",") if item.strip()}


def equal(actual: object, expected: object) -> bool:
    if isinstance(expected, float):
        return isinstance(actual, (int, float)) and math.isclose(
            float(actual), expected, rel_tol=0.0, abs_tol=1.0e-12
        )
    return actual == expected


def check_values(
    relative: str,
    values: dict[str, object],
    expected: dict[str, object],
    failures: list[str],
) -> None:
    for key, wanted in expected.items():
        if key not in values:
            failures.append(f"{relative}: missing explicit key {key}")
        elif not equal(values[key], wanted):
            failures.append(
                f"{relative}: {key}={values[key]!r}, expected {wanted!r}"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--index", type=Path, default=ROOT / "production_inputs" / "cases.tsv"
    )
    args = parser.parse_args()
    index_path = args.index.resolve()
    release_root = index_path.parents[1]
    failures: list[str] = []
    try:
        known_keys = declared_input_names(
            release_root / "src" / "nofdft_ks_spherical_skp_nlf_v8.f90"
        )
    except (OSError, ValueError) as exc:
        print(f"PRODUCTION INPUT VALIDATION FAILED\n  {exc}")
        return 1

    with index_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 8:
        failures.append(f"{index_path}: expected 8 case rows, found {len(rows)}")

    case_ids = [row["case_id"] for row in rows]
    if len(case_ids) != len(set(case_ids)):
        failures.append("cases.tsv: duplicate case_id")

    indexed_namelists: set[str] = set()
    for row in rows:
        case_id = row["case_id"]
        family = row["family_dir"]
        nucleus = row["nucleus"]
        expected_case = (
            f"per-nucleus/{nucleus}" if family == "per_nucleus" else f"shared/{nucleus}"
        )
        if case_id != expected_case:
            failures.append(f"cases.tsv: {case_id}: inconsistent family_dir/nucleus")

        for source_key in ("initial_density", "kernel"):
            source = release_root / row[source_key]
            if not source.is_file():
                failures.append(f"{case_id}: missing {source_key}: {row[source_key]}")

        common: dict[str, object] = {
            "run_ks": False,
            "run_teacher": False,
            "run_deriv_check": False,
            "run_minimize": False,
            "run_diagscf": False,
            "nlf_use_ks_density": False,
            "density_scale": 1.0,
            "a": int(row["A"]),
            "npart_q": int(row["npart_q"]),
            "kernel_mode": "dd",
            "kf_const": 0.8,
            "extra_kernel_file": "",
            "extra_g_file": "kernel.dat",
            "c_corr_shift": float(row["c_corr_shift"]),
            "s2_grad_in": float(row["s2_grad_in"]),
            "qmax": 60.0,
            "nq": 6000,
            "nr_tab": 4001,
            "dd_xqmax_in": 80.0,
            "dd_ny_base_in": 8001,
            "dd_kf_scale_in": 1.0,
            "kernel_lambda": 1.0,
            "use_semilocal": False,
        }

        for stage in STAGES:
            field = f"{stage}_nml"
            relative = row[field]
            expected_relative = f"production_inputs/{family}/{nucleus}/{stage}.nml"
            if relative != expected_relative:
                failures.append(
                    f"{case_id}: {field}={relative}, expected {expected_relative}"
                )
            indexed_namelists.add(relative)
            path = release_root / relative
            if not path.is_file():
                failures.append(f"{case_id}: missing {relative}")
                continue
            first_line = path.read_text(encoding="ascii").splitlines()[0]
            if not first_line.startswith(f"! Exact production input: {case_id}, stage "):
                failures.append(f"{relative}: missing standardized case/stage header")
            try:
                values, duplicates = parse_namelist(path)
            except (UnicodeDecodeError, ValueError) as exc:
                failures.append(f"{relative}: {exc}")
                continue
            if duplicates:
                failures.append(f"{relative}: duplicate assignments: {', '.join(duplicates)}")
            unknown = sorted(set(values) - known_keys)
            if unknown:
                failures.append(
                    f"{relative}: keys absent from Fortran /input/: {', '.join(unknown)}"
                )

            wanted = dict(common)
            wanted["nmu"] = int(row[f"{stage}_nmu"])
            if stage == "ite":
                wanted.update(
                    run_nlf=False,
                    run_derivative=False,
                    run_ite=True,
                    external_density_file="initial_density.dat",
                    ite_method="ITE_2A0",
                    ite_init_mode="file",
                    ite_initial_density_file="initial_density.dat",
                    ite_max_iter=8000,
                    ite_dt=1.0e-5,
                    ite_tol=0.0,
                    ite_monitor_every=20 if family == "per_nucleus" else 50,
                    ite_dt_adapt=True,
                    ite_dt_max=5.0e-5,
                    ite_dt_grow=1.1,
                    ite_dt_shrink=0.5,
                    ite_history_file="ite_history.dat",
                    ite_density_file="ite_density.dat",
                    ite_derivative_file="ite_derivative_radial.dat",
                    ite_summary_file="ite_summary.dat",
                    ite_kernel_file="ite_kernel_radial.dat",
                )
            elif stage == "energy":
                wanted.update(
                    run_nlf=False,
                    run_derivative=True,
                    run_ite=False,
                    external_density_file="ite_density.dat",
                    derivative_energy_file="energy.dat",
                    derivative_file="derivative_radial.dat",
                    derivative_kernel_file="derivative_kernel.dat",
                )
            else:
                wanted.update(
                    run_nlf=True,
                    run_derivative=False,
                    run_ite=False,
                    external_density_file="ite_density.dat",
                    nlf_file="nlf.dat",
                )
            check_values(relative, values, wanted, failures)

    discovered = {
        path.relative_to(release_root).as_posix()
        for path in (release_root / "production_inputs").glob("*/*/*.nml")
    }
    for missing in sorted(indexed_namelists - discovered):
        failures.append(f"cases.tsv indexes nonexistent namelist: {missing}")
    for extra in sorted(discovered - indexed_namelists):
        failures.append(f"unindexed production namelist: {extra}")

    if failures:
        print("PRODUCTION INPUT VALIDATION FAILED")
        for failure in failures:
            print(f"  {failure}")
        return 1
    print(f"PRODUCTION INPUTS VERIFIED: {len(rows)} cases, {len(discovered)} namelists")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
