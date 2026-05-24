#!/usr/bin/env python3
"""Run a small BRM-direct validation table and save CSV/JSON outputs."""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path

from benites_direct import estimate


def parse_points(text: str) -> list[tuple[float, float]]:
    out = []
    for item in text.split(","):
        rs, zeta = item.split(":")
        out.append((float(rs), float(zeta)))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--points", default="1:0,2:0,5:0,1:1,2:1,5:1",
                    help="comma-separated rs:zeta pairs")
    ap.add_argument("--m", type=int, default=17)
    ap.add_argument("--replicates", type=int, default=12)
    ap.add_argument("--q-scale", type=float, default=2.0)
    ap.add_argument("--seed", type=int, default=20261001)
    ap.add_argument("--out-prefix", type=Path,
                    default=Path("data/benites_direct_table"))
    args = ap.parse_args()

    results = []
    for i, (rs, zeta) in enumerate(parse_points(args.points)):
        est = estimate(rs, zeta, args.m, args.replicates,
                       args.seed + 1000 * i, args.q_scale)
        results.append(est)
        print(f"rs={rs:g} zeta={zeta:g}: "
              f"{est.full_kite_mry:.4f} +- {est.full_kite_err_mry:.4f} mRy")

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    json_path = args.out_prefix.with_suffix(".json")
    csv_path = args.out_prefix.with_suffix(".csv")
    json_path.write_text(json.dumps([asdict(x) for x in results], indent=2) + "\n")
    fields = [
        "rs", "zeta", "m", "replicates", "q_scale",
        "delta1_ry", "delta1_err_ry", "delta2_ry", "delta2_err_ry",
        "correction_ry", "correction_err_ry",
        "full_kite_mry", "full_kite_err_mry",
    ]
    with csv_path.open("w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=fields)
        wr.writeheader()
        for x in results:
            d = asdict(x)
            wr.writerow({k: d[k] for k in fields})
    print(f"wrote {csv_path} and {json_path}")


if __name__ == "__main__":
    main()
