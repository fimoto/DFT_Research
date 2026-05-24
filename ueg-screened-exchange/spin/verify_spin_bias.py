"""Independent checks of the spin-map interpolation-residual reanalysis.

This script does not import the figure generator. It checks the stored data,
covariance propagation, reference provenance, and exported figure numbers;
it does not rerun the Monte Carlo integral or establish confidence coverage.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'production'))
from common import check_source_compatibility


def run_checks():
    config = json.loads((ROOT / "production/config.json").read_text())
    lock = json.loads((ROOT / "production/reference_lock.json").read_text())["sha256"]
    raw = {}
    for rs in config["map"]["rs"]:
        for zeta in config["map"]["zeta"]:
            name = f"map/rs_{rs:g}_z_{zeta:g}.json"
            path = ROOT / "production/raw" / name
            assert hashlib.sha256(path.read_bytes()).hexdigest() == lock[name]
            row = json.loads(path.read_text())
            assert (row["rs"], row["zeta"]) == (rs, zeta)
            check_source_compatibility(row["metadata"]["sources"])
            raw[rs, zeta] = row
    assert len(raw) == 104
    sampling = [row["settings"] for row in raw.values()]
    assert all(setting == sampling[0] for setting in sampling)
    assert sampling[0]["replicates"] == 8 and sampling[0]["m"] == 17
    assert sampling[0]["seed"] == 20260713 and sampling[0]["q_scale"] == 2
    print("PROVENANCE PASS: all 104 raw records match reference hashes, exact or explicitly audited compatible sources, and common sampling settings.")

    # The saved replicate index is the common scramble index at every node.
    # Independent of the generator's difference formula, propagate each full
    # 3x3 cross-spin covariance matrix with weights [f-1, 1, -f].
    data = np.genfromtxt(ROOT / "figs/data/fig_spin_bias.csv", delimiter=",", names=True)
    assert len(data) == len(raw)
    seen = set()
    for point in data:
        rs, zeta = point["rs"], point["zeta"]
        assert (rs, zeta) not in seen
        seen.add((rs, zeta))
        fx = (((1 + zeta)**(4 / 3) + (1 - zeta)**(4 / 3)) / 2 - 1) / (2**(4 / 3) / 2 - 1)
        assert abs(point["fx"] - fx) < 1e-15
        weights = np.asarray([fx - 1, 1, -fx])
        nodes = [raw[rs, 0], raw[rs, zeta], raw[rs, 1]]
        direct_mean = np.dot(weights, [node["mean_mry"] for node in nodes])
        assert abs(direct_mean - point["residual_mry"]) < 6e-14
        term_mean = []
        variance_of_mean = 0.0
        for term in range(2):
            values = 1000 * np.column_stack([
                np.asarray(node["term_replicates_ry"])[:, term] for node in nodes])
            covariance = np.cov(values, rowvar=False, ddof=1)
            variance_of_mean += weights @ covariance @ weights / len(values)
            residuals = values @ weights
            term_mean.append(residuals.mean())
            assert abs(residuals.var(ddof=1) - weights @ covariance @ weights) < 3e-15
            for constant in (0.048179158, 1.0):
                assert np.max(np.abs((values + 1000 * constant) @ weights - residuals)) < 4e-13
        assert abs(point["se_mry"]**2 - variance_of_mean) < 1e-15
        assert abs(point["term1_residual_mry"] - term_mean[0]) < 3e-14
        assert abs(point["term2_residual_mry"] - term_mean[1]) < 3e-14
        if zeta in (0, 1):
            assert point["residual_mry"] == 0.0 and point["se_mry"] == 0.0
    assert seen == set(raw)
    print("NUM PASS: every exported residual agrees with direct subtraction of the three raw energy means.")
    print("NUM PASS: all 104 error bars agree with the independent 3x3 covariance-matrix propagation for each term.")
    print("NUM PASS: both endpoints vanish exactly and adding a common bare-exchange constant leaves every residual unchanged.")

    # Verify the cancellation algebra without choosing coefficient values.
    f, c, t, a1, a2, a3z, a30, a31, m4, m5, d4, d5 = s.symbols(
        "f c t a1 a2 a3z a30 a31 m4 m5 d4 d5")
    e0 = c + a1*t + a2*t**2 + a30*t**3
    ez = c + a1*m4*t + a2*m5*t**2 + a3z*t**3
    e1 = c + a1*(1+d4)*t + a2*(1+d5)*t**2 + a31*t**3
    residual = s.expand(ez - (1-f)*e0 - f*e1)
    assert s.diff(residual, c) == 0
    assert s.simplify(residual.coeff(t, 1).subs(m4, 1+f*d4)) == 0
    assert s.simplify(residual.coeff(t, 2) - a2*(m5-1-f*d5)) == 0
    assert s.simplify(residual.coeff(t, 3) - (a3z-a30-f*(a31-a30))) == 0
    print("EXACT PASS: common constants and the exchange-shaped leading term cancel; residual second and third coefficients assemble correctly.")

    macros = (ROOT / "spin/generated/spin_bias_numbers.tex").read_text()
    for rs, name in [(1, "One"), (10, "Ten"), (60, "Sixty")]:
        point, = data[(data["rs"] == rs) & (data["zeta"] == 0.8)]
        assert rf"\newcommand{{\BiasRs{name}}}{{{point['residual_mry']:+.3f}}}" in macros
        assert rf"\newcommand{{\BiasRs{name}SE}}{{{point['se_mry']:.3f}}}" in macros
        print(f"NUM PASS: zeta=0.8, rs={rs:g}: {point['residual_mry']:+.12f} +/- {point['se_mry']:.12f} mRy.")
    print("SCOPE: uncertainties are empirical replicate standard errors; neither confidence coverage nor sampler variance finiteness is established here.")
    print("SCOPE: stored-data reanalysis only; no Monte Carlo/HPC rerun, no fitted curve, and no claim about total correlation energy.")


if __name__ == "__main__":
    run_checks()
