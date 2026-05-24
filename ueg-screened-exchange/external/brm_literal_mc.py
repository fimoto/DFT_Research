#!/usr/bin/env python3
"""Literal Monte Carlo reconstruction of BRM Eqs. (28)--(32).

The published paper specifies the integral, the frequency map
``nu = a tan(u)``, and ``N_b = 10^9`` samples, but does not specify the
infinite-momentum proposal, random-number generator, batching, or error
estimator.  This module therefore reconstructs the most literal sampler that
can be inferred from the paper:

* rotate ``q1`` to the z axis and include the exact ``4*pi`` angular factor;
* sample ``k`` and ``q2`` uniformly in their unit Fermi balls;
* impose the two particle-hole Heaviside factors by rejection;
* retain the Coulomb denominators in the stochastic integrand;
* sample lambda and ``u in (0, pi/2)`` uniformly.

Two standard maps of ``q1 in (0,infinity)`` are supplied because the paper does
not state which one was used.  A third, explicitly truncated map is included
only to diagnose sensitivity to an undisclosed finite momentum cutoff.  For the
two infinite maps the target integral is independent of map and scale;
sensitivity at finite sample count is a sampling diagnostic.

This intentionally inefficient estimator is for comparison with
``benites_direct.py``, whose exit-cap and Coulomb-centred proposals evaluate the
same equations with much lower variance.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from scipy.stats import qmc

import benites_direct as bd


QMap = Literal["rational", "tangent", "finite"]
Sequence = Literal["sobol", "pseudo"]


def sample_unit_ball(u: np.ndarray) -> np.ndarray:
    """Map independent U(0,1) triples to a uniform unit-ball sample."""
    u = bd._clip_unit(u)
    radius = u[:, 0] ** (1.0 / 3.0)
    direction = bd._sample_direction(u[:, 1], u[:, 2])
    return radius[:, None] * direction


def q_map(t: np.ndarray, scale: float, kind: QMap) -> tuple[np.ndarray, np.ndarray]:
    """Return q and the complete radial measure 4*pi*q^2*dq/dt."""
    if not (scale > 0.0 and math.isfinite(scale)):
        raise ValueError("q scale must be finite and positive")
    t = bd._clip_unit(t)
    if kind == "rational":
        qv = scale * t / (1.0 - t)
        dqdt = scale / (1.0 - t) ** 2
    elif kind == "tangent":
        angle = 0.5 * np.pi * t
        qv = scale * np.tan(angle)
        dqdt = 0.5 * np.pi * scale / np.cos(angle) ** 2
    elif kind == "finite":
        # Diagnostic only: here ``scale`` is the hard upper cutoff q_max.
        qv = scale * t
        dqdt = np.full_like(t, scale)
    else:
        raise ValueError(f"unknown q map {kind!r}")
    return qv, 4.0 * np.pi * qv**2 * dqdt


def evaluate_term_literal(
    points: np.ndarray,
    rs: float,
    zeta: float,
    term: Literal[1, 2],
    q_scale: float,
    qmap: QMap,
) -> np.ndarray:
    """Literal axis-reduced estimator of BRM Eq. (28) or Eq. (29)."""
    if points.ndim != 2 or points.shape[1] != 9:
        raise ValueError("the axis-reduced literal sampler requires 9 coordinates")
    u = bd._clip_unit(points)
    lam = u[:, 0]
    qv, qweight = q_map(u[:, 1], q_scale, qmap)
    k = sample_unit_ball(u[:, 2:5])
    q2 = sample_unit_ball(u[:, 5:8])

    k_shift = k.copy()
    q2_shift = q2.copy()
    k_shift[:, 2] += qv
    q2_shift[:, 2] += qv
    valid = (
        np.einsum("ij,ij->i", k_shift, k_shift) > 1.0
    ) & (
        np.einsum("ij,ij->i", q2_shift, q2_shift) > 1.0
    )

    aq2 = 0.5 * qv**2 + qv * q2[:, 2]
    ak = 0.5 * qv**2 + qv * k[:, 2]
    tan_u = np.tan(0.5 * np.pi * u[:, 8])

    pref = 3.0 * bd.ALPHA * rs / (16.0 * np.pi**7)
    coupling_kernel = pref * lam**2 / np.maximum(qv, 1.0e-300) ** 4
    measure = qweight * bd.VBALL**2 * (0.5 * np.pi)
    common = measure * coupling_kernel

    out = np.zeros_like(qv)
    good = valid & (aq2 > 0.0) & (ak > 0.0) & np.isfinite(common)
    if not np.any(good):
        return out

    if term == 1:
        coulomb2 = np.einsum(
            "ij,ij->i", k[good] - q2[good], k[good] - q2[good]
        )
        dd = bd._divided_difference(
            qv[good], aq2[good], ak[good], tan_u[good], lam[good], rs, zeta
        )
        out[good] = common[good] * dd / np.maximum(coulomb2, 1.0e-300)
    elif term == 2:
        coulomb = k[good] + q2[good]
        coulomb[:, 2] += qv[good]
        coulomb2 = np.einsum("ij,ij->i", coulomb, coulomb)
        fa = bd.f_spin_sum(
            qv[good], aq2[good] * tan_u[good], lam[good], rs, zeta
        )
        fb = bd.f_spin_sum(
            qv[good], ak[good] * tan_u[good], lam[good], rs, zeta
        )
        out[good] = -common[good] * (fa + fb) / (
            np.maximum(coulomb2, 1.0e-300) * (aq2[good] + ak[good])
        )
    else:
        raise ValueError("term must be 1 or 2")
    return out


def evaluate_term_importance_map(
    points: np.ndarray,
    rs: float,
    zeta: float,
    term: Literal[1, 2],
    q_scale: float,
    qmap: QMap,
) -> np.ndarray:
    """BRM estimator with exit-cap and Coulomb-centred proposals.

    This is algebraically the same estimator as ``benites_direct.evaluate_term``
    but permits both published-infinite-domain diagnostics and an explicit hard
    cutoff through :func:`q_map`.
    """
    if points.ndim != 2 or points.shape[1] != 9:
        raise ValueError("the axis-reduced importance sampler requires 9 coordinates")
    u = bd._clip_unit(points)
    lam = u[:, 0]
    qv, qweight = q_map(u[:, 1], q_scale, qmap)
    k, kvol = bd.sample_exit_cap(qv, u[:, 2:5])
    direction = bd._sample_direction(u[:, 6], u[:, 7])

    qvec = np.zeros_like(k)
    qvec[:, 2] = qv
    if term == 1:
        center = k
        radius = np.full_like(qv, 2.0)
    elif term == 2:
        center = -k - qvec
        radius = np.linalg.norm(center, axis=1) + 1.0
    else:
        raise ValueError("term must be 1 or 2")
    rho = radius * u[:, 5]
    q2 = center + rho[:, None] * direction
    valid = bd._valid_exit_cap(q2, qv)

    aq2 = 0.5 * qv**2 + qv * q2[:, 2]
    ak = 0.5 * qv**2 + qv * k[:, 2]
    tan_u = np.tan(0.5 * np.pi * u[:, 8])
    pref = 3.0 * bd.ALPHA * rs / (16.0 * np.pi**7)
    coupling_kernel = pref * lam**2 / np.maximum(qv, 1.0e-300) ** 4
    common = (
        qweight * kvol * (4.0 * np.pi * radius) * (0.5 * np.pi) * coupling_kernel
    )

    out = np.zeros_like(qv)
    good = valid & (aq2 > 0.0) & (ak > 0.0) & np.isfinite(common)
    if not np.any(good):
        return out
    if term == 1:
        dd = bd._divided_difference(
            qv[good], aq2[good], ak[good], tan_u[good], lam[good], rs, zeta
        )
        out[good] = common[good] * dd
    else:
        fa = bd.f_spin_sum(
            qv[good], aq2[good] * tan_u[good], lam[good], rs, zeta
        )
        fb = bd.f_spin_sum(
            qv[good], ak[good] * tan_u[good], lam[good], rs, zeta
        )
        out[good] = -common[good] * (fa + fb) / (aq2[good] + ak[good])
    return out


def _points(m: int, sequence: Sequence, seed: int) -> np.ndarray:
    if sequence == "sobol":
        return qmc.Sobol(d=9, scramble=True, seed=seed).random_base2(m)
    if sequence == "pseudo":
        return np.random.default_rng(seed).random((2**m, 9))
    raise ValueError(f"unknown sequence {sequence!r}")


@dataclass
class LiteralEstimate:
    rs: float
    zeta: float
    m: int
    points_per_term_replicate: int
    replicates: int
    seed: int
    q_scale: float
    q_map: str
    sequence: str
    proposal: str
    delta1_ry: float
    delta1_err_ry: float
    delta2_ry: float
    delta2_err_ry: float
    full_kite_mry: float
    full_kite_err_mry: float
    replicate_values_delta1: list[float]
    replicate_values_delta2: list[float]


def estimate_literal(
    rs: float,
    zeta: float,
    m: int = 17,
    replicates: int = 8,
    seed: int = 20260807,
    q_scale: float = 2.0,
    qmap: QMap = "rational",
    sequence: Sequence = "sobol",
) -> LiteralEstimate:
    """Estimate the full kite from the literal BRM correction plus Onsager."""
    if m < 8:
        raise ValueError("m must be at least 8")
    if replicates < 2:
        raise ValueError("at least two replicates are needed for a scatter estimate")
    vals1: list[float] = []
    vals2: list[float] = []
    for rep in range(replicates):
        p1 = _points(m, sequence, seed + 2 * rep)
        p2 = _points(m, sequence, seed + 2 * rep + 1)
        vals1.append(float(np.mean(evaluate_term_literal(p1, rs, zeta, 1, q_scale, qmap))))
        vals2.append(float(np.mean(evaluate_term_literal(p2, rs, zeta, 2, q_scale, qmap))))
    a1 = np.asarray(vals1)
    a2 = np.asarray(vals2)
    mu1 = float(a1.mean())
    mu2 = float(a2.mean())
    err1 = float(a1.std(ddof=1) / np.sqrt(replicates))
    err2 = float(a2.std(ddof=1) / np.sqrt(replicates))
    full = bd.ONSAGER_RY + mu1 + mu2
    full_err = float(np.hypot(err1, err2))
    return LiteralEstimate(
        rs=rs,
        zeta=zeta,
        m=m,
        points_per_term_replicate=2**m,
        replicates=replicates,
        seed=seed,
        q_scale=q_scale,
        q_map=qmap,
        sequence=sequence,
        proposal="literal-uniform-balls",
        delta1_ry=mu1,
        delta1_err_ry=err1,
        delta2_ry=mu2,
        delta2_err_ry=err2,
        full_kite_mry=1000.0 * full,
        full_kite_err_mry=1000.0 * full_err,
        replicate_values_delta1=vals1,
        replicate_values_delta2=vals2,
    )


def estimate_importance_map(
    rs: float,
    zeta: float,
    m: int = 17,
    replicates: int = 8,
    seed: int = 20260807,
    q_scale: float = 2.0,
    qmap: QMap = "rational",
    sequence: Sequence = "sobol",
) -> LiteralEstimate:
    """Estimate the same BRM integral with the low-variance proposals."""
    if m < 8:
        raise ValueError("m must be at least 8")
    if replicates < 2:
        raise ValueError("at least two replicates are needed for a scatter estimate")
    vals1: list[float] = []
    vals2: list[float] = []
    for rep in range(replicates):
        p1 = _points(m, sequence, seed + 2 * rep)
        p2 = _points(m, sequence, seed + 2 * rep + 1)
        vals1.append(
            float(np.mean(evaluate_term_importance_map(p1, rs, zeta, 1, q_scale, qmap)))
        )
        vals2.append(
            float(np.mean(evaluate_term_importance_map(p2, rs, zeta, 2, q_scale, qmap)))
        )
    a1 = np.asarray(vals1)
    a2 = np.asarray(vals2)
    mu1 = float(a1.mean())
    mu2 = float(a2.mean())
    err1 = float(a1.std(ddof=1) / np.sqrt(replicates))
    err2 = float(a2.std(ddof=1) / np.sqrt(replicates))
    full = bd.ONSAGER_RY + mu1 + mu2
    full_err = float(np.hypot(err1, err2))
    return LiteralEstimate(
        rs=rs,
        zeta=zeta,
        m=m,
        points_per_term_replicate=2**m,
        replicates=replicates,
        seed=seed,
        q_scale=q_scale,
        q_map=qmap,
        sequence=sequence,
        proposal="exit-cap-coulomb-centred",
        delta1_ry=mu1,
        delta1_err_ry=err1,
        delta2_ry=mu2,
        delta2_err_ry=err2,
        full_kite_mry=1000.0 * full,
        full_kite_err_mry=1000.0 * full_err,
        replicate_values_delta1=vals1,
        replicate_values_delta2=vals2,
    )


def self_test() -> None:
    rng = np.random.default_rng(91)
    p = sample_unit_ball(rng.random((200000, 3)))
    assert np.all(np.einsum("ij,ij->i", p, p) <= 1.0 + 1.0e-14)
    assert abs(np.mean(np.einsum("ij,ij->i", p, p)) - 3.0 / 5.0) < 2.5e-3

    t = np.linspace(1.0e-5, 1.0 - 1.0e-5, 300000)
    for kind in ("rational", "tangent", "finite"):
        qv, weight = q_map(t, 1.7, kind)
        # Integral d^3q exp(-q), either complete or truncated at q=scale.
        got = float(np.trapezoid(weight * np.exp(-qv), t))
        if kind == "finite":
            exact = 8.0 * np.pi * (
                1.0 - np.exp(-1.7) * (1.0 + 1.7 + 0.5 * 1.7**2)
            )
        else:
            exact = 8.0 * np.pi
        assert abs(got - exact) < 2.0e-3, (kind, got, exact)

    tiny = qmc.Sobol(d=9, scramble=True, seed=123).random_base2(11)
    corr = evaluate_term_literal(tiny, 1.0e-12, 0.0, 1, 2.0, "rational")
    corr += evaluate_term_literal(tiny, 1.0e-12, 0.0, 2, 2.0, "rational")
    assert abs(float(np.mean(corr))) < 1.0e-10

    # With the rational map this must agree pointwise with the deposited
    # low-variance implementation.
    a = evaluate_term_importance_map(tiny, 3.0, 0.0, 1, 2.0, "rational")
    b = bd.evaluate_term(tiny, 3.0, 0.0, 1, 2.0)
    assert np.array_equal(a, b)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rs", type=float, default=5.0)
    ap.add_argument("--zeta", type=float, default=0.0)
    ap.add_argument("--m", type=int, default=17)
    ap.add_argument("--replicates", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260807)
    ap.add_argument("--q-scale", type=float, default=2.0)
    ap.add_argument(
        "--q-map", choices=("rational", "tangent", "finite"), default="rational"
    )
    ap.add_argument("--sequence", choices=("sobol", "pseudo"), default="sobol")
    ap.add_argument(
        "--proposal", choices=("literal", "importance"), default="literal"
    )
    ap.add_argument("--out", type=Path)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        print("[PASS] literal BRM sampler self-tests")
        return
    estimator = estimate_literal if args.proposal == "literal" else estimate_importance_map
    est = estimator(
        args.rs, args.zeta, args.m, args.replicates, args.seed,
        args.q_scale, args.q_map, args.sequence,
    )
    print(
        f"BRM Eqs. (28)-(32): rs={est.rs:g}, zeta={est.zeta:g}, "
        f"proposal={est.proposal}, map={est.q_map}, scale={est.q_scale:g}, "
        f"sequence={est.sequence}"
    )
    print(f"  Delta1 = {est.delta1_ry:+.8e} +- {est.delta1_err_ry:.2e} Ry")
    print(f"  Delta2 = {est.delta2_ry:+.8e} +- {est.delta2_err_ry:.2e} Ry")
    print(f"  full kite = {est.full_kite_mry:+.4f} +- {est.full_kite_err_mry:.4f} mRy")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(asdict(est), indent=2) + "\n")
        print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()
