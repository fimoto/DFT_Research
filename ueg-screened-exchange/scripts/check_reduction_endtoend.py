#!/usr/bin/env python3
"""End-to-end consistency check of the JPA Sec. 4.2 Lorentzian reduction.

This is a numerical consistency check, NOT a rigorous bound: the small-q tail
uses beta(Lambda) fitted from the lowest measured nodes, and the large-q tail
uses the leading term of JPA large-momentum expansion (app-eq:largeq-block) whose constant is not known.  Run with
several (qlo, qhi) to see the tail dependence.

K[b_Lambda] is evaluated directly from Eq. (KD) by Monte Carlo over p, k and a
Gauss-Legendre grid in q on [qlo, qhi]; the tails q < qlo (JPA fixed-width small-momentum expansion (app-eq:kappa):
kappa q log(1/q) + beta(inf) q) and q > qhi (JPA large-momentum expansion (app-eq:largeq-block): (6/pi^2)Lam/(q(q+2Lam)^2))
are added analytically INSIDE this script, and the statistical error is the
scatter over seeds.  Compared with the 1-D representation e2x*S(Lambda).

Usage: python scripts/check_reduction_endtoend.py [Lambda] [nseeds] [N]   (~8 min per seed at defaults)
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from ueg_exchange.constants import A0, E2X
from ueg_exchange import lorentzian as B

S_REF = {1.0: 2.33400, 10.0: 1.66894, 0.1: 1.03508, 0.01: 0.166379}   # 1-D formula, JPA Sec. 4.2


def J(w, a, b):
    aa, bb = np.abs(a), np.abs(b)
    same = a * b > 0
    return np.where(same, 1.0 / (2 * w * (aa + w) * (bb + w)),
                    -(aa + bb + 2 * w) / (2 * w * (aa + bb) * (aa + w) * (bb + w)))


def cap_volume(q):
    return np.where(q >= 2, 4 * np.pi / 3, np.pi * q * (12 - q * q) / 12)


def sample_cap(q, n, rng):
    out = np.empty((n, 3)); got = 0
    while got < n:
        m = int(1.5 * (n - got) / max(1e-9, cap_volume(q) / (4 * np.pi / 3))) + 10
        v = rng.normal(size=(m, 3)); v /= np.linalg.norm(v, axis=1)[:, None]; v *= rng.random(m)[:, None] ** (1 / 3)
        ok = v[:, 0] ** 2 + v[:, 1] ** 2 + (v[:, 2] + q) ** 2 > 1
        v = v[ok][: n - got]; out[got: got + len(v)] = v; got += len(v)
    return out


def K_block_grid(Lam, rng, N, nq, qlo, qhi):
    tq, wq = np.polynomial.legendre.leggauss(nq)
    lo, hi = np.log(qlo), np.log(qhi)
    t = 0.5 * (hi - lo) * tq + 0.5 * (hi + lo); wq = 0.5 * (hi - lo) * wq; qs = np.exp(t)
    total = 0.0
    per_node = []
    for q, w in zip(qs, wq):
        pp = sample_cap(q, N, rng); Vc = cap_volume(q)
        sp = np.where(rng.random(N) < 0.5, 1.0, -1.0)
        p = np.where(sp[:, None] > 0, pp, -pp - np.array([0, 0, q]))
        dp = p[:, 2] * q + q * q / 2
        Rmax = 2.0 + q; rho = Rmax * rng.random(N)
        n = rng.normal(size=(N, 3)); n /= np.linalg.norm(n, axis=1)[:, None]
        k = p + rho[:, None] * n
        kin = np.sum(k * k, 1) <= 1; kq = (k[:, 0] ** 2 + k[:, 1] ** 2 + (k[:, 2] + q) ** 2) > 1
        Dk = np.where(kin & kq, 1.0, np.where((~kin) & (~kq), -1.0, 0.0))
        dk = k[:, 2] * q + q * q / 2
        om = Lam * q
        val = np.where(Dk != 0, sp * Dk * om * om * J(om, dp, dk), 0.0)
        Xq = A0 * 4 * np.pi * (2 * Vc) * (4 * np.pi * Rmax) * val.mean()
        per_node.append(Xq)
        total += w * q * Xq
    return total, (qs, np.array(per_node))


def small_q_tail(Lam, qlo, grid_q, grid_X_list):
    """int_0^qlo e_Lambda dq using JPA fixed-width small-momentum expansion (app-eq:kappa), X = kappa(Lam) q log(1/q) + beta(Lam) q + o(q).

    beta(Lambda) is not available in closed form, so it is fitted from the two
    lowest measured grid nodes rather than replaced by beta(infinity): the two
    differ at finite Lambda, and using beta(infinity) here was an error.
    """
    kap = float(B.kappa(Lam))
    bs = [np.mean([X / q - kap * np.log(1 / q) for q, X in zip(grid_q[:2], gx[:2])])
          for gx in grid_X_list]
    b, b_sd = float(np.mean(bs)), float(np.std(bs, ddof=1)) if len(bs) > 1 else (float(np.mean(bs)), 0.0)
    if len(bs) > 1:
        b, b_sd = float(np.mean(bs)), float(np.std(bs, ddof=1))
    else:
        b, b_sd = float(bs[0]), float("nan")
    tail = kap * qlo ** 2 / 2 * (np.log(1 / qlo) + 0.5) + b * qlo ** 2 / 2
    return tail, b, b_sd


def large_q_tail(Lam, qhi):
    """int_qhi^inf e_Lambda dq from the leading term of JPA large-momentum expansion (app-eq:largeq-block).

    The remainder of D(ii) is C[1/(q^3(1+s_lambda)) + 1/q^4] with an unknown C, so no
    error is returned: vary qhi to see how much the tail moves.
    """
    from scipy.integrate import quad
    return quad(lambda q: 6 / np.pi ** 2 * Lam / (q * (q + 2 * Lam) ** 2), qhi, np.inf)[0]


if __name__ == "__main__":
    Lam = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    nseeds = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 100000
    qlo, qhi, nq = 0.02, 15.0, 32
    vals, node_vals, qs = [], [], None
    for s in range(nseeds):
        v, (qs, per_node) = K_block_grid(Lam, np.random.default_rng(100 + s), N, nq, qlo, qhi)
        vals.append(v); node_vals.append(per_node)
        print(f"  seed {s}: grid part = {1e3*v:8.3f} mRy")
    vals = np.array(vals); mean = vals.mean()
    err = vals.std(ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else float("nan")
    small, beta_fit, beta_sd = small_q_tail(Lam, qlo, qs, node_vals)
    large = large_q_tail(Lam, qhi)
    tot = mean + small + large
    ref = E2X * S_REF.get(Lam, float("nan"))
    print(f"Lambda={Lam}: grid {1e3*mean:.3f}(seed scatter {1e3*err:.3f})"
          f" + small-q tail {1e3*small:.3f} [beta({Lam}) fitted = {beta_fit:.4f} +- {beta_sd:.4f}]"
          f" + large-q tail {1e3*large:.3f} (leading term only)"
          f" = {1e3*tot:.3f} mRy   |   1-D formula e2x*S({Lam}) = {1e3*ref:.3f} mRy")
    print("  the two tails are model terms, not bounds; rerun with other qlo/qhi to see their size")
