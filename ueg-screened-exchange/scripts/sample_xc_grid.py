"""Regenerate figs/xc_grid.npy: the static geometric kernel X_c(q).

X_c is defined by Eq. (Xcdef) of the JPA paper.  Both loop momenta are drawn
uniformly in a box that covers the two crescents C(q).  Drawing both inside the
unit ball instead would not miss the support -- a point inside the ball can
still leave it after the shift -- but it restricts the pair to the sector in
which the two energy differences have the same sign, where the static frequency
kernel vanishes, so the estimate collapses to zero.  Settings below are those
used for the published grid.

    python scripts/sample_xc_grid.py            # 4e6 points per q, seed 3
"""
import numpy as np

N = 4_000_000
SEED = 3
QS = np.array([0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.55, 0.7, 0.85, 1.0, 1.2, 1.4,
               1.6, 1.8, 1.95, 2.05, 2.2, 2.5, 3.0, 3.5, 4.0, 5.0])


def static_frequency_kernel(a, b):
    return (np.sign(a * b) - 1.0) / (2.0 * (np.abs(a) + np.abs(b)))


def Xc(q, rng, n=N):
    L = 1.0 + q
    sc = np.array([1.0, 1.0, L])
    p = rng.uniform(-1, 1, size=(n, 3)) * sc
    k = rng.uniform(-1, 1, size=(n, 3)) * sc
    V = 8 * L
    qv = np.array([0.0, 0.0, q])
    Dp = ((np.linalg.norm(p, axis=1) < 1).astype(np.int8)
          - (np.linalg.norm(p + qv, axis=1) < 1).astype(np.int8))
    Dk = ((np.linalg.norm(k, axis=1) < 1).astype(np.int8)
          - (np.linalg.norm(k + qv, axis=1) < 1).astype(np.int8))
    D = (Dp * Dk).astype(float)
    m = D != 0
    dp = q * (p[:, 2] + q / 2)
    dk = q * (k[:, 2] + q / 2)
    r2 = np.sum((p[m] - k[m]) ** 2, axis=1)
    f = D[m] / np.maximum(r2, 1e-12) * static_frequency_kernel(dp[m], dk[m])
    # The estimator averages over all n draws, the n - M of them with Delta = 0
    # contributing zero.  The standard error must therefore use the variance of
    # the full sample, not that of the M nonzero values: dropping the zeros
    # ignores the fluctuation of M itself and understates the error (by about
    # 31 % at q = 5).
    norm = V * V / (16 * np.pi)
    s1 = f.sum()
    s2 = float(np.dot(f, f))
    val = norm * s1 / n
    var_mean = (s2 - s1 * s1 / n) / (n * (n - 1))
    err = norm * np.sqrt(max(var_mean, 0.0))
    return val, err


if __name__ == "__main__":
    rng = np.random.default_rng(SEED)
    out = [(q,) + Xc(q, rng) for q in QS]
    a = np.array(out)
    np.save("figs/xc_grid.npy", a.T)
    print("wrote figs/xc_grid.npy  (%d points per q, seed %d)" % (N, SEED))
