"""RPA ring sum, BRM Eq. (12), and its high-density expansion.

eps_r(rs) = gamma int dkappa kappa^3 int dx [log(1+c) - c],
    gamma = 3/(2 pi alpha_L^2 rs^2),  c = (q_tf_sq/kappa^2) lindhard_shape_spin(kappa, kappa x).
Validated against BRM Table VIII.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from .constants import ALPHA_L, Q_TF_SQ_COEFF
from .lindhard import lindhard_shape_spin


def eps_r(rs, zeta=0.0, nq=500, nx=500, lo=1e-8, hi=1e4):
    q_tf_sq = Q_TF_SQ_COEFF * rs
    tq, wq = leggauss(nq)
    llo, lhi = np.log(lo), np.log(hi)
    t = 0.5 * (lhi - llo) * tq + 0.5 * (lhi + llo)
    wq = 0.5 * (lhi - llo) * wq
    kap = np.exp(t)
    tx, wx = leggauss(nx)
    s = 0.5 * (tx + 1.0)
    x = s / (1.0 - s)
    wx = 0.5 * wx / (1.0 - s) ** 2
    K, X = np.meshgrid(kap, x, indexing="ij")
    c = q_tf_sq / K ** 2 * lindhard_shape_spin(K, K * X, zeta)
    F = np.log1p(c) - c
    inner = (F * wx).sum(1)
    gamma = 3.0 / (2 * np.pi * ALPHA_L ** 2 * rs ** 2)
    return gamma * (wq * kap ** 4 * inner).sum()   # kappa^3 * dkappa = kappa^4 dt
