"""Thin wrapper around the direct evaluation of BRM Eqs. (28)-(32).

The estimator itself lives in external/benites_direct.py (exit-cap +
Coulomb-centred importance sampling) and external/brm_literal_mc.py
(uniform-ball sampler, used only as a diagnostic).  Finiteness of the second
moment is not proved for either.  Empirically the production sampler behaves
as if it were finite (Hill tail index 3.8 and 2.7 for the two terms of BRM
Eq. (29), replicate scatter falling as N^-0.49), while the uniform-ball one
does not (1.3-1.6).  Every uncertainty reported by this module is a standard error of the
replicate mean, s_R/sqrt(R) with s_R the sample standard deviation over the R
independently randomized replicates; see LIMITATIONS.md.  This module locates them and adds the
q-resolved and tail diagnostics used in Secs. 8.3 and 10.
"""
import os
import sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))          # .../src/ueg_exchange
_ROOT = os.path.normpath(os.path.join(_HERE, os.pardir, os.pardir))   # repo root
_EXT = os.path.join(_ROOT, "external")
if _EXT not in sys.path:
    sys.path.insert(0, _EXT)


def _bd():
    import benites_direct as bd
    return bd


def kite(rs, zeta=0.0, m=19, replicates=8, q_scale=2.0):
    """Full kite in mRy per electron with replicate-to-replicate scatter (not a confidence interval)."""
    e = _bd().estimate(rs=rs, zeta=zeta, m=m, replicates=replicates, q_scale=q_scale)
    return e.full_kite_mry, e.full_kite_err_mry


def replicate_values(rs, zeta=0.0, m=19, replicates=8, q_scale=2.0, seed=1000):
    """Per-replicate values of Delta = kite - Onsager, in Ry.

    Returns the array of `replicates` independent estimates, so that a fit over
    several densities can be done replicate by replicate.  This is required
    whenever the densities share a seed: their errors are then correlated and a
    diagonal covariance is not valid.
    """
    e = _bd().estimate(rs=rs, zeta=zeta, m=m, replicates=replicates,
                       q_scale=q_scale, seed=seed)
    return (np.asarray(e.replicate_values_delta1)
            + np.asarray(e.replicate_values_delta2))


def correction(rs, zeta=0.0, m=19, replicates=8, q_scale=2.0, seed=None):
    """Delta = kite - Onsager, in Ry.

    A fixed seed makes the estimate at different rs share one random-number
    stream.  That is useful for differences, but it also correlates the errors
    across densities, so a fit over several rs must not treat them as independent;
    fit each replicate separately instead (scripts/check_ln2_window.py prints
    the measured correlation range for its own settings).
    """
    kw = {} if seed is None else {"seed": seed}
    e = _bd().estimate(rs=rs, zeta=zeta, m=m, replicates=replicates,
                       q_scale=q_scale, **kw)
    return e.correction_ry, e.correction_err_ry


def q_resolved(rs, edges, zeta=0.0, m=18, replicates=6, q_scale=2.0, seed=1000):
    """Contribution to Delta (Ry) binned in the transfer momentum q1."""
    from scipy.stats import qmc
    bd = _bd()
    tot = np.zeros(len(edges) - 1)
    for rep in range(replicates):
        for term in (1, 2):
            s = qmc.Sobol(d=9, scramble=True, seed=seed + 10 * rep + term)
            p = s.random_base2(m)
            q, _ = bd._q_map(bd._clip_unit(p[:, 1]), q_scale)
            v = bd.evaluate_term(p, rs, zeta, term, q_scale)
            idx = np.digitize(q, edges) - 1
            ok = (idx >= 0) & (idx < len(edges) - 1)
            tot += np.bincount(idx[ok], weights=v[ok], minlength=len(edges) - 1) / len(v)
    return tot / replicates


def tail_index(rs, zeta=0.0, term=1, n=250000, batches=32, k=1000, seed=2024):
    """Hill estimator of the tail index of the literal uniform-ball sampler."""
    import brm_literal_mc as lm
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(batches):
        p = rng.random((n, 9))
        vals.append(lm.evaluate_term_literal(p, rs, zeta, term, 2.0, "rational"))
    v = np.concatenate(vals)
    a = np.sort(np.abs(v))[::-1]
    return 1.0 / (np.mean(np.log(a[:k])) - np.log(a[k])), v
