#!/usr/bin/env python3
"""JPA Sec. 5.2: two numerical confirmations of the rs log^2 rs term.
(a) fit of the direct BRM evaluation at high density;
(b) the q-resolved window slope at rs = 1e-3.
Requires external/benites_direct.py.  Runtime ~20 min at default settings."""
import os, sys, warnings, numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from ueg_exchange import brm, constants as C, highdensity as H

M = int(os.environ.get("UEG_M", 19))
REP = int(os.environ.get("UEG_REP", 8))

print("(a) high-density fit of Delta = kite - Onsager")
print("    All densities share one random-number stream, so their errors are")
print("    correlated (the measured range is printed below).  The coefficients")
print("    are therefore fitted separately in each replicate and the scatter of")
print("    those fits gives the error; a diagonal covariance is not valid here")
print("    and neither is a chi^2 built from it.")
rss = np.array([0.002, 0.005, 0.01, 0.02, 0.05, 0.1])
L = np.log(rss)
# One call per density with the same seed: the estimator internally uses the
# seed pair (seed+2r, seed+2r+1) for the two terms of replicate r, so the
# replicates stay independent of each other while the densities share their
# points.  Building the replicates by hand from seed+r would make the second
# term of replicate r reuse the point set of the first term of replicate r+1.
per_rep = np.stack([brm.replicate_values(rs, m=M, replicates=REP, seed=1000)
                    for rs in rss], axis=1)
D = per_rep.mean(0)
E = per_rep.std(0, ddof=1)/np.sqrt(REP)
for rs, d, e in zip(rss, D, E):
    print(f"  rs={rs:<6} Delta={1e3*d:9.5f} +- {1e3*e:.5f} mRy")
corr = np.corrcoef(per_rep.T)
print(f"  density-density correlation: {corr[0,1]:.2f} to {corr[0,-1]:.2f}")
for cols, lab in [([L, np.ones_like(L)], "C1 L + C"),
                  ([L**2, L, np.ones_like(L)], "C2 L^2 + C1 L + C")]:
    A = np.stack(cols, 1)
    cs = np.array([np.linalg.lstsq(A, per_rep[r]/rss, rcond=None)[0]
                   for r in range(REP)])
    m = cs.mean(0)
    se = cs.std(0, ddof=1)/np.sqrt(REP)
    print(f"  {lab:20s} {np.round(m,6)} +- {np.round(se,6)}")
print(f"  analytic C_ln2 = {C.C_LN2:.6f} Ry")

print("\n(b) q-resolved window slope at rs = 1e-3")
rs = 1e-3
edges = np.exp(np.linspace(np.log(0.02), np.log(2.0), 13))
h = brm.q_resolved(rs, edges, m=M, replicates=6)
xs, ys = [], []
print("   q_mid    q dDelta/dq [mRy]")
for i in range(len(h)):
    qm = np.sqrt(edges[i]*edges[i+1]); val = 1e3*h[i]/np.log(edges[i+1]/edges[i])
    print(f"  {qm:7.4f}  {val:10.5f}")
    if 0.08 < qm < 0.7:
        xs.append(np.log(1/qm)); ys.append(val)
slope = np.polyfit(xs, ys, 1)[0]
pred = -(8*np.pi*C.A0*C.Q_TF_SQ_COEFF*rs/3)*H.calA()[1]*1e3
print(f"  fitted slope {slope:.5f} mRy   predicted {pred:.5f} mRy")
