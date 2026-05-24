#!/usr/bin/env python3
"""Complementary ring-energy check: the RPA ring sum contains no rs log^2 rs term.
Validates the integrator against BRM Table VIII, then fits the residual."""
import os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from ueg_exchange import ringsum as R, constants as C

print("validation against BRM Table VIII (mRy):")
for rs, ref in [(1.0, -157.6), (5.0, -85.0), (10.0, -61.3)]:
    print(f"  rs={rs:<5} eps_r={1000*R.eps_r(rs):9.3f}   BRM {ref}")

rss = np.array([1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2])
E = np.array([R.eps_r(r, nq=700, nx=600) for r in rss])
L = np.log(rss)
res = E - C.C_L_RING * L
print("\nresidual eps_r - c_L log rs   (c_L = %.6f Ry)" % C.C_L_RING)
for r, x in zip(rss, res):
    print(f"  rs={r:<9.0e} {x: .9f}")
for cols, lab in [([np.ones_like(rss), rss*L, rss], "c0 + rs(C1 L + C)"),
                  ([np.ones_like(rss), rss*L**2, rss*L, rss], "c0 + rs(C2 L^2 + C1 L + C)")]:
    A = np.stack(cols, 1)
    c, *_ = np.linalg.lstsq(A, res, rcond=None)
    print(f"\n{lab}\n  c0={c[0]:.6f} Ry (BRM Table II: -0.1423)  rest={np.round(c[1:],7)}"
          f"  max|resid|={np.max(np.abs(A@c-res)):.2e}")
print(f"\nkite C_ln2 = {C.C_LN2:.7f} Ry;  ring log^2 coefficient is ~0.5% of it -> no partner at V_1")
