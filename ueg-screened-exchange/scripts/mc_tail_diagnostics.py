#!/usr/bin/env python3
"""Sec. 10: heavy-tail diagnostics of the literal BRM sampler.
Reports the Hill tail index, the one-sidedness of the extreme tail, and the
share of the answer carried by rare samples."""
import os, sys, warnings, numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from ueg_exchange import brm

RS = float(os.environ.get("UEG_RS", 5.0))
N = int(os.environ.get("UEG_N", 250000))
B = int(os.environ.get("UEG_BATCHES", 32))

for term in (1, 2):
    alpha, v = brm.tail_index(RS, term=term, n=N, batches=B)
    a = np.abs(v); order = np.argsort(a)[::-1]
    print(f"term {term}: N={v.size}  mean={1e3*v.mean():.3f} mRy  Hill alpha={alpha:.3f}")
    for frac in (1e-6, 1e-5, 1e-4, 1e-3):
        k = max(1, int(v.size*frac)); sel = order[:k]
        print(f"   top {frac:g}: {1e3*v[sel].sum()/v.size:8.3f} mRy "
              f"({100*v[sel].sum()/v.sum():5.1f}% of term), negative fraction "
              f"{(v[sel] < 0).mean():.2f}")
print("\nalpha < 2 means the variance is infinite: sigma/sqrt(N) is not a "
      "confidence interval; the sample mean is unbiased in expectation but its finite-N median is displaced in a definite direction.")
