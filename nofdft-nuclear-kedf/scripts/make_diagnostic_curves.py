#!/usr/bin/env python3
"""
Fixed-density diagnostic NLF curves for the per-nucleus fitting figure.

For the released diagnostic kernels (kernels/diagnostics/):
  eg_NT14_{O16,Ca40} : per-nucleus fit WITHOUT the surface-tail band and
                       tail-margin penalty  -> "without tail-margin penalty"
  eg_VO14_O16        : NLF-value-only fit   -> optional fixed-density curve

each curve is NLF = 1/(1+chi^2) with chi evaluated on the FIXED KS density
using the same tau machinery as kernel_fit.py (this is the fit-level view; the
self-consistent endpoints come from the solver runs staged by
namelists/diagnostics/). Output: results/reference/diagnostics_v2/.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kernel_fit as K

ROOT = K.ROOT
OUT = os.path.join(ROOT, "results", "reference", "diagnostics_v2")
os.makedirs(OUT, exist_ok=True)

for tag, nuc in [("eg_NT14_O16", "O16"), ("eg_NT14_Ca40", "Ca40"),
                 ("eg_VO14_O16", "O16")]:
    kern = np.loadtxt(os.path.join(ROOT, "kernels", "diagnostics", tag + ".dat"))
    shift = float(open(os.path.join(ROOT, "kernels", "diagnostics",
                                    tag + "_shift.txt")).read())
    d = np.loadtxt(os.path.join(ROOT, "inputs", "ks", "ks_%s_nr250.dat" % nuc))
    r, rho = d[:, 0], d[:, 1]
    w = np.full_like(r, r[1] - r[0])
    k5 = np.loadtxt(os.path.join(ROOT, "results", "reference", "ks_nlf_%s.dat" % nuc))
    I = lambda x: np.interp(r, k5[:, 0], x)
    tauTF = np.maximum(I(k5[:, 2]), 1e-30)
    tauvW = I(k5[:, 3])
    t8 = np.loadtxt(os.path.join(K.GPD, "teacher_ddkf_%s.dat" % nuc))
    taudd = np.interp(r, t8[:, 0], t8[:, 6])
    fn, _ = K.spline_pair(np.interp(K.YTAB, kern[:, 0], kern[:, 1]))
    dtau = K.tau_response(r, w, rho, fn)
    tau = taudd + dtau - shift * np.maximum(rho, 0.0) ** (5.0 / 3.0)
    chi = (tau - tauvW) / tauTF
    nlf = 1.0 / (1.0 + chi ** 2)
    nlfKS = I(k5[:, 5])
    np.savetxt(os.path.join(OUT, "%s_%s_nlf_fixed.dat" % (nuc, tag)),
               np.column_stack([r, nlf, nlfKS]),
               header="r[fm]  NLF(fit, fixed KS density)  NLF_KS   [%s]" % tag)
    print("wrote %s_%s_nlf_fixed.dat" % (nuc, tag))
