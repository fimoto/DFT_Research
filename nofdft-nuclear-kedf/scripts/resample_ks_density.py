#!/usr/bin/env python3
"""
Resample the nr=500 KS species density onto the production nr=250 grid.

The production input densities inputs/ks/ks_<TAG>_nr250.dat are the LINEAR
interpolation of the nr=500 KS solution (namelists/ks/ks_<TAG>_nr500.nml) onto
the grid r_i = i * rbox/(nr+1), i=1..250, rbox=20 fm. Verified agreement with
the released files: max|drho| ~ 6e-10.

Usage: python3 resample_ks_density.py ks_O16_species.dat out_nr250.dat
"""
import sys
import numpy as np

src, dst = sys.argv[1], sys.argv[2]
d = np.loadtxt(src)
nr, rbox = 250, 20.0
r = np.arange(1, nr + 1) * rbox / (nr + 1)
rho = np.interp(r, d[:, 0], d[:, 1])
np.savetxt(dst, np.column_stack([r, rho]))
print("wrote %s (%d points)" % (dst, nr))
