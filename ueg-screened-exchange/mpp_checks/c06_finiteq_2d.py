"""SI Sec. S6 and SI Table S3: the first-order exchange polarization of the
two-dimensional gas at finite wave number, R(q,z) - R_2(z), at z = 0.3 and 1 for
q = 0.08, 0.04, 0.02 (main-text Sec. 3.1, property (iv)).

The integrator (finiteq2d.py) is described in its docstring.  Two node sets are
evaluated; the change between them is recorded.  The static comparison with
Table I of Das Sarma (z = 0, the four values Delta R/R(0) of SI Sec. S6) is NOT
reproduced here: at z = 0 the integrand is singular where delta_p or delta_k
vanishes, and this integrator needs nu = q z > 0.
"""
import numpy as np, mpmath as mp
from common import Recorder
from finiteq2d import R_finite, R2

R = Recorder("c06_finiteq_2d", "finite-q first-order exchange polarization in 2D (SI Sec. S6, Table S3)")
table = {('0.3', '0.08'): '-1.2186e-4', ('0.3', '0.04'): '-3.3626e-5', ('0.3', '0.02'): '-9.1990e-6',
         ('1', '0.08'): '-5.1759e-6', ('1', '0.04'): '0.9888e-6', ('1', '0.02'): '0.8176e-6'}
coarse, fine = (96, 48, 1024, 256), (144, 72, 1536, 384)
vals = {}
for (zs, qs), pr in table.items():
    z, q = float(zs), float(qs)
    v1 = R_finite(q, z, *coarse) - R2(z); v2 = R_finite(q, z, *fine) - R2(z)
    vals[(zs, qs)] = v2
    mant, ex = pr.split('e'); scale = 10 ** int(ex)
    R.digits("SI Table S3", f"R(q,z) - R_2(z) at z={zs}, q={qs} (in units of 1e{ex})", mp.mpf(v2) / scale, mant)
    R.note("SI Sec. S6", f"change between node sets at z={zs}, q={qs}", f"{abs(v2 - v1):.2e}")
r1 = vals[('0.3', '0.08')] / vals[('0.3', '0.04')]; r2 = vals[('0.3', '0.04')] / vals[('0.3', '0.02')]
R.close("SI Sec. S6", "z = 0.3: ratio of successive differences when q is halved, q = 0.08/0.04 ('about 3.6')", r1, 3.6, 0.1)
R.close("SI Sec. S6", "z = 0.3: ratio, q = 0.04/0.02 ('about 3.6')", r2, 3.6, 0.1)
R.exact("SI Sec. S6", "z = 1: the differences are below 1e-5 in magnitude and change sign", (max(abs(vals[('1', q)]) for q in ('0.08', '0.04', '0.02')) < 1e-5, vals[('1', '0.08')] < 0 < vals[('1', '0.04')]), (True, True))
R.note("SI Sec. S6", "static comparison with Das Sarma (z = 0)", "not reproduced by this script (see docstring)")
raise SystemExit(R.write(dict(node_sets=dict(coarse=coarse, fine=fine))))
