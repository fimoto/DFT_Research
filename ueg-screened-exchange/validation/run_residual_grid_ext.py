"""Estimator B interpolation residual, extended: r_s=10 with 2^15 points (n_q=48, converged grid),
and r_s=1 and 60 with 2^13 points; zeta=0,0.4,0.8,1 share the seed at each (r_s, n_q), so the
residual is formed replicate by replicate. Resumable."""
import os as _os, sys, json, time
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(_os.path.abspath(__file__)))))
sys.path.insert(0,'spin'); import reviewer_numerics as rn
out='validation/results/residual_grid_ext.json'; rec=json.load(open(out)) if _os.path.exists(out) else {}
jobs=[(10.0,15,48),(1.0,13,48),(60.0,13,48)]
for rs,power,nq in jobs:
    for z in (0.0,0.4,0.8,1.0):
        k=f'rs{rs:g}_p{power}_nq{nq}_z{z}'
        if k in rec: continue
        t=time.time(); v=rn.energy(rs=rs,zeta=z,power=power,reps=16,nq=nq,order=24,nk=48,seed=30301001)
        rec[k]=dict(rs=rs,zeta=z,power=power,nq=nq,mean_ry=float(v['mean']),se_ry=float(v['se']),replicates_ry=[float(x) for x in v['replicates']])
        json.dump(rec,open(out,'w'),indent=1); print(k,'%.4f +- %.4f mRy (%.0fs)'%(1000*v['mean'],1000*v['se'],time.time()-t),flush=True)
print('EXT_DONE')
