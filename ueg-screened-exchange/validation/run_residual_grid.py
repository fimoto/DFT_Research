"""Estimator B interpolation residual at rs=10 versus momentum grid (nq=24,48,96), same sample
size (2^13 points per replicate, 16 replicates) and same seed for all zeta, so the residual is
paired across zeta at each grid."""
import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # run from the repository root
import sys,json,time; from pathlib import Path
sys.path.insert(0,'spin'); import reviewer_numerics as rn
out=Path('validation/results/residual_grid_rs10.json'); rec=json.load(open(out)) if out.exists() else {}
for nq in (24,48,96):
    for z in (0.0,0.4,0.8,1.0):
        k=f'nq{nq}_z{z}'
        if k in rec: continue
        t=time.time(); v=rn.energy(rs=10.0,zeta=z,power=13,reps=16,nq=nq,order=24,nk=48,seed=30301001)
        rec[k]=dict(nq=nq,zeta=z,mean_ry=float(v['mean']),se_ry=float(v['se']),replicates_ry=[float(x) for x in v['replicates']],wall_s=time.time()-t)
        json.dump(rec,open(out,'w'),indent=1); print(k,'%.4f +- %.4f mRy (%.0fs)'%(1000*v['mean'],1000*v['se'],time.time()-t),flush=True)
print('GRID_RES_DONE')
