"""Checks requested in the 3 Oct 2026 reviews: (1) rs=2e4 seed dependence of the 16-replicate
validation runs; (2) Estimator B residual delta-e at rs=10 from common samples."""
import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # run from the repository root
import sys,json,time; from pathlib import Path
sys.path.insert(0,'production'); sys.path.insert(0,'spin')
from compute import direct; import reviewer_numerics as rn
V=Path('validation/results'); out=V/'revision_checks.json'
rec=json.load(open(out)) if out.exists() else {}
for seed in (20262001,20263001):
    k=f'rs2e4_z0_seed{seed}'
    if k in rec: continue
    t=time.time(); r=direct(2e4,0.0,dict(m=19,replicates=16,q_scale='qstar'),dict(seed=seed,chunk_power=15))
    rec[k]=dict(mean_mry=r['mean_mry'],se_mry=r['se_mry'],term_replicates_ry=r['term_replicates_ry']); json.dump(rec,open(out,'w'),indent=1)
    print(k,'%.3f +- %.3f (%.0fs)'%(r['mean_mry'],r['se_mry'],time.time()-t),flush=True)
for z in (0.0,1.0):
    k=f'axisB_rs10_z{z}'
    if k in rec: continue
    t=time.time(); v=rn.energy(rs=10.0,zeta=z,power=15,reps=16,nq=24,order=24,nk=48,seed=30301001)
    rec[k]=dict(mean_ry=float(v['mean']),se_ry=float(v['se']),replicates_ry=[float(x) for x in v['replicates']]); json.dump(rec,open(out,'w'),indent=1)
    print(k,'%.4f +- %.4f mRy (%.0fs)'%(1000*v['mean'],1000*v['se'],time.time()-t),flush=True)
print('REV_DONE')
