"""Independent imaginary-axis integration (spin/reviewer_numerics.energy) at intermediate polarization,
with refinement checks on the same geometry samples. Results in Ry inside JSON; printed in mRy."""
import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # run from the repository root
import sys,json,time,hashlib
from pathlib import Path
import numpy as np
sys.path.insert(0,'spin'); import reviewer_numerics as rn
V=Path('validation/results/axis'); V.mkdir(parents=True,exist_ok=True)
def save(name,v):
    d=dict(mean_ry=float(v['mean']),se_ry=float(v['se']),replicates_ry=[float(x) for x in v['replicates']],settings=v['settings'],Q=float(v['Q']),
           source_sha=hashlib.sha256(Path('spin/reviewer_numerics.py').read_bytes()).hexdigest())
    (V/f'{name}.json').write_text(json.dumps(d,indent=1)); return d
jobs=[]
for rs in [1.0,10.0,60.0]:
    for z in [0.4,0.8]:
        jobs.append((f'base_rs{rs:g}_z{z}',dict(rs=rs,zeta=z,power=15,reps=16,nq=24,order=24,nk=48,seed=30301001)))
CONV=[]
for rs in [1.0,10.0,60.0]:
    for z in [0.4,0.8]:
        CONV.append((f'nq96_rs{rs:g}_z{z}',dict(rs=rs,zeta=z,power=12,reps=16,nq=96,order=24,nk=48,seed=30301001)))
# refinements at rs=10, zeta=0.8 on the same samples (power=13, same seed) 
b=dict(rs=10.0,zeta=0.8,power=13,reps=16,nq=24,order=24,nk=48,seed=30301001)
jobs.append(('ref_base',b))
for tag,ch in [('nq48',dict(nq=48)),('order36',dict(order=36)),('nk96',dict(nk=96)),('qmin001',dict(qmin=0.001)),('nq96',dict(nq=96))]:
    jobs.append((f'ref_{tag}',{**b,**ch}))
# qmax doubled: default upper = 30*max(Q/x_sigma,1); use largest default (smallest x_sigma) times 2
Q=(4*rn.A_SCREEN*10.0/3)**.25; xmin=(1-0.8)**(1/3); jobs.append(('ref_qmax2',{**b,'qmax':2*30*max(Q/xmin,1)}))
# identity with the independent integrator: e(rs,1) vs e(2^-4/3 rs,0), power=13
jobs+=CONV
for rs in [1.0,10.0]:
    jobs.append((f'id_rs{rs:g}_z1',dict(rs=rs,zeta=1.0,power=13,reps=16,nq=24,order=24,nk=48,seed=30302001)))
    jobs.append((f'id_rs{rs:g}_z0',dict(rs=rs*2**(-4/3),zeta=0.0,power=13,reps=16,nq=24,order=24,nk=48,seed=30302001)))
for name,kw in jobs:
    if (V/f'{name}.json').exists(): print('exists',name); continue
    t=time.time(); v=rn.energy(**kw); d=save(name,v)
    print(f"{name}: {1000*d['mean_ry']:.4f} ± {1000*d['se_ry']:.4f} mRy  ({time.time()-t:.0f}s)",flush=True)
print('AXIS_DONE')
