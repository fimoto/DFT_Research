"""Validation driver (1 Oct 2026): direct BRM integration at arbitrary (rs, zeta) with
input checks, paired replicates stored per point. Uses production.direct unchanged."""
import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # run from the repository root
from pathlib import Path
import argparse, hashlib, json, platform, sys, time
import numpy as np, scipy
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'production'))
from compute import direct
def run(rs_list,zeta_list,m,reps,seed,scale,out,chunk_power=15):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    paths=['production/compute.py','production/common.py','production/runtime.py','external/benites_direct.py',
           'src/ueg_exchange/lindhard.py','src/ueg_exchange/constants.py']
    meta=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
              sources={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in paths},
              driver=hashlib.sha256(Path(__file__).resolve().read_bytes()).hexdigest(),date=time.strftime('%Y-%m-%d %H:%M'))
    settings=dict(m=m,replicates=reps,q_scale=scale); cfg=dict(seed=seed,chunk_power=min(chunk_power,m))
    for rs in rs_list:
        for zeta in zeta_list:
            dest=out/f"rs_{rs:.17g}_z_{zeta:.17g}.json"
            if dest.exists(): print('exists',dest); continue
            t=time.time(); row=direct(float(rs),float(zeta),settings,cfg)
            if not np.isfinite(row['term_replicates_ry']).all(): raise FloatingPointError('nonfinite replicate')
            row['validation_metadata']=meta; row['wall_s']=time.time()-t
            tmp=dest.with_suffix('.tmp'); tmp.write_text(json.dumps(row,indent=1,allow_nan=False)); tmp.replace(dest)
            print(f"rs={rs} zeta={zeta} m={m} R={reps} scale={scale}: {row['mean_mry']:.4f} +- {row['se_mry']:.4f} mRy ({row['wall_s']:.0f}s)",flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--rs',type=float,nargs='+',required=True); p.add_argument('--zeta',type=float,nargs='+',required=True)
    p.add_argument('--m',type=int,default=17); p.add_argument('--reps',type=int,default=16); p.add_argument('--seed',type=int,default=20261001)
    p.add_argument('--scale',default='2'); p.add_argument('--out',required=True)
    a=p.parse_args()
    if not 8<=a.m<=30 or a.reps<2: p.error('require 8<=m<=30, reps>=2')
    if not np.isfinite(a.rs).all() or min(a.rs)<=0: p.error('rs finite positive')
    if max(abs(z) for z in a.zeta)>1: p.error('|zeta|<=1')
    scale='qstar' if a.scale=='qstar' else float(a.scale)
    run(a.rs,a.zeta,a.m,a.reps,a.seed,scale,a.out)
