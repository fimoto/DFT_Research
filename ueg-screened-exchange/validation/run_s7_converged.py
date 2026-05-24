"""Recompute SM Sec. S7 (third implementation, hd_split.energy_split) on converged momentum grids.
Grid G0 = S7 as used: breaks (0, min(st,.5), max(2q*,2), inf), nodes (20,36,28), 40 frequency nodes.
Grid G1 = G0 plus a segment boundary at q=2 (kink of X at the ball-overlap threshold):
          breaks (0, min(st,.5), 2, max(2q*,3), inf), nodes (20,24,24,28), 56 frequency nodes.
Grid G2 = G1 with every node count doubled.  One JSON record per replicate (resumable)."""
import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # run from the repository root
import sys,json,time,hashlib
from pathlib import Path
from math import sqrt
sys.path.insert(0,'independent_audit/grid')
from hd_split import energy_split
from a7_mc import aTF
OUT=Path('validation/results/s7_converged.json')
rec=json.load(open(OUT)) if OUT.exists() else {'description':__doc__,'rows':[],
     'sources':{f:hashlib.sha256(Path('independent_audit/grid/'+f).read_bytes()).hexdigest() for f in ('hd_split.py','a7_mc.py')}}
def grid(rs,g):
    st=sqrt(aTF*rs); qs=max(1.0,(4*aTF*rs/3)**0.25)
    if g in ('G0','G0N250'): return (0.0,min(st,0.5),max(2*qs,2.0),None),(20,36,28),40
    br=(0.0,min(st,0.5),2.0,max(2*qs,3.0),None); ns=(20,24,24,28)
    if g=='G1': return br,ns,56
    if g=='G1nu80': return br,ns,80
    if g=='G2': return br,tuple(2*n for n in ns),56
N=100_000
jobs=[(30.0,g,i) for g in ('G0','G1','G2') for i in range(4)]+[(30.0,'G1nu80',i) for i in range(4)] \
    +[(5.0,g,i) for g in ('G1','G2') for i in range(4)] \
    +[(30.0,'G0N250',i) for i in range(2)] \
    +[(rs,'G1',i) for rs in (0.1,1.0,2.0,3.0,7.0,10.0,20.0,50.0) for i in range(4)]
for rs,g,i in jobs:
    if any(r['rs']==rs and r['grid']==g and r['seed']==3000+i for r in rec['rows']): continue
    br,ns,nnu=grid(rs,g); t=time.time(); NN=250_000 if g=='G0N250' else N
    e=1000*energy_split(rs,br,ns,nnu=nnu,N=NN,seed=3000+i)
    rec['rows'].append(dict(rs=rs,grid=g,seed=3000+i,N=NN,breaks=[b if b is not None else 'inf' for b in br],nodes=list(ns),nnu=nnu,energy_mry=e,wall_s=time.time()-t))
    tmp=OUT.with_suffix('.tmp'); tmp.write_text(json.dumps(rec,indent=1)); tmp.replace(OUT)
    print(f'rs={rs} {g} seed={3000+i}: {e:.4f} mRy ({time.time()-t:.0f}s)',flush=True)
print('S7_DONE')
