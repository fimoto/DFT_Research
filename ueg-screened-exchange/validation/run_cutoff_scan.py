"""Momentum-cutoff test of the difference from the published BRM values.
Estimator A (production settings m=19, 16 replicates, sampling scale 2, seed 20261004) with the
per-sample contributions of both terms accumulated in logarithmic bins of q, so that
e(q<q_c) is available for every cutoff q_c from one run per density."""
import os,sys,json,time
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(ROOT); sys.path.insert(0,'production'); sys.path.insert(0,'external')
import numpy as np
from scipy.stats import qmc
import benites_direct as bd
from compute import E2X
edges=np.concatenate([[0.],np.logspace(-3,4,281)]); nb=len(edges)-1
rs_list=[0.1,1,2,3,5,7,10,20,30,50,1e3,5e3,2e4,2e5,1e6]
m,reps,scale,seed,chunk=19,16,2.0,20261004,15
out={'edges':edges.tolist(),'settings':dict(m=m,replicates=reps,q_scale=scale,seed=seed),'e2x_ry':E2X,'rows':{}}
for rs in rs_list:
    t0=time.time(); hist=np.zeros((reps,2,nb))
    for rep in range(reps):
        for term in (1,2):
            sob=qmc.Sobol(9,scramble=True,seed=seed+2*rep+term-1)
            for _ in range(2**(m-chunk)):
                pts=sob.random(2**chunk); u=bd._clip_unit(pts); qv,_=bd._q_map(u[:,1],scale)
                vals=bd.evaluate_term(pts,rs,0.0,term,scale)
                hist[rep,term-1]+=np.histogram(qv,bins=edges,weights=vals)[0]
    hist/=2**m
    cum=np.cumsum(hist.sum(axis=1),axis=1)            # replicate x bin : correction integrated up to edge
    mean=cum.mean(axis=0); se=cum.std(axis=0,ddof=1)/np.sqrt(reps)
    out['rows'][str(rs)]=dict(cum_mean_ry=mean.tolist(),cum_se_ry=se.tolist(),total_mry=float(1000*(E2X+mean[-1])),total_se_mry=float(1000*se[-1]),
        term_totals_ry=hist.sum(axis=2).tolist())   # replicate x term totals, for the covariance across densities
    print(f'rs={rs:g}: total {1000*(E2X+mean[-1]):.4f} +- {1000*se[-1]:.4f} mRy ({time.time()-t0:.0f}s)',flush=True)
    json.dump(out,open('validation/results/cutoff_scan.json','w'))
print('CUTOFF_DONE')
