"""Skewness of plain pseudorandom Monte Carlo estimates of the BRM representation.
Stores every run value, so that the probability that a single run lies above the converged energy,
and the median shift, can be quantified with their errors."""
import os,sys,json,time
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(ROOT); sys.path.insert(0,'production'); sys.path.insert(0,'external')
import numpy as np
import benites_direct as bd
from compute import E2X
out={'settings':dict(runs=256,scale=2.0,generator='numpy PCG64',seed=20261006),'rows':{}}
for rs in (0.1,2,5,10,30,1e3):
    rng=np.random.default_rng([20261006,int(1000*rs)]); t0=time.time(); res={}
    for N in (2**10,2**12,2**14):
        vals=[]
        for r in range(256):
            vals.append(sum(bd.evaluate_term(rng.random((N,9)),rs,0.0,term,2.0).mean() for term in (1,2)))
        res[str(N)]=list(1000*(E2X+np.array(vals)))
    out['rows'][str(rs)]=res
    print(f'rs={rs:g} done ({time.time()-t0:.0f}s)',flush=True)
    json.dump(out,open('validation/results/skew_test.json','w'))
print('SKEW_DONE')
