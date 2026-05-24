"""Finite-sample behaviour of plain pseudorandom Monte Carlo for the BRM representation.
For each density, many independent runs of N samples per term are drawn with numpy's PCG64
(sampling scale 2, the production importance map). For each run we record the plain mean and two
robust estimates that discard rare large samples. The spread of these estimators relative to the
converged value measures how a finite-sample or outlier-rejecting estimator would bias the energy."""
import os,sys,json,time
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(ROOT); sys.path.insert(0,'production'); sys.path.insert(0,'external')
import numpy as np
import benites_direct as bd
from compute import E2X
rs_list=[0.1,1,2,5,10,30,1e3]
N_list=[2**12,2**16]; runs=64; scale=2.0
rng_seed=20261005
out={'settings':dict(N_list=N_list,runs=runs,scale=scale,generator='numpy PCG64',seed=rng_seed),'rows':{}}
for rs in rs_list:
    t0=time.time(); res={}
    rng=np.random.default_rng([rng_seed,int(1000*rs)])
    for N in N_list:
        plain=[];trim=[]
        for r in range(runs):
            tot_plain=0.0; tot_trim=np.zeros(2)
            for term in (1,2):
                v=bd.evaluate_term(rng.random((N,9)),rs,0.0,term,scale)
                tot_plain+=v.mean()
                a=np.sort(np.abs(v))
                for i,frac in enumerate((1e-3,1e-2)):            # discard the largest |v| fraction
                    k=max(1,int(round(frac*N))); thr=a[-k]
                    tot_trim[i]+=v[np.abs(v)<thr].sum()/N
            plain.append(tot_plain); trim.append(tot_trim)
        plain=np.array(plain); trim=np.array(trim)
        res[str(N)]=dict(mean_of_runs_mry=float(1000*(E2X+plain.mean())),se_of_mean_mry=float(1000*plain.std(ddof=1)/np.sqrt(runs)),
                         median_of_runs_mry=float(1000*(E2X+np.median(plain))),spread_mry=float(1000*plain.std(ddof=1)),
                         trimmed_1e3_mry=float(1000*(E2X+trim[:,0].mean())),trimmed_1e2_mry=float(1000*(E2X+trim[:,1].mean())))
    out['rows'][str(rs)]=res
    print(f'rs={rs:g} ({time.time()-t0:.0f}s): '+' | '.join(f"N=2^{int(np.log2(int(N)))}: mean {d['mean_of_runs_mry']:.2f} median {d['median_of_runs_mry']:.2f} trim0.1% {d['trimmed_1e3_mry']:.2f} trim1% {d['trimmed_1e2_mry']:.2f} spread {d['spread_mry']:.2f}" for N,d in res.items()),flush=True)
    json.dump(out,open('validation/results/tail_test.json','w'),indent=1)
print('TAIL_DONE')
