"""Estimator A at high density (main-text Table II, Fig. 6, SM Sec. S19): m=19, sampling scale 2.
compute.direct uses the 2R seeds seed+2r+t (replicate r=0..R-1, term t=0,1) at each call (128 seeds for R=64,
32 for R=16), so the master seeds of the densities are spaced by 10000 and checked for disjointness here. r_s<=0.1 use 64 replicates,
r_s=0.3 and 1 use 16. The term replicates are stored."""
import sys,json,time,os
ROOT=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); os.chdir(ROOT); sys.path.insert(0,'production')
from compute import direct
densities=(0.001,0.003,0.01,0.03,0.1,0.3,1.0); reps={0.3:16,1.0:16}
master=[61000000+10000*i for i in range(len(densities))]
used=[set(m+2*r+t for r in range(reps.get(rs,64)) for t in (0,1)) for m,rs in zip(master,densities)]
for i in range(len(used)):
    for j in range(i+1,len(used)):
        assert not (used[i]&used[j]), ('seed overlap between densities',densities[i],densities[j])
out={'settings':dict(m=19,q_scale=2.0,replicates={str(rs):reps.get(rs,64) for rs in densities},master_seeds=dict(zip(map(str,densities),master)),seed_rule='compute.direct draws seeds master+2*replicate+term; masters spaced by 10000; disjointness asserted'),'rows':{}}
for rs,m in zip(densities,master):
    t=time.time(); R=reps.get(rs,64); r=direct(rs,0.0,dict(m=19,replicates=R,q_scale=2.0),dict(seed=m,chunk_power=15))
    out['rows'][str(rs)]=dict(master_seed=m,replicates=R,mean_mry=r['mean_mry'],se_mry=r['se_mry'],term_replicates_ry=r['term_replicates_ry'])
    print(f'rs={rs}: {r["mean_mry"]:.4f} +- {r["se_mry"]:.4f} mRy ({R} replicates, master seed {m}, {time.time()-t:.0f}s)',flush=True)
    json.dump(out,open('validation/high_density/results/direct_high_density.json','w'),indent=1)
print('HD_DONE')
