"""High-density end with a q grid resolved on the screening window q <~ sqrt(q_tf_sq)."""
import numpy as np
from math import log,pi,sqrt
from a7_mc import Xq_freq, gl, lindhard_shape, coupling_screening_factor, A0, E2X, aTF
aL=(4/(9*pi))**(1/3); Clog2=-aL*(1-log(2))/(6*pi**3)
def energy_split(rs,breaks,nsegs,nnu=40,N=200_000,seed=1):
    rng=np.random.default_rng(seed); q_tf_sq=aTF*rs
    nodes=[]
    for i,n in enumerate(nsegs):
        lo,hi=breaks[i],breaks[i+1]
        if hi is None:
            t,w=gl(n,0.0,1.0); nodes+=list(zip(lo+2*t/(1-t),w*2/(1-t)**2))
        else:
            qs,w=gl(n,lo,hi); nodes+=list(zip(qs,w))
    tot=0.0
    for q,w in nodes:
        xs=q+q*q/2
        s,ws=gl(nnu,0.0,1.0); nus=xs*s/(1-s); jx=xs/(1-s)**2
        X,_=Xq_freq(q,nus,N,rng)
        c=np.array([q_tf_sq/q**2*lindhard_shape(q,x) for x in nus])
        tot+=w*4*np.pi*A0*(ws*jx*(coupling_screening_factor(c)-1.0)*X).sum()/np.pi
    return E2X+tot
if __name__=='__main__':
    for rs in [0.01,0.1]:
        q_tf_sq=aTF*rs; st=sqrt(q_tf_sq)
        print(f'  rs={rs}, sqrt(q_tf_sq)={st:.4f}   e2x+C2 rs log^2 rs = {1000*(E2X+Clog2*rs*log(rs)**2):.3f} mRy',flush=True)
        for tag,br,ns in [('coarse [0,inf)',(0.0,None),(48,)),
                          ('split at sqrt(q_tf_sq),1',(0.0,st,1.0,None),(24,24,24)),
                          ('finer',(0.0,st/3,st,3*st,1.0,None),(16,16,16,24,24))]:
            vals=[1000*energy_split(rs,br,ns,seed=900+i) for i in range(5)]
            m=np.mean(vals); se=np.std(vals,ddof=1)/np.sqrt(5)
            print(f'     {tag:<26} {m:8.3f} +- {se:5.3f} mRy',flush=True)
