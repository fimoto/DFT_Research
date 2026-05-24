"""Part 7: a completely independent Monte-Carlo evaluation of the kite energy.

Nothing is taken from the authors' code. The estimator is built from the
definition K[f] = A0 int d^3q d^3p d^3k (Delta_p Delta_q)/|p-k|^2
int dnu/pi f Re[1/((i nu - delta_p)(i nu - delta_k))].
"""
import numpy as np, sys, time
rng_master=np.random.default_rng(20260915)

VB=4*np.pi/3
def Voverlap(q): return np.pi*(4+q)*(2-q)**2/12 if q<2 else 0.0
A0=3/(32*np.pi**5)
E2X=np.log(2)/3-3*1.2020569031595943/(2*np.pi**2)
aL=(4/(9*np.pi))**(1/3); aTF=4*aL/np.pi

def sample_geometry(q,N,rng):
    """Sample p uniform on the Delta-support union and d=p-k with density
    proportional to 1/d^2 inside a ball of radius Rmax. Returns weights."""
    VU=2*VB-Voverlap(q)
    Rmax=2.0+q
    # p uniform in union of B and B-q*zhat : sample uniformly in the two balls
    which=rng.integers(0,2,N)
    u=rng.normal(size=(N,3)); u/=np.linalg.norm(u,axis=1)[:,None]
    r=rng.random(N)**(1/3)
    p=u*r[:,None]
    p[which==1,2]-=q                      # second ball centred at -q zhat
    # rejection to make it uniform on the union (overlap double counted)
    inB = (np.einsum('ij,ij->i',p,p)<=1.0)
    pq=p.copy(); pq[:,2]+=q
    inBq=(np.einsum('ij,ij->i',pq,pq)<=1.0)
    keep=~(inB&inBq)                      # Delta=0 on the intersection anyway
    Dp=np.where(inB,1.0,0.0)-np.where(inBq,1.0,0.0)
    # d with radial density uniform in [0,Rmax]
    ud=rng.normal(size=(N,3)); ud/=np.linalg.norm(ud,axis=1)[:,None]
    rd=rng.random(N)*Rmax
    d=ud*rd[:,None]
    k=p-d
    inBk=(np.einsum('ij,ij->i',k,k)<=1.0)
    kq=k.copy(); kq[:,2]+=q
    inBkq=(np.einsum('ij,ij->i',kq,kq)<=1.0)
    Dk=np.where(inBk,1.0,0.0)-np.where(inBkq,1.0,0.0)
    w=Dp*Dk*keep
    a=p[:,2]*q+q*q/2
    b=k[:,2]*q+q*q/2
    # p is drawn from density 1/(2 V_B) on the symmetric difference, so the
    # phase-space factor is 2 V_B (not the union volume) times 4 pi Rmax
    return w, a, b, 2*VB*4*np.pi*Rmax

def Xq_freq(q,nus,N,rng,weights=None):
    """returns X(q,nu) for every nu in nus, one shared geometry sample.
    If frequency weights are given, also returns the per-sample frequency-integrated
    values Y_j (audit F1, 1 Oct 2026): the samples are shared across frequencies, so the
    error of a frequency integral must be formed from Y_j, not from the per-frequency errors."""
    w,a,b,norm=sample_geometry(q,N,rng)
    nz=w!=0
    w=w[nz]; a=a[nz]; b=b[nz]
    out=np.empty(len(nus)); err=np.empty(len(nus))
    Y=np.zeros(len(w)) if weights is not None else None
    for i,nu in enumerate(nus):
        num=a*b-nu*nu
        val=w*num/(num*num+nu*nu*(a+b)**2)
        out[i]=norm*val.sum()/N
        err[i]=norm*np.sqrt(max(( (val**2).sum()/N-(val.sum()/N)**2 ),0)/N)
        if Y is not None: Y+=weights[i]*val
    if Y is None: return out,err
    return out,err,norm*Y

def gl(n,lo,hi):
    x,wt=np.polynomial.legendre.leggauss(n)
    return 0.5*(hi-lo)*x+0.5*(hi+lo), 0.5*(hi-lo)*wt

def lindhard_shape(q,nu):
    """imaginary-axis Lindhard response normalised to lindhard_shape(0,0)=1 (PRB Eq.24)"""
    kk,wk=gl(200,0.0,1.0)
    val=kk*np.log(((q*kk+q*q/2)**2+nu**2)/((q*kk-q*q/2)**2+nu**2))
    return (wk*val).sum()/(2*q)

def coupling_screening_factor(c):
    c=np.asarray(c,dtype=float)
    out=np.empty_like(c)
    small=c<1e-3
    cs=c[small]; out[small]=1-2*cs/3+cs**2/2-2*cs**3/5+cs**4/3
    cb=c[~small]; out[~small]=2*(cb-np.log1p(cb))/cb**2
    return out

def energy(rs,nq=48,nnu=40,N=200_000,qscale=None,bare=False,seed=1):
    rng=np.random.default_rng(seed)
    q_tf_sq=aTF*rs; qstar=(4*q_tf_sq/3)**0.25
    if qscale is None: qscale=max(2.0,qstar)
    t,wt=gl(nq,0.0,1.0)
    qs=qscale*t/(1-t); jq=qscale/(1-t)**2
    tot=0.0; var=0.0
    for q,wq,jj in zip(qs,wt,jq):
        xscale=q+q*q/2
        s,ws=gl(nnu,0.0,1.0)
        nus=xscale*s/(1-s); jx=xscale/(1-s)**2
        if bare:
            f=np.ones_like(nus)
        else:
            c=np.array([q_tf_sq/q**2*lindhard_shape(q,x) for x in nus])
            f=coupling_screening_factor(c)-1.0
        X,Xe,Y=Xq_freq(q,nus,N,rng,weights=ws*jx*f)
        contrib=4*np.pi*A0*(ws*jx*f*X).sum()/np.pi
        # per-node error from the sample variance of the frequency-integrated values Y_j
        # (the previous quadrature sum of per-frequency errors ignored their covariance)
        # Y has one entry per nonzero-weight sample; the N-len(Y) rejected samples contribute zeros
        econtrib=4*np.pi*A0*np.sqrt(max((Y**2).sum()/N-(Y.sum()/N)**2,0.0)/N)/np.pi
        tot+=wq*jj*contrib; var+=(wq*jj*econtrib)**2
    return tot, np.sqrt(var)

if __name__=='__main__':
    t0=time.time()
    v,e=energy(1.0,bare=True,qscale=2.0,N=400_000,seed=7)
    print(f'bare kite  K[v] = {v:.6f} +- {e:.6f} Ry   exact e2x = {E2X:.6f}   dev {(v-E2X)/e:+.2f} sigma  [{time.time()-t0:.0f}s]',flush=True)

def report(rs,**kw):
    import time; t0=time.time()
    d,e=energy(rs,**kw)
    tot=E2X+d
    print(f'rs={rs:<10g} eps_2b = {1000*tot:9.3f} +- {1000*e:6.3f} mRy   '
          f'(screening part {1000*d:9.3f})   [{time.time()-t0:.0f}s]',flush=True)
    return tot,e
