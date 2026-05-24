"""Independent imaginary-axis kite estimator. No asymptotic coefficients enter it.
Deterministic q/nu quadrature, independently scrambled Sobol replicates in p,k.
The standard error describes replicate variation, not quadrature or cutoff bias.
"""
import numpy as np
from scipy.integrate import simpson
from scipy.special import roots_legendre
from scipy.stats import qmc

A_NORMAL = 3/(32*np.pi**5)
A_SCREEN = 4/np.pi*(4/(9*np.pi))**(1/3)
V_BALL = 4*np.pi/3

def coupling(c):
    c=np.asarray(c,dtype=float)
    if np.any(c<0): raise ValueError('screening strength must be nonnegative')
    out=np.empty_like(c); small=c<1e-3
    v=c[small]
    out[small]=sum(2*(-v)**j/(j+2) for j in range(9))
    v=c[~small]; out[~small]=2*(v-np.log1p(v))/v**2
    return out

def gauss(a,b,n):
    t,w=roots_legendre(n)
    return (a+b)/2+(b-a)*t/2,(b-a)*w/2

def lindhard(q,nu,nk=48):
    """Eq. (20), direct k integral; log1p avoids high-frequency cancellation."""
    nu=np.atleast_1d(nu)
    edges=sorted(set([0., min(q/2,1.), 1.]))
    ans=np.zeros_like(nu)
    for lo,hi in zip(edges[:-1],edges[1:]):
        k,w=gauss(lo,hi,nk)
        den=(q*k[:,None]-q*q/2)**2+nu[None,:]**2
        ans+=np.sum((w*k)[:,None]*np.log1p(2*q**3*k[:,None]/den),axis=0)/(2*q)
    return ans

def spin_response(q,nu,zeta,nk=48):
    return sum(x*lindhard(q/x,nu/x**2,nk)/2
               for x in [(1+zeta)**(1/3),(1-zeta)**(1/3)] if x>0)

def frequency_grid(q,plasma,order=12):
    """Resolve low particle-hole and plasma scales; inverse-nu tail to infinity."""
    edge=np.unique([0., q*1e-4,q*1e-2,q,q*q,plasma,max(q*q,q,plasma)*10])
    xx=[];ww=[]
    for lo,hi in zip(edge[:-1],edge[1:]):
        if hi>lo:
            x,w=gauss(lo,hi,order);xx.extend(x);ww.extend(w)
    t,w=gauss(0,1,order);tail=edge[-1]
    xx.extend(tail/t);ww.extend(w*tail/t**2)
    return np.asarray(xx),np.asarray(ww)

def sphere(u):
    z=2*u[:,0]-1;ph=2*np.pi*u[:,1];rho=np.sqrt(np.maximum(0,1-z*z))
    return np.column_stack([rho*np.cos(ph),rho*np.sin(ph),z])

def ball(u):
    return sphere(u[:,:2])*u[:,2,None]**(1/3)

def kernel_integral(dp,dk,nu,pref):
    # Samples are chunked by the caller; arrays never span all samples and nodes.
    a=dp[:,None];b=dk[:,None];v=nu[None,:]
    return np.sum(((a*b-v*v)/((a*a+v*v)*(b*b+v*v)))*pref[None,:],axis=1)

def node_estimate(q,rs,u,order=12,nk=48,zeta=0.,spin_scale=1.,bare=False,chunk=2048):
    """Unit loop; c=c_zeta(x_sigma*q,x_sigma^2*nu). All occupation sectors."""
    if not bare:
        nu,w=frequency_grid(q,np.sqrt(A_SCREEN*rs/3)/spin_scale**2,order)
        c=A_SCREEN*rs*spin_response(spin_scale*q,spin_scale**2*nu,zeta,nk)/(spin_scale*q)**2
        pref=coupling(c)*w/np.pi
    def freq(a,b):
        if bare:
            den=2*(a-b)
            return np.divide(np.sign(b)-np.sign(a),den,out=np.zeros_like(a),where=den!=0)
        return kernel_integral(a,b,nu,pref)
    total=0.
    for start in range(0,len(u),chunk):
        h=u[start:start+chunk];p=ball(h[:,:3])
        if q>2:
            # Two same-ball terms, with 1/|p-k|^2 absorbed by uniform radius.
            k=p-2*h[:,5,None]*sphere(h[:,3:5]);inside=np.sum(k*k,axis=1)<1
            a=p[:,2]*q+q*q/2;b=k[:,2]*q+q*q/2
            same=2*V_BALL*8*np.pi*np.where(inside,freq(a,b),0.)
            # Two cross-ball terms; the balls are disjoint in this branch.
            k=ball(h[:,3:6]);k[:,2]-=q
            b=k[:,2]*q+q*q/2
            cross=-2*V_BALL**2*freq(a,b)/np.sum((p-k)**2,axis=1)
            total+=np.sum(same+cross)
        else:
            p[:,2]-=q*(h[:,6]<.5)
            ps=p.copy();ps[:,2]+=q
            np0=np.sum(p*p,axis=1)<1;np1=np.sum(ps*ps,axis=1)<1
            mult=np0.astype(float)+np1.astype(float)
            wp=2*V_BALL/mult
            k=p-(q+2)*h[:,5,None]*sphere(h[:,3:5]);ks=k.copy();ks[:,2]+=q
            dp=np0.astype(float)-np1.astype(float)
            dk=(np.sum(k*k,axis=1)<1).astype(float)-(np.sum(ks*ks,axis=1)<1).astype(float)
            a=p[:,2]*q+q*q/2;b=k[:,2]*q+q*q/2
            total+=np.sum(wp*4*np.pi*(q+2)*dp*dk*freq(a,b))
    return total/len(u)/(16*np.pi)

def energy(rs,zeta=0.,power=13,reps=8,nq=24,order=12,nk=48,qmin=.002,qmax=None,seed=714203,bare=False):
    """nq even subintervals on each side of q=2. Same scrambles at all q nodes.
    Propagate covariance by integrating each entire replicate first.
    """
    if nq%2: raise ValueError('nq must be even for nested Simpson grids')
    if reps<2: raise ValueError('at least two independent scrambles required')
    Q=(4*A_SCREEN*rs/3)**.25
    result=np.zeros(reps);profiles=[]
    for spin_scale in [(1+zeta)**(1/3),(1-zeta)**(1/3)]:
        if spin_scale==0:continue
        # At zeta=0 the two identical loops sum to one; do not compute twice.
        if zeta==0 and profiles:break
        upper=qmax if qmax is not None else 30*max(Q/spin_scale,1)
        logs=np.r_[np.linspace(np.log(qmin),np.log(2),nq+1),np.linspace(np.log(2),np.log(upper),nq+1)[1:]]
        qs=np.exp(logs);vals=np.empty((reps,len(qs)))
        for j in range(reps):
            u=qmc.Sobol(7,scramble=True,seed=seed+j).random_base2(power)
            for k,q in enumerate(qs):
                vals[j,k]=q*node_estimate(q,rs,u,order,nk,zeta,spin_scale,bare)
        integrals=(simpson(vals[:,:nq+1],x=logs[:nq+1],axis=1)+simpson(vals[:,nq:],x=logs[nq:],axis=1))*6/np.pi**3
        result+=(1 if zeta==0 else spin_scale**3/2)*integrals
        profiles.append({'q':qs,'qX':vals,'spin_scale':spin_scale})
    return {'mean':result.mean(),'se':result.std(ddof=1)/np.sqrt(reps),'replicates':result,'profiles':profiles,
            'Q':Q,'settings':dict(rs=rs,zeta=zeta,power=power,reps=reps,nq=nq,order=order,nk=nk,qmin=qmin,qmax=qmax,seed=seed,bare=bare)}
