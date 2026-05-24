"""High-density coefficients c_hd(zeta) and c_1(zeta) of the kite (main-text Sec. IV, SM Sec. S19)."""
import os as _os, sys as _sys
_os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # run from this directory
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import numpy as np, json
from scipy.integrate import quad
from math import pi, log, atan
import mpmath as mp
import hd_kernel as H
A0=3/(32*pi**5); aL=(4/(9*pi))**(1/3); aTF=4*aL/pi
h=H.h; a=H.a
Psi0=-2/9   # exact: lim[int_eps^inf (R-1)/c^2 dc - (2/3)log eps]
print('Psi0 =',Psi0)
# Gauss-Legendre nodes in z via z=t/(1-t)
tn,wn=np.polynomial.legendre.leggauss(160); tn=0.5*(tn+1); wn=0.5*wn
zn=tn/(1-tn); jac=1/(1-tn)**2
bn=np.array([H.b(z) for z in zn])
def Iz(fun): return float(np.sum(wn*jac*fun(zn)))
hv=np.vectorize(h); av=np.vectorize(a)
Ahat=Iz(lambda z: hv(z)*av(z)); Bhat=Iz(lambda z: hv(z)*bn/1.0) if False else float(np.sum(wn*jac*hv(zn)*bn))
print('A_hat GL =',Ahat,' exact',pi*(1-log(2))/6,'  B_hat GL =',Bhat)
def xs(z): return {+1:(1+z)**(1/3), -1:(1-z)**(1/3)}
def Hz(w,x): return 0.5*sum(x[t]*h(w/x[t]) for t in x if x[t]>0)
out={}
for zeta in (0.0,0.2,0.4,0.6,0.8,0.9,1.0):
    x=xs(zeta); A=B=al=0.
    for s in x:
        if x[s]==0: continue
        w=x[s]**3/2
        g=np.array([Hz(x[s]*z,x)/x[s]**2 for z in zn])
        As=float(np.sum(wn*jac*g*av(zn))); Bs=float(np.sum(wn*jac*g*bn))
        als=Psi0/2*As+float(np.sum(wn*jac*av(zn)*g*np.log(g)))/3
        A+=w*As; B+=w*Bs; al+=w*als
    chd=-(8*pi*pi*A0/3)*A*aTF; c1=2*chd*log(aTF)-2*A0*aTF*(8*pi*pi*al-(2/3)*B)
    out[zeta]=dict(c_hd=chd,c_1=c1,A=A,B=B,alpha=al); print(f'zeta={zeta}: c_hd={chd:.8f}  c_1={c1:.7f}  alpha={al:.6f}',flush=True)
chd0,c10=out[0.0]['c_hd'],out[0.0]['c_1']
print('rescaling check zeta=1: c_hd %.8f  c_1 %.7f'%(2**(-4/3)*chd0, 2**(-4/3)*(c10-(8/3)*log(2)*chd0)))
json.dump(dict(Psi0=Psi0,results=out),open('results/hd_coefficients.json','w'),indent=1)
