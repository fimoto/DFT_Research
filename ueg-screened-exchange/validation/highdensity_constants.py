"""Constants of the high-density expansion of the coupling-averaged kite energy (main text Sec. IV, SM S19):
e_kite = e_2x + c_hd r_s log^2 r_s + c_1 r_s log r_s + c_2 r_s + o(r_s).
c_hd and c_1 are computed here from convergent integrals; c_2 is fitted to the direct evaluations in
validation/results/highdensity/ by validation/highdensity_fit.py."""
import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # run from the repository root
import json, numpy as np
from scipy.integrate import quad
from math import pi, log, atan
h=lambda z: 1-z*atan(1/z) if z>0 else 1.0                      # long-wavelength Lindhard shape, nu = q z
def Rz(a,b,z):                                                   # Re[(iz-a)(iz-b)]^{-1}
    d=a*b-z*z; return d/(d*d+z*z*(a+b)**2)
A=lambda z: 8*pi*pi*(2*h(z)-1/(1+z*z))                           # coefficient of log(1/q) in X(q,qz)
def B(z):                                                        # O(1) term of X(q,qz): finite-part angular integral + diagonal terms
    F=lambda cp,ck: cp*ck*Rz(cp,ck,z)
    def inner(cp):
        f=lambda ck: (F(cp,ck)-F(cp,cp))/abs(cp-ck) if ck!=cp else 0.0
        return quad(f,-1,1,points=[cp],limit=200)[0]
    dbl=quad(inner,-1,1,limit=200)[0]
    diag=quad(lambda c: c*c*Rz(c,c,z)*(log(2)+1.5-log(abs(c))) if c!=0 else 0.0,-1,1,points=[0],limit=200)[0]
    return 2*pi*pi*dbl+4*pi*pi*diag
R=lambda u: 2*(u-log(1+u))/u**2 if u>1e-6 else 1-2*u/3+u*u/2   # coupling-averaged screening factor
a=quad(lambda z:h(z)*A(z),0,np.inf,limit=200)[0]
b=quad(lambda z:h(z)*B(z),0,1,limit=100)[0]+quad(lambda z:h(z)*B(z),1,np.inf,limit=100)[0]
J=quad(lambda z:A(z)*h(z)*log(h(z)),0,1,limit=200)[0]+quad(lambda z:A(z)*h(z)*log(h(z)),1,np.inf,limit=200)[0]
Psi1=quad(lambda u:(R(u)-1+2*u/3)/u**2,0,1,limit=200)[0]+quad(lambda u:(R(u)-1)/u**2,1,np.inf,limit=200)[0]
intB=quad(B,0,1,limit=100)[0]+quad(B,1,np.inf,limit=100)[0]    # check: 2 pi^3 (2 log 2 - 1)
alpha=-a/(3*pi); IA=(a*Psi1+2*J/3)/(2*pi); beta=-2*b/(3*pi)+IA
A0=3/(32*pi**5); aL=(4/(9*pi))**(1/3); aTF=4*aL/pi
chd=pi*A0*aTF*alpha; c1=2*pi*A0*aTF*(alpha*log(aTF)-beta)
out=dict(a=a,a_closed=4*pi**3/3*(1-log(2)),b=b,J=J,Psi1=Psi1,Psi1_closed=-2/9,intB=intB,intB_closed=2*pi**3*(2*log(2)-1),
         alpha=alpha,I_A=IA,beta_prime=beta,c_hd_Ry=chd,c_hd_closed_Ry=-aL*(1-log(2))/(6*pi**3),c1_Ry=c1,a_TF=aTF,A0=A0)
json.dump(out,open('validation/results/highdensity/constants.json','w'),indent=1)
for k,v in out.items(): print(f'{k:16s} {v:.12g}')
