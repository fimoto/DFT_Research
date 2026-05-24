"""Long-wavelength kernel of the kite: A(z)=8 pi^2 a(z) and the finite-part B(z)=b(z) of main-text Eq. (shell-expansion), SM Eq. (S-AB)."""
import numpy as np
from scipy.integrate import quad
from math import pi, log, atan, sqrt
h=lambda z: 1.0-z*atan(1.0/z) if z>0 else 1.0
a=lambda z: 2*h(z)-1/(1+z*z)
def Rc(c):
    if c<1e-6: return 1-2*c/3+c*c/2
    return 2*(c-log(1+c))/(c*c)
def Rdiag(c,z): return (c*c-z*z)/(c*c+z*z)**2
def Rpair(cp,ck,z):
    ab=cp*ck; return (ab-z*z)/((ab-z*z)**2+z*z*(cp+ck)**2)
def B2(z):   # 4 pi^2 int_{-1}^{1} c^2 R(c,c;z)[log2+3/2-log|c|] dc  = 8 pi^2 int_0^1 ...
    f=lambda c: c*c*Rdiag(c,z)*(log(2)+1.5-log(c))
    return 8*pi*pi*quad(f,0,1,points=[z] if z<1 else None,limit=200)[0]
def B1(z):   # 2 pi^2 int int [G(cp,ck)-G(cp,cp)]/|cp-ck|
    def inner(cp):
        G0=cp*cp*Rpair(cp,cp,z)
        g=lambda ck: (cp*ck*Rpair(cp,ck,z)-G0)/abs(cp-ck) if ck!=cp else 0.0
        return quad(g,-1,1,points=[cp],limit=200)[0]
    return 2*pi*pi*quad(inner,-1,1,limit=200)[0]
def b(z): return B1(z)+B2(z)
