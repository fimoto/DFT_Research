"""Direct two-dimensional quadrature of the kinematic coefficient A_{3,kin} (main-text Appendix A),
A_{3,kin} = 4 pi A_0 (104 pi^2/27) int_0^inf dx x^{-4} Phi_4(x^{-4}),
Phi_4(g) = int_0^inf du/pi R(g/(1+4u^2)) Re (iu-1/2)^{-4},
compared with the closed form -416/675 C_Gamma."""
import mpmath as mp
mp.mp.dps=15
pi=mp.pi; A0=3/(32*pi**5)
R=lambda c: 2*(c-mp.log(1+c))/c**2 if c>1e-6 else 1-2*c/3+c*c/2
def Phi4(g):
    return mp.quad(lambda u: R(g/(1+4*u*u))*mp.re((1j*u-mp.mpf(1)/2)**-4),[0,0.25,0.5,1,3,mp.inf])/pi
I=mp.quad(lambda x: x**-4*Phi4(x**-4),[0,0.3,0.6,1,2,5,mp.inf])
A3kin=4*pi*A0*(104*pi**2/27)*I
CG=mp.sqrt(2*pi)*mp.gamma(mp.mpf(3)/4)/(pi**2*mp.gamma(mp.mpf(1)/4)); exact=-mp.mpf(416)/675*CG
print('A_3,kin quadrature = %s   closed form = %s   relative difference = %s'%(mp.nstr(A3kin,12),mp.nstr(exact,12),mp.nstr((A3kin-exact)/exact,3)))
