#!/usr/bin/env python3
"""Reproduce every closed-form coefficient quoted in the manuscript.
Usage: python scripts/verify_analytics.py"""
import os, sys, numpy as np
from scipy.integrate import quad
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from ueg_exchange import constants as C, lorentzian as B, dispersion as T, highdensity as H, lowdensity as LD
from math import gamma

FAILED = []

def line(lab, num, ref, unit="", tol=1e-9):
    good = abs(num - ref) < tol * max(1.0, abs(ref))
    if not good:
        FAILED.append(lab)
    print(f"{'OK ' if good else '!! '}{lab:42s} {num: .12f}  ref {ref: .12f} {unit}")

print("== normalization ==")
line("A0 = 3/(32 pi^5)", C.A0, 2*2*(4*np.pi)**2*(2*np.pi)**-9*(3*np.pi**2)*0.25)
line("e_2x (Onsager)", C.E2X, np.log(2)/3 - 3*C.ZETA3/(2*np.pi**2), "Ry")

print("\n== geometry (JPA Appendix D) ==")
V = lambda r: np.pi*(4+r)*(2-r)**2/12
line("I_0 = int_{B^2} 1/|p-k|^2 = 4 pi^2", 4*np.pi*quad(V,0,2)[0], 4*np.pi**2)
line("int_0^inf C(q) dq = 8 pi^4/15",
     (np.pi/2)*2*np.pi**2*quad(lambda r: r*V(r),0,2)[0], 8*np.pi**4/15)

print("\n== Lorentzian kernel (JPA Lorentzian-kernel analysis) ==")
for Lam in (0.05,0.7,5.0):
    line(f"kappa({Lam})", B.kappa(Lam), B.kappa_quadrature(Lam))
for s_lambda in (0.3,0.5,1.0,3.0):
    ref = quad(lambda u: s_lambda**2*(0.25-u**2)/((u**2+s_lambda**2)*(u**2+0.25)**2),0,np.inf,limit=300)[0]
    line(f"p_loc({s_lambda})", B.p_loc(s_lambda), ref)
line("int_0^inf A(z) dz = 0", quad(B.A_shell,0,np.inf,limit=400)[0], 0.0)

print("\n== dispersion of the screening pole (JPA screening-pole analysis) ==")
line("A1 from the pole dispersion", T.A1_dispersion(), C.A1)
line("A1 closed form", C.A1, -4/(7*np.sqrt(2)*np.pi**1.5)*gamma(0.25)/gamma(0.75))
line("A2", C.A2, 8/15)
line("A3 (three profiles, reduced forms)", LD.A3_from_profiles(), C.A3, tol=1e-8)
cl = LD.A3_pieces_claimed()
line("  A3 cross  (Beta reduction)", LD.A3_cross_closed(), cl["cross"])
line("  A3 recoil (Beta reduction)", LD.A3_recoil_closed(), cl["recoil"])
line("  A3 Lindh. (Euler integral)", LD.A3_lindhard_closed(), cl["lindhard"], tol=1e-8)
line("  A3 cross  (raw profile quad)", LD.A3_cross_numeric(), cl["cross"], tol=1e-8)
line("  A3 recoil (raw profile quad)", LD.A3_recoil_numeric(), cl["recoil"], tol=2e-4)
line("  A3 Lindh. (raw profile quad)", LD.A3_lindhard_numeric(), cl["lindhard"], tol=2e-4)
k = C.q_star(1.0)
line("coef of rs^-1/4", C.A1/k, -0.2213782293, "Ry")
line("coef of rs^-1/2", C.A2/k**2, 0.5670608947, "Ry")
line("coef of rs^-3/4", C.A3/k**3, -0.4575317487, "Ry")

print("\n== high density (JPA Sec. 5.2 and Appendix G) ==")
num, closed = H.calA()
line("calA numeric vs closed", num, closed)
line("int h^2 dz", H.int_h2(), np.pi/3*(1-np.log(2)))
line("C_ln2", H.C_ln2_from_calA(), C.C_LN2, "Ry")
cm, cs = H.response_log_coefficients(0.3, 1.0)
if True:
    line("Pi_vtx^(1)^log + Pi_SE^(1)^log at z=0.3 (exchange cancellation)", cm + cs, 0.0)

D0_ = np.pi**3/18*np.log(2) - np.pi/4*1.2020569031595942854

print("\n== Mellin analysis of the Lorentzian kernel (JPA Sec. 4.4, Appendix F) ==")
from scipy.integrate import quad
from scipy.special import sici
import numpy as _np

def _osc_tail(f, a, period, nper=400):
    """Sum an oscillatory tail over successive half periods (alternating series)."""
    tot = 0.0
    x = a
    for _ in range(nper):
        y = x + period
        tot += quad(f, x, y, limit=200)[0]
        x = y
    return tot

for z in [0.3, 1.0, 5.0]:
    ref = quad(lambda t: z*_np.exp(-z*t)/(1+t*t), 0, _np.inf, limit=400)[0]
    line(f"Xi({z}) closed vs its integral", B.Xi_closed(z), ref, tol=1e-10)

# Xi -> 1 at infinity, so the defining Mellin integral converges only on
# -1 < Re s < 0; elsewhere K(s) is the analytic continuation.  The tail is
# subtracted against 1 to make the check well conditioned.
for s in [-0.9, -0.5, -0.2]:
    A = 60.0
    core = quad(lambda z: z**(s-1)*B.Xi_closed(z), 0, A, limit=800)[0]
    tail = quad(lambda z: z**(s-1)*(B.Xi_closed(z)-1.0), A, _np.inf, limit=800)[0] - A**s/s
    line(f"K({s}) closed vs Mellin integral", B.K_mellin(s), core+tail, tol=1e-3)

# j1(t/2)^2 = [1 - cos t - (t/2) sin t] * 2/t^2 + ... ; rather than fight the
# oscillatory tail, integrate against the exact envelope: the non-oscillatory
# part of j1(x)^2 is 1/(2x^2), so subtract it and integrate the remainder.
_j1 = lambda x: (_np.sin(x)/x - _np.cos(x))/x
def _b0_integrand(t):
    return 2*_j1(t/2)**2*(1-_np.log(t))/t
def _b0_smooth(t):                      # non-oscillatory envelope, 2*(2/t^2)*(1-log t)/t
    return 2*(2.0/t**2)*(1-_np.log(t))/t
_A = 40.0
_b0 = (quad(_b0_integrand, 0, 1, limit=400)[0]
       + quad(_b0_integrand, 1, _A, limit=1600)[0]
       + quad(lambda t: _b0_integrand(t)-_b0_smooth(t), _A, _np.inf, limit=1600)[0]
       + quad(_b0_smooth, _A, _np.inf, limit=800)[0])
line("b0 = gamma/2 - 3/8", B.b0_closed(), _b0, tol=1e-6)
line("C_0^Lambda = -C_1^Lambda/4", B.S_asymptote(50.0),
     1 + (_np.pi/(4*D0_))*(_np.log(50.0) - 0.25)/50.0, tol=1e-14)

for _rs, _k, _Z, _s in ((1.0,0.8372,1.0546,0.0120307),(2.0,1.2166,1.1078,0.0282590),(4.0,1.8142,1.2086,0.0689527)):
    _kk,_ZZ = H.yukawa_pole(_rs)
    line(f"Yukawa pole kappa_* at rs={_rs}", _kk, _k, tol=1e-4)
    line(f"Yukawa residue Z at rs={_rs}", _ZZ, _Z, tol=1e-4)
    line(f"static residual share at rs={_rs}", H.friedel_share(_rs), _s, tol=1e-5)


print("\n== self-energy insertion (exchange cancellation), independent of the identification ==")
for z in [0.0, 0.3, 1.0, 3.0]:
    h = 1 - z*_np.arctan(1/z) if z > 0 else 1.0
    hp = (-_np.arctan(1/z) + z/(1+z*z)) if z > 0 else -_np.pi/2
    line(f"Pi_SE^(1) log slope at z={z}", H.self_energy_response_log_slope(z),
         (h + z*hp)/4, tol=2e-6)

print("\n== spin ==")
for z in (0.0,0.4,1.0):
    line(f"M_3({z}) = 1", C.spin_moment(3,z), 1.0)
for rs in (2.0,5.0,1e4):
    line(f"e_kite({rs},1) = e_kite(2^-4/3 rs,0)",
         C.asymptotic_kite(rs,1.0), C.asymptotic_kite(rs*2**(-4/3),0.0), "Ry")

if FAILED:
    print(f"\n{len(FAILED)} mismatch(es): {FAILED}")
    sys.exit(1)
