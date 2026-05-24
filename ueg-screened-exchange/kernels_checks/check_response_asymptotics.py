#!/usr/bin/env python3
"""Independent quadratures for the self-energy and recoil asymptotics.

Run from any directory: python check_response_asymptotics.py
Dependencies: numpy, scipy, mpmath. Results are written beside this script.
The tests cover finite integrals and asymptotic coefficients. The uniform
remainder estimates are proved in the manuscript, not inferred from fits.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma
import mpmath as mp

mp.mp.dps = 45


def fx_float(x):
    x = float(x)
    if x == 0: return 1.0
    if x == 1: return 0.5
    if x < 1e-3:
        return 1 - sum(x**(2*n)/((2*n-1)*(2*n+1)) for n in range(1,7))
    if x > 100:
        return sum(x**(-2*n)/((2*n-1)*(2*n+1)) for n in range(1,7))
    return .5 + (1-x*x)/(4*x)*np.log(abs((1+x)/(1-x)))


def fx_mp(x):
    x = mp.mpf(x)
    if x == 0: return mp.mpf(1)
    if x == 1: return mp.mpf('.5')
    if x < mp.mpf('1e-5'):
        return 1-mp.fsum(x**(2*n)/((2*n-1)*(2*n+1)) for n in range(1,8))
    return mp.mpf('.5')+(1-x*x)/(4*x)*mp.log(abs((1+x)/(1-x)))


def sigma_1d(width):
    w=mp.mpf(width)
    def dd(a,b):
        return 1/a if a == b else mp.log1p((a-b)/b)/(a-b)
    def integrand(t):
        if t == 0: return mp.mpf(0)
        if t == 1:
            return -mp.mpf('.5') if w == 1 else -w*mp.log(w)/(w*w-1)
        H=dd(1+2*w-t,1+t)+dd(1+2*w+t,1-t)-2*dd(1+2*w+t,1+2*w-t)
        f=fx_mp(t)
        ell=2*mp.atanh(t)
        u=1+2*w*t
        return (1-f)/(t**3*(1+w*t))*(ell-u*2*mp.atanh(t/u))-2*t*t*f*H
    return 3/mp.pi**2*mp.quad(integrand,[0,mp.mpf('.25'),mp.mpf('.75'),1])


def radial_crosscheck():
    rows=[]
    for w in [.5,2.,10.]:
        def outer(r):
            fr=fx_float(r)
            value,_=quad(lambda s:s*(fr-fx_float(s))/((s*s-r*r)*((w+s/2)**2-r*r/4)),
                         1,np.inf,epsabs=2e-12,epsrel=2e-11,limit=200)
            return r*r*value
        value,error=quad(outer,0,1,epsabs=2e-12,epsrel=2e-11,limit=200)
        two=-6*w/np.pi**2*value
        one=float(sigma_1d(str(w)))
        rel=abs(two-one)/abs(one)
        assert rel < 2e-9,(w,one,two,rel)
        rows.append({'width':w,'one_dimensional':one,'double_radial_integral':two,
                     'relative_difference':rel,'outer_quadrature_error':6*w/np.pi**2*error})
    return rows


def finite_part_checks():
    # First integral is evaluated directly; the second is transformed s=1/t
    # to avoid an infinite tail and compared with its independent boundary value.
    i1=mp.quad(lambda r:r*r*fx_mp(r)*mp.log1p(-r*r),[0,.5,1])
    def j(t):
        if t==0:return mp.mpf(0)
        if t==1:return mp.mpf(0) # endpoint value does not affect quadrature
        if t < mp.mpf('1e-8'):
            # The product has a regular leading term t/9.
            return t/9+mp.mpf(4)*t**3/45
        return (1-fx_mp(t))/t**3*(mp.atanh(t)/t-1)
    i2=mp.quad(j,[0,.5,1])
    target1=mp.log(2)/2-mp.mpf(5)/8
    target2=mp.mpf(1)/8
    c=-i1/2-i2
    exactc=mp.mpf(3)/16-mp.log(2)/4
    assert abs(i1-target1)<mp.mpf('1e-35')
    assert abs(i2-target2)<mp.mpf('1e-30')
    # Assemble the nonlogarithmic coefficient without inserting its proposed value.
    coefficient=-6/mp.pi**2*((mp.log(2)-1)/4+c)
    exact=3/(8*mp.pi**2)
    assert abs(coefficient-exact)<mp.mpf('1e-30')
    return {k:mp.nstr(v,35) for k,v in {
        'logarithmic_moment':i1,'logarithmic_moment_error':abs(i1-target1),
        'exterior_moment':i2,'exterior_moment_error':abs(i2-target2),
        'finite_part':c,'finite_part_error':abs(c-exactc),
        'nonlogarithmic_coefficient':coefficient,'coefficient_error':abs(coefficient-exact)}.items()}


def large_width_checks():
    rows=[]
    for w in [100,1000,10000,100000]:
        value=sigma_1d(str(w))
        extracted=w*value+3/(2*mp.pi**2)*mp.log(w)
        error=extracted-3/(8*mp.pi**2)
        rows.append({'width':w,'energy':mp.nstr(value,25),
            'constant_estimate':mp.nstr(extracted,25),
            'constant_error':mp.nstr(error,20),
            'width_times_constant_error':mp.nstr(w*error,20)})
    assert abs(mp.mpf(rows[-1]['constant_error'])) < mp.mpf('2e-6')
    return rows


def local_recoil_checks():
    rows=[]
    for scaled_pole_frequency in [.5,2.,10.]:
        for q in [8.,32.,128.]:
            omega=scaled_pole_frequency*q*q
            def radial(r):
                def angular(t):
                    s=np.sqrt(r*r+q*q+2*r*q*t)
                    return (fx_float(r)-fx_float(s))/(omega+q*q/2+r*q*t)**2
                return r*r*quad(angular,-1,1,epsabs=1e-16,epsrel=2e-11)[0]
            exact=-3*omega/np.pi**2*quad(radial,0,1,epsabs=1e-16,epsrel=2e-11)[0]
            leading=-6*scaled_pole_frequency/(np.pi**2*q*q*(1+2*scaled_pole_frequency)**2)
            scaled=(exact-leading)*q**4
            rows.append({'scaled_pole_frequency':scaled_pole_frequency,'q':q,'local_self_energy':exact,
                         'leading_kernel':leading,'q4_times_error':scaled})
        tail=[r['q4_times_error'] for r in rows[-3:]]
        assert abs(tail[-1]-tail[-2]) < .01*max(abs(tail[-1]),1e-8)
    return rows


def recoil_checks():
    rows=[]
    for scaled_pole_frequency in [.5,1.,3.]:
        direct,_=quad(lambda u:scaled_pole_frequency*scaled_pole_frequency*(.25-u*u)/((u*u+scaled_pole_frequency*scaled_pole_frequency)*(u*u+.25)**2),0,np.inf,
                      epsabs=1e-12,epsrel=1e-12)
        closed=2*np.pi*scaled_pole_frequency/(1+2*scaled_pole_frequency)**2
        assert abs(direct-closed)<1e-11
        rows.append({'scaled_pole_frequency':scaled_pole_frequency,'frequency_integral':direct,'profile':closed,'absolute_difference':abs(direct-closed)})
    # This is the original two-dimensional (lambda,x) integral, without
    # substituting x=lambda**(1/4)*t or using the beta/gamma expression.
    def inner(lam):
        def integrand(x):
            scaled_pole_frequency=np.sqrt(1+lam/x**4)/2
            return lam/(x**4+lam)*scaled_pole_frequency/(1+2*scaled_pole_frequency)**2/x**2
        return quad(integrand,0,np.inf,epsabs=2e-11,epsrel=2e-11,limit=200)[0]
    direct=-12/np.pi**2*quad(lambda lam:lam*inner(lam),0,1,epsabs=2e-11,epsrel=2e-11,limit=200)[0]
    closed=-4/(7*np.sqrt(2)*np.pi**1.5)*gamma(.25)/gamma(.75)
    assert abs(direct-closed)<1e-9
    # A direct Fermi-ball integral supplies the independent normalization
    # of the self-energy large-q kernel; no kite coefficient is inserted here.
    moment=quad(lambda r:r*r*fx_float(r),0,1,epsabs=1e-13)[0]
    sigma_coefficient=-(4*moment)*direct
    assert abs(direct+sigma_coefficient)<1e-12
    return {'frequency_profile':rows,'A1_direct':direct,'A1_closed':closed,
            'A1_absolute_difference':abs(direct-closed),'fermi_ball_radial_moment':moment,
            'self_energy_A1':sigma_coefficient,'sum_of_leading_coefficients':direct+sigma_coefficient}


def run_all():
    return {'precision_decimal_digits':mp.mp.dps,'radial_crosscheck':radial_crosscheck(),
            'finite_part':finite_part_checks(),'large_width':large_width_checks(),
            'local_recoil':local_recoil_checks(),
            'recoil':recoil_checks(), 'status':'PASS'}

if __name__=='__main__':
    result=run_all()
    path=Path(__file__).with_name('response_asymptotics_results.json')
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
