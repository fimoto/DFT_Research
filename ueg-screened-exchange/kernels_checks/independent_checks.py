"""Independent checks written from manuscript definitions, no project-code imports."""
import json, math
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import quad
import mpmath as mp
mp.mp.dps=45
pi=math.pi
out={}
# 1. Direct frequency quadrature against each branch of J.
rows=[]
for omega,a,b in [(0.4,0.3,2),(3,-0.2,-1.5),(1.2,-.7,1.9),(.5,.8,.8)]:
    val=quad(lambda x:(a*b-x*x)/((a*a+x*x)*(b*b+x*x)*(omega*omega+x*x)),0,np.inf,epsabs=2e-12,epsrel=2e-12)[0]/pi
    al,be=abs(a),abs(b)
    expected=1/(2*omega*(al+omega)*(be+omega)) if a*b>0 else -(al+be+2*omega)/(2*omega*(al+be)*(al+omega)*(be+omega))
    rows.append(dict(omega=omega,a=a,b=b,quad=val,closed=expected,absdiff=abs(val-expected)))
out['frequency_J']=rows
# 2. Evaluate the low-density coefficient directly in t, independent beta expression.
f=lambda t: (1/(1+t**4))* (mp.sqrt(1+t**-4)/2)/(1+mp.sqrt(1+t**-4))**2/t**2
A1=-48/(7*mp.pi**2)*mp.quad(f,[0,.1,1,10,mp.inf])
closed=-4/(7*mp.sqrt(2)*mp.pi**mp.mpf('1.5'))*mp.gamma(mp.mpf('.25'))/mp.gamma(mp.mpf('.75'))
out['A1']={'direct':str(A1),'closed':str(closed),'absdiff':str(abs(A1-closed))}
# 3. Independent elementary finite-part identities.
def fxmp(r):
    if r==0:return mp.mpf(1)
    if r==1:return mp.mpf('.5')
    return mp.mpf('.5')+(1-r*r)/(4*r)*mp.log(abs((1+r)/(1-r)))
moment=mp.quad(lambda r:r*r*fxmp(r),[0,.5,1])
logmoment=mp.quad(lambda r:r*r*fxmp(r)*mp.log(1-r*r),[0,.5,1])
# Reciprocal variable avoids infinite-large-argument subtraction in f_x.
second=mp.quad(lambda t:((1-fxmp(t))/t**3)*(mp.atanh(t)/t-1),[0,.1,.5,1])
out['selfenergy_finite_part']={'moment':str(moment),'logmoment':str(logmoment),'logmoment_closed':str(mp.log(2)/2-mp.mpf(5)/8),'second':str(second),'finite_part':str(-logmoment/2-second),'finite_part_closed':str(mp.mpf(3)/16-mp.log(2)/4)}
# 4. Original crescent coordinates (q,p,cos theta), versus radial integrated expression.
def gauss(n,a=0,b=1):
    x,w=roots_legendre(n);return a+(b-a)*(x+1)/2,w*(b-a)/2
def fx(x):
    x=np.asarray(x);safe=np.where(x==0,1,x)
    val=.5+(1-x*x)/(4*safe)*np.log(np.abs((1+x)/(1-x)))
    # use asymptotic series where raw subtraction loses digits.
    t=1/np.maximum(x,1)
    tail=t*t/3+t**4/15+t**6/35+t**8/63+t**10/99+t**12/143+t**14/195
    return np.where(x>16,tail,np.where(x<1e-7,1-x*x/3-x**4/15,val))
def sigma_original(width,n):
    t,wt=gauss(n); c,wc=gauss(n,-1,1)
    total=0.
    # Split q at 1,2 and 6 for endpoint/control-domain changes.
    for qa,qb in [(0,1),(1,2),(2,6)]:
        qs,wqs=gauss(n,qa,qb)
        for q,wq in zip(qs,wqs):
            splits=[max(0,1-q),1]
            if 1<q<2:splits=[0,q-1,1]
            inner=0.
            for lo,hi in zip(splits[:-1],splits[1:]):
                ps=lo+(hi-lo)*t;wp=(hi-lo)*wt
                cmin=np.maximum(-1,(1-ps*ps-q*q)/(2*ps*q))
                cs=cmin[:,None]+(1-cmin[:,None])*(c[None,:]+1)/2
                ws=(1-cmin[:,None])*wc[None,:]/2
                ss=np.sqrt(ps[:,None]**2+q*q+2*ps[:,None]*q*cs)
                delta=q*q/2+q*ps[:,None]*cs;Omega=width*q
                inner+=2*pi*np.sum(wp[:,None]*ps[:,None]**2*ws*(fx(ps)[:,None]-fx(ss))/(Omega+delta)**2)
            total+=wq*(-3*width*q/(2*pi**3))*inner
    # q=6/t tail.
    for q,wq in zip(6/t,6*wt/t**2):
        ps=t;wp=wt;cs=c[None,:]
        ss=np.sqrt(ps[:,None]**2+q*q+2*ps[:,None]*q*cs)
        delta=q*q/2+q*ps[:,None]*cs;Omega=width*q
        inner=2*pi*np.sum(wp[:,None]*ps[:,None]**2*wc[None,:]*(fx(ps)[:,None]-fx(ss))/(Omega+delta)**2)
        total+=wq*(-3*Omega/(2*pi**3))*inner
    return total

def sigma_radial(width,n):
    r,wr=gauss(n);t,wt=gauss(n);s=1/t;ws=wt/t**2
    rr=r[:,None];ss=s[None,:]
    v=rr*rr*ss*(fx(r)[:,None]-fx(s)[None,:])/((ss*ss-rr*rr)*((width+ss/2)**2-rr*rr/4))
    return -6*width/pi**2*np.sum(wr[:,None]*ws[None,:]*v)
out['sigma_original_vs_radial']=[]
for width in [.5,2,10]:
    row={'width':width}
    for n in [48,80,128]:
        row[str(n)]={'original':sigma_original(width,n),'radial':sigma_radial(width,n)}
    out['sigma_original_vs_radial'].append(row)
    print('sigma',row,flush=True)
# 5. Fixed-domain static kernel quadrature, no series inserted.
def static_fixed(q,n):
    ts,ws=gauss(n,-1,1);a,wa=gauss(n);uu=a[:,None];rr=a[None,:]
    total=0.
    for t,wt in zip(ts,ws):
        Lt=q*t/2+math.sqrt(1-q*q*(1-t*t)/4)
        for s,w in zip(ts,ws):
            Ls=q*s/2+math.sqrt(1-q*q*(1-s*s)/4)
            den=Lt+rr*Ls
            H=1+q*uu*Lt*Ls*rr*(t+s)/den+q*q/4*((Lt*t-rr*Ls*s)**2/den**2-1)
            f=uu*Lt*Lt*Ls*Ls*rr/den**2/np.sqrt(H)
            total+=wt*w*np.sum(wa[:,None]*wa[None,:]*f)
    return pi/8*total
c0=pi/8*(2*math.log(2)-1);c2=pi/1152-pi*math.log(2)/96
out['static_endpoint']={'c0':c0,'c2':c2,'checks':[]}
for q in [.1,.05,.025]:
    val=static_fixed(q,32)
    out['static_endpoint']['checks'].append({'q':q,'Xc_over_q':val,'estimated_c2':(val-c0)/q**2,'fourth_order_residual':(val-c0-c2*q*q)/q**4})
# 6. Shell high-density integral independent direct h quadrature.
def h(z):
    if z>100:return 1/(3*z*z)-1/(5*z**4)+1/(7*z**6)
    return 1-z*math.atan(1/z)
I=quad(lambda z:h(z)**2,0,np.inf,epsabs=2e-13)[0]
out['shell_integral']={'integral_h2':I,'closed':pi/3*(1-math.log(2)),'coefficient_B':-(4/(9*pi))**(1/3)*(1-math.log(2))/(6*pi**3)}
# Numerical acceptance checks turn the recorded comparisons into executable
# verification. These tolerances test the sampled finite integrals; they do
# not certify the manuscript's uniform remainder estimates.
for row in out['frequency_J']:
    assert math.isfinite(row['absdiff']) and row['absdiff'] < 2e-12, row
assert mp.mpf(out['A1']['absdiff']) < mp.mpf('1e-35'), out['A1']
fp=out['selfenergy_finite_part']
assert abs(mp.mpf(fp['moment'])-mp.mpf(1)/4) < mp.mpf('1e-35'), fp
assert abs(mp.mpf(fp['logmoment'])-mp.mpf(fp['logmoment_closed'])) < mp.mpf('1e-35'), fp
assert abs(mp.mpf(fp['second'])-mp.mpf(1)/8) < mp.mpf('1e-24'), fp
assert abs(mp.mpf(fp['finite_part'])-mp.mpf(fp['finite_part_closed'])) < mp.mpf('1e-24'), fp
for row in out['sigma_original_vs_radial']:
    errors=[]
    for n in ['48','80','128']:
        vals=row[n]
        assert all(math.isfinite(v) and v < 0 for v in vals.values()), row
        errors.append(abs(vals['original']-vals['radial']))
    assert errors[2] < errors[1] < errors[0] and errors[2] < 2e-8, (row,errors)
endpoint_errors=[abs(row['estimated_c2']-c2) for row in out['static_endpoint']['checks']]
assert endpoint_errors[2] < endpoint_errors[1] < endpoint_errors[0], endpoint_errors
assert endpoint_errors[2] < 2e-6, endpoint_errors
assert all(abs(row['fourth_order_residual']) < .002 for row in out['static_endpoint']['checks'])
assert abs(out['shell_integral']['integral_h2']-out['shell_integral']['closed']) < 2e-12
out['status']='PASS'
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
