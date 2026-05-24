"""Check coupling-integral identities against direct finite rational kernels.
This is a scalar algebra check, not a new UEG energy calculation.
Requires Python 3 and SciPy. Run: python checks/check_series_transfer.py
"""
import json,math
from pathlib import Path
from scipy.integrate import quad
terms=[(.7,.3),(-.2,2.),(.5,5.)]
e0=sum(a for a,c in terms)
def F(r):return sum(a/(1+c*r) for a,c in terms)
def Fprime(r):return sum(-a*c/(1+c*r)**2 for a,c in terms)
def Q(fun):return quad(fun,0,1,epsabs=2e-12,epsrel=2e-12,limit=300)[0]
def e_kite_quad(r):return 2*Q(lambda lam:lam*F(lam*r))
def e_kite_closed(r):
 def mean(x):
  if x<.02:return sum(2*(-x)**m/(m+2) for m in range(20))
  return 2*(x-math.log1p(x))/x**2
 return sum(a*mean(c*r) for a,c in terms)
rows=[]
for r in [.001,.1,1.,10.,1000.]:
 b=e_kite_quad(r);d=2*Q(lambda lam:lam**2*Fprime(lam*r))
 c=2*F(r)-e0
 root=Q(lambda lam:2*lam*(2*F(lam*r)-e0)+2*lam**2*r*Fprime(lam*r))
 errors=[abs(b-e_kite_closed(r)),abs(b+r*d/2-F(r)),abs(2*b+r*d-e0-c),abs(root-c)]
 assert max(errors)<1e-10,(r,errors)
 rows.append({'r':r,'B_quadrature':b,'F_direct':F(r),'C_direct':c,'max_absolute_error':max(errors)})
moments=[]
for p in [.25,.5,.75]:
 moment=2*Q(lambda lam:lam**(1-p));logmoment=2*Q(lambda lam:lam**(1-p)*math.log(lam))
 assert abs(moment-2/(2-p))<1e-10
 assert abs(logmoment+2/(2-p)**2)<1e-10
 moments.append({'p':p,'moment':moment,'log_moment':logmoment})
result={'purpose':'coupling identities on finite scalar rational kernels, not UEG data','tests':rows,'endpoint_moments':moments}
out=Path(__file__).with_name('series_transfer_results.json');out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'passed':True,'max_error':max(row['max_absolute_error'] for row in rows),'output':str(out)}))
