"""Independent normalization bridge and derivative coefficients."""
import json, math
from pathlib import Path
from scipy.integrate import quad
import sympy as sp
pi=math.pi
rows=[]
for q,nu in [(.7,.4),(3.,2.)]:
    def radial(p):
        if p==0:return 0.
        cmin=max(-1.,(1-p*p-q*q)/(2*p*q))
        if cmin>=1:return 0.
        return p*p*quad(lambda c:(q*q/2+p*q*c)/(nu*nu+(q*q/2+p*q*c)**2),cmin,1,epsabs=1e-12,epsrel=1e-12)[0]
    splits=[max(0.,1-q),1.]
    if 1<q<2:splits=[0.,q-1,1.]
    crescent=sum(quad(radial,lo,hi,epsabs=1e-12,epsrel=1e-12)[0] for lo,hi in zip(splits[:-1],splits[1:]))
    q_tf_sq=.73
    Pi0=-q_tf_sq*crescent
    R=quad(lambda k:k*math.log(((q*k+q*q/2)**2+nu**2)/((q*k-q*q/2)**2+nu**2)),0,1,epsabs=1e-12,epsrel=1e-12)[0]/(2*q)
    q_tf_sq=.73
    rows.append({'q':q,'nu':nu,'Pi0_from_loop':Pi0,'Pi0_expected':-q_tf_sq*R,'relative_error':abs(Pi0+q_tf_sq*R)/abs(Pi0),'dielectric_product_error':abs(Pi0/q**2-(-q_tf_sq*R)/q**2)})
q_tf_sq=sp.symbols('q_tf_sq',positive=True)
normalization=sp.simplify(-3/(sp.pi**3*q_tf_sq**2)*(-q_tf_sq**2/(32*sp.pi**2))-3/(32*sp.pi**5))
assert normalization == 0, normalization
for row in rows:
    assert math.isfinite(row['relative_error']) and row['relative_error'] < 5e-11, row
    assert row['dielectric_product_error'] < 1e-12, row
l=sp.symbols('lambda',positive=True)
# Compare derivative marking λ² only versus the BRM integrand order by order.
mark=[]
for m in range(1,5):
    e_closed=2*l**(m+2)
    full=sp.simplify(sp.diff(e_closed,l)/l**(m+1))
    core=sp.simplify(sp.diff(l**2,l)*(2*l**m)/l**(m+1))
    coupling=sp.integrate(e_closed/l,(l,0,1))
    assert full == 2*(m+2) and core == 4
    assert coupling == sp.Rational(2,m+2)
    mark.append({'m':m,'e_closed_coefficient':2,'full_derivative_coefficient':int(full),'core_only_derivative_coefficient':int(core),'BRM_integrand_coefficient':2,'e_kite_coefficient':float(coupling),'e_kite_over_e_closed':float(coupling/2)})
out={'free_bubble_checks':rows,'K_Pi_vtx_minus_A0':str(normalization),'marking_checks':mark}
out['status']='PASS'
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
