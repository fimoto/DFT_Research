"""Independent mathematical audit of the supplied PRB manuscript.
Domain: real u; s=sqrt(1+m)>1; w in (0,1); q,Q,kappa>0.
Run: python analytic_audit.py. Requires SymPy only. All checks use exact arithmetic.
Contour closure, endpoint domination and integration interchange are justified
in mathematical_audit_ja.md; this script checks their algebraic certificates.
"""

def run_analytic_audit(record=None):
    import sympy as S
    j=S.I
    u,s,w,z,k,r,t,g=S.symbols('u s w z k r t g', positive=True)
    checks=[]
    def eq(name,a,b=0):
        diff=S.factor(S.cancel(a-b))
        if diff!=0: diff=S.simplify(diff)
        assert diff==0,(name,diff)
        checks.append(name)
        record(name, diff == 0) if record else print('EXACT PASS:',name)
    B=1+4*u*u; A=j*u-S.Rational(1,2); m=s*s-1
    us=j*s/2
    # Poles are determined from the unsimplified frequency integrands, not their proposed results.
    f2=B/(B+m)/A**2
    f4=B/(B+m)/A**4
    fc=4/(B+m)
    fL=-S.Rational(16,15)*(12*u*u-1)/(B*(B+m)**2*A**2)
    def simple(f,p): return S.cancel(f*(u-p)).subs(u,p)
    def double(f,p): return S.diff(S.cancel(f*(u-p)**2),u).subs(u,p)
    h2=S.factor(j*simple(f2,us));h4=S.factor(j*simple(f4,us))
    eq('screening residue: same sphere',h2,-m/(s*(s+1)**2))
    eq('screening residue: kinematic',h4,-4*m/(s*(s+1)**4))
    eq('screening residue: cross sphere',j*simple(fc,us),1/s)
    l0=S.factor(j*simple(fL,j/2)); ls=S.factor(j*double(fL,us))
    eq('Lindhard simple pole',l0,16/(15*m*m))
    eq('Lindhard double pole',ls,-8*(9*s**4-6*s**3+8*s*s-2*s-1)/(15*s**3*(s-1)**2*(s+1)**4))
    hL=S.factor(l0+ls)
    eq('Lindhard sum before coalescent limit',hL,8*(2*s**3-s*s+4*s+1)/(15*s**3*(s+1)**4))
    eq('coalescent m=0 limit',S.limit(hL,s,1),S.Rational(1,5))
    # Coupling average: t=sqrt(1+lambda*g). Integrate rational t functions.
    def avg(h,weight):
        return S.factor(S.integrate(4*t*(t*t-1)**weight*h.subs(s,t),(t,1,s)))
    F=S.factor(avg(h2,1)/(s*s-1)**2)
    JC=S.factor(avg(1/s,1)/(s*s-1)**2)
    Phi=S.factor(avg(h4,1)/(s*s-1)**2)
    Psi=S.factor(-avg(hL,2)/(s*s-1)**3)
    eq('F coupling integral',F,-S.Rational(4,3)*(s-1)/(s+1)**2)
    eq('JC coupling integral',JC,S.Rational(4,3)*(s+2)/(s+1)**2)
    print('Phi4(s)=',Phi,'\nPsi(s)=',Psi)
    # Endpoint checks on derived expressions, with s=sqrt(1+x^-4).
    x=S.symbols('x',positive=True)
    sx=S.sqrt(1+x**4)/x**2
    for name,expr,want in [('F origin',F.subs(s,sx)/x**2,-S.Rational(4,3)),('JC origin',JC.subs(s,sx)/x**2,S.Rational(4,3)),('Phi origin',Phi.subs(s,sx)/x**6,-16),('Psi origin',Psi.subs(s,sx)/x**8,-S.Rational(32,15))]:
        eq(name,S.limit(expr,x,0),want)
    # Independent logarithm expansion, including the next omitted coefficient.
    logratio=S.log((u*u+(S.Rational(1,2)+z)**2)/(u*u+(S.Rational(1,2)-z)**2))
    coef1=S.diff(logratio,z).subs(z,0)
    coef3=S.diff(logratio,z,3).subs(z,0)/6
    coef5=S.diff(logratio,z,5).subs(z,0)/120
    eq('Lindhard q^-2',coef1/6,4/(3*B))
    eq('Lindhard q^-4',coef3/10,16*(1-12*u*u)/(15*B**3))
    print('Lindhard q^-6 coefficient=',S.factor(coef5/14))
    # Geometric moments by explicit lens integration, independent of notebook constants.
    V=S.pi*(4+r)*(2-r)**2/12
    rho2=1-(z+r/2)**2
    a=1-r/2
    Jlens=2*S.pi*S.integrate(z*z*rho2+rho2**2/2,(z,0,a))
    eq('sphere Coulomb moment I0',4*S.pi*S.integrate(V,(r,0,2)),4*S.pi**2)
    eq('sphere quadratic moment I2',4*S.pi*S.integrate(Jlens+r*r*V/12,(r,0,2)),52*S.pi**2/27)
    eq('integrated overlap C',S.pi**3*S.integrate(r*V,(r,0,2)),8*S.pi**4/15)
    # Euler reduction certificate: differentiation of a rational function times a fixed weight.
    ld=-1/(2*w)-3/(4*(1-w))-13/(4*(1+w))
    r1=2*w*(1+w)*(11*w+21)/45
    r2=2*w*(1+w)*(7*w-3)/45
    r3=w*(1+w)*(15*w*w-8*w-3)/15
    for idx,(rr,pp,kk) in enumerate([(r1,1,S.Rational(16,15)),(r2,w,S.Rational(2,15)),(r3,4*w*w+w**3,S.Rational(1,5))]):
        eq('Euler derivative certificate '+str(idx+1),S.diff(rr,w)+ld*rr,pp-kk*(1+w)**2/(2*(1-w)))
    rr=S.factor(2*r1-r2+r3)
    eq('combined Euler certificate',S.diff(rr,w)+ld*rr,2-w+4*w*w+w**3-S.Rational(11,10)*(1+w)**2/(1-w))
    print('combined Euler r(w)=',rr)
    eq('M rational factor',S.Rational(16,15)*S.Rational(11,10)*2,S.Rational(176,75))
    # Finite-part integral: an integration-by-parts certificate gives a convergent integral.
    f=(s-1)**(-S.Rational(1,4))*(s+1)**(-S.Rational(1,4))
    h=(s-1)**S.Rational(3,4)*(s+1)**(-S.Rational(1,4))
    eq('finite part integration-by-parts certificate',S.diff(h,s),f/2+f/(2*(s+1)))
    eq('finite part transformed cutoff has no constant',S.limit(2*h.subs(s,S.sqrt(1+x**-4))-2/x,x,0),0)
    # Hence Fp integral f = -integral f/(s+1), all at s=1..infty;
    # convergent Beta B(3/4,1/2)/sqrt(2).
    G=S.gamma
    Cg=S.sqrt(2*S.pi)/S.pi**2*G(S.Rational(3,4))/G(S.Rational(1,4))
    def beta(a,b):return G(a)*G(b)/G(a+b)
    fp=S.Rational(2,3)*(-beta(S.Rational(3,4),S.Rational(1,2))/S.sqrt(2)-beta(S.Rational(3,4),S.Rational(3,2))/(2*S.sqrt(2)))
    eq('finite part coefficient via convergent Beta integrals',S.expand_func(fp),-S.Rational(8,5)*S.pi**2*Cg)
    A0=3/(32*S.pi**5)
    kin=4*S.pi*A0*104*S.pi**2/27*S.Rational(2,5)*(-2*S.sqrt(2)*beta(S.Rational(7,4),S.Rational(3,2)))
    eq('A3kin Beta to Gamma',S.expand_func(kin),-S.Rational(416,675)*Cg)
    eq('A3 decomposition',S.Rational(32,15)+S.Rational(416,675)+S.Rational(264,125),S.Rational(16408,3375))
    eq('S3 M6 weight',(S.Rational(32,15)+S.Rational(416,675))/S.Rational(16408,3375),S.Rational(1160,2051))
    eq('S3 M4M5 weight',S.Rational(264,125)/S.Rational(16408,3375),S.Rational(891,2051))
    # Remainder envelope integrated explicitly; Q>M.
    Q,M=S.symbols('Q M',positive=True)
    eq('integrated large-q error envelope',S.integrate(1/(Q**4*z**2),(z,M,Q))+S.integrate(z**-6,(z,Q,S.oo)),1/(M*Q**4)-4/(5*Q**5))
    # Cross energy shifts include the occupation sign. The linear term is odd in d_z.
    pz,kz,eps,kap,qq,lam,cc=S.symbols('p_z k_z eps kappa q lambda c',real=True)
    Ap=j*u+S.Rational(1,2)
    cross=-S.re(1/((A-eps*pz)*(Ap-eps*kz)))
    eq('occupation-weighted cross denominator',cross.subs(eps,0),4/B)
    eq('cross linear energy shift',S.diff(cross,eps).subs(eps,0),-8*(pz-kz)/B**2)
    # Derivative inequalities have nonnegative integrands for c>=0, lambda in [0,1].
    rl=2*lam/(1+lam*cc)
    eq('R-c|Rprime| positive integrand',rl+cc*S.diff(rl,cc),2*lam/(1+lam*cc)**2)
    eq('2R-c^2 Rsecond positive integrand',2*rl-cc**2*S.diff(rl,cc,2),4*lam*(1+2*lam*cc)/(1+lam*cc)**3)
    eq('Phi weak-screening coefficient',S.limit(Phi/(s*s-1),s,1),-S.Rational(1,6))
    eq('Psi zero-screening value',S.limit(Psi,s,1),-S.Rational(2,15))
    # Static q>2 integral: geometric expansion in d_z/q and d^2/q^2.
    dz,dsq=S.symbols('d_z d_squared',real=True)
    stat=1/((1+2*dz*eps+dsq*eps**2)*(1+dz*eps))
    eq('static leading integrand',stat.subs(eps,0),1)
    eq('static odd correction vanishes by inversion',S.diff(stat,eps).subs(eps,0),-3*dz)
    eq('static C_infty from two cross sectors',(4*S.pi/3)**2/(8*S.pi),2*S.pi/9)
    fi=S.Rational(1,2)+(kap/4+1/kap)*S.atan(kap/2)
    eq('static imaginary-pole uniqueness certificate',kap*S.diff(fi,kap)-2*fi,-S.Rational(1,2)-(kap/4+3/kap)*S.atan(kap/2))
    eq('static imaginary-pole asymptotic scale',S.limit(fi/kap,kap,S.oo),S.pi/8)
    eq('static imaginary-pole residue asymptote',S.limit(1/(1-kap*S.diff(fi,kap)/(2*fi)),kap,S.oo),2)
    print('ALL',len(checks),'EXACT CHECKS PASSED')
    return checks

if __name__ == "__main__":
    run_analytic_audit()
