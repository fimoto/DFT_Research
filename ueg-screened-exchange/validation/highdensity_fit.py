"""Compare the direct kite energies at small r_s with the high-density expansion and fit c_2."""
import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import json, glob, numpy as np
from math import log
C=json.load(open('validation/results/highdensity/constants.json')); chd=C['c_hd_Ry']; c1=C['c1_Ry']
E2X=0.0483583178362888
rows=[]
for f in sorted(glob.glob('validation/results/highdensity/rs*.json')):
    d=json.load(open(f)); rs=d['rs']; e=d['mean_mry']/1000; se=d['se_mry']/1000
    rows.append((rs,e,se))
rows.sort()
print(' r_s      direct-e2x (mRy)   c_hd term   c_1 term   residual/r_s (mRy)')
for rs,e,se in rows:
    de=e-E2X; t2=chd*rs*log(rs)**2; t1=c1*rs*log(rs)
    print(f'{rs:7.3f}   {1000*de:12.4f}   {1000*t2:9.4f}   {1000*t1:9.4f}   {1000*(de-t2-t1)/rs:9.3f} +- {1000*se/rs:.3f}')
sel=[(rs,e,se) for rs,e,se in rows if rs<=0.03]
rs=np.array([r[0] for r in sel]); de=np.array([r[1]-E2X for r in sel]); se=np.array([r[2] for r in sel])
y=(de-chd*rs*np.log(rs)**2-c1*rs*np.log(rs))/rs; w=1/(se/rs)**2
c2=np.sum(w*y)/np.sum(w); print('c_2 (weighted mean, r_s<=0.03): %.5f +- %.5f Ry'%(c2,1/np.sqrt(np.sum(w))))
X=np.vstack([np.log(rs),np.ones_like(rs)]).T; yy=de/rs-chd*np.log(rs)**2; W=np.diag(w)
coef=np.linalg.solve(X.T@W@X,X.T@W@yy); cov=np.linalg.inv(X.T@W@X)
print('free fit: c_1 = %.5f +- %.5f Ry (predicted %.5f), c_2 = %.5f +- %.5f Ry'%(coef[0],np.sqrt(cov[0,0]),c1,coef[1],np.sqrt(cov[1,1])))
json.dump(dict(c2_Ry=float(c2),c2_err_Ry=float(1/np.sqrt(np.sum(w))),c1_fit_Ry=float(coef[0]),c1_fit_err_Ry=float(np.sqrt(cov[0,0]))),open('validation/results/highdensity/fit.json','w'),indent=1)
