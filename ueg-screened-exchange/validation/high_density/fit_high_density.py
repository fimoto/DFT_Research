"""Fit of the direct high-density energies (main-text Sec. IV and V B, Table II, SM Sec. S19).
Reads results/direct_high_density.json (disjoint master seeds, so the densities are statistically independent;
term replicates stored) and, for the windows r_s<=0.1, 0.03 and 0.01, fits c_2 with the analytic c_1 and,
separately, both c_1 and c_2 by weighted least squares. The quoted errors are formal fit errors from the
replicate standard errors; the residual sums are diagnostics of the O(r_s) truncation, not exact chi-square
tests. Writes results/hd_fit.json, which make_fig_high_density.py reads (Fig. 6 uses the r_s<=0.1 window)."""
import os as _os
_os.chdir(_os.path.dirname(_os.path.abspath(__file__)))
import json, numpy as np
from math import log, pi
from scipy.special import zeta
J=json.load(open('results/direct_high_density.json')); D=J['rows']
C=json.load(open('results/hd_coefficients.json'))['results']['0.0']; chd,c1=C['c_hd'],C['c_1']
e2x=log(2)/3-3*zeta(3)/(2*pi**2)
rs_all=np.array(sorted(float(k) for k in D)); key=lambda r:[k for k in D if float(k)==r][0]
e=np.array([D[key(r)]['mean_mry'] for r in rs_all])/1000; s=np.array([D[key(r)]['se_mry'] for r in rs_all])/1000
# seed disjointness (the run script asserts it; re-check from the stored masters and replicate counts)
used=[set(D[key(r)]['master_seed']+2*i+t for i in range(D[key(r)]['replicates']) for t in (0,1)) for r in rs_all]
assert all(not (used[i]&used[j]) for i in range(len(used)) for j in range(i+1,len(used))), 'seed overlap'
out=dict(e2x_ry=e2x,c_hd=chd,c_1=c1,seed_rule=J['settings']['seed_rule'],rows=[],windows={})
for r,ei,si in zip(rs_all,e,s):
    out['rows'].append(dict(rs=float(r),replicates=D[key(r)]['replicates'],direct_mry=1000*ei,se_mry=1000*si,minus_e2x=1000*(ei-e2x),remainder_over_rs=(ei-e2x-chd*r*log(r)**2-c1*r*log(r))/r,rem_err=si/r))
print('seeds disjoint across densities; r_s, replicates, (e-e2x-logs)/r_s [Ry]:')
for row in out['rows']: print(f"  {row['rs']:<6g} {row['replicates']:>3d}  {row['remainder_over_rs']:+.5f} +- {row['rem_err']:.5f}")
for rmax in (0.1,0.03,0.01):
    m=rs_all<=rmax; rs=rs_all[m]; y=e[m]-e2x-chd*rs*np.log(rs)**2; w=1/s[m]**2
    yy=y-c1*rs*np.log(rs); c2=np.sum(w*yy*rs)/np.sum(w*rs**2); c2e=1/np.sqrt(np.sum(w*rs**2)); Q=np.sum(w*(yy-c2*rs)**2)
    A=np.vstack([rs*np.log(rs),rs]).T*np.sqrt(w)[:,None]; coef,_,_,_=np.linalg.lstsq(A,y*np.sqrt(w),rcond=None); covp=np.linalg.inv(A.T@A); Qf=float(np.sum((A@coef-y*np.sqrt(w))**2))
    out['windows'][str(rmax)]=dict(n_points=int(m.sum()),c2=float(c2),c2_err=float(c2e),Q_c1_fixed=float(Q),free_c1=float(coef[0]),free_c1_err=float(np.sqrt(covp[0,0])),free_c2=float(coef[1]),free_c2_err=float(np.sqrt(covp[1,1])),Q_free=Qf,free_c1_pull=float((coef[0]-c1)/np.sqrt(covp[0,0])))
    print(f'r_s<={rmax} ({m.sum()} pts): c_2 = {c2:+.5f} +- {c2e:.5f} Ry, Q = {Q:.2f} ({m.sum()-1} dof) | free c_1 = {coef[0]:.5f} +- {np.sqrt(covp[0,0]):.5f} ({(coef[0]-c1)/np.sqrt(covp[0,0]):+.1f} s.e. from analytic), c_2 = {coef[1]:+.5f} +- {np.sqrt(covp[1,1]):.5f}, Q = {Qf:.2f}')
w=out['windows']['0.1']; out.update(c2=w['c2'],c2_err=w['c2_err'],fit_window_rs_max=0.1)
for row in out['rows']:
    r=row['rs']; row['expansion_mry']=1000*(e2x+chd*r*log(r)**2+c1*r*log(r)+out['c2']*r)
json.dump(out,open('results/hd_fit.json','w'),indent=1); print('wrote results/hd_fit.json')
