"""Momentum-cutoff test of the difference from the published BRM values (main-text Sec. V C, SM Sec. S22).
All numbers quoted in the manuscript are produced here from results/cutoff_scan.json:
 - per density: published minus direct, its error (published and direct errors in quadrature), and the
   cutoff q_c that reproduces it;
 - a common cutoff fitted by least squares to all densities with published values;
 - standardized residual sums for stated subsets (conditional diagnostics: the published values are
   treated as independent, which their source does not establish);
 - the weighted mean difference at r_s<=50, with and without the covariance of the direct values across
   densities (they share random numbers), the published values again treated as independent;
 - a joint fit of an upper and a lower cutoff.
Writes results/cutoff_scan_summary.json."""
import os,json
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from math import pi
D=json.load(open('validation/results/cutoff_scan.json')); edges=np.array(D['edges'])[1:]
published={0.1:(43.53,0.18),1:(27.62,0.84),2:(20.31,0.33),3:(14.85,0.68),5:(8.04,0.34),7:(2.23,0.41),10:(-2.11,0.41),20:(-11.13,0.90),30:(-14.85,0.28),50:(-19.45,0.65),1e3:(-23.55,0.36),5e3:(-18.44,0.18),2e4:(-14.35,0.22),2e5:(-7.13,0.15),1e6:(-3.89,0.16)}
aL=(4/(9*pi))**(1/3); aTF=4*aL/pi; qstar=lambda rs:(4*aTF*rs/3)**0.25
rs=np.array(sorted(float(k) for k in D['rows'])); key=lambda r:[k for k in D['rows'] if float(k)==r][0]
cum={r:np.array(D['rows'][key(r)]['cum_mean_ry']) for r in rs}
diff=np.array([published[r][0]-D['rows'][key(r)]['total_mry'] for r in rs])
err_pub=np.array([published[r][1] for r in rs]); err_dir=np.array([D['rows'][key(r)]['total_se_mry'] for r in rs]); err=np.hypot(err_pub,err_dir)
def shift_upper(r,qc):   # e(q<q_c) - e  (mRy): what a cutoff q_c would add to the energy
    return 1000*(np.interp(np.log(qc),np.log(edges),cum[r])-cum[r][-1])
def shift_lower(r,qm):   # e(q>q_min) - e
    return -1000*np.interp(np.log(qm),np.log(edges),cum[r])
def implied_qc(r):
    f=np.array([shift_upper(r,q) for q in edges])-diff[list(rs).index(r)]
    idx=np.where(np.sign(f[:-1])!=np.sign(f[1:]))[0]
    if len(idx)==0: return None
    i=idx[-1]; x0,x1=np.log(edges[i]),np.log(edges[i+1]); return float(np.exp(x0+(x1-x0)*f[i]/(f[i]-f[i+1])))
grid=np.logspace(0.3,1.6,600)
chi=lambda pred: np.sum(((diff-pred)/err)**2)
best=min(((chi(np.array([shift_upper(r,q) for r in rs])),q) for q in grid))
qc=best[1]; pred=np.array([shift_upper(r,qc) for r in rs]); pulls=(diff-pred)/err
sub=lambda mask: float(np.sum(pulls[mask]**2))
# weighted mean at r_s<=50, published values independent; direct-side covariance from replicate totals
met=rs<=50; w=1/err[met]**2; mean=float(np.sum(w*diff[met])/np.sum(w)); mean_err=float(1/np.sqrt(np.sum(w)))
T={r:np.array(D['rows'][key(r)]['term_totals_ry']) for r in rs}; R=next(iter(T.values())).shape[0]
idx=[i for i,r in enumerate(rs) if met[i]]; covd=np.zeros((len(idx),)*2)
for a,i in enumerate(idx):
    for b,j in enumerate(idx):
        covd[a,b]=1e6*sum(np.cov(T[rs[i]][:,t],T[rs[j]][:,t],ddof=1)[0,1] for t in (0,1))/R
covtot=np.diag(err_pub[met]**2)+covd; wv=w/np.sum(w); mean_err_cov=float(np.sqrt(wv@covtot@wv))
# joint upper+lower cutoff
bestj=None
for q in np.logspace(0.9,1.3,50):
    for qm in np.logspace(-2.5,0.5,100):
        p=np.array([shift_upper(r,q)+shift_lower(r,qm) for r in rs]); c=chi(p)
        if bestj is None or c<bestj[0]: bestj=(c,q,qm)
rows=[]
print(f"{'rs':>8} {'pub-direct':>12} {'implied q_c':>12} {'pred(q_c=%.2f)'%qc:>16} {'pull':>6}")
for i,r in enumerate(rs):
    q1=implied_qc(r); rows.append(dict(rs=float(r),published_minus_direct_mry=float(diff[i]),err_mry=float(err[i]),implied_qc_kF=q1,implied_qc_over_qstar=None if q1 is None else q1/qstar(r),predicted_mry=float(pred[i]),pull=float(pulls[i])))
    print(f"{r:8g} {diff[i]:+8.2f}+-{err[i]:.2f} {'--' if q1 is None else f'{q1:12.2f}'} {pred[i]:+16.2f} {pulls[i]:+6.1f}")
summary=dict(common_cutoff_kF=float(qc),chi2_all=float(best[0]),n_all=int(len(rs)),
  chi2_rs_ge_1e3=sub(rs>=1e3),n_rs_ge_1e3=int(np.sum(rs>=1e3)),chi2_rs_ge_5e3=sub(rs>=5e3),n_rs_ge_5e3=int(np.sum(rs>=5e3)),chi2_rs_le_50=sub(met),n_rs_le_50=int(met.sum()),
  weighted_mean_rs_le_50=mean,weighted_mean_err_independent=mean_err,weighted_mean_err_with_direct_covariance=mean_err_cov,n_positive_rs_le_50=int(np.sum(diff[met]>0)),
  joint_fit=dict(qc_kF=float(bestj[1]),qmin_kF=float(bestj[2]),chi2=float(bestj[0])),rows=rows,
  note='The published errors are treated as independent and the standardized residual sums are conditional diagnostics, not exact chi-square tests.')
print(f"\ncommon cutoff q_c = {qc:.2f} k_F: sum of squared pulls = {best[0]:.2f} (15 points); subsets: rs>=1e3 {sub(rs>=1e3):.2f} ({np.sum(rs>=1e3)} pts), rs>=5e3 {sub(rs>=5e3):.2f} ({np.sum(rs>=5e3)} pts), rs<=50 {sub(met):.2f} ({met.sum()} pts)")
print(f"weighted mean (rs<=50, {met.sum()} points, {np.sum(diff[met]>0)} positive): {mean:+.3f} +- {mean_err:.3f} mRy (independent) / +- {mean_err_cov:.3f} mRy (direct covariance included)")
print(f"joint upper+lower cutoff: q_c = {bestj[1]:.2f}, q_min = {bestj[2]:.3f} k_F, sum of squared pulls = {bestj[0]:.2f}")
json.dump(summary,open('validation/results/cutoff_scan_summary.json','w'),indent=1)
