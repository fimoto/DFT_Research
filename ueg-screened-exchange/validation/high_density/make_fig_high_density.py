"""Figure of main-text Sec. IV: high-density kite energy, direct evaluations against the expansion.
Reads results/direct_high_density.json, results/hd_coefficients.json, results/hd_fit.json; writes
../../figs/fig_highdensity.pdf and the plotted series to results/fig_highdensity_*.dat."""
import os as _os
_os.chdir(_os.path.dirname(_os.path.abspath(__file__)))
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from math import log, pi
from scipy.special import zeta
e2x=log(2)/3-3*zeta(3)/(2*pi**2)
D=json.load(open('results/direct_high_density.json'))['rows']; C=json.load(open('results/hd_coefficients.json'))['results']['0.0']; F=json.load(open('results/hd_fit.json'))
chd,c1,c2=C['c_hd'],C['c_1'],F['c2']
key=lambda r:[k for k in D if float(k)==r][0]
rs=np.array(sorted(float(k) for k in D)); e=np.array([D[key(r)]['mean_mry'] for r in rs])/1000; s=np.array([D[key(r)]['se_mry'] for r in rs])/1000
two=lambda r: chd*r*np.log(r)**2+c1*r*np.log(r)
three=lambda r: two(r)+c2*r
rr=np.logspace(-3.3,0.1,300)
fit=rs<=0.1
fig,ax=plt.subplots(2,1,figsize=(3.4,4.3),sharex=True,gridspec_kw=dict(height_ratios=[1.25,1]))
a=ax[0]; a.set_xscale('log'); a.set_yscale('log')
a.plot(rr,-1000*three(rr),'k-',lw=1,label=r'two logarithms $+\,c_2 r_s$')
a.plot(rr,-1000*two(rr),'r--',lw=1,label='two logarithms')
a.errorbar(rs[fit],1000*(e2x-e[fit]),yerr=1000*s[fit],fmt='o',ms=3.5,mfc='none',color='C0',label='direct evaluation')
a.errorbar(rs[~fit],1000*(e2x-e[~fit]),yerr=1000*s[~fit],fmt='o',ms=3.5,color='C0',label='direct, outside the fit range')
a.set_ylabel(r'$e_{2x}-e_{\rm kite}$ (mRy)'); a.set_title('(a) High-density energy',fontsize=9); a.legend(fontsize=6.5,loc='lower right',frameon=False)
b=ax[1]; b.axhspan(-1,1,color='0.85',lw=0); b.axhline(0,color='k',lw=0.6)
b.plot(rs,(e-e2x-two(rs))/s,'r--',marker='s',ms=3.5,mfc='none',lw=1,label='vs two logarithms')
b.plot(rs,(e-e2x-three(rs))/s,'C0',marker='o',ms=3.5,mfc='none',lw=1,label=r'vs two logarithms $+\,c_2 r_s$')
b.set_yscale('symlog',linthresh=5); b.set_yticks([-1000,-100,-10,0,10]); b.set_xlabel(r'$r_s$'); b.set_ylabel(r'(direct $-$ model)/$\sigma$'); b.set_title('(b) Deviations in statistical-error units',fontsize=9); b.legend(fontsize=6.5,loc='upper right',frameon=False)
b.set_xticks([1e-3,1e-2,1e-1,1])
fig.tight_layout(); fig.savefig('../../figs/fig_highdensity.pdf'); print('wrote figs/fig_highdensity.pdf')
np.savetxt('results/fig_highdensity_points.dat',np.c_[rs,1000*(e2x-e),1000*s,(e-e2x-two(rs))/s,(e-e2x-three(rs))/s],header='rs  e2x-e_kite_mRy  error_mRy  (direct-two_logs)/sigma  (direct-two_logs_plus_c2)/sigma')
np.savetxt('results/fig_highdensity_curves.dat',np.c_[rr,-1000*two(rr),-1000*three(rr)],header='rs  e2x-two_logs_mRy  e2x-two_logs_plus_c2_mRy')
