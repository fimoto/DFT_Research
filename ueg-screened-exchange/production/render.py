#!/usr/bin/env python3
"""Generate every numerical PRB figure, table row and data-dependent claim.
No integral is evaluated here. All empirical inputs come from production/raw.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import runtime
if __name__ == "__main__": runtime.enter_cli()
import numpy as np
from scipy.interpolate import RectBivariateSpline, PchipInterpolator
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
matplotlib.rcdefaults()
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, ScalarFormatter, NullFormatter
from common import ROOT, config, source_hashes, save, stats, weights, key, q_star, A1, A2, A3, E2X, asymptotic_kite, check_reference, check_source_compatibility


def records(dataset):
    cfg = config()[dataset]; rows = []
    for rs in cfg['rs']:
        for z in cfg.get('zeta', [0.]):
            path = ROOT / 'production/raw' / dataset / (key(rs,z) + '.json')
            row = json.loads(path.read_text())
            check_source_compatibility(row['metadata']['sources'])
            check_reference(path.relative_to(ROOT / "production/raw"), row)
            rows.append(row)
    return rows


def direct_summary(row):
    a = np.asarray(row['term_replicates_ry'])
    return dict(rs=row['rs'], zeta=row['zeta'], mean_mry=float(1000*(E2X+a.mean(axis=0).sum())),
                se_mry=float(1000*np.sqrt(a.var(axis=0,ddof=1).sum()/len(a))))


def interpolate_zero(rs,z,E):
    spl=RectBivariateSpline(np.log(rs),z,E,kx=3,ky=3,s=0)
    f=lambda r: float(spl(np.log(r),1)[0,0]-spl(np.log(r),0)[0,0])
    changes=[(a,b) for a,b in zip(rs[:-1],rs[1:]) if f(a)*f(b)<0]
    if len(changes)!=1: raise ValueError('Expected exactly one resolved full-polarization crossing')
    return brentq(f,*changes[0]), spl


def summarize():
    cfg=config(); out={name:[direct_summary(r) for r in records(name)]
                          for name in ['tableI','tableII','tableIII_direct','tableIV']}
    out['reference']=json.loads((ROOT/'production/reference_brm.json').read_text())
    shells=[]
    for row in records('shell'):
        a=np.asarray(row['density_replicates_ry']); b=np.asarray(row['refined_density_replicates_ry'])
        mean=a.mean(axis=0);se=a.std(axis=0,ddof=1)/np.sqrt(len(a))
        shells.append(dict(rs=row['rs'],q=row['q'],x=row['x'],mean_ry=mean.tolist(),se_ry=se.tolist(),
                           frequency_shift_ry=(b.mean(axis=0)-mean).tolist(),
                           shown=(-mean>3*se).tolist()))
    out['shell']=shells
    grid=[]
    for row in records('grid'):
        q=np.asarray(row['q']);cor=np.asarray(row['correction_X_replicates']);fine=np.asarray(row['refined_correction_X_replicates']);bare=np.asarray(row['bare_X_replicates'])
        cols={}
        for name in ['published','extended']:
            nodes=np.loadtxt(ROOT/f'data/appendixB_grid_{name}.dat');idx=np.searchsorted(q,nodes)
            if not np.array_equal(q[idx],nodes): raise ValueError('Grid is not nested exactly')
            w=weights(nodes);C=1000*6/np.pi**3
            vals=1000*E2X+C*(cor[:,idx]@w); refined=1000*E2X+C*(fine[:,idx]@w)
            bare_error=C*(bare[:,idx]@w)-1000*E2X
            mu,se=stats(vals);bmu,bse=stats(bare_error)
            cols[name]=dict(mean_mry=mu,se_mry=se,replicates_mry=vals.tolist(),
                            frequency_shift_mry=float(np.mean(refined-vals)),
                            bare_discrepancy_mry=bmu,bare_se_mry=bse)
        grid.append(dict(rs=row['rs'],**cols))
    out['grid']=grid
    rows=records('map');rs=np.array(cfg['map']['rs']);z=np.array(cfg['map']['zeta']);R=cfg['map']['replicates']
    E=np.empty((len(rs),len(z)));se=np.empty_like(E);reps=np.empty((R,len(rs),len(z)))
    for row in rows:
        i=list(rs).index(row['rs']);j=list(z).index(row['zeta']);s=direct_summary(row)
        E[i,j]=s['mean_mry'];se[i,j]=s['se_mry'];reps[:,i,j]=1000*(E2X+np.asarray(row['term_replicates_ry']).sum(axis=1))
    crossing,spl=interpolate_zero(rs,z,E)
    jack=[interpolate_zero(rs,z,np.delete(reps,k,axis=0).mean(axis=0))[0] for k in range(R)]
    jse=float(np.sqrt((R-1)/R*np.sum((np.asarray(jack)-np.mean(jack))**2)))
    pp=PchipInterpolator(np.log(rs),E[:,-1]-E[:,0]);pzero=float(np.exp(brentq(pp,np.log(rs[0]),np.log(rs[-1]))))
    out['map']=dict(rs=rs.tolist(),zeta=z.tolist(),mean_mry=E.tolist(),se_mry=se.tolist(),crossing=crossing,
                    crossing_jackknife_se=jse,jackknife_crossings=jack,pchip_crossing=pzero)
    return out


def pm(mean,se,digits=3): return rf'${mean:.{digits}f}\pm{se:.{digits}f}$'


def density_label(rs):
    if rs<1000:return f'{rs:g}'
    exponent=int(np.floor(np.log10(rs)));mantissa=rs/10**exponent
    factor='' if mantissa==1 else f'{mantissa:g}'+r'\times'
    return '$'+factor+'10^{'+str(exponent)+'}$'


def write_tables(s,dest):
    p=dest/'spin/generated';p.mkdir(parents=True,exist_ok=True)
    refs={row[0]:row[1:] for row in s['reference']['tableI']}
    lines=[]
    for row in s['tableI']:
        rs=row['rs'];b,e=refs[rs]
        brm='---' if b is None else pm(b,e,2)
        lines.append(f"{density_label(rs)} & {pm(row['mean_mry'],row['se_mry'])} & ${1000*asymptotic_kite(rs):.3f}$ & {brm}\\\\")
    (p/'tableI.tex').write_text('\n'.join(lines)+'\n')
    refs={row[0]:row[1:] for row in s['reference']['tableII']};lines=[]
    for row in s['tableII']:
        b,e=refs[row['rs']]
        lines.append(f"{row['rs']:g} & {pm(row['mean_mry'],row['se_mry'])} & {pm(b,e,2)} & ${b-row['mean_mry']:+.2f}$\\\\")
    (p/'tableII.tex').write_text('\n'.join(lines)+'\n');lines=[]
    for row,direct in zip(s['grid'],s['tableIII_direct']):
        a,b=row['published'],row['extended'];brm,err=refs[row['rs']]
        lines.append(f"{row['rs']:g} & {pm(a['mean_mry'],a['se_mry'])} & {pm(b['mean_mry'],b['se_mry'])} & {pm(direct['mean_mry'],direct['se_mry'])} & {pm(brm,err,2)}\\\\")
    (p/'tableIII.tex').write_text('\n'.join(lines)+'\n')
    (p/'tableIV.tex').write_text('\n'.join(f"{r['rs']:g} & {pm(r['mean_mry'],r['se_mry'],2)}\\\\" for r in s['tableIV'])+'\n')
    headers={'tableI':('rccc',r'$r_s$ & direct BRM evaluation & three-term expansion & published BRM'),
             'tableII':('rrrr',r'$r_s$ & direct & published & difference'),
             'tableIII':('crrrr',r'$r_s$ & quad., $[0.006,10]$ & quad., $[10^{-4},60]$ & direct & published BRM'),
             'tableIV':('rr',r'$r_s$ & $e_{\rm kite}(r_s,0)$')}
    for name,(alignment,header) in headers.items():
        f=p/(name+'.tex');body=f.read_text()
        f.write_text('\\begin{tabular}{'+alignment+'}\n\\toprule\n'+header+'\\\\\n\\midrule\n'+body+'\\bottomrule\n\\end{tabular}\n')
    strong=s['tableI'];sr=[(r['mean_mry']-1000*asymptotic_kite(r['rs']))/r['se_mry'] for r in strong]
    residuals=[refs[r['rs']][0]-r['mean_mry'] for r in s['tableII']]
    m=s['map'];mac={'StrongFirstSigma':f'{sr[0]:.1f}','StrongSecondSigma':f'{sr[1]:.1f}',
      'StrongLateMaxSigma':f'{max(abs(v) for v in sr[2:]):.1f}',
      'StrongFirstDifference':f"{strong[0]['mean_mry']-1000*asymptotic_kite(strong[0]['rs']):.2f}",
      'PublishedResidualMean':f'{np.mean(residuals):+.2f}','PublishedResidualSD':f'{np.std(residuals,ddof=1):.2f}',
      'PositiveResidualCount':str(sum(v>0 for v in residuals)),
      'PolarizationCrossing':f"{m['crossing']:.1f}",'PolarizationCrossingSE':f"{m['crossing_jackknife_se']:.1f}",
      'PchipCrossing':f"{m['pchip_crossing']:.1f}",
      'MinimumGridRs':f"{min(s['tableIV'],key=lambda r:r['mean_mry'])['rs']:g}",
      'GridMaxDifference':f"{max(abs(r['extended']['mean_mry']-d['mean_mry']) for r,d in zip(s['grid'],s['tableIII_direct'])):.3f}",
      'GridFrequencyShift':f"{max(abs(r[g]['frequency_shift_mry']) for r in s['grid'] for g in ['published','extended']):.4f}"}
    for i,r in enumerate(strong[:3]):
        Q=q_star(r['rs']);val=Q**3*(r['mean_mry']/1000-A1/Q-A2/Q**2);err=Q**3*r['se_mry']/1000
        mac['ResidualThree'+['One','Two','Three'][i]]=f'{val:.3f}\\pm{err:.3f}'
    for name,suffix in [('published','Coarse'),('extended','Fine')]:
        v=s['grid'][0][name]
        mac['BareGrid'+suffix]=f"{v['bare_discrepancy_mry']:+.3f}\\pm{v['bare_se_mry']:.3f}"
    shift=max(abs(r[g]['frequency_shift_mry']) for r in s['grid'] for g in ['published','extended'])
    mantissa, exponent=f'{shift:.2e}'.split('e')
    mac['GridFrequencyShift']=mantissa+r'\times10^{'+str(int(exponent))+'}'
    (p/'numbers.tex').write_text('% Generated by production/render.py; do not edit.\n'+''.join(f'\\newcommand{{\\{k}}}{{{v}}}\n' for k,v in mac.items()))
    save(dest/'production/generated/summary.json',s)
    for name in ['tableI','tableII','tableIII_direct','tableIV']:
        (dest/f'production/generated/{name}.csv').write_text('rs,zeta,energy_mRy,standard_error_mRy\n'+''.join(f"{r['rs']:.17g},{r['zeta']:.17g},{r['mean_mry']:.17g},{r['se_mry']:.17g}\n" for r in s[name]))
    # Keep historical filenames synchronized for external plotting clients.
    data=dest/'figs/data';data.mkdir(parents=True,exist_ok=True)
    m=s['map'];matrix=[]
    for i,rs in enumerate(m['rs']):
        for j,z in enumerate(m['zeta']):matrix.append([rs,z,m['mean_mry'][i][j],m['se_mry'][i][j]])
    np.savetxt(data/'map.dat',matrix,fmt='%.17g',header='rs zeta energy_mRy standard_error_mRy; generated from production/raw/map')
    matrix=[]
    for row in s['shell']:
        Q=q_star(row['rs'])
        for x,q,y,e in zip(row['x'],row['q'],row['mean_ry'],row['se_ry']):
            u=np.sqrt(1+x**-4);inner=3/np.pi**2*(-4/3*(u-1)/(u+1)**2)/x**2
            matrix.append([row['rs'],x,q,Q**2*y,Q**2*e,inner])
    np.savetxt(data/'shell.dat',matrix,fmt='%.17g',header='rs x q qstar^2*density_Ry scaled_SE_Ry analytic_profile; generated from production/raw/shell')
    old=json.loads((ROOT/'data/tables.json').read_text())
    legacy={k:v for k,v in old.items() if not k.startswith(('companion_','prb_')) and k!='_note'}
    legacy['_note']='PRB entries are generated from production/raw. JPA entries belong to the separate companion project.'
    br1={r[0]:r[1:] for r in s['reference']['tableI']};br2={r[0]:r[1:] for r in s['reference']['tableII']}
    legacy['companion_tableI_strong_screening']={'columns':['rs','direct','direct_err','three_term','brm_tableIX','brm_err'],
      'rows':[[r['rs'],r['mean_mry'],r['se_mry'],1000*asymptotic_kite(r['rs']),*br1[r['rs']]] for r in s['tableI']]}
    legacy['companion_tableII_residuals']={'columns':['rs','direct','direct_err','brm','brm_err'],
      'rows':[[r['rs'],r['mean_mry'],r['se_mry'],*br2[r['rs']]] for r in s['tableII']]}
    legacy['companion_secVC_ferro']={'columns':['rs','e2b_unpolarized','err'],
      'rows':[[r['rs'],r['mean_mry'],r['se_mry']] for r in s['tableIV']]}
    legacy['companion_appendixB_direct']={'rows':[{'rs':g['rs'],'quadrature':g,'direct':d} for g,d in zip(s['grid'],s['tableIII_direct'])]}
    for k,v in legacy.items():
        if k.startswith('companion_'):v['provenance']='production/config.json, production/raw and production/render.py; energies mRy'
    save(dest/'data/tables.json',legacy)


def savefig(fig,path,date=None):
    date=date or datetime.datetime(2026,9,14,tzinfo=datetime.timezone.utc)
    fig.savefig(path,metadata={'CreationDate':date,'ModDate':date},bbox_inches='tight',pad_inches=.04)
    plt.close(fig)


def series(dest,name,columns,header,fmt='%.17g',blocks=None):
    """Write the plotted series of a PRB figure as whitespace-separated ASCII.

    columns: list of (label, 1-D array). The header names the figure panel and every
    column, so the figure can be redrawn in any plotting program."""
    data=dest/'figs/data';data.mkdir(parents=True,exist_ok=True)
    labels=[c[0] for c in columns];arr=np.column_stack([np.asarray(c[1],dtype=float) for c in columns])
    head=header+'\ncolumns: '+'  '.join(f'{i+1}:{l}' for i,l in enumerate(labels))+'\ngenerated by production/render.py from production/raw; nan marks an undefined value'
    if blocks is None:
        np.savetxt(data/name,arr,fmt=fmt,header=head);return
    # one block per value of `blocks`, separated by a blank line, so that plotting
    # programs that join points (gnuplot, xmgrace) start a new line for each series
    head+='\none block per series, separated by a blank line'
    lines=['# '+h for h in head.split('\n')];prev=None
    for key,row in zip(blocks,arr):
        if prev is not None and key!=prev: lines.append('')
        lines.append(' '.join(fmt%v for v in row));prev=key
    (data/name).write_text('\n'.join(lines)+'\n')


def plots(s,dest):
    folder=dest/'spin/figs';folder.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':9,'axes.labelsize':9,'legend.fontsize':7.5,'axes.linewidth':.7,'lines.linewidth':1.1,'figure.dpi':180})
    profile=lambda x:3/np.pi**2*(-4/3*(np.sqrt(1+x**-4)-1)/(np.sqrt(1+x**-4)+1)**2)/x**2
    fig,axs=plt.subplots(3,1,figsize=(3.4,6.1))
    for i,row in enumerate(s['shell']):
        Q=q_star(row['rs']);x=np.array(row['x']);q=np.array(row['q']);y=np.array(row['mean_ry']);e=np.array(row['se_ry']);good=np.array(row['shown'])
        label=rf"$r_s=10^{{{int(np.log10(row['rs']))}}}$";fmt=['o','s','^'][i];c=f'C{i}'
        axs[0].loglog(q[good],-y[good],fmt+'-',ms=3,mfc='w',color=c,label=label)
        axs[0].axvline(Q,color=c,ls=':',lw=.7)
        axs[1].loglog(x[good],-Q**2*y[good],fmt,ms=3,mfc='w',color=c,label=label)
        rel=Q**2*y/profile(x)-1;axs[2].errorbar(x[good],Q**2*rel[good],yerr=Q**4*e[good]/np.abs(profile(x[good])),fmt=fmt+'-',ms=3,mfc='w',color=c,lw=1,capsize=1.5,label=label)
    xx=np.geomspace(.25,6,300);axs[1].loglog(xx,-profile(xx),color='.2',label='leading profile')
    cols={k:[] for k in ('rs','qstar','x','q','neg_density','neg_scaled','resid','shown')}
    for row in s['shell']:
        Q=q_star(row['rs']);x=np.array(row['x']);y=np.array(row['mean_ry']);e=np.array(row['se_ry'])
        cols['rs']+=[row['rs']]*len(x);cols['qstar']+=[Q]*len(x);cols['x']+=list(x);cols['q']+=list(row['q'])
        cols['neg_density']+=list(-y);cols['neg_scaled']+=list(-Q**2*y);cols['resid']+=list(Q**2*(Q**2*y/profile(x)-1))
        cols['shown']+=[float(g) for g in row['shown']]
    series(dest,'spin_fig2_shell_points.dat',[('rs',cols['rs']),('q_star',cols['qstar']),('x=q/q_star',cols['x']),('q_kF',cols['q']),
        ('top:-e(q;rs)_Ry',cols['neg_density']),('middle:-q_star^2*e(q;rs)_Ry',cols['neg_scaled']),
        ('bottom:q_star^2*(data/profile-1)',cols['resid']),('shown_in_log_panels(1=yes)',cols['shown'])],
        'PRB Fig. 4 (fig_shell.pdf): shell-resolved kite density; shown=1 selects the points displayed in all three panels; vertical dotted lines at q=q_star',blocks=cols['rs'])
    series(dest,'spin_fig2_inner_profile.dat',[('x=q/q_star',xx),('middle:-profile(x)_Ry',-profile(xx))],
        'PRB Fig. 4 (fig_shell.pdf), middle panel: leading profile curve')
    axs[0].set(xlabel=r'$q$ (in units of $k_F$)',ylabel=r'$-e(q;r_s)$ (Ry)');axs[0].legend(frameon=False)
    axs[1].set(xlabel=r'$x=q/q_*$',ylabel=r'$-q_*^2 e(q;r_s)$ (Ry)');axs[1].legend(frameon=False)
    axs[2].axhline(0,color='.3',lw=.7);axs[2].set(xscale='log',yscale='symlog',xlabel=r'$x=q/q_*$',ylabel=r'$q_*^2\,[e(q;r_s)/e^{\rm lead}(q)-1]$');axs[2].set_ylim(-150,3);axs[2].set_yticks([0,-3,-10,-30,-100]);axs[2].set_yticklabels(['0','-3','-10','-30','-100']);axs[2].legend(frameon=False,loc='upper right')
    for a,title in zip(axs,['(a) Radial energy density','(b) Collapse onto the leading profile','(c) Relative deviation, scaled by $q_*^2$']):
        a.set_title(title,fontsize=9)
    fig.tight_layout();savefig(fig,folder/'fig_shell.pdf',datetime.datetime(2026,9,27,tzinfo=datetime.timezone.utc))
    rows=s['tableI'];rsd=np.array([r['rs'] for r in rows]);ed=np.array([r['mean_mry'] for r in rows]);sd=np.array([r['se_mry'] for r in rows])
    rs=np.geomspace(100,1e9,400);three=lambda r:1000*asymptotic_kite(r)
    interp=lambda r:1000*(.04836/(1+.10215*r)-.01382*r*np.log1p(1/(.46529*r+.00364*r**1.5)))
    fig,(ax,bx)=plt.subplots(2,1,figsize=(3.4,4.8),sharex=True)
    ax.loglog(rs,-three(rs),color='0.2',label='three-term expansion');ax.loglog(rs,-interp(rs),'--',color='C3',label='BRM interpolation')
    ax.errorbar(rsd,-ed,yerr=sd,fmt='o',ms=3,mfc='w',color='C0',label='direct evaluation')
    br=np.array([r for r in s['reference']['tableI'] if r[1] is not None]);ax.errorbar(br[:,0],-br[:,1],yerr=br[:,2],fmt='s',ms=3,mfc='w',color='C3',label='published BRM')
    ax.set(ylabel=r'$-e_{\rm kite}$ (mRy)',ylim=(.05,60));ax.legend(frameon=False,loc='lower left')
    bx.axhspan(-1,1,color='.88');bx.axhline(0,color='.3',lw=.7)
    bx.errorbar(rsd,(ed-three(rsd))/sd,yerr=1,fmt='o',ms=3,mfc='w',color='C0',label='vs three-term')
    bx.plot(rsd,(ed-interp(rsd))/sd,'s--',ms=3,mfc='w',color='C3',label='vs BRM interpolation')
    bx.set(xscale='log',xlabel=r'$r_s$',ylabel=r'(direct $-$ model)$/\sigma$');bx.legend(frameon=False)
    ax.set_title('(a) Low-density energy',fontsize=9)
    bx.set_title('(b) Deviations in statistical-error units',fontsize=9)
    fig.tight_layout();savefig(fig,folder/'fig_lowdensity.pdf')
    series(dest,'spin_fig3_lowdensity_points.dat',[('rs',rsd),('top:-e_kite_direct_mRy',-ed),('top:error_mRy',sd),
        ('bottom:(direct-three_term)/sigma',(ed-three(rsd))/sd),('bottom:(direct-BRM_interpolation)/sigma',(ed-interp(rsd))/sd)],
        'PRB Fig. 5 (fig_lowdensity.pdf): direct evaluation (Table I) in the top panel and its deviations in the bottom panel (error bar 1)')
    series(dest,'spin_fig3_lowdensity_curves.dat',[('rs',rs),('top:-three_term_mRy',-three(rs)),('top:-BRM_interpolation_mRy',-interp(rs))],
        'PRB Fig. 5 (fig_lowdensity.pdf), top panel: three-term expansion and BRM interpolation curves')
    series(dest,'spin_fig3_brm_published.dat',[('rs',br[:,0]),('top:-e_kite_BRM_mRy',-br[:,1]),('top:error_mRy',br[:,2])],
        'PRB Fig. 5 (fig_lowdensity.pdf), top panel: published BRM values (production/reference_brm.json)')
    Q=q_star(rsd);r1=Q*ed/1000;r2=Q**2*(ed/1000-A1/Q);r3=Q**3*(ed/1000-A1/Q-A2/Q**2)
    fig,ax=plt.subplots(figsize=(3.4,2.8))
    for power,vals,limit,fmt,label in [(1,r1,A1,'o',r'$q_*e_{\rm kite}\to A_1$'),(2,r2,A2,'s',r'$\mathcal{R}_2\to A_2$'),(3,r3,A3,'^',r'$\mathcal{R}_3\to A_3$')]:
        e=sd/1000*Q**power;good=np.ones(len(Q),dtype=bool) if power==1 else e<.25
        ax.errorbar(Q[good],vals[good],yerr=e[good],fmt=fmt,ms=3,mfc='w',capsize=1.5,color=f'C{power-1}',label=label)
        ax.axhline(limit,color=f'C{power-1}',ls='--',lw=.8)
    ax.set(xscale='log',xlim=(4,40),xlabel=r'$q_*$',ylabel='scaled energy (Ry)');ax.legend(frameon=False,loc='center right')
    ax.xaxis.set_major_locator(FixedLocator([5,10,20,40]))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    fig.tight_layout();savefig(fig,folder/'fig_residuals.pdf')
    e1,e2,e3=(sd/1000*Q**p for p in (1,2,3))
    series(dest,'spin_fig5_residuals.dat',[('rs',rsd),('q_star',Q),('q_star*e_kite_Ry',r1),('error1_Ry',e1),
        ('R2_Ry',r2),('error2_Ry',e2),('shown2(1=yes)',(e2<.25).astype(float)),('R3_Ry',r3),('error3_Ry',e3),('shown3(1=yes)',(e3<.25).astype(float))],
        f'PRB Fig. 9 (fig_residuals.pdf): scaled energies against q_star (log axis, 4..40); limits A1={A1:.17g} A2={A2:.17g} A3={A3:.17g} Ry drawn as dashed lines; R2, R3 are drawn only where shown=1')
    m=s['map'];r=np.asarray(m['rs']);z=np.asarray(m['zeta']);E=np.asarray(m['mean_mry']);cross,spl=interpolate_zero(r,z,E)
    lg=np.linspace(np.log(r[0]),np.log(r[-1]),400);zg=np.linspace(0,1,300);eg=spl(lg,zg);delta=eg-spl(lg,[0])
    RS,ZE=np.meshgrid(np.exp(lg),zg,indexing='ij')
    fig,ax=plt.subplots(1,2,figsize=(7.4,3.25),gridspec_kw={'width_ratios':[1,1.15]})
    vmax=np.abs(eg).max()
    im=ax[0].pcolormesh(RS,ZE,eg,cmap='RdBu_r',vmin=-vmax,vmax=vmax,
                        shading='gouraud',rasterized=True)
    cs=ax[0].contour(RS,ZE,eg,levels=[-20,-10,0,10,20],colors='k',linewidths=.5)
    ax[0].clabel(cs,fmt='%d',fontsize=7)
    # Mark the actual integration nodes so that interpolation is not mistaken for data.
    rr,zz=np.meshgrid(r,z,indexing='ij')
    ax[0].plot(rr.ravel(),zz.ravel(),'.',color='0.15',ms=1.8,alpha=.55,clip_on=False)
    ax[0].set(xscale='log',xlabel=r'$r_s$',ylabel=r'Spin polarization $\zeta$',
              xlim=(r[0],r[-1]),ylim=(0,1))
    cb=fig.colorbar(im,ax=ax[0],pad=.02)
    cb.set_label(r'$e_{\rm kite}$ (mRy per electron)',fontsize=8)
    ax[0].set_title('(a) Kite energy',fontsize=10)

    # Form differences before estimating variances; retain common random numbers
    # across polarization nodes, and combine the two independent terms in quadrature.
    raw={(row['rs'],row['zeta']):np.asarray(row['term_replicates_ry']) for row in records('map')}
    point_rows=[];curve_cols=[('rs',np.exp(lg))]
    for zz,color,marker,ls in [(0.4,'#0072B2','s','-'),(0.8,'#D55E00','o','--'),(1.0,'#009E73','^','-.')]:
        paired=np.array([1000*(raw[rs,zz]-raw[rs,0.0]) for rs in r])
        means=paired.mean(axis=1).sum(axis=1)
        errors=np.sqrt(paired.var(axis=1,ddof=1).sum(axis=1)/paired.shape[1])
        np.testing.assert_allclose(means,E[:,list(z).index(zz)]-E[:,0],rtol=0,atol=3e-14)
        curve=(spl(lg,[zz])-spl(lg,[0])).ravel()
        ax[1].plot(np.exp(lg),curve,ls,color=color,lw=1.1)
        ax[1].errorbar(r,means,yerr=errors,fmt=marker,color=color,mfc='white',
                       ms=3.5,capsize=1.7,elinewidth=.8,label=rf'$\zeta={zz:g}$')
        point_rows.extend(zip(r,[zz]*len(r),means,errors))
        curve_cols.append((f'delta_e_zeta={zz:g}_mRy',curve))
    ax[1].axhline(0,color='0.35',lw=.8)
    ax[1].plot(cross,0,'k*',ms=7,zorder=5)
    ax[1].annotate(rf'$\zeta=1$: $r_s\simeq{cross:.0f}$',xy=(cross,0),xytext=(22,-3.8),
                    fontsize=8,arrowprops=dict(arrowstyle='->',lw=.7))
    ax[1].set(xscale='log',xlabel=r'$r_s$',
              ylabel=r'$\Delta e_{\rm kite}=e_{\rm kite}(r_s,\zeta)-e_{\rm kite}(r_s,0)$ (mRy)',
              xlim=(r[0],r[-1]),ylim=(-5,15))
    ax[1].set_title('(b) Polarization energy from the kite',fontsize=10)
    ax[1].legend(frameon=False,loc='upper right')
    fig.tight_layout(w_pad=1.5);savefig(fig,folder/'fig_map.pdf')
    series(dest,'spin_fig4_map_surface.dat',[('rs',RS.ravel()),('zeta',ZE.ravel()),('panel_a:e_kite_mRy',eg.ravel()),('auxiliary:delta_e_mRy',delta.ravel())],
        f'PRB Fig. 6(a): interpolated energy surface, {len(lg)} x {len(zg)} points; '
        'contours -20,-10,0,10,20 mRy; dots are integration nodes in figs/data/map.dat. Column 4 is retained as auxiliary data, not a panel.',fmt='%.10g')
    points=np.array(point_rows)
    series(dest,'spin_fig4_polarization_points.dat',[(name,points[:,i]) for i,name in enumerate(['rs','zeta','delta_e_mRy','paired_standard_error_mRy'])],
        'PRB Fig. 6(b): direct polarization differences; within-term paired subtraction before combining independent term variances',blocks=points[:,1])
    series(dest,'spin_fig4_polarization_curves.dat',curve_cols,
        'PRB Fig. 6(b): cubic-spline guides through the direct spin-map means; no additional integrals')
    series(dest,'spin_fig4_map_crossing.dat',[('rs_crossing_at_zeta=1',[cross]),('delta_e_mRy',[0.0])],
        'PRB Fig. 6(b): star at the zero of the full-polarization difference; ordinate is now energy difference, not polarization')
    # Compatibility copies: the old PRB entry points now use these same products.
    (dest/'figs').mkdir(exist_ok=True)
    for f in folder.glob('*.pdf'): shutil.copyfile(f,dest/'figs'/f.name)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT);a=p.parse_args()
    # Finish all products before publishing them; a failed plot cannot truncate
    # the currently released PDF. Each completed file is replaced atomically.
    with tempfile.TemporaryDirectory(prefix='prb-render-stage-') as tmp:
        stage=Path(tmp);s=summarize();write_tables(s,stage);plots(s,stage)
        expected=json.loads((ROOT/'production/product_lock.json').read_text())['sha256']
        actual={str(f.relative_to(stage)):hashlib.sha256(f.read_bytes()).hexdigest()
                for f in stage.rglob('*') if f.is_file()}
        if actual != expected:
            raise RuntimeError('Generated products differ from product_lock.json; no output was replaced')
        for source in sorted(stage.rglob('*')):
            if not source.is_file():continue
            target=a.output/source.relative_to(stage);target.parent.mkdir(parents=True,exist_ok=True)
            pending=target.with_name(target.name+'.new')
            shutil.copyfile(source,pending);pending.replace(target)
    print('Generated 4 numerical figures, 4 table bodies, numerical text macros and CSV/JSON data.')


if __name__=='__main__':main()
