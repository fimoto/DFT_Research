"""Generate LaTeX for new Supplemental Material sections S8-S11 from validation/results results."""
import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # run from the repository root
import json,glob,math
from pathlib import Path
import numpy as np
import io,contextlib
with contextlib.redirect_stdout(io.StringIO()): exec(open("analyze_validation.py").read())
V=Path('validation/results')
L=[]
def tab(hdr,rows,caption,label,cols):
    L.append(r'\begin{table}[h]\centering\caption{%s}\label{%s}'%(caption,label))
    wide=len(cols)>=7
    if wide: L.append(r'\resizebox{\linewidth}{!}{%')
    L.append(r'\begin{tabular}{%s}\hline'%cols); L.append(' & '.join(hdr)+r'\\\hline')
    for r in rows: L.append(' & '.join(r)+r'\\')
    L.append(r'\hline\end{tabular}'+('}' if wide else '')+r'\end{table}')
# ---- S8
L.append(r'\section{Spin-dependent checks of the direct integration}\label{sm:spin-checks}')
L.append(r"""All results in this and the following three sections were obtained with Estimator~A
(\texttt{production.direct}) at $m=19$ (i.e.\ $2^{19}$ Sobol points per term and
randomization), 16 independent Sobol randomizations, master seed 20261001, and the
momentum-scale parameter 2 unless stated otherwise. Energies are in mRy per electron.
Errors are $\sqrt{(s_1^2+s_2^2)/R}$ as in Sec.~IV of the main text. The replicate
values, settings, and source hashes of every run are stored in
\texttt{validation\_20261001/}.""")
L.append(r"""\paragraph{Full-polarization identity.} Equation~(49) of the main text,
$e_{\rm kite}(r_s,1)=e_{\rm kite}(2^{-4/3}r_s,0)$, was tested by two separate direct
integrations: the fully polarized response at $r_s$ and the unpolarized response at
$2^{-4/3}r_s$, computed in the code as \texttt{rs*2**(-4/3)}. Both used the same
momentum-scale parameter and the same Sobol points, so that the comparison can be made
replicate by replicate. Table~\ref{sm:tab:identity} lists the energies, the paired
difference of the replicate means, and the largest difference between the per-replicate
term values. The agreement is at the level of floating-point round-off: the
spin-polarized response and the rescaled exchange loop reduce, in the implementation,
to the same integrand.""")
I=load('identity'); rows=[]
for (rs,z),r in sorted(I.items()):
    if z!=1: continue
    rs2=round(rs*2**(-4/3),12); b=[v for k,v in I.items() if k[1]==0 and abs(k[0]-rs2)<1e-9*max(1,rs2)][0]
    D=T(r)-T(b); m,s=paired(D)
    rows.append([f'{rs:g}',f'${r["mean_mry"]:.4f}\\pm{r["se_mry"]:.4f}$',f'${b["mean_mry"]:.4f}\\pm{b["se_mry"]:.4f}$',f'${m:.1e}$',f'${np.abs(D).max():.1e}$'])
tab(['$r_s$','$e(r_s,1)$','$e(2^{-4/3}r_s,0)$','paired difference','max $|$term difference$|$ (Ry)'],rows,
 'Full-polarization identity, Eq.~(49), from two separate direct integrations ($m=19$, $R=16$, scale 2, seed 20261001). The paired difference is the mean over replicates of the difference of the two full-energy replicate values; its statistical error is of the same order as the quoted value and both are round-off.','sm:tab:identity','ccccc')
N=load('nearfull')
if N:
    L.append(r"""\paragraph{Near-full polarization.} Table~\ref{sm:tab:nearfull} gives the energies at
$\zeta=0.95$ and $0.99$, where the minority Fermi radius $x_-=(1-\zeta)^{1/3}$ is $0.37$ and
$0.22$. The values interpolate smoothly between the $\zeta=0.8$ and $\zeta=1$ entries of
Table~\ref{sm:tab:fig8}, and the statistical errors do not grow as the minority population
empties.""")
    rows=[[f'{rs:g}',f'{z}',f'${r["mean_mry"]:.4f}\\pm{r["se_mry"]:.4f}$'] for (rs,z),r in sorted(N.items())]
    tab(['$r_s$','$\\zeta$','$e_{\\rm kite}$'],rows,'Kite energy near full polarization ($m=19$, $R=16$, scale 2).','sm:tab:nearfull','ccc')
# ---- S9
levels={m:load(f'fig8_m{m}') for m in (17,19,21)}; B=load('fig8_m19_seedB')
L.append(r'\section{Convergence of the interpolation residuals of Fig.~8}\label{sm:fig8-convergence}')
L.append(r"""The residual $\delta e_{\rm kite}=e(r_s,\zeta)-(1-f_x)e(r_s,0)-f_xe(r_s,1)$ is formed
replicate by replicate from runs that share the same Sobol randomizations across the three
polarizations (common random numbers), and its error is the standard error of these paired
residuals. Table~\ref{sm:tab:fig8} repeats this at three nested sample sizes,
$m=17,19,21$ ($2^{17}$--$2^{21}$ points per term and randomization) with the same seed, so
that the change between sample sizes can also be evaluated replicate by replicate, and at a
second seed block with twice as many randomizations. The paired changes between sample
sizes are at most two standard errors ($-0.0023\pm0.0011$~mRy at $r_s=1$, $\zeta=0.8$ from
$2^{19}$ to $2^{21}$ points), and the two seed blocks differ by at most $2.1$ standard errors
($0.0064\pm0.0030$~mRy at $r_s=10$, $\zeta=0.8$). The residuals are therefore established to
about $0.005$~mRy: the sign and magnitude shown in Fig.~8 are confirmed, and the error bars
of Fig.~8 should be read as statistical standard errors at the stated sample size. The
energies themselves are listed in Table~\ref{sm:tab:fig8-energies}.""")
rows=[]
for rs in (1.0,10.0,60.0):
    for z in (0.4,0.8):
        Ds={m:resid(levels[m],rs,z) for m in (17,19,21)}
        c=[('$%+.4f\\pm%.4f$'%paired(Ds[m])) if Ds[m] is not None else '--' for m in (17,19,21)]
        d1=('$%+.4f\\pm%.4f$'%paired(Ds[19]-Ds[17])) if Ds[17] is not None and Ds[19] is not None else '--'
        d2=('$%+.4f\\pm%.4f$'%paired(Ds[21]-Ds[19])) if Ds[21] is not None and Ds[19] is not None else '--'
        bb=('$%+.4f\\pm%.4f$'%paired(resid(B,rs,z))) if resid(B,rs,z) is not None else '--'
        rows.append([f'{rs:g}',f'{z}']+c+[d1,d2,bb])
tab(['$r_s$','$\\zeta$','$m=17$','$m=19$','$m=21$','$(19)-(17)$','$(21)-(19)$','seed B, $R=32$, $m=19$'],rows,
 'Interpolation residual $\\delta e_{\\rm kite}$ (mRy) at nested sample sizes (seed 20261001, $R=16$) and at an independent seed block (20262001, $R=32$). Differences between sample sizes are paired by replicate.','sm:tab:fig8','cccccccc')
rows=[]
for (rs,z),r in sorted(levels[21].items() or levels[19].items()):
    r19=levels[19][(rs,z)]; r21=levels[21].get((rs,z))
    rows.append([f'{rs:g}',f'{z}',f'${r19["mean_mry"]:.4f}\\pm{r19["se_mry"]:.4f}$',(f'${r21["mean_mry"]:.4f}\\pm{r21["se_mry"]:.4f}$' if r21 else '--')])
tab(['$r_s$','$\\zeta$','$e_{\\rm kite}$, $m=19$','$e_{\\rm kite}$, $m=21$'],rows,'Direct energies (mRy) underlying Table~\\ref{sm:tab:fig8}.','sm:tab:fig8-energies','cccc')
# ---- S10 axis
A={Path(f).stem:json.load(open(f)) for f in glob.glob(str(V/'axis/*.json'))}
if A:
    L.append(r'\section{Independent imaginary-axis integration at intermediate polarization}\label{sm:axis}')
    L.append(r"""Estimator~B (\texttt{spin/reviewer\_numerics.energy}) integrates Eq.~(10) of the main
text on the imaginary-frequency axis with deterministic frequency and momentum quadrature and
randomized sampling of the internal momenta; it shares no code with the BRM representation
used by Estimator~A. Before comparing it with Estimator~A at $\zeta>0$, we tested its
convergence. Table~\ref{sm:tab:axis-bare} evaluates the bare limit, whose exact value is
$e_{2x}=48.358$~mRy: with the default momentum grid ($n_q=24$ log-spaced Simpson nodes per
side) Estimator~B is low by $0.4$~mRy, and the estimate also drifts with the sample size;
with $n_q\\ge64$ and $2^{17}$ points per replicate it agrees with the exact value within one
standard error. The quadrature columns of Table~I of the main text were computed with
$n_q=32$ at $r_s\\ge10^3$, where the energy density is spread over the moving scale; the
refinement to $n_q=64$ recorded in the verification notebook changes the $r_s=10^3$ value by
$0.002\\pm0.005$~mRy, so those columns are unaffected. The grid sensitivity found here
concerns the intermediate densities $r_s\\lesssim60$, where the bare contribution with its
kink at $q=2$ dominates. Table~\ref{sm:tab:axis-refine} varies the quadrature settings one at a time
on the same geometry samples at $r_s=10$, $\zeta=0.8$: the frequency order and the Lindhard
grid $n_k$ are converged (the paired changes vanish), while the momentum grid is not at
$n_q=24$. The row labelled $q_{\min}$ halves the lower cutoff at fixed $n_q$ and therefore also
coarsens the grid; it is not a clean cutoff test. These two tables establish that the
default settings of Estimator~B carry a systematic error of $0.2$--$0.5$~mRy, so that
agreement with Estimator~A can only be expected at $n_q\gtrsim96$.""")
    bt=json.load(open(V/'axis/bare_tests.json'))
    tab(['$n_q$','points per replicate','Estimator B (mRy)','deviation from $e_{2x}$ ($\\sigma$)'],
        [[str(r['nq']),f"$2^{{{r['power']}}}$",f"${r['mean_mry']:.3f}\\pm{r['se_mry']:.3f}$",f"${(r['mean_mry']-bt['exact_e2x_mry'])/r['se_mry']:+.1f}$"] for r in bt['rows']],
        'Bare limit of Estimator B ($R=16$, order 24, $n_k=48$, seed 7) against the exact second-order exchange energy $e_{2x}=48.358$ mRy.','sm:tab:axis-bare','cccc')
    L.append(r"""Table~\ref{sm:tab:axis} compares the two estimators at $\zeta=0.4$ and $0.8$. The
columns with $n_q=24$ were computed with $2^{15}$ points per replicate and show the systematic
offset of the default grid; the columns with $n_q=96$ ($2^{12}$ points per replicate, the largest sample size
affordable here at that grid) agree with Estimator~A within $1.3$ standard errors at all six
points. Table~\\ref{sm:tab:axis-identity} repeats the full-polarization identity with
Estimator~B: the two integrations use the same geometry samples and agree exactly, as for
Estimator~A. The statistical error of Estimator~B is large
because its geometry sample is shared across the frequency and momentum grids; the error of
the difference is the quadrature sum of the two errors. These comparisons confirm Estimator~A
at the level of $0.3$--$0.5$~mRy at intermediate polarization; the $0.01$~mRy level of the
interpolation residuals in Fig.~8 is established by the common-random-number convergence of
Estimator~A itself (Sec.~\ref{sm:fig8-convergence}), not by Estimator~B.""")
    rows=[]
    for rs in (1.0,10.0,60.0):
        for z in (0.4,0.8):
            a=A.get(f'base_rs{rs:g}_z{z}'); d=levels[19].get((rs,z)); c=A.get(f'nq96_rs{rs:g}_z{z}')
            if not a: continue
            am,ase=1000*a['mean_ry'],1000*a['se_ry']; dm,dse=d['mean_mry'],d['se_mry']
            cc=(f"${1000*c['mean_ry']:.3f}\\pm{1000*c['se_ry']:.3f}$",f"${1000*c['mean_ry']-dm:+.3f}\\pm{math.hypot(1000*c['se_ry'],dse):.3f}$") if c else ('--','--')
            rows.append([f'{rs:g}',f'{z}',f'${dm:.3f}\\pm{dse:.3f}$',f'${am:.3f}\\pm{ase:.3f}$',f'${(am-dm):+.3f}\\pm{math.hypot(ase,dse):.3f}$',cc[0],cc[1]])
    tab(['$r_s$','$\\zeta$','A ($m=19$)','B, $n_q=24$, $2^{15}$','B$-$A','B, $n_q=96$, $2^{12}$','B$-$A'],rows,
     'Direct BRM integration (Estimator A) and independent imaginary-axis integration (Estimator B: 16 replicates, order 24, $n_k=48$, seed 30301001) at intermediate polarization (mRy).','sm:tab:axis','ccccccc')
    base=A.get('ref_base')
    if base:
        rows=[]
        for tag,desc in [('nq48','$n_q$: 24 $\\to$ 48'),('nq96','$n_q$: 24 $\\to$ 96'),('order36','frequency order: 24 $\\to$ 36'),('nk96','$n_k$: 48 $\\to$ 96'),('qmin001','$q_{\\min}$: 0.002 $\\to$ 0.001 (same $n_q$)'),('qmax2','$q_{\\max}$: default $\\to$ 2$\\times$default')]:
            r=A.get(f'ref_{tag}')
            if not r: continue
            D=1000*(np.array(r['replicates_ry'])-np.array(base['replicates_ry']))
            rows.append([desc,f'${1000*r["mean_ry"]:.3f}\\pm{1000*r["se_ry"]:.3f}$',f'${D.mean():+.5f}\\pm{D.std(ddof=1)/math.sqrt(len(D)):.5f}$'])
        tab(['change','energy (mRy)','paired change (mRy)'],rows,'Estimator B at $r_s=10$, $\\zeta=0.8$ (power 13, $R=16$, seed 30301001): one quadrature setting changed at a time on the same geometry samples. Baseline: $%.3f\\pm%.3f$ mRy.'%(1000*base['mean_ry'],1000*base['se_ry']),'sm:tab:axis-refine','lcc')
    rows=[]
    for rs in (1.0,10.0):
        a1,a0=A.get(f'id_rs{rs:g}_z1'),A.get(f'id_rs{rs:g}_z0')
        if a1 and a0:
            D=1000*(np.array(a1['replicates_ry'])-np.array(a0['replicates_ry']))
            rows.append([f'{rs:g}',f'${1000*a1["mean_ry"]:.3f}\\pm{1000*a1["se_ry"]:.3f}$',f'${1000*a0["mean_ry"]:.3f}\\pm{1000*a0["se_ry"]:.3f}$',f'${D.mean():+.4f}\\pm{D.std(ddof=1)/math.sqrt(len(D)):.4f}$'])
    if rows: tab(['$r_s$','$e(r_s,1)$','$e(2^{-4/3}r_s,0)$','paired difference'],rows,'Full-polarization identity with Estimator B (power 13, $R=16$, seed 30302001, same geometry samples).','sm:tab:axis-identity','cccc')
# ---- S11 low density
Lo=load('lowdensity')
if Lo:
    L.append(r'\section{Low-density spin dependence}\label{sm:lowdensity-spin}')
    L.append(r"""Table~\ref{sm:tab:lowdensity} compares direct integrations at $\zeta=0,0.4,0.8,1$ with
the three-term expansion Eq.~(46), using the momentum-scale parameter $q_*$ as in Table~I.
The last three columns scale the data as in Appendix~B: $q_*e$ against $A_1M_4(\zeta)$,
$q_*^2(e-A_1M_4/q_*)$ against $A_2M_5(\zeta)$, and the scaled interpolation residual
$q_*^2\delta e$ against $A_2\mathcal B(\zeta)$ of Eq.~(57), with $\delta e$ formed
replicate by replicate as in Sec.~\ref{sm:fig8-convergence}. Figure~\ref{sm:fig:lowdensity}
shows the same quantities. All three scaled quantities approach their limits from below as
$q_*$ increases, at every polarization. The approach is consistent with the next term of the
expansion: for the second column the predicted value including $A_3S_3(\zeta)/q_*$ is
0.457, against the measured $0.471\pm0.005$ at $r_s=10^3$, $\zeta=0$. The $r_s^{-3/4}$ coefficient itself is
not separately resolved by these data, as stated in Appendix~B. The four polarizations at a
given density share the same Sobol points, so their deviations from Eq.~(46) are strongly
correlated (for example the common $-2.5$ to $-2.9$ standard-error offset at $r_s=2\times10^4$),
and the entries are 16-replicate runs, less precise than the 64-replicate entries of Table~I.
""")
    rows=[]
    for (rs,z),r in sorted(Lo.items()):
        q=Q(rs); e=r['mean_mry']/1000; se=r['se_mry']/1000
        three=A1*M(4,z)/q+A2*M(5,z)/q**2+A3*S3(z)/q**3
        D=resid(Lo,rs,z) if 0<z<1 else None
        c3=('$%.4f\\pm%.4f$ (%.4f)'%(q*q*paired(D)[0]/1000,q*q*paired(D)[1]/1000,A2*Bz(z))) if D is not None else '--'
        rows.append([f'{rs:g}',f'{z}',f'${r["mean_mry"]:.3f}\\pm{r["se_mry"]:.3f}$',f'{1000*three:.3f}',f'${(e-three)/se:+.1f}$',
                     '$%.4f\\pm%.4f$ (%.4f)'%(q*e,q*se,A1*M(4,z)),'$%.3f\\pm%.3f$ (%.3f)'%(q*q*(e-A1*M(4,z)/q),q*q*se,A2*M(5,z)),c3])
    tab(['$r_s$','$\\zeta$','direct (mRy)','Eq.~(46) (mRy)','$(\\mathrm{direct}-\\mathrm{Eq.~46})/\\sigma$','$q_*e$ ($A_1M_4$)','$q_*^2(e-A_1M_4/q_*)$ ($A_2M_5$)','$q_*^2\\delta e$ ($A_2\\mathcal B$)'],rows,
     'Low-density direct integrations at fixed spin polarization ($m=19$, $R=16$, scale $q_*$, seed 20261001) compared with the three-term expansion. Scaled quantities are in Ry; the predicted limits are in parentheses.','sm:tab:lowdensity','cccccccc')
    L.append(r"""\begin{figure}[h]\centering\includegraphics[width=0.95\linewidth]{fig_sm_lowdensity_spin.pdf}
\caption{Scaled low-density data of Table~\ref{sm:tab:lowdensity}. (a) $q_*e_{\rm kite}$ versus $\zeta$ at four densities, with the limit $A_1M_4(\zeta)$ (line). (b) $q_*^2(e_{\rm kite}-A_1M_4/q_*)$ with the limit $A_2M_5(\zeta)$. (c) $q_*^2\delta e_{\rm kite}$ with the limit $A_2\mathcal B(\zeta)$. Bars are one statistical standard error; in (c) the residual is formed replicate by replicate.}\label{sm:fig:lowdensity}\end{figure}""")
open(V/'supplement_additions_20261001.tex','w').write('\n'.join(L)); print('written',len(L),'lines')
