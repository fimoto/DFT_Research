"""Aggregate validation/results results into tables (paired replicates, common random numbers)."""
import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # run from the repository root
import json,glob,math,sys
from pathlib import Path
import numpy as np
V=Path('validation/results')
def load(d):
    out={}
    for f in sorted(glob.glob(str(V/d/'*.json'))):
        r=json.load(open(f)); out[(round(r['rs'],12),round(r['zeta'],12))]=r
    return out
def T(r): return np.asarray(r['term_replicates_ry'])   # (R,2) Ry
def fx(z): return (((1+z)**(4/3)+(1-z)**(4/3))/2-1)/(2**(1/3)-1)
def M(n,z): return ((1+z)**(n/3)+(1-z)**(n/3))/2
AL=(4/(9*math.pi))**(1/3); ATF=4*AL/math.pi; A1=-0.2146937609891; A2=8/15; A3=-0.4173253626244
def Q(rs): return (4*ATF*rs/3)**.25
def S3(z): return (1160*M(6,z)+891*M(4,z)*M(5,z))/2051
def Bz(z): return M(5,z)-1-(2**(2/3)-1)*fx(z)
E2X=math.log(2)/3-3*1.2020569031595943/(2*math.pi**2)
def paired(D):   # D (R,2) Ry -> mean,se in mRy
    return 1000*D.mean(axis=0).sum(), 1000*math.sqrt(D.var(axis=0,ddof=1).sum()/len(D))
def check_same(a,b):
    assert a['settings']['m']==b['settings']['m'] and a['settings']['seed']==b['settings']['seed'] and len(T(a))==len(T(b))
out=[]
# 1 identity
I=load('identity')
if I:
    out.append('## Full-polarization identity e(rs,1) = e(2^{-4/3} rs, 0)  [m=19, R=16, scale 2, seed 20261001, same Sobol points]')
    out.append('| rs | e(rs,1) mRy | e(2^-4/3 rs,0) mRy | paired diff mRy | max |term diff| per replicate (Ry) |'); out.append('|---|---|---|---|---|')
    for (rs,z),r in sorted(I.items()):
        if z!=1: continue
        rs2=round(rs*2**(-4/3),12); b=[v for k,v in I.items() if k[1]==0 and abs(k[0]-rs2)<1e-9*max(1,rs2)]
        if not b: continue
        b=b[0]; check_same(r,b); D=T(r)-T(b); m,s=paired(D)
        out.append(f'| {rs:g} | {r["mean_mry"]:.4f} ± {r["se_mry"]:.4f} | {b["mean_mry"]:.4f} ± {b["se_mry"]:.4f} | {m:+.2e} ± {s:.1e} | {np.abs(D).max():.1e} |')
# 2 fig8 convergence
def resid(d,rs,z):
    r0,rz,r1=d.get((rs,0.0)),d.get((rs,z)),d.get((rs,1.0))
    if not (r0 and rz and r1): return None
    check_same(r0,rz); check_same(r0,r1); f=fx(z)
    D=T(rz)-(1-f)*T(r0)-f*T(r1); return D
levels={m:load(f'fig8_m{m}') for m in (17,19,21)}
if any(levels.values()):
    out.append('\n## Interpolation residual delta e = e(rs,z) - (1-f_x)e(rs,0) - f_x e(rs,1)  [scale 2, seed 20261001, R=16; paired replicates]')
    out.append('| rs | zeta | m=17 | m=19 | m=21 | (m19-m17) paired | (m21-m19) paired |'); out.append('|---|---|---|---|---|---|---|')
    for rs in (1.0,10.0,60.0):
        for z in (0.4,0.8):
            Ds={m:resid(levels[m],rs,z) for m in (17,19,21)}
            cells=[]
            for m in (17,19,21):
                cells.append('%+.4f ± %.4f'%paired(Ds[m]) if Ds[m] is not None else '—')
            d1='%+.4f ± %.4f'%paired(Ds[19]-Ds[17]) if Ds[17] is not None and Ds[19] is not None else '—'
            d2='%+.4f ± %.4f'%paired(Ds[21]-Ds[19]) if Ds[21] is not None and Ds[19] is not None else '—'
            out.append(f'| {rs:g} | {z} | {cells[0]} | {cells[1]} | {cells[2]} | {d1} | {d2} |')
    out.append('\nEnergies (mRy), m=19:')
    out.append('| rs | zeta | e direct |'); out.append('|---|---|---|')
    for (rs,z),r in sorted(levels[19].items()): out.append(f'| {rs:g} | {z} | {r["mean_mry"]:.4f} ± {r["se_mry"]:.4f} |')
B=load('fig8_m19_seedB')
if B:
    out.append('\n## Independent seed block (seed 20262001, R=32, m=19) at rs=10')
    out.append('| zeta | delta e block B | delta e block A (R=16) |'); out.append('|---|---|---|')
    for z in (0.4,0.8):
        Db=resid(B,10.0,z); Da=resid(levels[19],10.0,z)
        out.append(f'| {z} | %+.4f ± %.4f | %s |'%(*paired(Db),'%+.4f ± %.4f'%paired(Da) if Da is not None else '—'))
# 5 low density
L=load('lowdensity')
if L:
    out.append('\n## Low-density spin dependence [scale q*, m=19, R=16, seed 20261001]')
    out.append('| rs | zeta | e direct mRy | three-term Eq.(46) mRy | (direct-3term)/σ | Q e vs A1 M4 | Q^2(e-A1M4/Q) vs A2 M5 | Q^2 δe vs A2 B(ζ) |'); out.append('|---|---|---|---|---|---|---|---|')
    for (rs,z),r in sorted(L.items()):
        q=Q(rs); e=r['mean_mry']/1000; se=r['se_mry']/1000
        three=A1*M(4,z)/q+A2*M(5,z)/q**2+A3*S3(z)/q**3
        c1='%.4f (%.4f)'%(q*e,A1*M(4,z)); c2='%.3f ± %.3f (%.3f)'%(q*q*(e-A1*M(4,z)/q),q*q*se,A2*M(5,z))
        D=resid(L,rs,z) if 0<z<1 else None
        c3='%.4f ± %.4f (%.4f)'%(q*q*paired(D)[0]/1000,q*q*paired(D)[1]/1000,A2*Bz(z)) if D is not None else '—'
        out.append(f'| {rs:g} | {z} | {r["mean_mry"]:.4f} ± {r["se_mry"]:.4f} | {1000*three:.4f} | {(e-three)/se:+.2f} | {c1} | {c2} | {c3} |')
N=load('nearfull')
if N:
    out.append('\n## Near-full polarization [scale 2, m=19, R=16]')
    out.append('| rs | zeta | e direct mRy |'); out.append('|---|---|---|')
    for (rs,z),r in sorted(N.items()): out.append(f'| {rs:g} | {z} | {r["mean_mry"]:.4f} ± {r["se_mry"]:.4f} |')
txt='\n'.join(out); open(V/'summary.md','w').write(txt); print(txt)
