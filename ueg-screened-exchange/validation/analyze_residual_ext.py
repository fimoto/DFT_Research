"""Estimator B interpolation residuals (paired across polarizations) from results/residual_grid_ext.json,
compared with Estimator A (main-text Table IV; m=19 validation values at r_s=10)."""
import os,json,math
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
G=json.load(open('validation/results/residual_grid_ext.json'))
fx=lambda z:(((1+z)**(4/3)+(1-z)**(4/3))/2-1)/(2**(1/3)-1)
A={(1.0,0.4):(-0.2176,0.0023),(1.0,0.8):(-0.6298,0.0057),(10.0,0.4):(0.0944,0.0015),(10.0,0.8):(0.2234,0.0026),(60.0,0.4):(0.2704,0.0039),(60.0,0.8):(0.6384,0.0076)}
out={}
for rs,power,nq in [(10.0,15,48),(1.0,13,48),(60.0,13,48)]:
    try: r={z:np.array(G[f'rs{rs:g}_p{power}_nq{nq}_z{z}']['replicates_ry']) for z in (0.0,0.4,0.8,1.0)}
    except KeyError: continue
    for z in (0.4,0.8):
        D=1000*(r[z]-(1-fx(z))*r[0.0]-fx(z)*r[1.0]); m,s=D.mean(),D.std(ddof=1)/math.sqrt(len(D)); a,sa=A[(rs,z)]
        out[f'rs{rs:g}_z{z}']=dict(B=m,B_err=s,A=a,A_err=sa,pull=(m-a)/math.hypot(s,sa))
        print(f'rs={rs:<4g} zeta={z}: B = {m:+.4f} +- {s:.4f}   A = {a:+.4f} +- {sa:.4f}   difference = {m-a:+.4f} +- {math.hypot(s,sa):.4f}')
json.dump(out,open('validation/results/residual_ext_summary.json','w'),indent=1)
