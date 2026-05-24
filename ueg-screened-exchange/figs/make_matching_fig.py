"""PRB Fig. (matching): the moving momentum scale q_* and the fixed-momentum region.

Left: the leading moving-scale kernel of the coupling-averaged kite,
    e^lead(q) = (3/pi^2) F((q_*/q)^4) / q^2,   F(g) = -(4/3)(sqrt(1+g)-1)/(sqrt(1+g)+1)^2,
(PRB Eqs. (Fprofile) and (A1-integral)), for r_s = 10^3 and 10^6 (q_* = (4 a r_s/3)^(1/4)). The area under each curve is the
magnitude -A_1/q_* of the leading energy (logarithmic measure); the peak moves with q_* while the fixed-momentum
region q ~ 1 (Fermi scale, shaded) stays put.
Right: the same two kernels in the scaled variable x = q/q_*, multiplied by q_* so that the
curves collapse onto one profile whose logarithmic area is -A_1. On this axis the fixed-momentum region
q ~ 1 (shaded, one band per r_s) moves toward x -> 0 as the density decreases, and the plotted profile vanishes there linearly in x: the fixed region enters only at the next order, through the
matching that gives A_2. The limiting slope is shown for comparison.
Output: figs/fig_matching.pdf (drawn from the leading kernel formula; not a production data figure).
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
aL = (4/(9*np.pi))**(1/3); a = 4*aL/np.pi
def qstar(rs): return (4*a*rs/3)**0.25
def frequency_profile(g):                       # coupling-averaged leading frequency profile, PRB Eq. (Fprofile)
    r = np.sqrt(1+g); return -4/3*(r-1)/(r+1)**2
def leading(q, qs):               # leading kernel e^lead(q) = (3/pi^2) F((q_*/q)^4)/q^2, whose q integral is A_1/q_*
    return 3/np.pi**2*frequency_profile((qs/q)**4)/q**2
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.9))
q = np.logspace(-0.5, 2.5, 600)
for rs, ls, col in ((1e3, '-', 'C0'), (1e6, '--', 'C3')):
    qs = qstar(rs); y = leading(q, qs)
    ax1.plot(q, -y*q, ls, color=col, label=r'$r_s=10^{%d}$, $q_*=%.1f$' % (np.log10(rs), qs))
    ax1.axvline(qs, color=col, lw=0.6, ls=ls, alpha=0.7)
ax1.axvspan(0.3, 2.0, color='0.85', zorder=0)
ax1.set_xscale('log'); ax1.set_xlabel(r'$q$ (in units of $k_F$)'); ax1.set_ylabel(r'$-q\,e^{\rm lead}(q)$ (Ry)')
ax1.legend(fontsize=7, loc='upper right', frameon=False); ax1.set_title('(a) Leading radial energy density', fontsize=9)
# right: the same two kernels in the scaled variable x = q/q_*, showing the collapse onto one
# profile and where the fixed-momentum region (q ~ 1) sits on this axis
x = np.logspace(-2.2, 1.6, 600)
for rs, ls, col, lw, al in ((1e3, '-', 'C0', 3.2, 0.45), (1e6, '--', 'C3', 1.2, 1.0)):
    qs = qstar(rs); qq = qs*x; y = -leading(qq, qs)*qq*qs         # -q_* q e^lead(q): the same area, scaled by q_*
    ax2.plot(x, y, ls, color=col, lw=lw, alpha=al, label=r'$r_s=10^{%d}$' % np.log10(rs))
    ax2.axvspan(0.3/qs, 2.0/qs, color='0.85', zorder=0)         # fixed-momentum region q in [0.3, 2]
ax2.plot(x[x < 0.25], 4*x[x < 0.25]/np.pi**2, ':', color='k', lw=1.3, label=r'small-$x$ limit $4x/\pi^2$')
ax2.set_xscale('log'); ax2.set_xlabel(r'$x=q/q_*$'); ax2.set_ylabel(r'$-q_*\,q\,e^{\rm lead}(q)$')
ax2.set_title('(b) Scaled form', fontsize=9); ax2.set_ylim(0, 0.135); ax2.legend(fontsize=7, loc='upper right', frameon=False, borderaxespad=0.6)
fig.tight_layout()
out = Path(__file__).resolve().parent / 'fig_matching.pdf'
fig.savefig(out, metadata={'CreationDate': None, 'ModDate': None}); print('wrote', out)
