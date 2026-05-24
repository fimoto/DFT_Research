"""Figure of SM S11: scaled low-density energies versus spin polarization, from the stored replicates."""
import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # run from the repository root
import io, contextlib, math, sys
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
with contextlib.redirect_stdout(io.StringIO()):
    exec(open('validation/analyze_validation.py').read())   # loaders, constants, residual construction
out = sys.argv[1] if len(sys.argv) > 1 else 'validation/results/fig_sm_lowdensity_spin.pdf'
Lo = load('lowdensity'); zs = np.linspace(0, 1, 201)
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.2))
cols = {1000: 'C0', 5000: 'C1', 20000: 'C2', 100000: 'C3'}
for rs in (1000.0, 5000.0, 20000.0, 100000.0):
    pts = [(z, Lo[(rs, z)]) for z in (0.0, 0.4, 0.8, 1.0) if (rs, z) in Lo]
    q = Q(rs); zz = [p[0] for p in pts]
    e = np.array([p[1]['mean_mry']/1000 for p in pts]); se = np.array([p[1]['se_mry']/1000 for p in pts])
    ax[0].errorbar(zz, q*e, yerr=q*se, fmt='o', color=cols[rs], ms=4,
                   label=f'$r_s=10^{{{math.log10(rs):.2g}}}$' if rs in (1000.0, 100000.0) else f'$r_s={rs:g}$')
    ax[1].errorbar(zz, q*q*(e - A1*np.array([M(4, z) for z in zz])/q), yerr=q*q*se, fmt='o', color=cols[rs], ms=4)
    for z in (0.4, 0.8):
        D = resid(Lo, rs, z)
        if D is not None:
            m, s = paired(D); ax[2].errorbar([z], [q*q*m/1000], yerr=[q*q*s/1000], fmt='o', color=cols[rs], ms=4)
ax[0].plot(zs, [A1*M(4, z) for z in zs], 'k-', lw=1, label='$A_1M_4(\\zeta)$'); ax[0].set_ylabel('$q_*\\,e_{\\rm kite}$ (Ry)')
ax[1].plot(zs, [A2*M(5, z) for z in zs], 'k-', lw=1, label='$A_2M_5(\\zeta)$'); ax[1].set_ylabel('$q_*^2(e_{\\rm kite}-A_1M_4/q_*)$ (Ry)')
ax[2].plot(zs, [A2*Bz(z) for z in zs], 'k-', lw=1, label='$A_2\\mathcal{B}(\\zeta)$'); ax[2].set_ylabel('$q_*^2\\,\\delta e_{\\rm kite}$ (Ry)')
for a, t in zip(ax, ['(a)', '(b)', '(c)']): a.set_xlabel('$\\zeta$'); a.text(0.03, 0.9, t, transform=a.transAxes)
ax[0].legend(fontsize=7, loc='lower left'); ax[1].legend(fontsize=7); ax[2].legend(fontsize=7)
plt.tight_layout(); plt.savefig(out); print('wrote', out)
