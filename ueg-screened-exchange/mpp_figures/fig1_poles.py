"""Main-text Fig. 1: positions s of the Mellin poles as functions of d.

Reads ../mpp_checks/figure_data/fig1_poles.csv and fig1_coincidences.csv
(written by ../mpp_checks/fig_data.py).
"""
import csv, os
import numpy as np
import matplotlib.pyplot as plt
import figstyle as st

NAME = 'fig1_poles'
XLIM = (1.45, 5.0)
YLIM = (0.0, 8.2)


def read(name):
    with open(os.path.join(st.DATA, name)) as f:
        return list(csv.DictReader(f))


def main():
    st.apply()
    poles = read('fig1_poles.csv')
    d = np.array([float(r['d']) for r in poles])
    ncol = sorted(k for k in poles[0] if k.startswith('s_N'))
    col = {k: np.array([float(r[k]) for r in poles]) for k in poles[0]}

    fig = plt.figure(figsize=st.figsize(NAME, 0.91))
    ax = fig.add_axes([0.125, 0.135, 0.845, 0.845])

    # screening poles s = N(d-1), N = 1..8
    for k in ncol:
        ax.plot(d, col[k], color=st.GRAY_LIGHT, lw=st.LW['thin'], zorder=1)
    for N, (x, y) in {1: (4.62, 3.22), 2: (4.36, 6.30), 3: (3.70, 7.50)}.items():
        ax.text(x, y, rf'$N\!=\!{N}$', color='0.55', fontsize=st.FONT['annotation'] - 0.5,
                ha='left', va='center', zorder=2)

    # kernel poles
    ax.plot(d, col['s_insertion'], color=st.BLUE, lw=st.LW['main'], zorder=3,
            label=r'polarization insertion, $s=2$')
    ax.plot(d, col['s_rings'], color=st.ORANGE, lw=st.LW['main'], zorder=3,
            label=r'ring kernel, $s=d+1$')
    ax.plot(d, col['s_q2_insertion'], color=st.BLUE, lw=st.LW['secondary'], ls=(0, (3.7, 1.6)),
            zorder=3, label=r'insertion of order $q^2$, $s=4$')

    # coincidences: open circles at d = 2, 3; filled dots at the noninteger
    # resonant dimensions of the solid lines inside the plotted range
    for r in read('fig1_coincidences.csv'):
        x, s = float(r['d']), float(r['s'])
        if not (XLIM[0] <= x < XLIM[1]):
            continue
        if r['integer_d'] == 'True':
            ax.plot(x, s, 'o', ms=7.0, mfc='white', mec='black', mew=1.1, zorder=5)
        elif r['kernel_pole'] != 'q^2 insertion s=4':
            ax.plot(x, s, 'o', ms=4.2, mfc='0.35', mec='0.35', zorder=5)

    box = dict(fc='white', ec='none', pad=0.6)
    fa = st.FONT['annotation']
    ax.text(3.10, 2.12, r'$C$, exchange part', fontsize=fa, ha='left', va='bottom', bbox=box, zorder=4)
    ax.text(3.11, 3.84, r'$A\,\log r_s$, $C$ (Lindhard part)', fontsize=fa, ha='left', va='top', bbox=box, zorder=4)
    ax.text(2.10, 1.86, r'$a_{2\ell}\,r_s^2\log r_s$', fontsize=fa, ha='left', va='top', bbox=box, zorder=4)
    ax.text(2.10, 2.86, r'$a_\ell\,r_s\log r_s$', fontsize=fa, ha='left', va='top', bbox=box, zorder=4)
    ax.text(2.10, 4.14, r'$r_s^3\log r_s$', fontsize=fa, ha='left', va='bottom', bbox=box, zorder=4)

    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_xticks(np.arange(1.5, 5.01, 0.5))
    ax.set_yticks(range(0, 9))
    ax.set_xlabel(r'dimension $d$')
    ax.set_ylabel(r'pole position $s$')
    leg = ax.legend(loc='lower right', frameon=True, framealpha=1.0, edgecolor='none',
                    borderaxespad=0.3, handlelength=2.2)
    leg.set_zorder(6)
    return st.save(fig, NAME)


if __name__ == '__main__':
    print(main())
