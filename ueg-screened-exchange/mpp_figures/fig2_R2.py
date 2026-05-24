"""Main-text Fig. 2: R_2(z) and h_2(z), with R_3 and h_3 for comparison (a),
and R_2(z) on log-log axes with its static value and high-frequency tail (b).

Reads ../mpp_checks/figure_data/fig2_R2.csv and fig2_R2_loglog.csv
(written by ../mpp_checks/fig_data.py).
"""
import csv, os
import numpy as np
import matplotlib.pyplot as plt
import figstyle as st

NAME = 'fig2_R2'


def read(name):
    with open(os.path.join(st.DATA, name)) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def main():
    st.apply()
    a = read('fig2_R2.csv')
    b = read('fig2_R2_loglog.csv')

    fig = plt.figure(figsize=st.figsize(NAME, 0.40))
    axa = fig.add_axes([0.065, 0.175, 0.415, 0.795])
    axb = fig.add_axes([0.575, 0.175, 0.415, 0.795])

    # (a)
    m = a['z'] <= 3.0
    z = a['z'][m]
    axa.plot(z, a['h2'][m], color=st.GRAY_DARK, lw=st.LW['secondary'], label=r'$h_2(z)$')
    axa.plot(z, a['R2'][m], color=st.BLUE, lw=st.LW['main'], label=r'$R_2(z)$')
    axa.plot(z, a['h3'][m], color=st.GRAY_DARK, lw=1.0, ls=(0, (3.7, 1.6)), label=r'$h_3(z)=h(z)$')
    axa.plot(z, a['R3'][m], color=st.ORANGE, lw=st.LW['secondary'], ls=(0, (3.7, 1.6)), label=r'$R_3(z)$')
    r20 = 2 / np.pi
    axa.plot(0, r20, 'o', ms=5, mfc='white', mec=st.BLUE, mew=1.0, clip_on=False, zorder=5)
    axa.annotate(r'$R_2(0)=2/\pi$', xy=(0.04, r20 + 0.005), xytext=(0.42, 0.66),
                 fontsize=st.FONT['annotation'], ha='left', va='center',
                 arrowprops=dict(arrowstyle='-', color='0.3', lw=0.6, shrinkA=1, shrinkB=2))
    axa.set_xlim(0, 3)
    axa.set_ylim(0, 1.02)
    axa.set_xlabel(r'$z=\nu/q$')
    axa.legend(loc='upper right', handlelength=2.4)
    axa.text(0.97, 0.58, '(a)', transform=axa.transAxes, ha='right', va='center', fontsize=st.FONT['panel'])

    # (b)
    axb.loglog(b['z'], b['R2'], color=st.BLUE, lw=st.LW['main'], label=r'$R_2(z)$')
    axb.loglog(b['z'], b['tail_5_over_12pi_z4'], color='black', lw=st.LW['thin'], ls=(0, (1, 1.65)),
               label=r'$5/(12\pi z^4)$')
    axb.axhline(b['static_2_over_pi'][0], color='black', lw=st.LW['thin'], ls=(0, (3.7, 1.6)),
                label=r'$2/\pi$')
    axb.set_xlim(10 ** -2.2, 10 ** 2.2)
    axb.set_ylim(1e-9, 2.0)
    axb.set_yticks([1e-8, 1e-6, 1e-4, 1e-2, 1e0])
    axb.set_xlabel(r'$z=\nu/q$')
    axb.legend(loc='lower left', handlelength=2.4)
    axb.text(0.97, 0.58, '(b)', transform=axb.transAxes, ha='right', va='center', fontsize=st.FONT['panel'])
    return st.save(fig, NAME)


if __name__ == '__main__':
    print(main())
