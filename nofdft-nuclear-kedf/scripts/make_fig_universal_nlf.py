#!/usr/bin/env python3
"""
Regenerate the shared-multi-nucleus-fit NLF / density comparison figure
(manuscript Fig. 4) from public raw data only:

  results/reference/ks_nlf_<tag>.dat   (KS: r, rho, tau_TF, tau_vW, tau_KS, NLF_KS)
  results/raw/v5_<tag>_d.dat           (final OFDFT species density)
  results/raw/v5nlf_<tag>.dat          (species-resolved OFDFT NLF; col 7 = NLF)

Top row: species-resolved NLF(r) for the shared multi-nucleus fit vs the KS NLF.
Bottom row: total density rho_tot(r)=2 rho_q(r) for the final OFDFT density vs KS.
A=140 (rightmost) was not used in fitting or model selection. The displayed
d_NLF is the density-weighted distance used in paper Table III, recomputed here
so the figure and table use the same number.
"""
import numpy as np, os, argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TAGS = ['O16', 'Ca40', 'A80', 'A140']
TITLE = {'O16': r'$^{16}$O', 'Ca40': r'$^{40}$Ca', 'A80': r'$A=80$',
         'A140': r'$A=140$ (not used in fit)'}
RCUT = {'O16': 7.0, 'Ca40': 7.0, 'A80': 9.0, 'A140': 9.0}


def dnlf(r, model, ks_r, ks_nlf, ks_rho, rcut):
    ksn = np.interp(r, ks_r, ks_nlf)
    rho = np.interp(r, ks_r, ks_rho)
    w = r**2 * rho
    m = r < rcut
    return float(np.sqrt(np.sum(w[m] * (model[m] - ksn[m])**2) / np.sum(w[m])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--raw', required=True)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--out', default='fig4.png')
    args = ap.parse_args()

    fig, ax = plt.subplots(2, 4, figsize=(13.5, 6.0), sharex='col')
    for j, t in enumerate(TAGS):
        k = np.loadtxt(os.path.join(args.ref, 'ks_nlf_%s.dat' % t))
        d = np.loadtxt(os.path.join(args.raw, 'v5_%s_d.dat' % t))
        n = np.loadtxt(os.path.join(args.raw, 'v5nlf_%s.dat' % t))
        r, model = n[:, 0], n[:, 6]
        dd = dnlf(r, model, k[:, 0], k[:, 5], k[:, 1], RCUT[t])

        a = ax[0, j]
        ks_line, = a.plot(k[:, 0], k[:, 5], 'k-', lw=1.4, label='KS reference')
        shared_line, = a.plot(r, model, 'm-', lw=1.8, label='shared multi-nucleus fit')
        a.set_xlim(0, RCUT[t]); a.set_ylim(0, 1.05)
        a.set_title(TITLE[t])
        a.text(0.95, 0.06, r'$d_{\rm NLF}=%.3f$' % dd, transform=a.transAxes,
               ha='right', va='bottom', fontsize=10)
        if j == 0:
            a.set_ylabel('species NLF')

        b = ax[1, j]
        b.plot(k[:, 0], 2 * k[:, 1], 'k-', lw=1.4)          # total density = 2*species (N=Z)
        b.plot(d[:, 0], 2 * d[:, 1], 'm-', lw=1.8)
        b.set_xlim(0, RCUT[t]); b.set_xlabel('r [fm]')
        if j == 0:
            b.set_ylabel(r'$\rho_{\rm tot}$ [fm$^{-3}$]')

    fig.legend(
        [ks_line, shared_line], ['KS reference', 'shared multi-nucleus fit'],
        loc='lower center', bbox_to_anchor=(0.5, 0.005), ncol=2,
        fontsize=9, frameon=False,
    )
    fig.tight_layout(rect=(0, 0.065, 1, 1))
    fig.savefig(args.out, dpi=150)
    print('wrote', args.out)


if __name__ == '__main__':
    main()
