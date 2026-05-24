#!/usr/bin/env python3
"""
Deterministic regeneration of manuscript Tables II and III from public raw data.

Reads ONLY:
  results/reference/ks_reference_energies.dat   (KS-SkP binding energies)
  results/reference/ks_nlf_<tag>.dat            (KS tau/NLF reference)
  <raw>/dens_pernuc_<tag>.dat, nlf_pernuc_<tag>.dat, pnE_<tag>.dat  (per-nucleus)
  <raw>/hist_pernuc_<tag>.dat                    (per-nucleus ITE history)
  <raw>/v5_<tag>_d.dat, v5nlf_<tag>.dat, v5e64_<tag>.dat  (shared multi-nucleus)

Writes:
  tables/table1.tex, tables/table1.json   (per-nucleus)
  tables/table2.tex, tables/table2.json   (shared fit; A140 = unseen test)

No manuscript numbers are hard-coded here; every value is computed from the raw
files, so `make_tables.py` + `validate_tables.py` prove the manuscript is
reproducible from the released data.
"""
import numpy as np, json, os, re, argparse

TAGS = ['O16', 'Ca40', 'A80', 'A140']
PRETTY = {'O16': r'$^{16}$O', 'Ca40': r'$^{40}$Ca', 'A80': 'A80', 'A140': 'A140'}
RCUT = {'O16': 7.0, 'Ca40': 7.0, 'A80': 9.0, 'A140': 9.0}


def rrms(r, rho):
    return float(np.sqrt(np.sum(r**4 * rho) / np.sum(r**2 * rho)))


def read_e64(path):
    m = re.search(r'E_total_NeqZ_MeV\s+(-?\d+\.\d+E[+-]\d+)', open(path).read())
    return float(m.group(1))


def read_egrad(path):
    # E_grad_species_MeV is per-species; for N=Z the total (both species) is 2x,
    # consistent with E^{(64)} being the total energy and with the released
    # E_grad checksums (quoted as two-species totals).
    m = re.search(r'E_grad_species_MeV\s+(-?\d+\.\d+E[+-]\d+)', open(path).read())
    return 2.0 * float(m.group(1)) if m else 0.0


def nlf_distance(nlf_file, ks_file, rcut):
    d = np.loadtxt(nlf_file)          # col 2 = rho_species, col 7 = NLF
    k = np.loadtxt(ks_file)           # col 2 = rho, col 5 = tau_KS, col 6 = NLF_KS
    r = d[:, 0]; model = d[:, 6]
    ksn = np.interp(r, k[:, 0], k[:, 5])
    rho = np.interp(r, k[:, 0], k[:, 1])
    w = r**2 * rho
    m = r < rcut
    return float(np.sqrt(np.sum(w[m] * (model[m] - ksn[m])**2) / np.sum(w[m])))


def load_ks_ref(ref_dir):
    e = {}
    for line in open(os.path.join(ref_dir, 'ks_reference_energies.dat')):
        if line.startswith('#') or not line.strip():
            continue
        p = line.split()
        e[p[0]] = float(p[4])
    return e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--raw', required=True, help='directory with the raw production files')
    ap.add_argument('--ref', required=True, help='results/reference directory')
    ap.add_argument('--out', required=True, help='output tables directory')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    Eks = load_ks_ref(args.ref)
    ksnlf = {t: os.path.join(args.ref, 'ks_nlf_%s.dat' % t) for t in TAGS}

    # ---------- Paper Table II: per-nucleus ----------
    t1 = {}
    for t in TAGS:
        d = np.loadtxt(os.path.join(args.raw, 'dens_pernuc_%s.dat' % t))
        k = np.loadtxt(ksnlf[t])
        h = np.loadtxt(os.path.join(args.raw, 'hist_pernuc_%s.dat' % t))
        # Use the last ACCEPTED monitor row (col index 14 = accepted flag). A
        # rolled-back (rejected) trial is still written to the history with the
        # pre-rollback energy/residual, so the final row can be a rejected step
        # for the heavy nuclei; taking it would report a rejected E/R_EL.
        if h.ndim == 2 and h.shape[1] > 14:
            # a rolled-back trial is written to the history but is neither an
            # accepted nor a *monitored* point; require both flags so the row
            # is a scheduled monitor point that passed the energy-acceptance test.
            acc = h[(h[:, 13] == 1.0) & (h[:, 14] == 1.0)]
            hend = acc[-1] if len(acc) else h[-1]
            hbeg = acc[0] if len(acc) else h[0]
        else:
            hend, hbeg = h[-1], h[0]
        E64 = read_e64(os.path.join(args.raw, 'pnE_%s.dat' % t))
        Egrad = read_egrad(os.path.join(args.raw, 'pnE_%s.dat' % t))
        dd = nlf_distance(os.path.join(args.raw, 'nlf_pernuc_%s.dat' % t), ksnlf[t], RCUT[t])
        t1[t] = dict(
            E_KS=Eks[t], E64=E64, dE_abs=abs(hend[4] - hbeg[4]), R_EL=float(hend[11]),
            rho_c=2 * d[0, 1], rho_c_KS=2 * k[0, 1],
            rrms=rrms(d[:, 0], d[:, 1]), rrms_KS=rrms(k[:, 0], k[:, 1]),
            dNLF=dd)

    # ---------- Paper Table III: shared fit (A140 was not used in fitting) ----------
    t2 = {}
    for t in TAGS:
        d = np.loadtxt(os.path.join(args.raw, 'v5_%s_d.dat' % t))
        k = np.loadtxt(ksnlf[t])
        E64 = read_e64(os.path.join(args.raw, 'v5e64_%s.dat' % t))
        Egrad = read_egrad(os.path.join(args.raw, 'v5e64_%s.dat' % t))
        dd = nlf_distance(os.path.join(args.raw, 'v5nlf_%s.dat' % t), ksnlf[t], RCUT[t])
        t2[t] = dict(
            E_KS=Eks[t], E64=E64, dE_ref=E64 - Eks[t], E_grad=Egrad,
            rho_c=2 * d[0, 1], rho_c_KS=2 * k[0, 1],
            rrms=rrms(d[:, 0], d[:, 1]), rrms_KS=rrms(k[:, 0], k[:, 1]),
            dNLF=dd, holdout=(t == 'A140'))

    json.dump(t1, open(os.path.join(args.out, 'table1.json'), 'w'), indent=2)
    json.dump(t2, open(os.path.join(args.out, 'table2.json'), 'w'), indent=2)

    # ---------- TeX (rows only; the manuscript wraps them) ----------
    with open(os.path.join(args.out, 'table1.tex'), 'w') as f:
        f.write('% Auto-generated by make_tables.py -- do not edit by hand.\n')
        f.write('% Per-nucleus dd-OFDFT: E^{(64)}, |Delta E|, R_EL, rho_c, r_rms, d_NLF.\n')
        for t in TAGS:
            x = t1[t]
            f.write('%s & $%.3f$ & $%.3f$ & $%.3f$ & $%.4f$ & $%.4f\\,(%.4f)$ & $%.3f\\,(%.3f)$ & $%.4f$ \\\\\n'
                    % (PRETTY[t], x['E_KS'], x['E64'], x['dE_abs'], x['R_EL'],
                       x['rho_c'], x['rho_c_KS'], x['rrms'], x['rrms_KS'], x['dNLF']))
    with open(os.path.join(args.out, 'table2.tex'), 'w') as f:
        f.write('% Auto-generated by make_tables.py -- do not edit by hand.\n')
        f.write('% Shared multi-nucleus dd-OFDFT (A140 was not used in fitting): E^{(64)}, Delta E_ref, E_grad, rho_c, r_rms, d_NLF.\n')
        for t in TAGS:
            x = t2[t]
            tag = PRETTY[t] + (r'$^{\dagger}$' if x['holdout'] else '')
            f.write('%s & $%.2f$ & $%+.1f$ & $%.2f$ & $%.4f\\,(%.4f)$ & $%.3f\\,(%.3f)$ & $%.4f$ \\\\\n'
                    % (tag, x['E64'], x['dE_ref'], x['E_grad'],
                       x['rho_c'], x['rho_c_KS'], x['rrms'], x['rrms_KS'], x['dNLF']))

    print('Paper Table II  (per-nucleus):')
    print('  %-6s %10s %10s %8s %8s %8s' % ('nuc', 'E_KS', 'E64', 'R_EL', 'rho_c', 'dNLF'))
    for t in TAGS:
        x = t1[t]
        print('  %-6s %10.3f %10.3f %8.4f %8.4f %8.4f' % (t, x['E_KS'], x['E64'], x['R_EL'], x['rho_c'], x['dNLF']))
    print('Paper Table III (shared multi-nucleus fit; A140 was not used in fitting):')
    print('  %-6s %10s %10s %8s %8s %8s' % ('nuc', 'E64', 'dE_ref', 'E_grad', 'rho_c', 'dNLF'))
    for t in TAGS:
        x = t2[t]
        print('  %-6s %10.2f %+10.1f %8.2f %8.4f %8.4f' % (t, x['E64'], x['dE_ref'], x['E_grad'], x['rho_c'], x['dNLF']))
    print('\nWrote table1.{tex,json}, table2.{tex,json} to', args.out)


if __name__ == '__main__':
    main()
