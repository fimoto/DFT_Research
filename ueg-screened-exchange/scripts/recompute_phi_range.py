"""Recompute the supplied Phi grid on one range of abscissae, resumably.

For every grid abscissa with YMIN < y < YMAX (or YMIN <= y for --ymax inf),
Phi(y) = sample_phi_grid.phi(y, nu, nv, add_tail=True) with the given u cut and,
for y <= pi, the given geometric growth of the u panels. Every --check-every-th
point, and the last one, is repeated at (nu_hi, nv_hi) and, if --check-umax-scale
is given, with that u cut; the run stops if a check differs by more than --tol.
Results are written point by point to --output, so an interrupted run resumes.
scripts/build_phi_grid_corrected.py merges the ranges into the corrected grid.

The y < 4 range distributed in figs/phi_recomputed_y_lt_4.json was produced by

    python scripts/recompute_phi_range.py --ymin 0 --ymax 4 --umax-scale 40 \\
        --check-umax-scale 80 --output figs/phi_recomputed_y_lt_4.json

(about an hour on one core; the settings are stored in the file and a resumed run
with other settings is refused). The 78 points above y = 60, in
figs/phi_recomputed_y_gt_60.json, were produced by

    python scripts/recompute_phi_range.py --ymin 60 --ymax inf --nu 24 --nv 16 \\
        --nu-hi 32 --nv-hi 20 --umax-scale 1 --check-umax-scale 2 --check-every 20 \\
        --tol 5e-9 --output figs/phi_recomputed_y_gt_60.json

(about four hours on one core). A u cut of 200 with the tail correction agrees with
a cut of 400 to 3e-10 at y = 64.8, and every 20th point and the last are repeated at
higher order and with the doubled cut.
"""
import argparse, json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sample_phi_grid as S


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--ymin', type=float, required=True)
    p.add_argument('--ymax', type=float, required=True)
    p.add_argument('--nu', type=int, default=32); p.add_argument('--nv', type=int, default=20)
    p.add_argument('--nu-hi', type=int, default=48); p.add_argument('--nv-hi', type=int, default=24)
    p.add_argument('--umax-scale', type=float, default=10.0)
    p.add_argument('--check-umax-scale', type=float, default=None)
    p.add_argument('--growth', type=float, default=0.1)
    p.add_argument('--v-panels-per-period', type=float, default=8.0)
    p.add_argument('--check-every', type=int, default=20)
    p.add_argument('--tol', type=float, default=1e-9)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    ys, ph = np.load(os.path.join(HERE, '..', 'figs', 'phi_grid.npy'))
    idx = [int(i) for i in np.nonzero((ys > a.ymin) & (ys < a.ymax))[0]]
    settings = dict(nu=a.nu, nv=a.nv, nu_hi=a.nu_hi, nv_hi=a.nv_hi, umax_scale=a.umax_scale,
                    check_umax_scale=a.check_umax_scale, growth=a.growth,
                    v_panels_per_period=a.v_panels_per_period,
                    check_every=a.check_every, tol=a.tol, ymin=a.ymin, ymax=a.ymax)
    rec = json.load(open(a.output)) if os.path.exists(a.output) else dict(settings=settings, points={})
    if rec['settings'].get('v_panels_per_period') is None and settings['v_panels_per_period'] == 8.0:
        rec['settings']['v_panels_per_period'] = 8.0   # records written before the option existed
    if rec['settings'] != settings:
        raise SystemExit(f'{a.output} was produced with other settings: {rec["settings"]}')
    S.GROWTH = a.growth
    S.V_PANELS_PER_PERIOD = a.v_panels_per_period
    for n, i in enumerate(idx):
        if str(i) in rec['points']:
            continue
        t0 = time.perf_counter()
        S.UMAX_SCALE = a.umax_scale
        v = S.phi(ys[i], nu=a.nu, nv=a.nv, add_tail=True)
        row = dict(y=float(ys[i]), shipped=float(ph[i]), phi=float(v))
        if n % a.check_every == 0 or i == idx[-1]:
            row['phi_order_check'] = float(S.phi(ys[i], nu=a.nu_hi, nv=a.nv_hi, add_tail=True))
            checks = [row['phi_order_check']]
            if a.check_umax_scale:
                S.UMAX_SCALE = a.check_umax_scale
                row['phi_cut_check'] = float(S.phi(ys[i], nu=a.nu, nv=a.nv, add_tail=True))
                checks.append(row['phi_cut_check'])
            worst = max(abs(c - v) for c in checks)
            if worst > a.tol:
                raise SystemExit(f'not converged at y={ys[i]}: {v} vs {checks}')
        row['seconds'] = round(time.perf_counter() - t0, 1)
        rec['points'][str(i)] = row
        with open(a.output + '.tmp', 'w') as f:
            json.dump(rec, f, indent=1, sort_keys=True); f.write('\n')
        os.replace(a.output + '.tmp', a.output)
        print(f'{len(rec["points"])}/{len(idx)} y={ys[i]:.6g} phi={v:+.12e} shipped={ph[i]:+.12e} '
              f'diff={v - ph[i]:+.2e} [{row["seconds"]}s]', flush=True)
    print('ALLDONE', flush=True)


if __name__ == '__main__':
    main()
