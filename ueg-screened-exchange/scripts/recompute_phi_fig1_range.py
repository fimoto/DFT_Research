"""Recompute the shipped Phi grid on 4 < y < 60, the range of the right panel of the JPA Phi figure (current Fig. 5; the output filename retains its historical prefix).

Why: that panel plots y^2 Phi(y), which multiplies grid errors by y^2. With
the fully split panel rule of scripts/sample_phi_grid.py the shipped grid was
found to be off by 3.3e-5 (2.8 percent of Phi, 0.11 in y^2 Phi) at y = 57.3.

What: for every grid abscissa with 4 < y < 60, Phi(y) from
sample_phi_grid.phi(y, nu=NU_LO, nv=NV_LO, add_tail=True). Every CHECK_EVERY-th
point, and the largest y, is also evaluated at (NU_HI, NV_HI) as a convergence
check; the run fails if the two differ by more than TOL. Results are appended
point by point to figs/phi_fig1_range_recomputed.json, so an interrupted run
resumes where it stopped. scripts/build_phi_grid_corrected.py merges this range
and scripts/recompute_phi_range.py's y < 4 range into figs/phi_grid_corrected.npy.
(Until 2026-09-16 this script wrote that file itself, with only these points.)

    python scripts/recompute_phi_fig1_range.py      # about an hour on one core
"""
import json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sample_phi_grid as S

YMIN, YMAX = 4.0, 60.0
NU_LO, NV_LO = 32, 20
NU_HI, NV_HI = 48, 24
CHECK_EVERY = 10
TOL = 1e-9
FIGS = os.path.join(HERE, '..', 'figs')
OUT = os.path.join(FIGS, 'phi_fig1_range_recomputed.json')


def main():
    ys, ph = np.load(os.path.join(FIGS, 'phi_grid.npy'))
    idx = [int(i) for i in np.nonzero((ys > YMIN) & (ys < YMAX))[0]]
    rec = json.load(open(OUT)) if os.path.exists(OUT) else {
        'rule': f'sample_phi_grid.phi(y, nu={NU_LO}, nv={NV_LO}, add_tail=True), UMAX_SCALE={S.UMAX_SCALE}',
        'check_rule': f'nu={NU_HI}, nv={NV_HI} at every {CHECK_EVERY}th point and the largest y',
        'tolerance': TOL, 'points': {}}
    for n, i in enumerate(idx):
        if str(i) in rec['points']:
            continue
        t0 = time.perf_counter()
        v = S.phi(ys[i], nu=NU_LO, nv=NV_LO, add_tail=True)
        row = dict(y=float(ys[i]), shipped=float(ph[i]), phi=float(v))
        if n % CHECK_EVERY == 0 or i == idx[-1]:
            row['phi_check'] = float(S.phi(ys[i], nu=NU_HI, nv=NV_HI, add_tail=True))
            if abs(row['phi_check'] - v) > TOL:
                raise SystemExit(f'not converged at y={ys[i]}: {v} vs {row["phi_check"]}')
        row['seconds'] = round(time.perf_counter() - t0, 1)
        rec['points'][str(i)] = row
        tmp = OUT + '.tmp'
        with open(tmp, 'w') as f:
            json.dump(rec, f, indent=1, sort_keys=True)
            f.write('\n')
        os.replace(tmp, OUT)
        print(f'{len(rec["points"])}/{len(idx)} y={ys[i]:.4f} phi={v:+.12e} shipped={ph[i]:+.12e} '
              f'diff={v - ph[i]:+.2e} y2diff={ys[i]**2 * (v - ph[i]):+.2e} [{row["seconds"]}s]', flush=True)
    d = np.array([rec['points'][str(i)]['phi'] - ph[i] for i in idx])
    y = ys[idx]
    print(f'all {len(idx)} points present: max |diff| {np.max(abs(d)):.2e}, '
          f'max |y^2 diff| {np.max(abs(y**2 * d)):.2e} at y={y[np.argmax(abs(y**2*d))]:.2f}; '
          'run scripts/build_phi_grid_corrected.py to write the corrected grid')

if __name__ == '__main__':
    main()
