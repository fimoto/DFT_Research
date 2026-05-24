"""Build figs/phi_grid_corrected.npy from the supplied grid and the recomputed ranges.

The supplied grid figs/phi_grid.npy (900 points) is kept unchanged. Every point that
has a recomputed value is replaced:
  * 0 < y < 4    figs/phi_recomputed_y_lt_4.json      (scripts/recompute_phi_range.py)
  * 4 < y < 60   figs/phi_fig1_range_recomputed.json  (scripts/recompute_phi_fig1_range.py)
  * y > 60       figs/phi_recomputed_y_gt_60.json     (scripts/recompute_phi_range.py)
Every point is recomputed; the abscissae are not changed.
A range file must be complete; a partial one is refused.

    python scripts/build_phi_grid_corrected.py
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, '..', 'figs')
RANGES = [('phi_recomputed_y_lt_4.json', 0.0, 4.0), ('phi_fig1_range_recomputed.json', 4.0, 60.0),
          ('phi_recomputed_y_gt_60.json', 60.0, np.inf)]


def main():
    ys, ph = np.load(os.path.join(FIGS, 'phi_grid.npy'))
    new = ph.copy(); replaced = np.zeros(len(ys), bool)
    for name, lo, hi in RANGES:
        pts = json.load(open(os.path.join(FIGS, name)))['points']
        idx = [int(i) for i in np.nonzero((ys > lo) & (ys < hi))[0]]
        missing = [i for i in idx if str(i) not in pts]
        if missing:
            raise SystemExit(f'{name}: {len(missing)} of {len(idx)} points missing; finish that run first')
        for i in idx:
            p = pts[str(i)]
            assert p['y'] == float(ys[i]) and p['shipped'] == float(ph[i]), (name, i)
            new[i] = p['phi']; replaced[i] = True
        d = np.abs(new[idx] - ph[idx])
        print(f'{name}: {len(idx)} points, max |recomputed - supplied| {d.max():.3e}, '
              f'max relative {np.max(d / np.maximum(np.abs(new[idx]), 1e-300)):.3e}')
    np.save(os.path.join(FIGS, 'phi_grid_corrected.npy'), np.stack([ys, new]))
    print(f'wrote figs/phi_grid_corrected.npy: {replaced.sum()} of {len(ys)} points recomputed, '
          f'{(~replaced).sum()} supplied')


if __name__ == '__main__':
    main()
