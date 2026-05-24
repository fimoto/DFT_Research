"""Supplemental Material Table S-grid: kite energy at r_s = 1e8 under three momentum grids.

Third implementation (a7_mc.py, hd_split.py); imports nothing from src/, production/
or pysecdec_kite/. Every replicate value is stored at full precision in
table_S_grid.json, point by point, so an interrupted run resumes and the table can be
audited replicate by replicate.

  row 1  single rational map q = s t/(1-t), s = q_*, 48 nodes      a7_mc.energy
  row 2  split [0,1] u [1,4] u [4,inf), 16/16/40 nodes            hd_split.energy_split
  row 3  same split, 24/24/64 nodes, frequency order 56 (vs 40)   hd_split.energy_split

8 replicates per row, N = 1000000 geometry samples per momentum node, seeds 5000+i,
6000+i, 7000+i. The error is the sample standard deviation of the replicates over
sqrt(8). The three-term expansion A1/q_* + A2/q_*^2 + A3/q_*^3 is printed for reference.

Why N = 1e6: at N = 250000 the replicate distribution has a heavy tail at this
density (one replicate of row 2 at -1.37 against a median of -2.10), and the
8-replicate standard error understates the scatter; row 3 then gave
-2.026 +- 0.034, while N = 1e6 gives values consistent with the expansion.

Two checks, recorded in the same file but not in the table, test the grid and the
frequency order of row 3 (N = 250000, 8 replicates): a finer momentum split
[0,1] u [1,4] u [4,20] u [20,100] u [100,inf) with 24/24/32/32/48 nodes, and
frequency order 80.

    python3 table_S_grid.py        # about 40 minutes on one core
"""
import json, os, time
from math import gamma, sqrt, pi
import numpy as np
from a7_mc import energy, aTF, E2X
from hd_split import energy_split

RS, REPS, N = 1e8, 8, 1_000_000
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'table_S_grid.json')
ROWS = [
    dict(label='single rational map, 48 nodes', kind='map', nq=48, nnu=40, seed0=5000),
    dict(label='split [0,1] u [1,4] u [4,inf), 16/16/40 nodes', kind='split',
         breaks=[0.0, 1.0, 4.0, None], nsegs=[16, 16, 40], nnu=40, seed0=6000),
    dict(label='split, 24/24/64 nodes, refined frequency order', kind='split',
         breaks=[0.0, 1.0, 4.0, None], nsegs=[24, 24, 64], nnu=56, seed0=7000),
    dict(label='check: finer split 0,1,4,20,100,inf, 24/24/32/32/48 nodes, order 56', kind='split',
         breaks=[0.0, 1.0, 4.0, 20.0, 100.0, None], nsegs=[24, 24, 32, 32, 48], nnu=56,
         seed0=9000, N=250_000, in_table=False),
    dict(label='check: split 24/24/64 nodes, frequency order 80', kind='split',
         breaks=[0.0, 1.0, 4.0, None], nsegs=[24, 24, 64], nnu=80, seed0=10000,
         N=250_000, in_table=False),
]


def one(row, i):
    seed = row['seed0'] + i
    n = row.get('N', N)
    # Previously recorded, resumable grids retain their historical nxi key.
    frequency_order = row['nnu'] if 'nnu' in row else row['nxi']
    if row['kind'] == 'map':
        d, _ = energy(RS, nq=row['nq'], nnu=frequency_order, N=n, seed=seed)   # scale defaults to q_*
        return 1000 * (E2X + d)
    return 1000 * energy_split(RS, tuple(row['breaks']), tuple(row['nsegs']),
                               nnu=frequency_order, N=n, seed=seed)


def main():
    qstar = (4 * aTF * RS / 3) ** 0.25
    A1 = -4 / (7 * sqrt(2) * pi ** 1.5) * gamma(0.25) / gamma(0.75)
    A2 = 8 / 15
    A3 = -16408 / 3375 * sqrt(2 * pi) / pi ** 2 * gamma(0.75) / gamma(0.25)
    three = 1000 * (A1 / qstar + A2 / qstar ** 2 + A3 / qstar ** 3)
    rec = json.load(open(OUT)) if os.path.exists(OUT) else dict(
        rs=RS, replicates=REPS, N=N, q_star=qstar, three_term_mry=three,
        error_definition='sample SD of replicate values / sqrt(replicates)',
        rows=[dict(r, values_mry=[]) for r in ROWS])
    # Resume check (audit F2, 1 Oct 2026): refuse to continue a record that was
    # produced under different physical or sampling settings.
    stored_rows = [{k: v for k, v in r.items() if k not in ('values_mry', 'mean_mry', 'se_mry')} for r in rec['rows']]
    if rec.get('rs') != RS or rec.get('replicates') != REPS or rec.get('N') != N or stored_rows != ROWS:
        raise SystemExit(f'{OUT} was produced with different settings; use a new output file')
    for row in rec['rows']:
        while len(row['values_mry']) < REPS:
            i = len(row['values_mry']); t0 = time.time()
            row['values_mry'].append(one(row, i))
            with open(OUT + '.tmp', 'w') as f:
                json.dump(rec, f, indent=1); f.write('\n')
            os.replace(OUT + '.tmp', OUT)
            print(f"{row['label']}: replicate {i} {row['values_mry'][-1]:.6f} mRy [{time.time()-t0:.0f}s]", flush=True)
        v = np.array(row['values_mry'])
        row['mean_mry'] = float(v.mean()); row['se_mry'] = float(v.std(ddof=1) / sqrt(len(v)))
    with open(OUT + '.tmp', 'w') as f:
        json.dump(rec, f, indent=1); f.write('\n')
    os.replace(OUT + '.tmp', OUT)
    print(f'\nr_s = 1e8, q_* = {qstar:.2f}; three-term expansion {three:.3f} mRy')
    for row in rec['rows']:
        tag = '' if row.get('in_table', True) else '   (check, N=%d)' % row['N']
        print(f"  {row['label']:<70} {row['mean_mry']:8.3f} +- {row['se_mry']:.3f} mRy"
              f"  {(row['mean_mry'] - three) / row['se_mry']:+.1f} sigma vs expansion{tag}")


if __name__ == '__main__':
    main()
