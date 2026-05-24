"""Run every mpp_checks script and write results/summary.json.

    python mpp_checks/run_all.py            # all checks, Monte Carlo 4 x 2.5e7 samples (about 1.5 h)
    python mpp_checks/run_all.py --quick    # skips c06, c07 and c12 (several minutes each) and uses 4 x 2e6 samples
    python mpp_checks/run_all.py --full-mc  # Monte Carlo with 4 x 1e8 samples, as quoted in SI Sec. S19
    python mpp_checks/run_all.py --skip=c12_static_dassarma   # leave out a script (repeatable)
"""
import subprocess, sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
quick = '--quick' in sys.argv; full_mc = '--full-mc' in sys.argv
scripts = ['c01_lindhard_d', 'c02_ring_d', 'c03_response_3d', 'c04_R2', 'c05_a2l', 'c06_finiteq_2d', 'c07_screened_line',
           'c08_spin_2d', 'c09_kite_3d', 'c10_hd_fit', 'c11_lowdensity_kite', 'c12_static_dassarma']
skip = [a.split('=', 1)[1] for a in sys.argv if a.startswith('--skip=')]
scripts = [s for s in scripts if s not in skip]
if quick:
    scripts = [s for s in scripts if s not in ('c06_finiteq_2d', 'c07_screened_line', 'c12_static_dassarma')]
summary = {}
t0 = time.time()
for s in scripts:
    args = [sys.executable, os.path.join(HERE, s + '.py')]
    if s == 'c09_kite_3d':
        args += ['--mc-samples', '1e8' if full_mc else ('2e6' if quick else '2.5e7'), '--mc-runs', '4']
    t = time.time(); rc = subprocess.call(args, cwd=HERE)
    res = json.load(open(os.path.join(HERE, 'results', s + '.json')))
    summary[s] = dict(exit_code=rc, n_checks=res['n_checks'], n_fail=res['n_fail'], seconds=round(time.time() - t, 1))
tot = sum(v['n_checks'] for v in summary.values()); bad = sum(v['n_fail'] for v in summary.values())
json.dump(dict(scripts=summary, n_checks=tot, n_fail=bad, mode='quick' if quick else ('full-mc' if full_mc else 'default'),
               seconds=round(time.time() - t0, 1)), open(os.path.join(HERE, 'results', 'summary.json'), 'w'), indent=1)
print(f'\n{tot} checks, {bad} failed')
sys.exit(1 if bad else 0)
