#!/usr/bin/env python3
"""Recompute an explicitly selected PRB dataset, with bounded sample memory.

Run each dataset separately; see REPRODUCE_SPIN.md for the sequence.
Existing current records are skipped; --force recomputes them from the integral.
--output permits an independent rerun without touching the released records.
Verified recomputation into production/raw preserves the historical locked file;
use a separate --output to save a fresh record with current source provenance.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import runtime
if __name__ == "__main__": runtime.enter_cli()
import numpy as np
from scipy.stats import qmc
from common import ROOT, config, metadata, save, stats, q_star, E2X, key, check_reference, check_source_compatibility


def direct(rs, zeta, settings, cfg):
    sys.path.insert(0, str(ROOT / 'external'))
    import benites_direct as bd
    m, reps = settings['m'], settings['replicates']
    scale = float(q_star(rs)) if settings['q_scale'] == 'qstar' else settings['q_scale']
    chunk_power = min(cfg['chunk_power'], m)
    terms = []
    for rep in range(reps):
        row = []
        for term in (1, 2):
            sob = qmc.Sobol(9, scramble=True, seed=cfg['seed'] + 2 * rep + term - 1)
            partial = []
            for _ in range(2 ** (m - chunk_power)):
                points = sob.random(2 ** chunk_power)
                partial.append(float(np.sum(bd.evaluate_term(points, rs, zeta, term, scale))))
            row.append(float(np.sum(partial) / 2 ** m))
        terms.append(row)
    a = np.asarray(terms)
    mean = E2X + a.mean(axis=0).sum()
    se = np.sqrt(np.sum(a.var(axis=0, ddof=1)) / reps)
    return {'rs': rs, 'zeta': zeta, 'settings': dict(m=m, replicates=reps, q_scale=scale,
            seed=cfg['seed'], chunk_power=chunk_power, generator='scipy scrambled Sobol, d=9'),
            'term_replicates_ry': terms, 'mean_mry': float(1000 * mean), 'se_mry': float(1000 * se),
            'error_definition': 'sqrt((sample_var(term1)+sample_var(term2))/R); independent term scrambles'}


def shell(rs, settings):
    sys.path.insert(0, str(ROOT / 'spin'))
    import reviewer_numerics as rn
    x = np.geomspace(settings['x_min'], settings['x_max'], settings['nx'])
    q = q_star(rs) * x
    values = np.empty((settings['replicates'], len(q)))
    refined = np.empty_like(values)
    for rep in range(settings['replicates']):
        u = qmc.Sobol(7, scramble=True, seed=settings['seed'] + rep).random_base2(settings['power'])
        for i, v in enumerate(q):
            for dest, order in [(values, settings['order']), (refined, settings['refined_order'])]:
                dest[rep, i] = 6 / np.pi**3 * rn.node_estimate(v, rs, u, order=order, nk=settings['nk'])
    return {'rs': rs, 'settings': settings, 'x': x.tolist(), 'q': q.tolist(),
            'density_replicates_ry': values.tolist(), 'refined_density_replicates_ry': refined.tolist(),
            'error_definition': 'sample SD of full density replicate means / sqrt(R); frequency refinement separate'}


def grid(rs, settings):
    sys.path.insert(0, str(ROOT / 'external/kite-numerics/mc'))
    import xc_mc as xc
    q = np.loadtxt(ROOT / 'data/appendixB_grid_extended.dat')
    shape = (settings['replicates'], len(q))
    bare, correction, refined = np.empty(shape), np.empty(shape), np.empty(shape)
    n = settings['samples_per_node_per_replicate']
    for i, v in enumerate(q):
        # External literature-reference kernel retains its original nxi keyword.
        kernels = [xc.make_dyn_kernel(v, rs, nxi=order, nseg=settings['nseg'])
                   for order in (settings['nnu'], settings['refined_nnu'])]
        for rep in range(settings['replicates']):
            seed = settings['seed'] + 10000 * rep + i
            bare[rep, i] = xc.xc_samples(v, n, np.random.default_rng(seed), kernel='static').mean()
            for dest, kernel in [(correction, kernels[0]), (refined, kernels[1])]:
                cache = {(v, rs, 0., 1.): kernel}
                dest[rep, i] = xc.xc_samples(v, n, np.random.default_rng(seed), kernel='diff', rs=rs,
                                            _kern_cache=cache).mean()
    return {'rs': rs, 'settings': settings, 'q': q.tolist(), 'bare_X_replicates': bare.tolist(),
            'correction_X_replicates': correction.tolist(), 'refined_correction_X_replicates': refined.tolist(),
            'error_definition': 'integrate each replicate on each q grid, then sample SD/sqrt(R)'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dataset', required=True, choices=['tableI','tableII','tableIII_direct','tableIV','map','shell','grid'])
    p.add_argument('--rs', type=float, nargs='+', help='Optional subset of configured density nodes')
    p.add_argument('--zeta', type=float, nargs='+', help='Optional subset of map polarizations')
    p.add_argument('--output', type=Path, default=ROOT / 'production/raw')
    p.add_argument('--force', action='store_true')
    p.add_argument('--prepare-release', action='store_true', help='Maintainer: write unsealed candidates only to a separate output directory')
    args = p.parse_args(); cfg = config(); settings = cfg[args.dataset]
    rsv = args.rs if args.rs is not None else settings['rs']
    zv = (args.zeta if args.zeta is not None else cfg['map']['zeta']) if args.dataset == 'map' else [0.]
    if not set(rsv).issubset(settings['rs']): p.error('rs must be a configured node')
    if args.dataset == 'map' and not set(zv).issubset(cfg['map']['zeta']): p.error('zeta must be a configured node')
    if args.prepare_release and args.output.resolve() == (ROOT/'production/raw').resolve():
        p.error('--prepare-release requires a separate --output; released data cannot be overwritten')
    meta = metadata()
    for rs_input in rsv:
        # Use the exact configured representation even when CLI parsing gives 5.0.
        rs = next(v for v in settings['rs'] if v == rs_input)
        for zeta_input in zv:
            zeta = next(v for v in cfg['map']['zeta'] if v == zeta_input) if args.dataset == 'map' else 0.0
            path = args.output / args.dataset / (key(rs,zeta) + '.json')
            if path.exists() and not args.force:
                old = json.loads(path.read_text())
                check_source_compatibility(old['metadata']['sources'])
                check_reference(Path(args.dataset)/path.name, old)
                print('CURRENT', args.dataset, rs, zeta, flush=True); continue
            start = time.perf_counter()
            if args.dataset == 'shell': result = shell(rs, settings)
            elif args.dataset == 'grid': result = grid(rs, settings)
            else: result = direct(rs, zeta, settings, cfg)
            result.update(dataset=args.dataset, metadata=meta)
            if not args.prepare_release: check_reference(Path(args.dataset)/path.name, result)
            preserve = (args.output.resolve() == (ROOT/'production/raw').resolve()
                        and path.exists() and not args.prepare_release)
            if not preserve: save(path, result)
            print('VERIFIED (historical record preserved)' if preserve else 'DONE',
                  args.dataset, rs, zeta, f"{time.perf_counter()-start:.1f}s",
                  result.get('mean_mry',''), result.get('se_mry',''), flush=True)


if __name__ == '__main__': main()
