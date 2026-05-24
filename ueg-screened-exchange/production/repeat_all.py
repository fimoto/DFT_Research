#!/usr/bin/env python3
"""Fresh integral evaluation of every PRB record, with exact numerical checks.

Historical comparisons permit only audited source-metadata revisions. Separate
fresh passes must still agree in complete bytes, including current provenance.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import runtime
if __name__ == '__main__': runtime.enter_cli()
from common import ROOT, config, key, save, check_reference


def jobs():
    for dataset, settings in config().items():
        if not isinstance(settings, dict) or 'rs' not in settings: continue
        for rs in settings['rs']:
            for z in settings.get('zeta', [0.]):
                yield dataset, rs, z


def evaluate(job, destination, prepare=False):
    dataset, rs, z = job
    cmd=[sys.executable, str(ROOT/'production/compute.py'), '--dataset', dataset,
         '--rs', str(rs), '--output', str(destination), '--force']
    if dataset=='map': cmd+=['--zeta',str(z)]
    if prepare: cmd+=['--prepare-release']
    r=subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode: raise RuntimeError(' '.join(cmd)+'\n'+r.stdout+r.stderr)
    return r.stdout.strip()


def run(destination, workers, prepare=False):
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(evaluate, job, destination, prepare) for job in jobs()]
        for i,f in enumerate(as_completed(futures),1):
            print(f'[{i}/{len(futures)}] {f.result()}',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--jobs',type=int,default=2)
    p.add_argument('--passes',type=int,default=2)
    p.add_argument('--report',type=Path)
    a=p.parse_args()
    if a.jobs<1 or a.passes<1:p.error('jobs and passes must be positive')
    hashes={}
    for n in range(a.passes):
        with tempfile.TemporaryDirectory(prefix='prb-fresh-') as td:
            destination=Path(td);run(destination,a.jobs)
            current={}
            for job in jobs():
                dataset,rs,z=job;rel=Path(dataset)/(key(rs,z)+'.json')
                content=(destination/rel).read_bytes()
                check_reference(rel,json.loads(content))
                current[str(rel)]=hashlib.sha256(content).hexdigest()
            if hashes and current!=hashes:raise RuntimeError('Passes differ')
            hashes=current
        print(f'PASS {n+1}: {len(hashes)} exact historical numerical payloads; complete fresh-file bytes agree across passes.',flush=True)
    report={'passes':a.passes,'records_per_pass':len(hashes),'comparison':'exact numerical payload, settings and non-source metadata against historical records; only registered source revisions allowed; complete bytes across fresh passes',
            'runtime':runtime.contract(),'sha256':hashes,'status':'PASS'}
    if a.report:save(a.report,report)


if __name__=='__main__':main()
