#!/usr/bin/env python3
"""Inspect current Appendix B measurements. --recompute reruns both integrals.
No historical direct-column or grid-error constants are printed.
"""
import argparse, subprocess, sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--recompute',action='store_true');args=p.parse_args()
if args.recompute:
    for name in ['grid','tableIII_direct']:
        subprocess.run([sys.executable,str(root/'production/compute.py'),'--dataset',name,'--force'],check=True)
sys.path.insert(0,str(root/'production'))
from render import summarize
s=summarize()
print('Current recorded integral results' if not args.recompute else 'Freshly recomputed integral results')
for g,d in zip(s['grid'],s['tableIII_direct']):
    print('rs =',g['rs'])
    for name in ['published','extended']:
        v=g[name]
        print(name, v['mean_mry'], '+/-',v['se_mry'],'mRy; bare discrepancy',v['bare_discrepancy_mry'],'+/-',v['bare_se_mry'],'mRy; frequency shift',v['frequency_shift_mry'])
    print('direct',d['mean_mry'],'+/-',d['se_mry'],'mRy')
