"""Shared configuration, units and traceable numerical records."""
import hashlib
import json
import platform
import sys
from pathlib import Path
import runtime
runtime.configure()
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from ueg_exchange.constants import Q_TF_SQ_COEFF, A1, A2, A3, E2X, q_star, asymptotic_kite


def config():
    return json.loads((ROOT / 'production/config.json').read_text())


def source_hashes():
    paths = ['production/config.json', 'production/common.py', 'production/compute.py',
             'production/runtime.py', 'production/runtime_lock.json',
             'external/benites_direct.py', 'src/ueg_exchange/lindhard.py',
             'src/ueg_exchange/constants.py', 'spin/reviewer_numerics.py',
             'external/kite-numerics/mc/xc_mc.py', 'external/kite-numerics/analytic/lindhard.py',
             'data/appendixB_grid_published.dat', 'data/appendixB_grid_extended.dat']
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}


def metadata():
    runtime.validate()
    return {'sources': source_hashes(), 'runtime': runtime.contract(),
            'units': 'energies Ry unless key ends in mry; momenta k_F'}


def encoded(payload):
    return (json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + '\n').encode('utf-8')


def check_source_compatibility(recorded_sources):
    """Accept exact sources or an explicitly audited, hash-pinned revision.

    Historical raw records retain the hashes of the code that produced them.
    Path keys use the public package name; their recorded hash values still
    identify the originally executed sources. Compatibility never permits an
    unknown source revision. The register documents each accepted revision.
    """
    current = source_hashes()
    if recorded_sources == current:
        return
    register = ROOT / 'production/source_compatibility.json'
    if register.exists():
        for revision in json.loads(register.read_text()).get('revisions', []):
            if (recorded_sources == revision['recorded_sources']
                    and current == revision['current_sources']):
                return
    raise RuntimeError('Numerical source hashes differ without an audited compatibility entry')


def check_reference(relative_path, payload):
    lock = ROOT / 'production/reference_lock.json'
    if not lock.exists():
        raise RuntimeError('Missing reference_lock.json; refusing to publish an unverified result')
    expected = json.loads(lock.read_text())['sha256'].get(str(relative_path))
    actual = hashlib.sha256(encoded(payload)).hexdigest()
    if expected == actual:
        return
    # A fresh run uses current, truthful source hashes. If the only difference
    # from a locked historical record is an audited source revision, retain
    # those fresh hashes while requiring every numerical value, setting and
    # other metadata field to be exactly equal to the historical reference.
    historical_path = ROOT / 'production/raw' / relative_path
    if historical_path.is_file():
        historical = json.loads(historical_path.read_text())
        if hashlib.sha256(encoded(historical)).hexdigest() == expected:
            check_source_compatibility(historical['metadata']['sources'])
            if payload.get('metadata', {}).get('sources') == source_hashes():
                normalized = json.loads(json.dumps(payload))
                normalized['metadata']['sources'] = historical['metadata']['sources']
                # Translate only the two renamed frequency-grid keys at the
                # archived-data boundary; settings values remain strictly locked.
                historical_settings = historical.get('settings', {})
                settings = normalized.get('settings', {})
                for new, old in (('nnu', 'nxi'), ('refined_nnu', 'refined_nxi')):
                    if old in historical_settings and new in settings and old not in settings:
                        settings[old] = settings.pop(new)
                if hashlib.sha256(encoded(normalized)).hexdigest() == expected:
                    return
    raise RuntimeError(f'Recomputed result differs from the released reference: {relative_path}. '
                       'Existing results were not replaced.')


def save(path, payload):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_bytes(encoded(payload))
    tmp.replace(path)


def stats(values):
    a = np.asarray(values, dtype=float)
    return float(a.mean()), float(a.std(ddof=1) / np.sqrt(len(a)))


def weights(q):
    q = np.asarray(q); w = np.zeros_like(q)
    w[1:] += np.diff(q) / 2; w[:-1] += np.diff(q) / 2
    w[0] += q[0] / 2  # linear segment from X(0)=0
    return w


def key(rs, zeta=0):
    return f'rs_{rs:g}_z_{zeta:g}'
