"""Fixed runtime settings for the PRB entry points, before NumPy is imported."""
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
FIXED_ENV = {'OPENBLAS_NUM_THREADS':'1', 'OMP_NUM_THREADS':'1',
             'MKL_NUM_THREADS':'1', 'BLIS_NUM_THREADS':'1',
             'VECLIB_MAXIMUM_THREADS':'1', 'NUMEXPR_NUM_THREADS':'1',
             'PYTHONHASHSEED':'0', 'SOURCE_DATE_EPOCH':'1789344000',
             'TZ':'UTC'}

_NEEDS_REEXEC = any(os.environ.get(k) != v for k, v in FIXED_ENV.items())

def configure():
    for name,value in FIXED_ENV.items(): os.environ[name]=value


def contract():
    return json.loads((ROOT/'production/runtime_lock.json').read_text())


def validate():
    expected=contract()
    actual={'python':platform.python_version(),
            'packages':{name:importlib.metadata.version(name) for name in expected['packages']}}
    if actual!=expected:
        raise RuntimeError('Runtime differs from production/runtime_lock.json; no result was written. '
                           f'Expected {expected}; found {actual}')
    # A notebook may have loaded BLAS before importing the production module.
    # Keep the controller alive for the lifetime of this module.
    from threadpoolctl import threadpool_limits
    global _controller
    _controller=threadpool_limits(limits=1)


def enter_cli():
    # Python's hash seed is only read at interpreter startup.
    restart=_NEEDS_REEXEC
    configure()
    if restart: os.execv(sys.executable,[sys.executable,*sys.argv])
    validate()


configure()
