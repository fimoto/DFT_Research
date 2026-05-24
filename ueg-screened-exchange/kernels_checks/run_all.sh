#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
python check_response_asymptotics.py > response_asymptotics.log
python check_exchange_wick_counts.py > wick_counts.txt
python verify_exchange_chain_coefficients.py --max-order 4 --json exchange_chain_coefficients.json > exchange_chain_coefficients.log
python check_series_transfer.py > series_transfer.log
python normalization_checks.py > normalization_checks.log
python independent_checks.py > independent_checks.log
python - <<'PY'
import importlib.metadata
import json
import platform
from pathlib import Path
record = {'python': platform.python_version(),
          'packages': {p: importlib.metadata.version(p)
                       for p in ['numpy', 'scipy', 'mpmath', 'sympy']}}
Path('environment.json').write_text(json.dumps(record, indent=2)+'\n')
PY
printf '%s\n' 'All verification scripts completed.'
