#!/usr/bin/env python3
"""Compatibility entry point: redraw the current PRB release from production/raw.
Use production/compute.py explicitly to regenerate the underlying integrals.
"""
from pathlib import Path
import subprocess, sys
root=Path(__file__).resolve().parents[1]
subprocess.run([sys.executable,str(root/'production/render.py')],check=True)
