"""Thin shim: the single source of the direct BRM evaluator is external/benites_direct.py."""
import os, sys, importlib.util
_p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benites_direct.py")
_spec = importlib.util.spec_from_file_location("_benites_direct_src", _p)
_m = importlib.util.module_from_spec(_spec)
sys.modules["_benites_direct_src"] = _m
_spec.loader.exec_module(_m)
globals().update({k: v for k, v in vars(_m).items() if not k.startswith("__")})
del _p, _spec, _m
