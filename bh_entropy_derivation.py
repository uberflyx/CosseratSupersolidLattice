"""
This file has moved to gravity/bh_entropy_derivation.py

All imports are re-exported for backwards compatibility.
Update references to: gravity/bh_entropy_derivation.py
"""
import importlib.util, os, sys
_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gravity', 'bh_entropy_derivation.py')
_spec = importlib.util.spec_from_file_location('bh_entropy_derivation', _path,
           submodule_search_locations=[os.path.dirname(_path)])
_mod = importlib.util.module_from_spec(_spec)
sys.modules['bh_entropy_derivation'] = _mod
_spec.loader.exec_module(_mod)
from bh_entropy_derivation import *
