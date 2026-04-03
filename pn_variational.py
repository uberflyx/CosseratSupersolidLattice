"""
This file has moved to foundations/pn_variational.py

All imports are re-exported for backwards compatibility.
Update references to: foundations/pn_variational.py
"""
import importlib.util, os, sys
_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'foundations', 'pn_variational.py')
_spec = importlib.util.spec_from_file_location('pn_variational', _path,
           submodule_search_locations=[os.path.dirname(_path)])
_mod = importlib.util.module_from_spec(_spec)
sys.modules['pn_variational'] = _mod
_spec.loader.exec_module(_mod)
from pn_variational import *
