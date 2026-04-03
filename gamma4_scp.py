"""
This file has moved to foundations/gamma4_scp.py

All imports are re-exported for backwards compatibility.
Update references to: foundations/gamma4_scp.py
"""
import importlib.util, os, sys
_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'foundations', 'gamma4_scp.py')
_spec = importlib.util.spec_from_file_location('gamma4_scp', _path,
           submodule_search_locations=[os.path.dirname(_path)])
_mod = importlib.util.module_from_spec(_spec)
sys.modules['gamma4_scp'] = _mod
_spec.loader.exec_module(_mod)
from gamma4_scp import *
