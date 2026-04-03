"""
This file has moved to foundations/cosserat_transfer_matrix.py

All imports are re-exported for backwards compatibility.
Update references to: foundations/cosserat_transfer_matrix.py
"""
import importlib.util, os, sys
_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'foundations', 'cosserat_transfer_matrix.py')
_spec = importlib.util.spec_from_file_location('cosserat_transfer_matrix', _path,
           submodule_search_locations=[os.path.dirname(_path)])
_mod = importlib.util.module_from_spec(_spec)
sys.modules['cosserat_transfer_matrix'] = _mod
_spec.loader.exec_module(_mod)
from cosserat_transfer_matrix import *
