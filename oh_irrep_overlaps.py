"""
This file has moved to foundations/oh_irrep_overlaps.py

All imports are re-exported for backwards compatibility.
Update references to: foundations/oh_irrep_overlaps.py
"""
import importlib.util, os, sys
_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'foundations', 'oh_irrep_overlaps.py')
_spec = importlib.util.spec_from_file_location('oh_irrep_overlaps', _path,
           submodule_search_locations=[os.path.dirname(_path)])
_mod = importlib.util.module_from_spec(_spec)
sys.modules['oh_irrep_overlaps'] = _mod
_spec.loader.exec_module(_mod)
from oh_irrep_overlaps import *
