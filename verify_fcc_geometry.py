"""
This file has moved to hadrons/verify_fcc_geometry.py

All imports are re-exported for backwards compatibility.
Update references to: hadrons/verify_fcc_geometry.py
"""
import importlib.util, os, sys
_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hadrons', 'verify_fcc_geometry.py')
_spec = importlib.util.spec_from_file_location('verify_fcc_geometry', _path,
           submodule_search_locations=[os.path.dirname(_path)])
_mod = importlib.util.module_from_spec(_spec)
sys.modules['verify_fcc_geometry'] = _mod
_spec.loader.exec_module(_mod)
from verify_fcc_geometry import *
