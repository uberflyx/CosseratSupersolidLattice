"""
This file has moved to hadrons/three_mechanisms.py

All imports are re-exported for backwards compatibility.
Update references to: hadrons/three_mechanisms.py
"""
import importlib.util, os, sys
_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hadrons', 'three_mechanisms.py')
_spec = importlib.util.spec_from_file_location('three_mechanisms', _path,
           submodule_search_locations=[os.path.dirname(_path)])
_mod = importlib.util.module_from_spec(_spec)
sys.modules['three_mechanisms'] = _mod
_spec.loader.exec_module(_mod)
from three_mechanisms import *
