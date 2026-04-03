"""
This file has moved to foundations/theta_ch_ab_initio_v2.py

All imports are re-exported for backwards compatibility.
Update references to: foundations/theta_ch_ab_initio_v2.py
"""
import importlib.util, os, sys
_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'foundations', 'theta_ch_ab_initio_v2.py')
_spec = importlib.util.spec_from_file_location('theta_ch_ab_initio_v2', _path,
           submodule_search_locations=[os.path.dirname(_path)])
_mod = importlib.util.module_from_spec(_spec)
sys.modules['theta_ch_ab_initio_v2'] = _mod
_spec.loader.exec_module(_mod)
from theta_ch_ab_initio_v2 import *
