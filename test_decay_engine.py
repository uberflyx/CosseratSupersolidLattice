"""
This file has moved to decays/test_decay_engine.py

All imports are re-exported for backwards compatibility.
Update references to: decays/test_decay_engine.py
"""
import importlib.util, os, sys
_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'decays', 'test_decay_engine.py')
_spec = importlib.util.spec_from_file_location('test_decay_engine', _path,
           submodule_search_locations=[os.path.dirname(_path)])
_mod = importlib.util.module_from_spec(_spec)
sys.modules['test_decay_engine'] = _mod
_spec.loader.exec_module(_mod)
from test_decay_engine import *
