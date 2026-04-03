"""
This file has moved to neutrinos/neutrino_predictions.py

All imports are re-exported for backwards compatibility.
Update references to: neutrinos/neutrino_predictions.py
"""
import importlib.util, os, sys
_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'neutrinos', 'neutrino_predictions.py')
_spec = importlib.util.spec_from_file_location('neutrino_predictions', _path,
           submodule_search_locations=[os.path.dirname(_path)])
_mod = importlib.util.module_from_spec(_spec)
sys.modules['neutrino_predictions'] = _mod
_spec.loader.exec_module(_mod)
from neutrino_predictions import *
