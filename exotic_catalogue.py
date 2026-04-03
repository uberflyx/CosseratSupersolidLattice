"""
This file has moved to hadrons/exotic_catalogue.py

All imports are re-exported for backwards compatibility.
Update references to: hadrons/exotic_catalogue.py
"""
import importlib.util, os, sys
_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hadrons', 'exotic_catalogue.py')
_spec = importlib.util.spec_from_file_location('exotic_catalogue', _path,
           submodule_search_locations=[os.path.dirname(_path)])
_mod = importlib.util.module_from_spec(_spec)
sys.modules['exotic_catalogue'] = _mod
_spec.loader.exec_module(_mod)
from exotic_catalogue import *
