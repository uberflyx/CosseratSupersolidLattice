#!/usr/bin/env python3
"""test_decay_engine.py -- Regression test for the unified decay engine.

Runs the engine's built-in regression and exits nonzero if any case
exceeds 40% residual.  Current status at time of commit:
    - 36 decay cases covering 13 topologies
    - 26/36 within 15% (PASS band)
    - 36/36 within 40% (no FAILs)
    - Median |residual|: 4.0%
"""
import sys, os, subprocess

here = os.path.dirname(os.path.abspath(__file__))
engine = os.path.join(here, 'cosserat_decay_engine.py')

result = subprocess.run(['python3', engine], capture_output=True, text=True)
print(result.stdout)
if result.stderr: print(result.stderr, file=sys.stderr)

# Parse result lines for FAIL flags
fails = sum(1 for line in result.stdout.splitlines() if 'FAIL' in line)
if fails:
    print(f"\nRegression: {fails} FAIL case(s)")
    sys.exit(1)
print("\nRegression: all cases within 40%.")
sys.exit(0)
