#!/usr/bin/env python3
"""test_decay_engine.py — Regression tests for the decay engine."""
import sys, numpy as np
from cosserat_decay_engine import build_table

table = build_table()
n = 0; n15 = 0; fails = 0
for mo, nm, mp, mobs, dp, dobs, u in table:
    if not dobs: continue
    n += 1
    r = abs((dp - dobs) / dobs * 100)
    if r < 15: n15 += 1
    if r > 40:
        print(f"  FAIL: {nm}: {r:.1f}% > 40%")
        fails += 1

print(f"\n  {n15}/{n} within 15% | {n}/{n} within 40%: "
      f"{'PASS' if fails == 0 else 'FAIL'}")
sys.exit(fails)
