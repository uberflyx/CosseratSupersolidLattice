"""Compare Fiedler value lambda_2 of the strange decuplet under two
constructions, using the repository's own FCC geometry as ground truth.

  old (decay chapter):  shell + ALL 4 voids + n_strange hex extensions
  new (spectral, void-swap): shell + (4 - n_strange) voids + n_strange hex ext

Width scales as sqrt(lambda_2 / lambda_2^Delta) * shadow in the decay engine,
so a change in lambda_2 directly changes the predicted strong width.
"""
# ============================================================================
# SUPERSEDED / EXPLORATORY -- not the authoritative decay path.
# The production engine is decays/cosserat_decay_engine.py, which derives
# g_piNN = N_H = 13 and the decuplet widths (Delta -3.7%, Sigma* -2.4%) from
# graph invariants with NO fitted coupling.  These scripts are session
# explorations kept for the record; their absolute scale is less accurate.
# See decays/README.md.
# ============================================================================

import sys as _sys, os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_os.path.dirname(_HERE))
_sys.path[:0] = [_ROOT, _os.path.join(_ROOT, 'spectral_mass')]
import numpy as np
from cosserat_graph import FCCLattice, Defect

lat = FCCLattice()

def lam2(defect):
    pn, A, L = defect.graph_matrices()
    ev = np.sort(np.linalg.eigvalsh(L))
    return len(pn), len(defect.nn_edge_list()), ev[1]

def build(n_voids_keep, n_strange):
    """shell + a chosen number of voids + n_strange hex extensions."""
    d = Defect(lat)
    d.add_shell()                 # 13 nodes
    # add only n_voids_keep of the 4 voids
    for vi in range(n_voids_keep):
        n = f'v{vi}'; d.nodes.add(n); d.roles[n] = 'void'
        d.pos[n] = lat.void_positions[vi].copy()
    d.has_voids = n_voids_keep > 0
    for pi in range(n_strange):   # one strange arm per plane
        d.add_strange_ext(pi)
    return d

print(f"ext nodes per plane = {lat.n_ext_per_plane}\n")

# Delta validation
d_delta = build(4, 0)
print("Delta  (shell+4 voids):           N=%2d  E=%2d  lambda_2=%.4f   [target 2.722]" % lam2(d_delta))

print("\n--- Sigma* (|S|=1) ---")
old_s = build(4, 1)      # decay chapter: keep all 4 voids, add 1 hex arm
new_s = build(3, 1)      # spectral void-swap: 3 voids + 1 hex arm
print("OLD add-hex (4 voids+1 arm):       N=%2d  E=%2d  lambda_2=%.4f" % lam2(old_s))
print("NEW void-swap (3 voids+1 arm):     N=%2d  E=%2d  lambda_2=%.4f" % lam2(new_s))

print("\n--- Xi* (|S|=2) ---")
old_x = build(4, 2)
new_x = build(2, 2)
print("OLD add-hex (4 voids+2 arm):       N=%2d  E=%2d  lambda_2=%.4f" % lam2(old_x))
print("NEW void-swap (2 voids+2 arm):     N=%2d  E=%2d  lambda_2=%.4f" % lam2(new_x))

# width impact: Gamma ~ sqrt(lambda_2/lambda_2_Delta) * shadow ; shadow same in both
L2D = 2.722
import math
def wfac(d): 
    _,_,l2 = lam2(d); return math.sqrt(l2/L2D)
print("\n--- width factor sqrt(lambda_2/lambda_2^Delta), shadow held fixed ---")
print("Sigma*: OLD %.4f  NEW %.4f  ratio NEW/OLD = %.3f" % (wfac(old_s), wfac(new_s), wfac(new_s)/wfac(old_s)))
print("Xi*   : OLD %.4f  NEW %.4f  ratio NEW/OLD = %.3f" % (wfac(old_x), wfac(new_x), wfac(new_x)/wfac(old_x)))
