#!/usr/bin/env python3
"""
fault_patch_far_field.py
========================
How far a Shockley fault patch on one facet of the coordination shell reaches.

A strange partial on a {111} facet bounds a fault patch, and the patch's two
faces sit in adjacent Peierls valleys with mutual disregistry b = ell/sqrt3.
This script asks how much of that reaches the sites OFF the patch, in
particular the <200> triple in the facet's own layer, which is closer to the
shell centre than the patch's far face (the extension triangle at sqrt3 ell)
but has no bond across the slipped interface.

Model: a large FCC block with the contact law (central k_n = 1, tangential
k_t = r_t k_n, r_t = 1/(3 pi - 5)), an eigen-slip b on the six bonds joining
the facet to the layer-4 triangle, outer boundary fixed, Phi u = f solved
sparsely.  Being linear, the model relaxes the patch to a fraction of b; the
number that matters is the RATIO of off-patch to on-patch displacement,
which is what the elastic field of the bounding loop sets and which does
not depend on how much of b the harmonic bonds retain.

Result: the <200> triple moves by 0.23 of the patch faces, unchanged
between 935 and 1700 sites.  Scaled to a patch held at the full Shockley
step, that is 0.07 to 0.13 ell against the zero-point envelope 0.404 ell.

Author: M. Cox, with Claude (Anthropic).  License: MIT.
"""

import numpy as np
from itertools import product
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

R_T = 1.0 / (3 * np.pi - 5)
A = 2.0                              # cubic constant in a/2 units -> sites at integer coords
ELL = np.sqrt(2.0)                   # nearest-neighbour separation
SIGMA = 0.404 * ELL                  # zero-point envelope, per axis, in a/2 units
RMAX = 7.5                           # block radius in a/2 units
RFIX = 6.0                           # sites beyond this are held fixed

sites = [np.array(c, float) for c in product(range(-9, 10), repeat=3)
         if sum(c) % 2 == 0 and np.dot(c, c) <= RMAX ** 2 + 1e-9]
X = np.array(sites)
n = len(X)
idx = {tuple(map(int, p)): i for i, p in enumerate(X)}
free = np.linalg.norm(X, axis=1) <= RFIX
print(f"sites {n}, free {free.sum()}")

# nearest-neighbour bonds
NN = [np.array(v) for v in product((-1, 0, 1), repeat=3) if sum(abs(x) for x in v) == 2]
bonds = []
for i, p in enumerate(X):
    for v in NN:
        q = tuple(map(int, p + v))
        if q in idx and idx[q] > i:
            bonds.append((i, idx[q]))
print(f"bonds {len(bonds)}")

# the fault patch: facet (layer 2) <-> layer-4 triangle
facet = [(1, 1, 0), (1, 0, 1), (0, 1, 1)]
tri = [(2, 1, 1), (1, 2, 1), (1, 1, 2)]
patch = set()
for f in facet:
    for t in tri:
        if np.dot(np.subtract(t, f), np.subtract(t, f)) == 2:
            patch.add((idx[f], idx[t]))
print(f"patch bonds {len(patch)} (expect 6)")

# Shockley partial in the (111) plane: b = (a/6)[11-2] -> in a/2 units (1/3, 1/3, -2/3)
b_vec = np.array([1.0, 1.0, -2.0]) / 3.0
print(f"|b| = {np.linalg.norm(b_vec):.4f} = {np.linalg.norm(b_vec)/ELL:.4f} ell  (Shockley: ell/sqrt3 = {1/np.sqrt(3):.4f} ell)")

# stiffness matrix and eigen-slip forces
K = lil_matrix((3 * n, 3 * n))
f = np.zeros(3 * n)
I3 = np.eye(3)
for (i, j) in bonds:
    d = X[j] - X[i]
    u = d / np.linalg.norm(d)
    Kb = np.outer(u, u) + R_T * (I3 - np.outer(u, u))
    for (a, c, s) in ((i, i, 1), (j, j, 1), (i, j, -1), (j, i, -1)):
        K[3 * a:3 * a + 3, 3 * c:3 * c + 3] += s * Kb
    if (i, j) in patch:
        # the bond wants relative displacement u_j - u_i = +b (layer-4 side slips by b)
        f[3 * j:3 * j + 3] += Kb @ b_vec
        f[3 * i:3 * i + 3] -= Kb @ b_vec
K = csr_matrix(K)

# fixed boundary: solve on free dofs only
fdof = np.repeat(free, 3)
u = np.zeros(3 * n)
u[fdof] = spsolve(K[fdof][:, fdof], f[fdof])
U = u.reshape(n, 3)

def disp(p):
    return U[idx[tuple(p)]]

print("\nDisplacement magnitude of the sites in question (units of ell), threshold sigma = 0.404")
groups = {
    "shell centre (0,0,0)": [(0, 0, 0)],
    "facet triangle, layer 2": facet,
    "extension triangle, layer 4": tri,
    "<200> triple, layer 2": [(2, 0, 0), (0, 2, 0), (0, 0, 2)],
    "<220> triple, layer 4": [(2, 2, 0), (2, 0, 2), (0, 2, 2)],
    "other shell nodes, layer 0": [(1, -1, 0), (-1, 1, 0), (1, 0, -1), (-1, 0, 1), (0, 1, -1), (0, -1, 1)],
    "other shell nodes, layer -2": [(-1, -1, 0), (-1, 0, -1), (0, -1, -1)],
    "layer-4 second ring (3,1,0)-type": [(3, 1, 0), (3, 0, 1), (1, 3, 0), (0, 3, 1), (1, 0, 3), (0, 1, 3)],
}
for name, pts in groups.items():
    mags = [np.linalg.norm(disp(p)) / ELL for p in pts]
    print(f"  {name:<36} |u|/ell = {np.mean(mags):.3f}  (min {min(mags):.3f}, max {max(mags):.3f})"
          f"   {'OUT of registry' if np.mean(mags) > 0.404 else 'in registry'}")

print("\nCross-interface disregistry (relative displacement along b) per bond, in units of |b|:")
def disreg(p, q):
    return np.dot(disp(q) - disp(p), b_vec) / np.dot(b_vec, b_vec)
pairs = {
    "patch bonds (facet -> extension triangle)": [(fp, t) for fp in facet for t in tri if np.dot(np.subtract(t, fp), np.subtract(t, fp)) == 2],
    "facet -> <220> triple (outside the patch, same interface)": [(fp, t) for fp in facet for t in [(2, 2, 0), (2, 0, 2), (0, 2, 2)] if np.dot(np.subtract(t, fp), np.subtract(t, fp)) == 2],
    "<200> triple -> layer 4 (outside the patch)": [(p, t) for p in [(2, 0, 0), (0, 2, 0), (0, 0, 2)] for t in tri + [(2, 2, 0), (2, 0, 2), (0, 2, 2)] if np.dot(np.subtract(t, p), np.subtract(t, p)) == 2],
    "extension triangle -> layer 6 (next interface out)": [(t, s) for t in tri for s in [(3, 2, 1), (3, 1, 2), (2, 3, 1), (1, 3, 2), (2, 1, 3), (1, 2, 3), (2, 2, 2)] if np.dot(np.subtract(s, t), np.subtract(s, t)) == 2],
}
for name, pp in pairs.items():
    vals = [disreg(p, q) for p, q in pp]
    print(f"  {name:<58} n={len(vals):<2} mean {np.mean(vals):+.3f}  (in ell: {np.mean(vals)*np.linalg.norm(b_vec)/ELL:+.3f})")
print(f"\nsigma in units of |b|: {SIGMA/np.linalg.norm(b_vec):.3f}")
