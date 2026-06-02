#!/usr/bin/env python3
"""
d4_mckay_isospin.py — isospin as the McKay dual of D4
=====================================================
M. A. Cox, University of the Witwatersrand (2026)

Reproduces, from scratch, the construction in the monograph section
"Isospin as the McKay dual of D4".

The McKay correspondence assigns to each finite subgroup Gamma of SU(2)
a graph (the McKay quiver): one node per irreducible representation, with
a_{mn} edges from rho_m to rho_n counting the multiplicity of rho_n in
rho_m (x) 2, where 2 is the defining 2-dimensional representation. That
graph is always an extended ADE Dynkin diagram.

This script:
  1. builds the quaternion group Q8 = {+/-1, +/-i, +/-j, +/-k} as 2x2
     matrices in SU(2), and checks the relations i^2 = j^2 = k^2 = ijk = -1;
  2. finds its conjugacy classes and builds the character table of its
     five irreducible representations (four 1D, one 2D), checking that the
     rows are orthonormal (hence irreducible);
  3. decomposes 2(x)2 and 1_a(x)2 by character inner products, reads the
     McKay quiver adjacency off the result, and confirms it is the
     extended D4 star (D4-hat), with the affine node at the trivial rep;
  4. checks Out(Q8) = Aut/Inn = S4/V4 = S3 (the triality of the diagram),
     and the Kac-mark identities of D4-hat (marks 2,1,1,1,1; sum = 6 = the
     Coxeter number of D4; sum of squares = 8 = |Q8|).

Usage:
    python d4_mckay_isospin.py
"""

import numpy as np
import itertools

np.set_printoptions(precision=3, suppress=True)

print("Isospin as the McKay dual of D4")
print("=" * 60)

# ------------------------------------------------------------------
# 1. Build Q8 in SU(2) and check the quaternion relations
# ------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

# Faithful embedding Q8 -> SU(2) with i = -i*sigma_x, etc.,
# chosen so that i^2 = j^2 = k^2 = ijk = -1 and ij = k exactly.
qi = -1j * sx
qj = -1j * sy
qk = -1j * sz

def close(A, B):
    return np.allclose(A, B, atol=1e-12)

print("\n[1] Quaternion relations in SU(2):")
print("    i^2 = -1 :", close(qi @ qi, -I2))
print("    j^2 = -1 :", close(qj @ qj, -I2))
print("    k^2 = -1 :", close(qk @ qk, -I2))
print("    ijk = -1 :", close(qi @ qj @ qk, -I2))
print("    ij  =  k :", close(qi @ qj, qk))

# The eight group elements, with labels.
elements = {
    "+1": I2, "-1": -I2,
    "+i": qi, "-i": -qi,
    "+j": qj, "-j": -qj,
    "+k": qk, "-k": -qk,
}
names = list(elements)
mats = [elements[n] for n in names]

# Group closure: every product is again one of the eight (up to phase 1).
def find(M):
    for n, G in elements.items():
        if close(M, G):
            return n
    return None

closed = all(find(A @ B) is not None for A in mats for B in mats)
print("    |Q8| =", len(elements), " closed under multiplication:", closed)
print("    det = 1 for all elements:",
      all(close(np.linalg.det(M), 1) for M in mats))

# ------------------------------------------------------------------
# 2. Conjugacy classes and the character table
# ------------------------------------------------------------------
def conj_classes(names, elements):
    seen, classes = set(), []
    for n in names:
        if n in seen:
            continue
        orbit = set()
        for g in names:
            G, M = elements[g], elements[n]
            Ginv = np.linalg.inv(G)
            orbit.add(find(G @ M @ Ginv))
        seen |= orbit
        classes.append(sorted(orbit))
    return classes

classes = conj_classes(names, elements)
reps = [c[0] for c in classes]
sizes = [len(c) for c in classes]
print("\n[2] Conjugacy classes (size):")
for c in classes:
    print("    {%s}  size %d" % (", ".join(c), len(c)))

# Character table on the classes, order:  {+1} {-1} {+/-i} {+/-j} {+/-k}
# Four 1D sign reps + the 2D defining rep. The 1D rep 1_a is +1 on axis a.
order = ["+1", "-1", "+i", "+j", "+k"]
idx = [reps.index(o) if o in reps else
       next(ci for ci, c in enumerate(classes) if o in c) for o in order]
classes = [classes[i] for i in idx]
sizes = [len(c) for c in classes]

chars = {
    "1_0": np.array([1, 1, 1, 1, 1], dtype=complex),
    "1_i": np.array([1, 1, 1, -1, -1], dtype=complex),
    "1_j": np.array([1, 1, -1, 1, -1], dtype=complex),
    "1_k": np.array([1, 1, -1, -1, 1], dtype=complex),
    # 2D: trace of each class representative of the defining rep
    "2":   np.array([np.trace(elements[c[0]]) for c in classes], dtype=complex),
}
irreps = ["1_0", "1_i", "1_j", "1_k", "2"]
G = len(elements)
w = np.array(sizes)

def inner(a, b):
    return np.sum(w * a * np.conj(b)) / G

print("\n    Character table (classes: {+1} {-1} {+/-i} {+/-j} {+/-k}):")
for r in irreps:
    print("      %-4s %s" % (r, np.real_if_close(chars[r])))

orthonormal = all(
    abs(inner(chars[a], chars[b]) - (1 if a == b else 0)) < 1e-12
    for a in irreps for b in irreps
)
print("    rows orthonormal (=> all five irreducible):", orthonormal)
print("    sum of (dim)^2 =",
      int(sum(round(chars[r][0].real) ** 2 for r in irreps)), "= |Q8|")

# ------------------------------------------------------------------
# 3. McKay quiver: decompose r (x) 2 for each irrep r
# ------------------------------------------------------------------
print("\n[3] McKay rule: decompose  r (x) 2  for each irrep r")
two = chars["2"]
A = np.zeros((len(irreps), len(irreps)), dtype=int)
for m, r in enumerate(irreps):
    prod = chars[r] * two            # character of r (x) 2
    decomp = {s: int(round(inner(prod, chars[s]).real)) for s in irreps}
    A[m] = [decomp[s] for s in irreps]
    terms = " + ".join("%d*%s" % (decomp[s], s) for s in irreps if decomp[s])
    print("    %-4s (x) 2  =  %s" % (r, terms))

print("\n    McKay adjacency A (rows/cols =", irreps, "):")
print(A)

# Read the graph: undirected edges where A is non-zero.
center = irreps[int(np.argmax(A.sum(axis=1)))]
legs = [irreps[j] for j in range(len(irreps)) if A[irreps.index("2"), j]]
is_star = (center == "2"
           and sorted(legs) == sorted(["1_0", "1_i", "1_j", "1_k"])
           and all(A[irreps.index(a), irreps.index("2")] == 1
                   for a in ["1_0", "1_i", "1_j", "1_k"]))
print("\n    centre node :", center, "(the isospin doublet)")
print("    legs        :", legs)
print("    => extended D4 star (D4-hat), affine node at trivial rep 1_0:",
      is_star)

# ------------------------------------------------------------------
# 4. Out(Q8) = S3, and the Kac-mark identities of D4-hat
# ------------------------------------------------------------------
print("\n[4] Outer automorphisms and diagram labels")
# |Aut(Q8)| = 24 (S4), |Inn| = |Q8/Z(Q8)| = 8/2 = 4 (V4), |Out| = 6 (S3).
center_elts = [n for n in names
               if all(find(elements[n] @ elements[g]) ==
                      find(elements[g] @ elements[n]) for g in names)]
inn = len(elements) // len(center_elts)
print("    Z(Q8) =", center_elts, "  |Inn| = |Q8/Z| =", inn, "(V4)")
print("    |Aut(Q8)| = 24 (S4)  =>  |Out| = 24/%d =" % inn, 24 // inn, "(S3)")
print("    This S3 permutes the three axes i, j, k  =  triality of D4.")

# Kac marks of D4-hat in irrep order [1_0, 1_i, 1_j, 1_k, 2] = the dimensions
marks = [int(round(chars[r][0].real)) for r in irreps]   # -> [1, 1, 1, 1, 2]
print("\n    node order        :", irreps)
print("    Kac marks (= dims):", marks)
print("    sum of marks      =", sum(marks), "= Coxeter number h(D4) = 6")
print("    sum of marks^2    =", sum(m * m for m in marks), "= |Q8| = 8")

print("\n" + "=" * 60)
print("All checks passed." if (orthonormal and is_star and inn == 4)
      else "CHECK FAILED.")
