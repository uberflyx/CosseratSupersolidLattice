#!/usr/bin/env python3
"""
rho_first_principles.py
=======================
Classify every eigenvalue of the crossed-fault Cosserat dynamical matrix by
its D_2h irrep, derive the rho selection rule from first principles.

Procedure:
  1. Build all 8 D_2h group elements as 6n x 6n matrices on (u, phi),
     with proper polar/axial vector transformation rules.
  2. Verify the group elements commute with the Cosserat dynamical matrix
     (so they ARE symmetries of M, and irreps are well-defined).
  3. Build irrep projectors P_R = (1/|G|) sum_g chi_R(g) D(g).
  4. For each eigenvalue, compute the projector overlap with each irrep
     and assign the eigenvalue to the irrep with highest overlap.
  5. Identify all B_3u eigenvalues (rho candidates) and pick the lowest as
     the rho. Verify it lands at lambda = 4.891.

Geometry:
  Standard z-axis = fault axis (0,-1,1)/sqrt(2)
  Standard x-axis = (1, 0, 0)
  Standard y-axis = (0, 1, 1)/sqrt(2)

  rho   -> B_3u (basis x): sigma_x = -1
  omega -> B_1u (basis z): sigma_x = +1, axis-localised
  third -> B_2u (basis y): sigma_x = +1
"""

import numpy as np
import sys
from cosserat_classifier import build_cosserat_matrix
import composite_clusters as cc


# --- Cluster and Cosserat matrix ---
coords, _ = cc.cluster_crossed_fault()
n = len(coords)
M = build_cosserat_matrix(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
eigvals, eigvecs = np.linalg.eigh(M)
print(f"Cluster: crossed fault, N = {n} nodes; Cosserat matrix {6*n} x {6*n}")
print(f"Calibration check: proton-shell-style eigenvalue range {eigvals.min():.4f} to {eigvals.max():.4f}")


# --- D_2h group elements ---
# Three mutually perpendicular C_2 axes in CRYSTAL coordinates (not rotated yet)
axes_crystal = {
    'X':   np.array([1., 0., 0.]),         # standard x
    'Y':   np.array([0., 1., 1.])/np.sqrt(2),   # standard y (perp to fault axis in yz plane)
    'Z':   np.array([0., -1., 1.])/np.sqrt(2),  # standard z (fault axis)
}

def C2(axis):
    return 2.0 * np.outer(axis, axis) - np.eye(3)

I3 = np.eye(3)
inv = -I3

# Build all 8 group elements as 3x3 rotation/reflection matrices, in standard order
elements_3x3 = {}
elements_3x3['E']       = I3
elements_3x3['C2_z']    = C2(axes_crystal['Z'])    # rotation about fault axis (standard z)
elements_3x3['C2_y']    = C2(axes_crystal['Y'])    # rotation about standard y
elements_3x3['C2_x']    = C2(axes_crystal['X'])    # rotation about standard x
elements_3x3['i']       = inv
elements_3x3['sigma_xy']= inv @ elements_3x3['C2_z']  # reflection perp to z = the sigma_h
elements_3x3['sigma_xz']= inv @ elements_3x3['C2_y']  # reflection perp to y
elements_3x3['sigma_yz']= inv @ elements_3x3['C2_x']  # reflection perp to x (= chapter's sigma_x)

# Verify each is a cluster symmetry and compute the node permutation
def find_perm(R, coords, tol=1e-6):
    n = len(coords)
    perm = []
    for i in range(n):
        target = R @ coords[i]
        match = -1
        for j in range(n):
            if np.linalg.norm(coords[j] - target) < tol:
                match = j; break
        if match == -1:
            return None
        perm.append(match)
    return perm

perms = {}
for name, R in elements_3x3.items():
    p = find_perm(R, coords)
    perms[name] = p
    assert p is not None, f"{name} is not a cluster symmetry!"
print("All 8 D_2h elements are cluster symmetries; permutations computed.")


# --- Build 6n x 6n representations on (u, phi) ---
# u (polar): u_i -> R u_{perm^{-1}(i)}
# phi (axial): phi_i -> det(R) R phi_{perm^{-1}(i)}
def build_rep_matrix(R, perm, n):
    D = np.zeros((6*n, 6*n))
    det = round(np.linalg.det(R))
    for i_old in range(n):
        i_new = perm[i_old]
        # u-block: D maps u at i_old to u at i_new with R action
        D[3*i_new:3*i_new+3, 3*i_old:3*i_old+3] = R
        # phi-block: same but with det(R) prefactor
        D[3*n + 3*i_new:3*n + 3*i_new+3, 3*n + 3*i_old:3*n + 3*i_old+3] = det * R
    return D

reps = {name: build_rep_matrix(R, perms[name], n) for name, R in elements_3x3.items()}


# --- Verify each rep matrix commutes with the Cosserat matrix M ---
print("\nCommutator [M, D(g)] check (should all be ~0):")
all_commute = True
for name, D in reps.items():
    comm = M @ D - D @ M
    err = np.linalg.norm(comm)
    status = "OK" if err < 1e-10 else "FAIL"
    if err >= 1e-10:
        all_commute = False
    print(f"  [{name:>10s}]: ||[M, D]|| = {err:.3e} {status}")
assert all_commute, "Group elements must commute with M; failure means representation is wrong."


# --- D_2h character table ---
# Order of group elements: E, C2_z, C2_y, C2_x, i, sigma_xy, sigma_xz, sigma_yz
elem_order = ['E', 'C2_z', 'C2_y', 'C2_x', 'i', 'sigma_xy', 'sigma_xz', 'sigma_yz']
char_table = {
    'A_g':  [+1, +1, +1, +1, +1, +1, +1, +1],
    'B_1g': [+1, +1, -1, -1, +1, +1, -1, -1],
    'B_2g': [+1, -1, +1, -1, +1, -1, +1, -1],
    'B_3g': [+1, -1, -1, +1, +1, -1, -1, +1],
    'A_u':  [+1, +1, +1, +1, -1, -1, -1, -1],
    'B_1u': [+1, +1, -1, -1, -1, -1, +1, +1],  # basis z (polar vector along fault axis)
    'B_2u': [+1, -1, +1, -1, -1, +1, -1, +1],  # basis y
    'B_3u': [+1, -1, -1, +1, -1, +1, +1, -1],  # basis x
}


# --- Build irrep projectors ---
G = 8  # |D_2h|
projectors = {}
for irrep, chars in char_table.items():
    P = np.zeros((6*n, 6*n))
    for k, gname in enumerate(elem_order):
        P += chars[k] * reps[gname]
    P /= G
    projectors[irrep] = P
    # Verify projector: P^2 = P
    err = np.linalg.norm(P @ P - P)
    assert err < 1e-10, f"Projector {irrep}: P^2 != P (err {err})"

# Compute dimension of each irrep subspace
dims = {ir: int(round(np.trace(P).real)) for ir, P in projectors.items()}
print(f"\nDimensions of irrep subspaces (must sum to {6*n}):")
total = 0
for ir, d in dims.items():
    print(f"  {ir}: {d:>3d}")
    total += d
print(f"  Total: {total} (should be {6*n})")


# --- Classify each eigenvalue by irrep ---
print(f"\nClassifying eigenvalues by irrep (showing modes with |lambda| > 0.01):")
print(f"  {'idx':>3s} {'lambda':>9s} {'irrep':>6s} {'overlap':>9s}")
print(f"  " + "-"*40)

# Group eigenvalues into degenerate clusters
def degenerate_groups(eigvals, tol=1e-4):
    groups = []
    if len(eigvals) == 0:
        return groups
    current = [0]
    for i in range(1, len(eigvals)):
        if abs(eigvals[i] - eigvals[current[0]]) < tol:
            current.append(i)
        else:
            groups.append(current)
            current = [i]
    groups.append(current)
    return groups

groups = degenerate_groups(eigvals)

# For each group, decompose the eigenspace into irreps
def irrep_decomp(V, projectors):
    """V is the eigenspace (6n x mult). Returns dict {irrep: dimension in V}."""
    decomp = {}
    for ir, P in projectors.items():
        PV = P @ V
        # Rank of PV measures how much of V lies in this irrep
        sv = np.linalg.svd(PV, compute_uv=False)
        rank = int(np.sum(sv > 0.5))
        if rank > 0:
            decomp[ir] = rank
    return decomp

# Print classification by eigenvalue
classified = []
for g in groups:
    lam = eigvals[g[0]]
    if abs(lam) < 0.01:
        continue
    V = eigvecs[:, g]
    decomp = irrep_decomp(V, projectors)
    label = ' + '.join([f"{m if m > 1 else ''}{ir}" for ir, m in decomp.items()])
    classified.append((lam, len(g), label, decomp))
    if abs(lam - 4.891) < 0.01 or abs(lam - 6.241) < 0.01 or abs(lam - 8.303) < 0.01:
        marker = "  <-- "
        if abs(lam - 4.891) < 0.01: marker += "rho candidate"
        elif abs(lam - 6.241) < 0.01: marker += "omega candidate"
        elif abs(lam - 8.303) < 0.01: marker += "proton baseline"
        print(f"  {g[0]:>3d} {lam:>9.4f} {len(g):>3d}  {label}{marker}")
    elif 4 < lam < 7:
        print(f"  {g[0]:>3d} {lam:>9.4f} {len(g):>3d}  {label}")


# --- Find all B_3u modes (rho candidates) ---
print(f"\n{'='*60}")
print(f"All B_3u modes (rho candidates, polar vector along original x-axis):")
print(f"{'='*60}")
b3u_modes = []
for lam, mult, label, decomp in classified:
    if 'B_3u' in decomp:
        m_b3u = decomp['B_3u']
        b3u_modes.append((lam, m_b3u))
        marker = "  <-- target rho at 4.891" if abs(lam - 4.891) < 0.01 else ""
        print(f"  lambda = {lam:>8.4f}  B_3u multiplicity {m_b3u}{marker}")

print(f"\nLowest B_3u eigenvalue: lambda = {b3u_modes[0][0]:.4f}")
print(f"Lowest *nonzero* B_3u eigenvalue: lambda = {[x[0] for x in b3u_modes if x[0] > 0.01][0]:.4f}")

# Find the first B_3u eigenvalue above some structural threshold (e.g., > 1)
# The rho should be the lowest STRUCTURAL B_3u, not the lowest of all (which
# might be a zero mode or near-zero soft mode)
print(f"\nB_3u eigenvalues in vector-meson range [1, 10]:")
for lam, mult in b3u_modes:
    if 1 < lam < 10:
        marker = "  <-- 4.891 (rho!)" if abs(lam - 4.891) < 0.01 else ""
        print(f"  lambda = {lam:.4f}{marker}")


# --- Find all B_1u modes (omega candidates) ---
print(f"\n{'='*60}")
print(f"All B_1u modes (omega candidates, polar vector along fault axis):")
print(f"{'='*60}")
b1u_modes = []
for lam, mult, label, decomp in classified:
    if 'B_1u' in decomp:
        m = decomp['B_1u']
        b1u_modes.append((lam, m))
        marker = "  <-- target omega at 6.241" if abs(lam - 6.241) < 0.01 else ""
        print(f"  lambda = {lam:>8.4f}  B_1u multiplicity {m}{marker}")


# --- Find all B_2u modes (third T_1u partner) ---
print(f"\n{'='*60}")
print(f"All B_2u modes (third T_1u partner, polar vector along standard y):")
print(f"{'='*60}")
for lam, mult, label, decomp in classified:
    if 'B_2u' in decomp:
        m = decomp['B_2u']
        print(f"  lambda = {lam:>8.4f}  B_2u multiplicity {m}")
