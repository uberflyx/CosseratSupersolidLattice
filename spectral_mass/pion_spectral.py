#!/usr/bin/env python3
"""
pion_spectral.py
================
First-principles spectral derivation of the pion mass on the smallest
possible FCC cluster: the cell pair (N = 2).

The cell pair is two bonded FCC nearest-neighbour nodes; it is the structural
unit of the up and down quarks, the host for a single Shockley partial
dislocation (Sec. sec:cosserat_micro), and the smallest defect-supporting
cluster in the framework.

If the spectral mass formula is universal, it must give the pion mass
correctly at N = 2.  The naive leading order is striking:

  N * m_0 = 2 * 70.0253 MeV = 140.05 MeV

against the charged pion PDG mass of 139.57 MeV.  The bare cell-pair mass is
within 0.34% of the charged pion, suggesting the cell pair IS the pion's
structural cluster, with a small spectral correction.

This script builds the cell-pair 12 x 12 Cosserat dynamical matrix, computes
the spectrum, identifies the pseudoscalar (J^P = 0^-) mode, and applies the
master formula.

Master formula: m = N * m_0 - N * (4 - lambda) * m_e
For the pion at PDG, the required eigenvalue is:
  pi+: lambda = 4 + (139.57 - 140.05)/(2 * 0.511) = 3.530
  pi0: lambda = 4 + (134.98 - 140.05)/(2 * 0.511) = -0.96
"""

import numpy as np

m_e = 0.51099895069
alpha = 7.2973525643e-3
m_0 = m_e / alpha

PDG_pi_plus = 139.57039  # MeV
PDG_pi_zero = 134.9768   # MeV
PDG_pi_iso = (2 * PDG_pi_plus + PDG_pi_zero) / 3  # isospin average

print(f"m_e = {m_e:.6f} MeV")
print(f"m_0 = m_e/alpha = {m_0:.4f} MeV")
print(f"\nPDG: m(pi+) = {PDG_pi_plus} MeV")
print(f"PDG: m(pi0) = {PDG_pi_zero} MeV")
print(f"PDG: isospin-averaged pion mass = {PDG_pi_iso:.4f} MeV")
print(f"\nN * m_0 at N=2: {2*m_0:.4f} MeV (bare cell pair)")
print(f"  Against pi+: dev = {(2*m_0 - PDG_pi_plus)/PDG_pi_plus*100:+.4f}%")
print(f"  Against pi0: dev = {(2*m_0 - PDG_pi_zero)/PDG_pi_zero*100:+.4f}%")
print(f"  Against isoavg: dev = {(2*m_0 - PDG_pi_iso)/PDG_pi_iso*100:+.4f}%")

print(f"\nRequired spectral eigenvalues:")
print(f"  pi+: lambda = {4 + (PDG_pi_plus - 2*m_0)/(2*m_e):.4f}")
print(f"  pi0: lambda = {4 + (PDG_pi_zero - 2*m_0)/(2*m_e):.4f}")
print(f"  isoavg: lambda = {4 + (PDG_pi_iso - 2*m_0)/(2*m_e):.4f}")


# --- Build the cell pair Cosserat dynamical matrix ---
# Two FCC nodes at (0,0,0) and (1,1,0) in a/2 units, separated by NN distance sqrt(2)
positions = np.array([
    [0., 0., 0.],
    [1., 1., 0.],
])
n = len(positions)

# Cosserat matrix: 6n x 6n = 12 x 12
# Conventions matching the rest of the framework's spectral analysis:
#   u-u: K * (rhat rhat^T) at diagonal, -K * rhat rhat^T at off-diagonal
#   phi-phi: gamma * I at diagonal, -gamma * I at off-diagonal
#   u-phi: alpha_Cos coupling per bond
#
# Per the eta-family derivation (Sec. subsec:eta_family_masses) the central NN
# coupling at alpha_Cos = 1 reproduces the framework's universal sum rules.

def build_dm(positions, alpha_cos=1.0, K=1.0, gamma=1.0):
    n = len(positions)
    D = np.zeros((6*n, 6*n))
    bonds = []
    NN_DIST = np.sqrt(2.0)  # in a/2 units
    for i in range(n):
        for j in range(i+1, n):
            d = np.linalg.norm(positions[j] - positions[i])
            if abs(d - NN_DIST) < 1e-6:
                bonds.append((i, j))
    for (i, j) in bonds:
        r = positions[j] - positions[i]
        rh = r / np.linalg.norm(r)
        rr = np.outer(rh, rh)
        # u-u block
        D[3*i:3*i+3, 3*j:3*j+3] -= K * rr
        D[3*j:3*j+3, 3*i:3*i+3] -= K * rr
        D[3*i:3*i+3, 3*i:3*i+3] += K * rr
        D[3*j:3*j+3, 3*j:3*j+3] += K * rr
        # phi-phi block
        bi, bj = 3*n + 3*i, 3*n + 3*j
        D[bi:bi+3, bj:bj+3] -= gamma * np.eye(3)
        D[bj:bj+3, bi:bi+3] -= gamma * np.eye(3)
        D[bi:bi+3, bi:bi+3] += gamma * np.eye(3)
        D[bj:bj+3, bj:bj+3] += gamma * np.eye(3)
        # Cosserat coupling
        for (a, b) in [(i, j), (j, i)]:
            ua, pb = 3*a, 3*n + 3*b
            D[ua:ua+3, pb:pb+3] += alpha_cos * np.eye(3)
            D[pb:pb+3, ua:ua+3] += alpha_cos * np.eye(3)
        for a in (i, j):
            ua, pa = 3*a, 3*n + 3*a
            D[ua:ua+3, pa:pa+3] += alpha_cos * np.eye(3)
            D[pa:pa+3, ua:ua+3] += alpha_cos * np.eye(3)
    return 0.5*(D + D.T), bonds


D, bonds = build_dm(positions)
print(f"\nCell-pair dynamical matrix: {D.shape}")
print(f"Number of NN bonds: {len(bonds)}")

eigvals, eigvecs = np.linalg.eigh(D)
print(f"\nFull spectrum (eigenvalue, multiplicity):")
rounded = np.round(eigvals, 4)
unique = np.unique(rounded)
for u in unique:
    mult = int(np.sum(rounded == u))
    if abs(u) > 1e-6:
        m_pred = 2 * m_0 - 2 * (4 - u) * m_e
        print(f"  lambda = {u:>8.4f}  (mult {mult})  ->  m = {m_pred:.3f} MeV")
    else:
        print(f"  lambda = {u:>8.4f}  (mult {mult})  ->  zero mode (rigid translation/rotation)")

# Identify the pseudoscalar mode: J^P = 0^-
# Under the cell pair's symmetry (D_inf_h along the bond axis), the irreps are:
#   Sigma_g+ (scalar even): bond longitudinal symmetric stretch
#   Sigma_u+ (vector along bond): one node moves +, other -
#   Sigma_g- (pseudoscalar): pure rotation parallel to bond? actually subtle
#   Pi_g, Pi_u (transverse vector/pseudovector)
# 
# The pion (J^P = 0^-) is a pseudoscalar. In the cell pair, the candidate is
# the mode that's antisymmetric under spatial inversion (P-odd) but scalar
# under rotations.
# 
# For 2 nodes along a bond:
# - The microrotation phi along the bond direction is parity-odd (axial vector,
#   axial=g under O_h but along the bond it's the parity-odd mode of the
#   cell-pair's reduced symmetry).
# - Inversion swaps the two nodes; the symmetric combination (phi_1 + phi_2)
#   along the bond is parity-even (g); the antisymmetric (phi_1 - phi_2)
#   along the bond is parity-odd (u).
# 
# The pseudoscalar J^P = 0^- on the cell pair is the antisymmetric microrotation
# along the bond axis.

# Project each eigenvector onto the candidate pseudoscalar mode
bond_axis = (positions[1] - positions[0]) / np.linalg.norm(positions[1] - positions[0])

# Pseudoscalar mode: phi_1 along bond axis - phi_2 along bond axis (normalized)
# In the 6n = 12 vector layout: [u_0_x, u_0_y, u_0_z, u_1_x, u_1_y, u_1_z, 
#                                phi_0_x, phi_0_y, phi_0_z, phi_1_x, phi_1_y, phi_1_z]
pseudo_mode = np.zeros(12)
pseudo_mode[6:9] = bond_axis    # phi at node 0 along bond
pseudo_mode[9:12] = -bond_axis  # phi at node 1 along bond, opposite sign
pseudo_mode /= np.linalg.norm(pseudo_mode)

print(f"\nProjection of each eigenvector onto pseudoscalar mode (antisymm phi along bond):")
print(f"  {'idx':>3s} {'lambda':>8s} {'overlap':>10s} {'character':>20s}")
for i in range(12):
    overlap = pseudo_mode @ eigvecs[:, i]
    if abs(overlap) > 0.1:
        marker = " <-- pseudoscalar candidate"
    else:
        marker = ""
    print(f"  {i:>3d} {eigvals[i]:>8.4f} {abs(overlap):>10.4f}{marker}")

# Find the eigenvalue with maximum overlap with the pseudoscalar mode
overlaps = np.abs(eigvecs.T @ pseudo_mode)
best_idx = np.argmax(overlaps)
lam_pion = eigvals[best_idx]
m_pion_pred = 2 * m_0 - 2 * (4 - lam_pion) * m_e

print(f"\nMost pseudoscalar-like eigenvalue: lambda = {lam_pion:.4f}")
print(f"  Master formula: m = 2*m_0 - 2*(4 - lambda)*m_e = {m_pion_pred:.4f} MeV")
print(f"  PDG pi+ = {PDG_pi_plus} MeV (dev: {(m_pion_pred - PDG_pi_plus)/PDG_pi_plus*100:+.3f}%)")
print(f"  PDG pi0 = {PDG_pi_zero} MeV (dev: {(m_pion_pred - PDG_pi_zero)/PDG_pi_zero*100:+.3f}%)")
print(f"  PDG iso = {PDG_pi_iso:.4f} MeV (dev: {(m_pion_pred - PDG_pi_iso)/PDG_pi_iso*100:+.3f}%)")
