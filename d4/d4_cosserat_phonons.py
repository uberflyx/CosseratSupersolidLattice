#!/usr/bin/env python3
"""
d4_cosserat_phonons.py — 4D Cosserat phonon structure
=====================================================
Explores the phonon branch structure of the D4 Cosserat lattice.

Key question: how does the 4D Cosserat theory reduce to 3D?

In 3D:
  - Displacement u_i: 3 components (vector)
  - Microrotation phi_i: 3 components (pseudovector, Hodge dual of SO(3))
  - Curl coupling: eps_{ijk} d_j phi_k  (vector -> vector)
  - Total: 6 dofs/node, 12 with dipolar cell

In 4D:
  - Displacement u_mu: 4 components (vector)
  - Microrotation phi_{mu,nu}: 6 components (antisymmetric 2-form, SO(4))
  - Coupling: d^nu phi_{mu,nu}  (codifferential, 2-form -> 1-form)
  - Total: 10 dofs/node, 20 with dipolar cell

The crucial test: does the k4=0 sector of the 4D theory reproduce
the 3D Cosserat equations exactly?

M. A. Cox, University of the Witwatersrand (2026)
"""

import numpy as np
from itertools import combinations

# ================================================================
# PART 1: THE 3+1 DECOMPOSITION OF SO(4) ROTATIONS
# ================================================================

print("=" * 72)
print("  PART 1: SO(4) DECOMPOSITION")
print("=" * 72)

# SO(4) generators: 6 independent rotation planes
planes_4d = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]  # (μ,ν) pairs
labels_4d = ['xy', 'xz', 'xw', 'yz', 'yw', 'zw']

# Split into spatial (SO(3)) and mixed (boost-like)
spatial_rotations = [(1,2), (0,2), (0,1)]  # yz, xz, xy → phi_1, phi_2, phi_3
mixed_rotations  = [(0,3), (1,3), (2,3)]   # xw, yw, zw → phi_14, phi_24, phi_34

print(f"""
  SO(4) has {len(planes_4d)} generators (rotation planes):
    Spatial SO(3): {[labels_4d[planes_4d.index(p)] for p in spatial_rotations]}
      → 3D Cosserat microrotation phi_i
    Mixed (spatial-compact): {[labels_4d[planes_4d.index(p)] for p in mixed_rotations]}
      → New vector field phi_{{i4}} in 3D projection

  Dimension count per node:
    3D: 3 (u_i) + 3 (phi_i) = 6
    4D: 4 (u_mu) + 6 (phi_{{mu,nu}}) = 10
    
  With dipolar 2-site cell:
    3D: 12 branches
    4D: 20 branches
""")


# ================================================================
# PART 2: SELF-DUAL / ANTI-SELF-DUAL DECOMPOSITION
# ================================================================

print("=" * 72)
print("  PART 2: CHIRAL DECOMPOSITION  SO(4) ≅ SU(2)_L × SU(2)_R")
print("=" * 72)

# In 4D Euclidean space, the Hodge star maps 2-forms to 2-forms.
# Self-dual: *phi = +phi  (3 components)
# Anti-self-dual: *phi = -phi  (3 components)
#
# Explicit basis (using eps_{1234} = +1, coordinates x,y,z,w):
#   phi^+_1 = (phi_{xy} + phi_{zw}) / sqrt(2)
#   phi^+_2 = (phi_{xz} - phi_{yw}) / sqrt(2) 
#   phi^+_3 = (phi_{xw} + phi_{yz}) / sqrt(2)
#
#   phi^-_1 = (phi_{xy} - phi_{zw}) / sqrt(2)
#   phi^-_2 = (phi_{xz} + phi_{yw}) / sqrt(2)
#   phi^-_3 = (phi_{xw} - phi_{yz}) / sqrt(2)

# Verify: the Hodge dual of phi_{mu,nu} in 4D is
# (*phi)_{mu,nu} = (1/2) eps_{mu,nu,rho,sigma} phi^{rho,sigma}

eps4 = np.zeros((4,4,4,4))
for perm in [(0,1,2,3), (0,1,3,2), (0,2,1,3), (0,2,3,1), (0,3,1,2), (0,3,2,1),
             (1,0,2,3), (1,0,3,2), (1,2,0,3), (1,2,3,0), (1,3,0,2), (1,3,2,0),
             (2,0,1,3), (2,0,3,1), (2,1,0,3), (2,1,3,0), (2,3,0,1), (2,3,1,0),
             (3,0,1,2), (3,0,2,1), (3,1,0,2), (3,1,2,0), (3,2,0,1), (3,2,1,0)]:
    # Compute sign from number of transpositions
    a = list(perm)
    sign = 1
    for i in range(4):
        while a[i] != i:
            j = a[i]
            a[i], a[j] = a[j], a[i]
            sign *= -1
    eps4[perm] = sign

# Verify eps_{0123} = +1
assert eps4[0,1,2,3] == 1.0

# Build the 6x6 Hodge star matrix in the basis of independent 2-forms
# Ordered: (01), (02), (03), (12), (13), (23)
basis_pairs = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]

hodge_matrix = np.zeros((6,6))
for a, (mu, nu) in enumerate(basis_pairs):
    for b, (rho, sigma) in enumerate(basis_pairs):
        hodge_matrix[a, b] = 0.5 * eps4[mu, nu, rho, sigma]

print("\n  Hodge star matrix (* on 2-forms in 4D):")
print("  Basis: (01)=xy, (02)=xz, (03)=xw, (12)=yz, (13)=yw, (23)=zw")
for i in range(6):
    row = "  " + " ".join(f"{hodge_matrix[i,j]:+5.1f}" for j in range(6))
    print(f"    {row}")

# Eigenvalues of the Hodge star
evals_hodge = np.linalg.eigvalsh(hodge_matrix)
print(f"\n  Eigenvalues of *: {sorted(evals_hodge)}")
print(f"  → 3 eigenvalues +1 (self-dual), 3 eigenvalues -1 (anti-self-dual)")

# Find eigenvectors
evals_h, evecs_h = np.linalg.eigh(hodge_matrix)
sd_vecs = evecs_h[:, evals_h > 0.5].T   # self-dual
asd_vecs = evecs_h[:, evals_h < -0.5].T  # anti-self-dual

print(f"\n  Self-dual eigenvectors (each row):")
for i, v in enumerate(sd_vecs):
    components = " + ".join(f"{v[j]:+.3f}·{labels_4d[j]}" for j in range(6) if abs(v[j]) > 0.01)
    print(f"    phi^+_{i+1} = {components}")

print(f"\n  Anti-self-dual eigenvectors:")
for i, v in enumerate(asd_vecs):
    components = " + ".join(f"{v[j]:+.3f}·{labels_4d[j]}" for j in range(6) if abs(v[j]) > 0.01)
    print(f"    phi^-_{i+1} = {components}")

# KEY OBSERVATION: each self-dual / anti-self-dual component
# MIXES a spatial rotation with a compact-direction rotation
print(f"""
  KEY RESULT: Each chiral component mixes spatial and compact rotations.
  
  The 3+1 split (spatial vs compact) and the chiral split (SD vs ASD)
  are INDEPENDENT decompositions of the same 6D space.
  
  Neither is a refinement of the other.
""")


# ================================================================
# PART 3: THE k4=0 FACTORISATION
# ================================================================

print("=" * 72)
print("  PART 3: COUPLING STRUCTURE AT k4 = 0")
print("=" * 72)

# The 4D Cosserat coupling in the u_mu equation:
#   F_mu = kappa_c * d^nu phi_{mu,nu}    (codifferential of phi)
#
# In 3D (mu = i, summing nu over spatial only at k4=0):
#   F_i = kappa_c * d_j phi_{ij}
#       = kappa_c * (curl of phi)_i   [by Hodge duality in 3D]
#
# The mixed rotation phi_{i4} does NOT contribute at k4=0
# because its coupling requires d_4 which vanishes.
#
# The phi equation for spatial rotations (at k4=0):
#   J phi_ij_tt = gamma nabla^2 phi_{ij} + kappa_c (d_i u_j - d_j u_i) - 2 kappa_c phi_{ij}
#   → IDENTICAL to 3D Cosserat equation for phi_i
#
# The phi equation for mixed rotations (at k4=0):
#   J phi_{i4}_tt = gamma nabla^2 phi_{i4} + kappa_c d_i u_4 - 2 kappa_c phi_{i4}
#   → couples phi_{i4} to u_4 via GRADIENT (not curl!)
#   → independent massive system

print(f"""
  4D COSSERAT EQUATIONS (Eringen tensorial form):
  
  u equation:
    rho u_mu_tt = (lam+mu) d_mu (div u) + (mu+kappa) nabla^2 u_mu 
                  + kappa d^nu phi_{{mu,nu}}
  
  phi equation:
    J phi_{{mu,nu}}_tt = gamma nabla^2 phi_{{mu,nu}} 
                        + kappa (d_mu u_nu - d_nu u_mu) 
                        - 2 kappa phi_{{mu,nu}}
  
  AT k4 = 0 (zero mode in compact direction, d_4 = 0):
  ──────────────────────────────────────────────────────
  
  SYSTEM A: spatial sector (u_i, phi_i)
    F_i = kappa * d_j phi_{{ij}} = kappa * (curl phi)_i     ← 3D CURL
    phi_{{ij}} driven by (d_i u_j - d_j u_i)               ← 3D CURL
    
    → EXACTLY the 3D Cosserat equations.
    → 6 dofs/node × 2 (dipolar) = 12 branches ← THE EXISTING TABLE
    
  SYSTEM B: compact sector (u_4, phi_{{i4}})
    F_4 = kappa * d_j phi_{{4j}} = -kappa * div(phi_4)     ← DIVERGENCE
    phi_{{i4}} driven by kappa * d_i u_4                    ← GRADIENT
    
    → A SEPARATE massive scalar-vector system.
    → 4 dofs/node × 2 (dipolar) = 8 branches ← NEW, ALL MASSIVE
    
  The factorisation is EXACT at k4 = 0.
  Systems A and B are completely independent.
""")


# ================================================================
# PART 4: MIXING AT k4 ≠ 0
# ================================================================

print("=" * 72)
print("  PART 4: CHIRAL MIXING AT k4 ≠ 0")
print("=" * 72)

print(f"""
  At k4 ≠ 0, the coupling acquires a new term:
  
    F_i = kappa * (d_j phi_{{ij}} + i*k4 * phi_{{i4}})
                   ^^^^^^^^^^^      ^^^^^^^^^^^^^^^^^
                   System A           System B leaks in!
  
  And in the phi equation:
    
    kappa * (d_i u_j - d_j u_i)   → spatial curl (System A)
    kappa * (d_i u_4 - i*k4*u_i)  → MIXES u_i into System B!
    
  The two systems HYBRIDISE at k4 ≠ 0.
  The mixing strength is proportional to k4.
  
  For the first KK level: k4 = 2pi/(3 ell), so k4*ell = 2pi/3 ≈ 2.09.
  This is an ORDER-ONE mixing — not a small perturbation!
  
  THE CHIRAL STRUCTURE BECOMES RELEVANT HERE:
  
  The self-dual combination phi^+ = (phi_spatial + phi_compact)/sqrt(2)
  and anti-self-dual phi^- = (phi_spatial - phi_compact)/sqrt(2)
  mix spatial and compact rotations in equal measure.
  
  At k4 ≠ 0, the hybridised system naturally diagonalises in the
  chiral basis, not the spatial/compact basis.
  
  This means: ABOVE THE KK SCALE (~121 MeV), the physical eigenmodes
  are not "spatial rotations" or "compact rotations" but LEFT-HANDED
  and RIGHT-HANDED rotations in the SO(4) sense.
""")


# ================================================================
# PART 5: PHYSICAL CONSEQUENCES
# ================================================================

print("=" * 72)
print("  PART 5: PHYSICAL CONSEQUENCES")
print("=" * 72)

# Compute the System B mass gap
N2 = 1.0 / np.pi  # Cosserat coupling number
kappa_over_mu = 2 * N2 / (1 - N2)  # kappa_c / mu
omega0_sq = 2 * kappa_over_mu  # 2*kappa_c/(J) in units where mu=1, J~rho*ell^2

m0_MeV = 0.51099895069 / (1/137.035999177)  # node mass
alpha = 1/137.035999177
ell_fm = alpha * 0.197326980 / (0.51099895069 / 1000)  # ell in fm... 
# Actually: ell = r_e = alpha * hbar_c / (m_e c) in natural units
# r_e = alpha * hbar/(m_e c) = 2.8179 fm
ell_fm = 2.8179

# KK mass from D4 dynamical matrix
m_KK_T = np.sqrt(3) * m0_MeV  # transverse KK
m_KK_L = 3 * m0_MeV            # longitudinal KK

# System B gap mass (from the -2*kappa*phi term)
# The gap frequency is omega_0 = sqrt(2*kappa_c/J)
# In the monograph: this is the same gap that appears in the
# Cosserat optical branch. For the f2(1270), the optical gap
# is at ~ 18*m0 (from the Born cluster).
# But for System B, the gap is from the BARE kappa_c coupling,
# not the Born-cluster-enhanced version.

print(f"""
  SYSTEM B MASS SPECTRUM (compact-direction modes):
  
  The 8 branches of System B split as:
  
    u_4 sector (scalar, 2 branches with dipolar cell):
      Acoustic: compact-direction Goldstone mode
                (massless before gap, but acquires gap from Cosserat coupling)
                Gap ~ sqrt(2*kappa_c/J) × m_0 ≈ {np.sqrt(2*kappa_over_mu):.2f} × m_0
      Optical:  zone-boundary mode
    
    phi_{{i4}} sector (3-component vector, 6 branches with dipolar cell):
      3 acoustic: massive vector field
      3 optical: massive vector field at higher gap
    
  KK masses (from D4 dynamical matrix):
    Transverse: m_T = sqrt(3) × m_0 = {m_KK_T:.1f} MeV
    Longitudinal: m_L = 3 × m_0 = {m_KK_L:.1f} MeV
  
  
  WHY THE CHIRAL STRUCTURE MATTERS:
  ─────────────────────────────────
  
  1. Below KK scale (~121 MeV): Systems A and B are independent.
     All familiar physics (photon, graviton, Dirac equation, 
     12-branch table) is UNCHANGED. The 3D Cosserat theory is exact.
     
  2. At the KK scale: Systems A and B hybridise. The physical
     eigenmodes become chiral (self-dual / anti-self-dual).
     The SO(4) ≅ SU(2)_L × SU(2)_R structure is dynamically relevant.
     
  3. The D4 lattice has a BUILT-IN chirality from the Z3 stacking:
     A→B→C→A has a handedness. This handedness distinguishes
     self-dual from anti-self-dual at the KK scale.
     
  4. The monograph's weak force arises from the evanescent T_{{1g}}
     channel, which is a CROSS-COUPLING between displacement and
     microrotation (the Sigma_12 off-diagonal self-energy). In the
     4D picture, this cross-coupling naturally splits into chiral
     components above the KK scale.

  OPEN QUESTIONS (flagged for development):
  ──────────────────────────────────────────
  
  • Is the electroweak chirality (W couples only to left-handed 
    fermions) a consequence of the D4 chiral structure?
    
  • The Z3 stacking acts differently on phi^+ and phi^- because
    it mixes spatial and compact rotations. Does this select one
    chirality for the weak interaction?
    
  • The Dirac equation emerges from the curl coupling at the BZ 
    boundary. In 4D, the codifferential coupling d^nu phi_{{mu,nu}}
    has a richer structure. Does the BZ-boundary analysis in 4D
    produce the Standard Model's chiral fermion content?
    
  • System B contains a massive vector field (phi_{{i4}}, 3 components,
    spin-1 from 3D perspective). What are its quantum numbers?
    Could it be identified with a known particle or channel?
""")

# ================================================================
# PART 6: BRANCH COUNTING SUMMARY
# ================================================================

print("=" * 72)
print("  SUMMARY: COMPLETE D4 PHONON BRANCH COUNT")
print("=" * 72)

print(f"""
  ┌─────────────────────────────────────────────────────────────────┐
  │                                                                 │
  │  4D COSSERAT THEORY: 10 dofs/node × 2 (dipolar) = 20 branches │
  │                                                                 │
  │  3+1 DECOMPOSITION AT k4 = 0:                                  │
  │                                                                 │
  │  SYSTEM A: spatial sector (u_i, phi_i)  ← 12 branches          │
  │    Curl coupling: d_j phi_{{ij}} = (curl phi)_i                  │
  │    EXACTLY reproduces 3D Cosserat / existing Table              │
  │    Contains: Goldstone, photon, graviton, f2, Higgs, ...       │
  │                                                                 │
  │  SYSTEM B: compact sector (u_4, phi_{{i4}})  ← 8 branches       │
  │    Gradient coupling: d_i u_4 (NOT curl!)                       │
  │    All branches massive (gap from kappa_c)                      │
  │    u_4: scalar field (compact-direction displacement)           │
  │    phi_{{i4}}: vector field (spatial-compact rotation)            │
  │                                                                 │
  │  CHIRAL STRUCTURE:                                              │
  │    SO(4) ≅ SU(2)_L × SU(2)_R                                  │
  │    Self-dual phi^+ = (phi_spatial + phi_compact)/sqrt(2)       │
  │    Anti-self-dual phi^- = (phi_spatial - phi_compact)/sqrt(2)  │
  │    k4=0: chiral decomposition hidden (A and B independent)     │
  │    k4≠0: chiral decomposition is the natural eigenbasis        │
  │                                                                 │
  │  THE 12-BRANCH TABLE IS THE k4=0 PROJECTION OF SYSTEM A.      │
  │  IT IS EXACT AND UNCHANGED.                                    │
  │                                                                 │
  └─────────────────────────────────────────────────────────────────┘

  Why "4D × 2 × 2 = 16" is wrong:
  ─────────────────────────────────
  In 3D, dim(SO(3)) = 3 = dim(R³) — a coincidence that makes the
  "×2 Cosserat" counting work. In 4D, dim(SO(4)) = 6 ≠ 4 = dim(R⁴).
  The rotation sector has MORE modes than the translation sector.
  Correct count: (4 + 6) × 2 = 20, not (4 + 4) × 2 = 16.
""")

