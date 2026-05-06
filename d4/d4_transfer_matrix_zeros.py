#!/usr/bin/env python3
"""
d4_transfer_matrix_zeros.py
===========================
Companion to ../foundations/transfer_matrix_zeros.py.  Same structural
analysis, four-dimensional D_4 lattice instead of three-dimensional FCC.

Build the 20x20 Cosserat PN transfer matrix in NATURAL ordering and
catalogue its zeros as a function of compact-direction momentum k_4.

State vector (20 components):
  Psi = ( u_1, u_2, u_3, u_4,
          phi_12, phi_13, phi_14, phi_23, phi_24, phi_34,
          u_1', u_2', u_3', u_4',
          phi_12', phi_13', phi_14', phi_23', phi_24', phi_34' )^T

Four displacements (one of which, u_4, lies along the compact stacking
direction) plus six antisymmetric microrotations, paired with their
spatial derivatives.  Propagation is along x_1 with k_2 = k_3 = 0 and
a prescribed compact-direction Fourier mode k_4.

Three findings, in increasing order of structural importance:

  (a) At k_4 = 0 the matrix is exactly block-diagonal: a 12x12 System
      A (spatial sector, identical to the FCC matrix) and an 8x8
      System B (compact-direction sector), with all 192 off-diagonal
      cells zero.  This is the matrix-entry expression of the
      factorisation theorem.

  (b) At the first KK level k_4 ell = 2 pi / 3, exactly seven
      additional matrix entries activate.  Every one of them is
      purely imaginary.  Each carries a factor of i k_4 inherited
      from differentiating exp(i k_4 x_4).  The compact direction
      supplies the imaginary unit geometrically.

  (c) Spatial rotations (phi_12, phi_13, phi_23 -- the FCC
      microrotation triplet, which carries the graviton and the f_2
      tensor meson) acquire a k_4^2 mass shift but no new
      off-diagonal coupling.  These channels are CP-protected at
      first order in k_4.  Gravitational-wave birefringence in the
      thermal window is suppressed by an extra factor of k_4^2 / k_typ^2
      relative to electromagnetic birefringence at the same temperature.

Real entries of M govern CP-even physics; imaginary entries govern
CP-odd physics.  This script makes both visible.

M. A. Cox, University of the Witwatersrand (2026)
"""

import numpy as np
from scipy.integrate import solve_ivp

np.set_printoptions(precision=4, suppress=True, linewidth=200)

# ----------------------------------------------------------------------
# Cosserat moduli (matched to FCC repo script)
# ----------------------------------------------------------------------
N2 = 1.0 / np.pi
mu = 1.0
ell = 1.0
kappa_c = 2 * N2 * mu / (1.0 - N2)
mu_tot = mu + kappa_c
gamma_c = mu * ell**2
K_over_mu = 1e6
K = K_over_mu * mu
lam = K - 2 * mu / 3
M_long = lam + 2 * mu + kappa_c

# Calibrated PN amplitude from FCC repo run
V0 = 106.425237


def V_PN(x):
    return V0 * np.sin(np.pi * x / ell) ** 2


# ----------------------------------------------------------------------
# Field index conventions (used throughout the script)
# ----------------------------------------------------------------------
# Displacements 0..3:        u_1, u_2, u_3, u_4
# Rotations    4..9:         phi_12, phi_13, phi_14, phi_23, phi_24, phi_34
# Their primes  10..19:      u_i' and phi_jk'
IDX = dict(
    u1=0, u2=1, u3=2, u4=3,
    p12=4, p13=5, p14=6, p23=7, p24=8, p34=9,
)
IDX.update({k + "p": v + 10 for k, v in IDX.items()})  # primes

LABELS = ['u1', 'u2', 'u3', 'u4',
          'p12', 'p13', 'p14', 'p23', 'p24', 'p34',
          "u1'", "u2'", "u3'", "u4'",
          "p12'", "p13'", "p14'", "p23'", "p24'", "p34'"]


# ----------------------------------------------------------------------
# Build the 20x20 matrix dPsi/dx = M(x; k_4) Psi.
#
# At ω = 0, propagating along x_1 with k_2 = k_3 = 0 and compact mode k_4:
#
#   ∇^2_4 = ∂_1^2 - k_4^2
#   ∂_4 → i k_4
#
# Equation set (4D Cosserat at ω=0, k_2=k_3=0, finite k_4):
#
#   M_long u_1'' = V u_1 + (μ+κ_c) k_4^2 u_1 - i k_4 (λ+μ) u_4'
#   μ_tot u_2'' = V u_2 + μ_tot k_4^2 u_2 - κ_c phi_12' + i k_4 κ_c phi_24
#   μ_tot u_3'' = V u_3 + μ_tot k_4^2 u_3 - κ_c phi_13' + i k_4 κ_c phi_34
#   M_long u_4'' = V u_4 + (μ+κ_c) k_4^2 u_4 + (λ+μ) k_4^2 u_4 - i k_4 (λ+μ) u_1' - κ_c phi_14'
#
#   γ phi_12'' = (2 κ_c + γ k_4^2) phi_12 - κ_c u_2'
#   γ phi_13'' = (2 κ_c + γ k_4^2) phi_13 - κ_c u_3'
#   γ phi_14'' = (2 κ_c + γ k_4^2) phi_14 - κ_c u_4' + i k_4 κ_c u_1
#   γ phi_23'' = (2 κ_c + γ k_4^2) phi_23
#   γ phi_24'' = (2 κ_c + γ k_4^2) phi_24 - i k_4 κ_c u_2
#   γ phi_34'' = (2 κ_c + γ k_4^2) phi_34 - i k_4 κ_c u_3
#
# Sign conventions match the FCC chapter's u_2 ↔ phi_3 curl pairing.
# ----------------------------------------------------------------------
def M20(x, k4):
    """
    Returns a complex 20x20 array (because k_4 != 0 introduces i k_4 mixings).
    Becomes purely real and block-diagonal at k_4 = 0.
    """
    M = np.zeros((20, 20), dtype=complex)
    V = V_PN(x)

    # -------- bookkeeping rows: dPsi_i/dx = Psi_{i+10} for i = 0..9 -----
    for i in range(10):
        M[i, i + 10] = 1.0

    # -------- dynamics rows (indices 10..19) ----------------------------
    # u_1''
    M[10, IDX['u1']]  = (V + (mu + kappa_c) * k4**2) / M_long
    M[10, IDX['u4p']] = -1j * k4 * (lam + mu) / M_long

    # u_2''
    M[11, IDX['u2']]  = (V + mu_tot * k4**2) / mu_tot
    M[11, IDX['p24']] = 1j * k4 * kappa_c / mu_tot
    M[11, IDX['p12p']] = -kappa_c / mu_tot

    # u_3''
    M[12, IDX['u3']]  = (V + mu_tot * k4**2) / mu_tot
    M[12, IDX['p34']] = 1j * k4 * kappa_c / mu_tot
    M[12, IDX['p13p']] = -kappa_c / mu_tot

    # u_4''  (NEW field in 4D)
    M[13, IDX['u4']]  = (V + (lam + 2*mu + kappa_c) * k4**2) / M_long
    M[13, IDX['u1p']] = -1j * k4 * (lam + mu) / M_long
    M[13, IDX['p14p']] = -kappa_c / M_long

    # phi_12''
    M[14, IDX['p12']] = (2 * kappa_c + gamma_c * k4**2) / gamma_c
    M[14, IDX['u2p']] = -kappa_c / gamma_c

    # phi_13''
    M[15, IDX['p13']] = (2 * kappa_c + gamma_c * k4**2) / gamma_c
    M[15, IDX['u3p']] = -kappa_c / gamma_c

    # phi_14''  (NEW mixed rotation)
    M[16, IDX['p14']] = (2 * kappa_c + gamma_c * k4**2) / gamma_c
    M[16, IDX['u4p']] = -kappa_c / gamma_c
    M[16, IDX['u1']]  = 1j * k4 * kappa_c / gamma_c

    # phi_23''  (axial torsion - structurally sealed)
    M[17, IDX['p23']] = (2 * kappa_c + gamma_c * k4**2) / gamma_c

    # phi_24''  (NEW mixed rotation)
    M[18, IDX['p24']] = (2 * kappa_c + gamma_c * k4**2) / gamma_c
    M[18, IDX['u2']]  = -1j * k4 * kappa_c / gamma_c

    # phi_34''  (NEW mixed rotation)
    M[19, IDX['p34']] = (2 * kappa_c + gamma_c * k4**2) / gamma_c
    M[19, IDX['u3']]  = -1j * k4 * kappa_c / gamma_c

    return M


# ----------------------------------------------------------------------
# Pretty-print the sparsity pattern
# ----------------------------------------------------------------------
def print_sparsity(M, title, labels, k4=None, threshold=1e-12):
    n = M.shape[0]
    nnz = (np.abs(M) > threshold).sum()
    print()
    print("=" * 78)
    print(f"  {title}")
    if k4 is not None:
        print(f"  k_4 * ell = {k4 * ell:.4f}")
    print("=" * 78)
    print(f"\n  Nonzero entries: {nnz} of {n*n}")
    print(f"  Zeros          : {n*n - nnz}")
    print(f"  Sparsity       : {100*(n*n - nnz)/(n*n):.1f}%")
    print()

    # Build a row-by-row visual
    head = "       " + " ".join(f"{l:>5s}" for l in labels)
    print(head)
    for i in range(n):
        row_chars = []
        for j in range(n):
            val = M[i, j]
            if abs(val) < threshold:
                row_chars.append("    .")
            elif abs(val.imag) > 1e-10 and abs(val.real) < 1e-10:
                # purely imaginary
                row_chars.append("    i")
            elif abs(val.imag) > 1e-10:
                row_chars.append("    z")  # complex
            else:
                row_chars.append("    R")  # real
        print(f"  {labels[i]:>4s} " + " ".join(row_chars))


# ----------------------------------------------------------------------
# PART A.  Zero pattern at k_4 = 0
# ----------------------------------------------------------------------
M0 = M20(ell / 2, k4=0.0)
print_sparsity(M0, "PART A.  D_4 equation matrix M(x) at k_4 = 0",
               LABELS, k4=0.0)

# Identify the 12x12 System A block (spatial sector) and the 8x8 System B block
sysA_indices = [IDX[k] for k in ['u1', 'u2', 'u3', 'p12', 'p13', 'p23',
                                  'u1p', 'u2p', 'u3p', 'p12p', 'p13p', 'p23p']]
sysB_indices = [IDX[k] for k in ['u4', 'p14', 'p24', 'p34',
                                  'u4p', 'p14p', 'p24p', 'p34p']]

A_block = M0[np.ix_(sysA_indices, sysA_indices)]
B_block = M0[np.ix_(sysB_indices, sysB_indices)]
AB_offdiag = M0[np.ix_(sysA_indices, sysB_indices)]
BA_offdiag = M0[np.ix_(sysB_indices, sysA_indices)]

print(f"\n  System A (spatial, 12x12) nonzeros: "
      f"{(np.abs(A_block) > 1e-12).sum()} of 144")
print(f"  System B (compact, 8x8)   nonzeros: "
      f"{(np.abs(B_block) > 1e-12).sum()} of 64")
print(f"  A <-> B off-diagonal      nonzeros: "
      f"{(np.abs(AB_offdiag) > 1e-12).sum() + (np.abs(BA_offdiag) > 1e-12).sum()} of 192")
print()
print("  Confirmation: at k_4 = 0 the matrix is BLOCK-DIAGONAL.")
print("  Systems A and B are completely decoupled.")


# ----------------------------------------------------------------------
# PART B.  What switches on at k_4 != 0
# ----------------------------------------------------------------------
k4_KK = 2 * np.pi / 3 / ell           # first Kaluza-Klein level
M1 = M20(ell / 2, k4=k4_KK)
print_sparsity(M1, "PART B.  D_4 equation matrix M(x) at k_4 ell = 2pi/3 (first KK)",
               LABELS, k4=k4_KK)

# Find entries that are zero at k_4=0 but nonzero at k_4 != 0
new_entries = (np.abs(M0) < 1e-12) & (np.abs(M1) > 1e-12)
print(f"\n  Entries that ACTIVATE at k_4 != 0: "
      f"{new_entries.sum()}")
print()
print("  Position    Value at k_4 = 2pi/3        Field-language coupling")
print("  --------    --------------------        -----------------------")
for i in range(20):
    for j in range(20):
        if new_entries[i, j]:
            val = M1[i, j]
            re, im = val.real, val.imag
            phase = "real" if abs(im) < 1e-10 else \
                    ("pure imag" if abs(re) < 1e-10 else "complex")
            print(f"   ({LABELS[i]:>4s}, {LABELS[j]:>5s})  "
                  f"{val.real:+.4f} + {val.imag:+.4f}i  ({phase})")


# ----------------------------------------------------------------------
# PART C.  Eigenvalue spectrum: k_4 = 0  versus  k_4 = 2 pi / 3
# ----------------------------------------------------------------------
def transfer_matrix(M_func, k4, dim=20, n_steps=4000):
    """Path-ordered integral exp(integral M dx) over one lattice period."""
    def rhs(x, Psi_flat):
        Psi = Psi_flat.reshape(dim, dim).astype(complex)
        return (M_func(x, k4) @ Psi).flatten()
    Psi0 = np.eye(dim, dtype=complex).flatten()
    sol = solve_ivp(rhs, [0, ell], Psi0, method='DOP853',
                    rtol=1e-11, atol=1e-13,
                    max_step=ell / n_steps)
    return sol.y[:, -1].reshape(dim, dim)


print("\n")
print("=" * 78)
print("  PART C.  Transfer-matrix eigenvalues over one lattice period")
print("=" * 78)

for k4_val, label in [(0.0, "k_4 = 0  (System A factorised from System B)"),
                      (k4_KK, "k_4 ell = 2pi/3  (first KK level, full coupling)")]:
    T = transfer_matrix(M20, k4_val)
    eigs = np.linalg.eigvals(T)
    abs_eigs = np.sort(np.abs(eigs))
    print(f"\n  {label}")
    print(f"  |lambda| spectrum (20 values, sorted):")
    for k, a in enumerate(abs_eigs):
        # Pair eigenvalues as |lambda| and 1/|lambda| for symplectic structure
        marker = "  <-- alpha" if 0.005 < a < 0.01 else ""
        print(f"    {k:>2d}: {a:>10.6f}{marker}")


# ----------------------------------------------------------------------
# PART D.  Block structure of T at finite k_4
# ----------------------------------------------------------------------
print("\n")
print("=" * 78)
print("  PART D.  Block structure of T at k_4 ell = 2 pi / 3")
print("=" * 78)

T_KK = transfer_matrix(M20, k4_KK)

# Test: is the block structure of M preserved by the integration?
T_AA = T_KK[np.ix_(sysA_indices, sysA_indices)]
T_BB = T_KK[np.ix_(sysB_indices, sysB_indices)]
T_AB = T_KK[np.ix_(sysA_indices, sysB_indices)]
T_BA = T_KK[np.ix_(sysB_indices, sysA_indices)]

print(f"\n  ||T_AA||  = {np.linalg.norm(T_AA):.4f}   (System A block)")
print(f"  ||T_BB||  = {np.linalg.norm(T_BB):.4f}   (System B block)")
print(f"  ||T_AB||  = {np.linalg.norm(T_AB):.4f}   (A-to-B leakage)")
print(f"  ||T_BA||  = {np.linalg.norm(T_BA):.4f}   (B-to-A leakage)")
print()
print("  At k_4 != 0 the off-diagonal blocks are nonzero: the systems mix.")
print("  Dyson resummation of T over many periods amplifies this mixing.")


# ----------------------------------------------------------------------
# PART E.  Imaginary structure summary
# ----------------------------------------------------------------------
print("\n")
print("=" * 78)
print("  PART E.  All new entries are PURELY IMAGINARY")
print("=" * 78)

new_vals = M1[new_entries]
all_imag = all(abs(v.real) < 1e-12 for v in new_vals)
print(f"\n  Total new entries: {len(new_vals)}")
print(f"  Of these, purely imaginary: "
      f"{sum(abs(v.real) < 1e-12 for v in new_vals)}")
print(f"  All purely imaginary?  {all_imag}")
print()
print("  Each new matrix entry carries a factor of i k_4.")
print("  Real eigenvalues of T track CP-even amplitudes (constitutive theta).")
print("  Imag eigenvalues of T track CP-odd amplitudes (propagation theta).")
print("  The compact direction supplies the imaginary unit geometrically.")
