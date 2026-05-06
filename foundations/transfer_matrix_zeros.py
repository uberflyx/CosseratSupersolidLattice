#!/usr/bin/env python3
"""
transfer_matrix_zeros.py
========================
Companion to cosserat_transfer_matrix.py.

Build the 12x12 Cosserat PN transfer matrix on the FCC lattice in
NATURAL ordering
  Psi = (u_1, u_2, u_3, phi_1, phi_2, phi_3,
         u_1', u_2', u_3', phi_1', phi_2', phi_3')
and read off the zero structure of the equation matrix M(x).

This complements cosserat_transfer_matrix.py, which reorders the basis
into block-diagonal form (sectors A + B + C + D).  That ordering hides
the structural zeros visible in the natural ordering, where the matrix
has 16 nonzeros out of 144 entries and six categories of zero entries
each carry definite physical meaning (PN potential never enters
microrotation rows; compression and axial torsion are sealed; the two
photon polarisations don't mix; matched curl couplings have opposite
signs; five irreps are absent at the leading reciprocal-lattice shells).

M. A. Cox, University of the Witwatersrand (2026)
"""

import numpy as np
from scipy.integrate import solve_ivp

np.set_printoptions(precision=4, suppress=True, linewidth=200)

# ------- Cosserat moduli (same as repo script) -------
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
gamma_torsion = gamma_c

# Calibrated PN amplitude from repo run: V0 ~ 106.4 gives |lambda_C| = alpha
V0 = 106.425237

def V_PN(x):
    return V0 * np.sin(np.pi * x / ell)**2


# ----------------------------------------------------------------------
# Build the 12x12 dPsi/dx = M(x) Psi matrix in NATURAL ordering
# ----------------------------------------------------------------------
# Psi = (u_1, u_2, u_3, phi_1, phi_2, phi_3,
#        u_1', u_2', u_3', phi_1', phi_2', phi_3')
# index: 0     1     2     3       4       5
#        6     7     8     9       10      11
# ----------------------------------------------------------------------
def M12_natural(x):
    M = np.zeros((12, 12))
    V = V_PN(x)

    # First block: dPsi_i/dx = Psi_{i+6}  for i = 0..5  (definitions)
    for i in range(6):
        M[i, i + 6] = 1.0

    # Second block: dPsi_{i+6}/dx for i = 0..5 (the dynamics)

    # ---- u_1 (compression sector A): M_long * u_1'' = V(x) u_1 ----
    # d(u_1')/dx = V/M_long * u_1
    M[6, 0] = V / M_long          # u_1'  <- u_1

    # ---- u_2 (in-plane shear, sector D): mu_tot u_2'' + kappa_c phi_3' = V u_2 ----
    # d(u_2')/dx = (V u_2 + kappa_c phi_3') / mu_tot
    M[7, 1] = V / mu_tot           # u_2' <- u_2  (from V(x))
    M[7, 11] = kappa_c / mu_tot    # u_2' <- phi_3'  (curl coupling)

    # ---- u_3 (anti-plane shear, sector C): mu_tot u_3'' - kappa_c phi_2' = V u_3 ----
    # d(u_3')/dx = (V u_3 - kappa_c phi_2') / mu_tot
    M[8, 2] = V / mu_tot
    M[8, 10] = -kappa_c / mu_tot

    # ---- phi_1 (axial torsion, sector B): gamma phi_1'' = 2 kappa_c phi_1 ----
    M[9, 3] = 2 * kappa_c / gamma_torsion

    # ---- phi_2 (anti-plane microrotation): gamma phi_2'' + kappa_c u_3' = 2 kappa_c phi_2 ----
    M[10, 4] = 2 * kappa_c / gamma_c
    M[10, 8] = kappa_c / gamma_c

    # ---- phi_3 (in-plane microrotation): gamma phi_3'' - kappa_c u_2' = 2 kappa_c phi_3 ----
    M[11, 5] = 2 * kappa_c / gamma_c
    M[11, 7] = -kappa_c / gamma_c

    return M


# ----------------------------------------------------------------------
# Integrate through one period to get the transfer matrix
# ----------------------------------------------------------------------
def transfer_matrix(M_func, dim=12):
    def rhs(x, Psi_flat):
        Psi = Psi_flat.reshape(dim, dim)
        return (M_func(x) @ Psi).flatten()
    Psi0 = np.eye(dim).flatten()
    sol = solve_ivp(rhs, [0, ell], Psi0, method='DOP853',
                    rtol=1e-12, atol=1e-14, max_step=ell / 2000)
    return sol.y[:, -1].reshape(dim, dim)


# ----------------------------------------------------------------------
# Analyse the zero pattern of the equation matrix M(x) (mid-period)
# ----------------------------------------------------------------------
print("=" * 78)
print("PART A.  Zero pattern of the EQUATION matrix M(x) at mid-period")
print("(this is where the dynamics lives -- T = exp(integral M dx))")
print("=" * 78)

M_mid = M12_natural(ell / 2)
labels = ['u1', 'u2', 'u3', 'p1', 'p2', 'p3',
          "u1'", "u2'", "u3'", "p1'", "p2'", "p3'"]
print(f"\n  Pattern (X = nonzero, . = zero):")
print(f"        {' '.join(f'{l:>6s}' for l in labels)}")
for i in range(12):
    row = ' '.join('     X' if abs(M_mid[i, j]) > 1e-10 else '     .'
                    for j in range(12))
    print(f"  {labels[i]:>4s} {row}")

print(f"\n  Nonzero entries: {(np.abs(M_mid) > 1e-10).sum()} of 144")
print(f"  Zeros          : {144 - (np.abs(M_mid) > 1e-10).sum()}")


# ----------------------------------------------------------------------
# Analyse the zero pattern of the transfer matrix T itself
# ----------------------------------------------------------------------
print()
print("=" * 78)
print("PART B.  Zero pattern of the FULL TRANSFER MATRIX T")
print("(the propagator over one lattice period)")
print("=" * 78)

T = transfer_matrix(M12_natural, 12)
print(f"\n  Pattern (X = nonzero with |T_ij| > 1e-6, . = zero):")
print(f"        {' '.join(f'{l:>6s}' for l in labels)}")
for i in range(12):
    row = ' '.join('     X' if abs(T[i, j]) > 1e-6 else '     .'
                    for j in range(12))
    print(f"  {labels[i]:>4s} {row}")

print(f"\n  Nonzero entries: {(np.abs(T) > 1e-6).sum()} of 144")

# Eigenvalues
evals = np.linalg.eigvals(T)
print(f"\n  Eigenvalues (sorted by |lambda|):")
for ev in sorted(evals, key=abs):
    print(f"    lambda = {ev.real:+12.6f} {ev.imag:+8.4f}i    |lambda| = {abs(ev):.6f}")


# ----------------------------------------------------------------------
# Identify the four blocks explicitly
# ----------------------------------------------------------------------
print()
print("=" * 78)
print("PART C.  Block decomposition: which fields couple to which")
print("=" * 78)

# Identify components in each sector
sectors = {
    'A (compression, u_1)':       [0, 6],         # u_1, u_1'
    'B (axial torsion, phi_1)':   [3, 9],         # phi_1, phi_1'
    'C (anti-plane: u_3, phi_2)': [2, 4, 8, 10],  # u_3, phi_2, u_3', phi_2'
    'D (in-plane: u_2, phi_3)':   [1, 5, 7, 11],  # u_2, phi_3, u_2', phi_3'
}

# Check for cross-block coupling in T
print(f"\n  Cross-block coupling matrix (max |T_ij| between sectors):")
print(f"  {'':<32s} {'A':>10s} {'B':>10s} {'C':>10s} {'D':>10s}")

sector_names = list(sectors.keys())
for sname1 in sector_names:
    inds1 = sectors[sname1]
    row = f"  {sname1:<32s}"
    for sname2 in sector_names:
        inds2 = sectors[sname2]
        block = T[np.ix_(inds1, inds2)]
        max_abs = np.max(np.abs(block))
        row += f"{max_abs:>10.2e} "
    print(row)

print(f"\n  Off-diagonal blocks should all be ~ 0  (block-diagonal structure)")
print(f"  This is the C_4v symmetry of propagation along x_1.")
print(f"  No cross-channel manipulation at LINEAR order.")


# ----------------------------------------------------------------------
# Cross-coupling at SECOND order: the mechanism the framework uses
# ----------------------------------------------------------------------
print()
print("=" * 78)
print("PART D.  Where the cross-coupling lives: the self-energy matrix")
print("=" * 78)

print(f"""
  At linear order, T is block-diagonal: each force has its own channel,
  no manipulation between them.

  At second order, the self-energy matrix Sigma couples sectors that have
  the same O_h irrep but different field content:

    Sigma = alpha * v v^T,    v = (1, 1/sqrt(pi))

  This is a RANK-1 matrix.  Its eigenvalues are:
    Sigma_+ = alpha (1 + 1/pi)  (the displacement-microrotation mode)
    Sigma_- = 0                  (the null mode)

  The off-diagonal element Sigma_12 = alpha/sqrt(pi) is what makes the
  weak force exist:
    (Sigma_12)^2 / (Sigma_11 * 2) = (alpha^2/pi) / (2 alpha) = alpha/(2 pi)
    chirality theta_ch = alpha^2/(2 pi)   (with the alpha already in Sigma)

  THIS IS THE MANIPULATION PATHWAY THE FRAMEWORK ALREADY USES:
    photon (T_1u, sector C+D)  --[Sigma_12]-->  weak (T_1g)

  The question for the engineer is: are there OTHER such pathways?
""")

# Check: what other pairs of sectors share an irrep?
print(f"  O_h irreps and their sector content:")
print(f"  ---------------------------------------------------")
print(f"  T_1u  (3-dim, polar vector, EM):    u_1, u_2, u_3   -> sectors A, D, C")
print(f"  T_1g  (3-dim, axial vector, weak):  phi_1, phi_2, phi_3 -> sectors B, C, D")
print(f"  A_1g  (1-dim, scalar, gravity):     div(u) -- couples to compression A")
print(f"  T_2g, E_g, A_2u, T_2u, A_2g, A_1u, E_u: angular overlaps zero on {{111}}+{{200}}")

print(f"""
  Three coupling pathways exist by simple representation theory:

   (1) T_1u <-> T_1g :  photon channel <-> microrotation channel
       This is the realised weak-force pathway.  alpha -> alpha^2/(2 pi).

   (2) T_1u <-> A_1g :  photon channel <-> compression channel
       Coupling element: divergence of u in the curl term.
       For pure transverse photon, div(u) = 0, so this coupling is
       killed at linear order: nabla x (nabla phi) = 0.
       This is the FIVE STAR ZERO -- it would be the EM-gravity coupling.

   (3) T_2g <-> T_2g (strong sector self-coupling):
       Stacking-fault energy.  The strong force lives here but is
       commensurate (no PN barrier), so the transfer matrix doesn't
       see it -- it lives outside this 12x12 block.
""")

# ----------------------------------------------------------------------
# The interesting zero: T_1u <-> A_1g
# ----------------------------------------------------------------------
print()
print("=" * 78)
print("PART E.  The interesting zero:  EM <-> Gravity coupling")
print("=" * 78)

print(f"""
  The framework says:
    grad x grad(Phi) = 0   =>   transverse photon CANNOT excite compression
    div(curl A) = 0        =>   compression CANNOT generate transverse photon

  Both identities are vector calculus, exact on a continuum.  On the
  discrete lattice they are NOT exact.  The errors are:
    grad x grad on lattice ~ a / lambda      (a = lattice spacing)
    div(curl) on lattice    ~ a / lambda
  i.e. the "zeros" leak by O(ka) per cycle.

  At wavelength lambda = lambda_C (Compton wavelength of electron),
  ka = 2 pi a / lambda_C.  With a = r_e = lambda_C * alpha,
  ka = 2 pi alpha ~ 0.046.

  So the EM-gravity leakage per period is ~ alpha (in amplitude).
  The probability of a photon producing compression: ~ alpha^2 ~ 5e-5.

  This is the order of magnitude of the Suchard-Pugh-type bounds on
  EM-gravity coupling.  The framework predicts this naturally because
  the lattice spacing breaks the continuum identity at exactly O(alpha).
""")
