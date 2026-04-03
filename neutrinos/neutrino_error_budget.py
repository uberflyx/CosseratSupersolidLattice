#!/usr/bin/env python3
"""
neutrino_error_budget.py
========================
Companion script for: "The Cosserat Supersolid" monograph (M. A. Cox, 2026).
Repository: https://github.com/uberflyx/CosseratSupersolidLattice

Computes zero-parameter predictions for all neutrino-sector observables
and propagates theoretical uncertainties via Monte Carlo.

No experimental neutrino data are used as input.  The sole measured
quantity is the electron mass m_e.  The fine-structure constant alpha is
derived from FCC Cosserat geometry (monograph Ch. 4).

See monograph Secs. 4.7 (chirality derivation and error budget)
and Ch. 13 (Z3 cyclic Hamiltonian, mixing angles) for full derivations.
"""

import numpy as np

# ══════════════════════════════════════════════════════════════
# INPUTS
# ══════════════════════════════════════════════════════════════

# Fine-structure constant (derived from FCC Cosserat geometry, monograph Ch. 4)
alpha = 1.0 / 137.035999177

# Electron mass (CODATA 2022) — the framework's sole dimensional input
m_e_MeV = 0.51099895069   # MeV/c^2

# Unit conversion
hbar_c = 197.3269804       # hbar*c in MeV·fm


# ══════════════════════════════════════════════════════════════
# DERIVED LATTICE PARAMETERS
# ══════════════════════════════════════════════════════════════

m_0_MeV = m_e_MeV / alpha              # node mass (= m_e / alpha)
m_0_eV  = m_0_MeV * 1e6                # node mass in eV
r_e_fm  = alpha * hbar_c / m_e_MeV     # lattice spacing = classical electron radius
N2      = 1.0 / np.pi                  # Cosserat coupling number (rolling contact)


# ══════════════════════════════════════════════════════════════
# NEUTRINO-SECTOR PARAMETERS (all derived, monograph Sec. 4.7)
# ══════════════════════════════════════════════════════════════

# Lattice chirality: theta_ch = alpha^2 / (2*pi)
# Centrosymmetry of FCC forces chirality to enter at 2nd order in alpha.
theta_ch = alpha**2 / (2 * np.pi)

# Lightest neutrino mass: m1 = theta_ch^2 * m0 = alpha^3 * m_e / (4*pi^2)
# N^2 cancels between theta_ch and the constitutive normalisation.
m1_eV = alpha**3 * m_e_MeV * 1e6 / (4 * np.pi**2)

# Chiral phase shift per stacking hop (NLO, per-leg EM vertex screening)
# Bare: delta = N^4 = 1/pi^2.  Dressed: delta = (1-alpha)^2 / pi^2.
delta = (1 - alpha)**2 / np.pi**2

# Hopping ratio h/A from strain + curvature overlap at 60 degrees
# h/A = cos(60) + N^2 * cos(60) = (1 + N^2) / 2 = (pi+1) / (2*pi)
rho = (1 + N2) / 2


# ══════════════════════════════════════════════════════════════
# Z3 EIGENVALUE SOLVER
# ══════════════════════════════════════════════════════════════

def compute_spectrum(m1, delta_val, rho_val):
    """
    Solve the Z3 cyclic Hamiltonian for three neutrino mass eigenvalues.

    The mass matrix M = A*I + h*(omega*P + conj(omega)*P^T) has exact
    eigenvalues E_k = A + 2*h*cos(phi + 2*pi*k/3) for k = 0, 1, 2,
    with phi = 2*pi/3 + delta, h = rho*A, and A fixed by requiring
    the lightest eigenvalue to equal m1.

    Returns sorted eigenvalues [eV], on-site energy A, hopping h.
    """
    phi = 2 * np.pi / 3 + delta_val
    A   = m1 / (1 + 2 * rho_val * np.cos(phi))
    h   = rho_val * A
    E   = np.array([A + 2*h*np.cos(phi + 2*np.pi*k/3) for k in range(3)])
    E.sort()
    return E, A, h


# ══════════════════════════════════════════════════════════════
# CENTRAL PREDICTIONS
# ══════════════════════════════════════════════════════════════

masses, A, h = compute_spectrum(m1_eV, delta, rho)
m1, m2, m3   = masses
sum_nu       = m1 + m2 + m3
Dm21_sq      = m2**2 - m1**2
Dm31_sq      = m3**2 - m1**2

# Experimental values for comparison
Dm21_exp, Dm21_err = 7.50e-5, 0.12e-5     # JUNO 2025
Dm31_exp, Dm31_err = 2.534e-3, 0.023e-3   # NuFit 6.0 (normal ordering)

pull_21 = (Dm21_sq - Dm21_exp) / Dm21_err
pull_31 = (Dm31_sq - Dm31_exp) / Dm31_err


# ══════════════════════════════════════════════════════════════
# ERROR PROPAGATION (Monte Carlo)
#
# The eigenvalue formula is nonlinear (cosines), so we propagate
# the theoretical uncertainties numerically.  Each input (m1, delta,
# rho) is Gaussian-perturbed by its theoretical error, the eigenvalues
# are recomputed, and the spread of 200k output samples gives the
# propagated ± on each observable.
#
# The uncertainties are estimates of uncalculated higher-order
# corrections (monograph Sec. 4.7), not experimental errors.
# ══════════════════════════════════════════════════════════════

N_samples  = 200_000
dm1_rel    = np.sqrt(2) * 2 * alpha   # ~2.1%, sources (i)+(ii) in quadrature
ddelta_rel = alpha                      # ~0.73%, NLO screening ambiguity
drho_rel   = 0.0005                     # ~0.05%, h/A derivation precision

rng = np.random.default_rng(seed=42)

m1_s    = m1_eV * (1 + dm1_rel    * rng.standard_normal(N_samples))
delta_s = delta  * (1 + ddelta_rel * rng.standard_normal(N_samples))
rho_s   = rho    * (1 + drho_rel   * rng.standard_normal(N_samples))

mc = np.zeros((N_samples, 3))
for i in range(N_samples):
    mc[i], _, _ = compute_spectrum(m1_s[i], delta_s[i], rho_s[i])

mc_sums = mc.sum(axis=1)
mc_Dm21 = mc[:, 1]**2 - mc[:, 0]**2
mc_Dm31 = mc[:, 2]**2 - mc[:, 0]**2

err = {
    'm1':    np.std(mc[:, 0]),
    'm2':    np.std(mc[:, 1]),
    'm3':    np.std(mc[:, 2]),
    'sum':   np.std(mc_sums),
    'Dm21':  np.std(mc_Dm21),
    'Dm31':  np.std(mc_Dm31),
}


# ══════════════════════════════════════════════════════════════
# OUTPUT
# ══════════════════════════════════════════════════════════════

print("=" * 70)
print("COSSERAT SUPERSOLID — NEUTRINO SECTOR PREDICTIONS")
print("M. A. Cox, University of the Witwatersrand (2026)")
print("=" * 70)

print(f"\n--- Inputs ---")
print(f"  alpha       = 1/{1/alpha:.9f}")
print(f"  m_e         = {m_e_MeV} MeV")
print(f"  m_0         = {m_0_MeV:.4f} MeV")
print(f"  ell = r_e   = {r_e_fm:.4f} fm")
print(f"  N^2         = 1/pi = {N2:.6f}")

print(f"\n--- Derived neutrino parameters ---")
print(f"  theta_ch    = alpha^2/(2pi) = {theta_ch:.6e}")
print(f"  delta (NLO) = (1-alpha)^2/pi^2 = {delta:.6f}")
print(f"  h/A         = (pi+1)/(2pi) = {rho:.6f}")
print(f"  A           = {A*1e3:.4f} meV")
print(f"  h           = {h*1e3:.4f} meV")

print(f"\n{'=' * 70}")
print(f"ZERO-PARAMETER PREDICTIONS")
print(f"{'=' * 70}")
print(f"  m_1       = {m1*1e3:8.2f} +/- {err['m1']*1e3:.2f} meV")
print(f"  m_2       = {m2*1e3:8.2f} +/- {err['m2']*1e3:.2f} meV")
print(f"  m_3       = {m3*1e3:8.2f} +/- {err['m3']*1e3:.2f} meV")
print(f"  sum(m_nu) = {sum_nu*1e3:8.1f}  +/- {err['sum']*1e3:.1f}  meV")
print(f"  3A        = {3*A*1e3:8.1f}  meV  [trace theorem: sum = 3A, exact]")
print(f"")
print(f"  Dm^2_21   = ({Dm21_sq:.4e} +/- {err['Dm21']:.2e}) eV^2")
print(f"  Dm^2_31   = ({Dm31_sq:.4e} +/- {err['Dm31']:.2e}) eV^2")
print(f"")
print(f"  Mass ordering:    Normal (m1 < m2 < m3)")
print(f"  Dirac/Majorana:   Dirac (m_bb = 0)")

print(f"\n{'=' * 70}")
print(f"COMPARISON WITH EXPERIMENT")
print(f"{'=' * 70}")
print(f"  Dm^2_21:  pred {Dm21_sq:.4e},  JUNO {Dm21_exp:.4e} +/- {Dm21_err:.2e}")
print(f"            pull = {pull_21:+.2f} sigma")
print(f"  Dm^2_31:  pred {Dm31_sq:.4e},  NuFit {Dm31_exp:.4e} +/- {Dm31_err:.2e}")
print(f"            pull = {pull_31:+.2f} sigma")
print(f"  Combined chi^2 = {pull_21**2 + pull_31**2:.3f}  (2 dof)")

print(f"\n--- Mixing angles (monograph Ch. 13-14) ---")
print(f"  sin^2(theta_12) = 0.310   (pull +0.1 sigma vs JUNO)")
print(f"  sin^2(theta_23) = 0.556   (pull -0.3 sigma vs NuFit 6.0)")
print(f"  sin^2(theta_13) = 0.0222  (pull +0.0 sigma vs NuFit 6.0)")
print(f"  delta_CP         = 183 deg (pull +0.3 sigma vs NuFit 6.0)")

print(f"\n--- Error budget (monograph Sec. 4.7) ---")
print(f"  dm1/m1 = {dm1_rel*100:.1f}%  [O(alpha) corrections to theta_ch]")
print(f"  ddelta/delta = {ddelta_rel*100:.2f}%  [NLO screening ambiguity]")
print(f"  drho/rho = {drho_rel*100:.2f}%  [h/A derivation precision]")
print(f"  MC samples = {N_samples}")
