#!/usr/bin/env python3
r"""
fluid_crystal_defect_energy.py
================================
How the SAME defect winding stores energy in the TWO orders of the supersolid.

The electron is one topological object read two ways (monograph Sec. on the
electron, "the Burgers vector and the circulation are the same winding read in
the two orders"):

    * CRYSTAL reading  -> a screw dislocation, Burgers vector b = ell,
      elastic self-energy per length  E_cr/L = (mu b^2 / 4pi) ln(R/r_c).

    * CONDENSATE reading -> a quantised vortex, circulation kappa = h/m0,
      superflow kinetic energy per length  E_sf/L = (rho_s kappa^2 / 4pi) ln(R/xi).

Same winding number (charge = circulation, one integer). But the ENERGY stored
in that winding is NOT the same in the two orders. This script computes the
ratio, using only quantities the framework already fixes:

    ell = hbar/(m0 c)      (Mott bootstrap)      -> hbar = m0 c ell
    mu  = rho c^2,  rho = m0/ell^3               (shear modulus from c)
    kappa = h/m0 = 2 pi hbar / m0
    rho_s = f_s rho,  f_s = 4/5                   (superfluid fraction bootstrap)

Taking the two logarithms equal (both cutoffs are the single core scale ell,
r_c ~ xi ~ ell), the log cancels and the ratio is a pure number.

This is NOT the electron's inertial mass. The inertial (rest) mass is the
Peierls-Nabarro KINK mass, m_e = alpha m0, the barrier to moving the defect one
lattice valley. The condensate is a superfluid and puts up no such barrier, so
it adds (almost) nothing to the kink mass. What the ratio below shows is which
order carries the defect's stored FIELD energy -- and it is overwhelmingly the
condensate. Crystal sets the inertia; condensate carries the field. That is the
division of labour, made quantitative.
"""

import numpy as np

f_s = 4.0 / 5.0          # superfluid fraction (bootstrap)
alpha = 1.0 / 137.035999 # fine-structure constant (for the mass-scale sanity line)

# ---------------------------------------------------------------------------
# Energy-per-length prefactors (the ln(R/core) is common and cancels).
#
# Crystal screw (Burgers b = ell), mu = rho c^2, rho = m0/ell^3:
#   E_cr/L = (mu b^2 / 4pi) ln = (m0 c^2 / ell^3)(ell^2)/(4pi) ln
#          = (m0 c^2 / ell) * (1/(4pi)) * ln
#   -> prefactor in units of (m0 c^2 / ell):  P_cr = 1/(4pi)
#
# Condensate vortex, kappa = 2 pi hbar/m0, hbar = m0 c ell  => kappa = 2 pi c ell,
#   kappa^2 = 4 pi^2 c^2 ell^2,  rho_s = f_s m0/ell^3:
#   E_sf/L = (rho_s kappa^2 / 4pi) ln = (f_s m0/ell^3)(4 pi^2 c^2 ell^2)/(4pi) ln
#          = (m0 c^2 / ell) * (f_s * pi) * ln
#   -> prefactor in units of (m0 c^2 / ell):  P_sf = f_s * pi
# ---------------------------------------------------------------------------
P_cr = 1.0 / (4.0 * np.pi)
P_sf = f_s * np.pi

ratio = P_sf / P_cr           # E_sf / E_cr at equal logs
ratio_symbolic = 4.0 * np.pi**2 * f_s

print("=" * 70)
print("FLUID vs CRYSTAL: energy stored in the SAME electron winding")
print("=" * 70)
print(f"  superfluid fraction         f_s            = {f_s:.4f}")
print()
print(f"  crystal prefactor  P_cr  [m0 c^2/ell]      = 1/(4pi)  = {P_cr:.5f}")
print(f"  condensate prefactor P_sf [m0 c^2/ell]     = f_s*pi   = {P_sf:.5f}")
print("-" * 70)
print(f"  ENERGY RATIO  E_sf / E_cr  (equal logs)    = {ratio:.3f}")
print(f"  closed form   4 pi^2 f_s                   = {ratio_symbolic:.3f}")
assert abs(ratio - ratio_symbolic) < 1e-9
print()
print("  => the condensate carries ~%.0fx the field energy the crystal does," % ratio)
print("     for one and the same winding. Same topology, different energy.")
print()

# --- sanity: neither field energy is the inertial mass -------------------
# Bare field energy (crystal, ln ~ 4pi so the prefactor 1/(4pi) gives ~ m0 c^2):
# the inertial mass is instead the PN kink mass alpha*m0, ~137x lighter.
print("-" * 70)
print("  Sanity on the mass scale (why the field energy is not the mass):")
print(f"    PN kink (inertial) mass   m_e = alpha m0  = {alpha:.5f} m0")
print(f"    both field energies are O(m0 c^2), i.e. ~1/alpha = {1/alpha:.1f}x heavier")
print("    than the rest mass. They are renormalised self-energy, as the")
print("    classical electromagnetic self-energy of a charge is, not inertia.")
print("    The superfluid puts up no Peierls barrier, so it moves the defect")
print("    for free and contributes ~nothing to the kink mass. The crystal's")
print("    periodic potential is what sets the inertia at alpha m0.")
print()

# --- the two faces of the coupling, for context --------------------------
print("=" * 70)
print("CONTEXT: the coupling has two faces, only one is small")
print("=" * 70)
g_over_Ksf = 1.0e-40      # g ~ mu, K_sf ~ 1e40 mu  (monograph eq:g_value scale)
print("  RADIATIVE face (mode mixing, second-sound shift, bremsstrahlung):")
print(f"    correction ~ g^2/(K_c K_sf) ~ mu/K_sf ~ {g_over_Ksf:.0e}   (tiny)")
print("  REACTIVE face (near-field Magnus grip on a charged vortex):")
print("    exceeds gravity by ~1e45   (the framework's 'full strength')")
print("  The 125-order gap between them is why 'insignificant' is only ever")
print("  true of the radiative channel. The reactive channel dominates.")
print("=" * 70)
