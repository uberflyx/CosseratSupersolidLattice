#!/usr/bin/env python3
r"""
stellar_cooling_secondsound.py
================================
THE GATEKEEPER. The monograph's superfluid sector names one uncomputed
number as deciding the whole engineering programme: the stellar-cooling
bound on second-sound emission. This script scopes it with the framework's
own monopole rule and the measured solar neutrino luminosity, turning a
flagged worry into a quantitative pass/fail with an explicit margin.

THE RULE (Sec. secondsound_emission/constraints): second sound couples to a
nuclear process only through its MONOPOLE CONTENT epsilon = |delta-rho/rho|,
the net fractional change in defect (node) count. The coupling constant is
order unity; the smallness is all in epsilon. A process with epsilon ~ 1
(annihilation) drives at full strength; a pure multipole rearrangement
(gamma, most reactions) has epsilon ~ 0.

THE TEST: stars burn at the MeV scale, below the ~180 MeV second-sound gap,
so they CAN radiate the mode. If solar hydrogen burning carried order-unity
monopole content, second sound would drain the Sun at a rate rivalling its
neutrino luminosity, which helioseismology and neutrino fits forbid. The
framework survives only if the per-reaction monopole content is small.

WHAT SETS epsilon FOR REAL BURNING. The net reaction 4p -> 4He + 2e+ + 2nu
changes the DEFECT LEDGER only through the beta vertices (p<->n flips the
screw/edge partial character; Sec. open_em_correction), NOT through the
strong rearrangement that merely re-docks existing nodes into the alpha.
Each beta flip carries a monopole projection set by the chiral mixing that
ties the weak vertex to compression: epsilon_beta ~ theta_ch per flip, the
same second-order-protected amplitude that gives the neutrino its mass.
So epsilon_solar ~ theta_ch, NOT ~1. This is the physical input; the rest
is bookkeeping.

NORMALISATION: the second-sound luminosity is scaled from the neutrino
luminosity, L_2nd/L_nu = epsilon^2, with the order-unity coupling and the
phase-space/velocity factor set to one, and the exclusion is taken as
L_2nd < 0.1 L_nu (a tenth of the neutrino luminosity).

RESULT: with epsilon ~ theta_ch, the second-sound luminosity is
theta_ch^2 = 7.2e-11 of the solar neutrino luminosity, so the margin against
the bound is 0.1/7.2e-11 = 1.4e9: the bound is cleared by ~9 orders. Even
the conservative epsilon ~ theta_ch (not theta_ch^2) reading passes
handily, while an epsilon ~ 1 reading gives L_2nd = L_nu and overshoots the
bound tenfold, which is why the question was worth computing rather than
asserting.
"""

import numpy as np

# ---- constants -----------------------------------------------------------
alpha  = 1/137.035999
theta  = alpha**2/(2*np.pi)          # chirality, 8.5e-6
c      = 2.99792458e8
N2     = 1/np.pi                     # Cosserat coupling number
v2     = np.sqrt((3 - 4*N2)/(1 - N2))*c   # second-sound speed, 1.592 c; it does not
                                     # enter the ratio L_2nd/L_nu ~ epsilon^2 below

# ---- solar energy budget (measured) --------------------------------------
L_sun      = 3.828e26                # W, photon luminosity
# solar neutrino luminosity ~ 2% of photon luminosity (pp-chain, standard)
L_nu_sun   = 0.02*L_sun              # W carried by neutrinos
# the neutrino channel is the tight one: any exotic drain rivalling it is
# excluded by combined helioseismology + neutrino-flux fits to the
# ten-percent level. The bound used is L_exotic/L_nu < 0.1.
bound_frac = 0.1

print("="*72)
print("STELLAR-COOLING GATEKEEPER: second sound from solar burning")
print("="*72)
print(f"  solar photon luminosity   L_sun    = {L_sun:.3e} W")
print(f"  solar neutrino luminosity L_nu     = {L_nu_sun:.3e} W (~2% of L_sun)")
print(f"  exclusion (conservative)  L_2nd/L_nu < {bound_frac}")
print()

# ---- second-sound luminosity, scaled from the neutrino luminosity --------
# Both channels ride the SAME MeV nuclear power source. Second sound differs
# from the (already order-unity-coupled) neutrino channel by (a) its monopole
# content epsilon per reaction, entering the RATE as epsilon^2, and (b) an
# O(1) phase-space/velocity factor we set to 1 for a scoping estimate.
# So L_2nd / L_nu ~ epsilon^2  (both sourced by the same beta vertices).
for label, eps in [("epsilon ~ theta_ch (protected beta monopole)", theta),
                   ("epsilon ~ theta_ch^2 (fully second-order)",     theta**2),
                   ("epsilon ~ 1 (naive full-strength, the worry)",  1.0)]:
    ratio = eps**2
    L2 = ratio*L_nu_sun
    verdict = "PASS" if ratio < bound_frac else "FAIL"
    print(f"  {label}")
    print(f"     L_2nd/L_nu = epsilon^2 = {ratio:.2e}   L_2nd = {L2:.2e} W   -> {verdict}"
          f" (margin {bound_frac/ratio:.1e})" if ratio<bound_frac
          else f"     L_2nd/L_nu = epsilon^2 = {ratio:.2e}   L_2nd = {L2:.2e} W   -> {verdict}"
               f" (over by {ratio/bound_frac:.1e})")
print()
print("="*72)
print("WHY epsilon ~ theta_ch, NOT 1")
print("="*72)
print("  The net solar reaction 4p -> He4 + 2e+ + 2nu re-docks existing nodes")
print("  into an alpha (a MULTIPOLE strong rearrangement, epsilon~0) plus two")
print("  beta flips (p<->n). Only the beta flips change the screw/edge partial")
print("  ledger, and each carries the chiral monopole projection theta_ch that")
print("  ties the weak vertex to compression (the same amplitude that makes the")
print("  neutrino massive). So the star's monopole content is chiral-suppressed,")
print("  which is exactly what lets the framework clear its own gatekeeper.")
print()
margin = bound_frac/theta**2
print(f"  BASELINE FLIPPED (L_2nd/L_nu = epsilon^2, bound L_2nd < {bound_frac} L_nu):")
print(f"  the worry (epsilon~1) overshoots the bound {1.0/bound_frac:.0f}-fold; the physics")
print(f"  (epsilon~theta_ch) clears it with margin {bound_frac}/{theta**2:.1e} = {margin:.1e},")
print(f"  ~{int(np.floor(np.log10(margin)))} orders. The gate turns on whether burning is")
print("  chiral-suppressed, and the framework says it is.")
print("="*72)
