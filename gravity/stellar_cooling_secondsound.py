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
monopole content, second sound would drain the Sun at about its neutrino
luminosity, L_nu = 0.0242 L_sun. The measured 8B neutrino flux caps any
novel energy loss at 0.1 L_sun (Gondolo & Raffelt, PRD 79, 107301 (2009),
Eq. 3), so an order-unity drain sits at about a quarter of the bound: the
Sun alone does not decide it, and the per-reaction monopole content does.

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
phase-space/velocity factor set to one, and the exclusion is the 8B bound
L_2nd < 0.1 L_sun. The neutrino luminosity is computed from the luminosity
constraint, L_nu/L_sun = sum_i <E_nu>_i Phi_i / sum_i alpha_i Phi_i, with
<E_nu>_i and alpha_i from Bahcall, PRC 65, 025801 (2002), Table I, and the
GS98 reference fluxes Phi_i of Bergstrom et al., JHEP 03 (2016) 132, Table 1;
the denominator reproduces the solar constant to 0.5%.

RESULT: with epsilon ~ theta_ch, the second-sound luminosity is
theta_ch^2 = 7.2e-11 of the solar neutrino luminosity, i.e. 1.74e-12 L_sun,
so the margin against the bound is 0.1/1.74e-12 = 5.7e10: the bound is
cleared by nearly 11 orders. An epsilon ~ 1 reading gives L_2nd = L_nu, a
quarter of the bound, close enough that the unknown O(1) coupling could
carry it across; the chiral-suppressed reading is where no O(1) factor can
matter, which is why the question was worth computing rather than asserting.
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
L_sun      = 3.828e26                # W, photon luminosity (IAU 2015 nominal)
# Neutrino luminosity from the luminosity constraint:
#   L_sun/(4 pi AU^2) = sum_i alpha_i Phi_i,  neutrino energy flux = sum_i <E>_i Phi_i.
# alpha_i, <E>_i [MeV]: Bahcall, PRC 65, 025801 (2002), Table I (17F <E> not
# tabulated there; its flux is negligible and it is bracketed out below).
# Phi_i [cm^-2 s^-1]: GS98 reference fluxes, Bergstrom et al., JHEP 03 (2016) 132, Table 1.
phi   = dict(pp=5.98e10, Be7=5.00e9, pep=1.44e8, N13=2.96e8, O15=2.23e8,
             F17=5.52e6, B8=5.58e6, hep=8.04e3)
alpha_i = dict(pp=13.0987, Be7=12.6008, pep=11.9193, N13=3.4577, O15=21.5706,
               F17=2.3630, B8=6.6305, hep=3.7370)
E_nu  = dict(pp=0.2668, Be7=0.814, pep=1.445, N13=0.706, O15=0.996, F17=0.0,
             B8=6.735, hep=9.628)
nu_flux  = sum(phi[k]*E_nu[k] for k in phi)        # MeV cm^-2 s^-1 in neutrinos
lum_flux = sum(phi[k]*alpha_i[k] for k in phi)     # MeV cm^-2 s^-1 into the star
S_const  = 1361.0*1e-4/1.602176634e-13             # solar constant in MeV cm^-2 s^-1
nu_frac  = nu_flux/lum_flux                         # L_nu/L_sun = 0.0242
L_nu_sun = nu_frac*L_sun
# Exclusion: the all-flavour 8B flux caps any novel energy loss at
# L_x < 0.1 L_sun (Gondolo & Raffelt, PRD 79, 107301 (2009), Eq. 3); the older
# helioseismic bound is 0.2 L_sun. The bound is on the PHOTON luminosity.
bound_frac = 0.1                                    # L_2nd/L_sun < 0.1

print("="*72)
print("STELLAR-COOLING GATEKEEPER: second sound from solar burning")
print("="*72)
print(f"  luminosity-constraint check: sum alpha_i Phi_i = {lum_flux:.4e} MeV/cm2/s"
      f" vs solar constant {S_const:.4e} ({(lum_flux/S_const - 1)*100:+.1f}%)")
print(f"  solar photon luminosity   L_sun    = {L_sun:.3e} W")
print(f"  solar neutrino luminosity L_nu     = {L_nu_sun:.3e} W ({nu_frac*100:.2f}% of L_sun)")
print(f"  exclusion (8B flux)       L_2nd/L_sun < {bound_frac}")
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
    ratio = eps**2                       # L_2nd/L_nu
    frac  = ratio*nu_frac                # L_2nd/L_sun
    L2 = ratio*L_nu_sun
    verdict = "PASS" if frac < bound_frac else "FAIL"
    print(f"  {label}")
    print(f"     L_2nd/L_nu = {ratio:.2e}, L_2nd/L_sun = {frac:.2e}, L_2nd = {L2:.2e} W"
          f"   -> {verdict} (margin {bound_frac/frac:.1e})" if frac < bound_frac
          else f"     L_2nd/L_nu = {ratio:.2e}, L_2nd/L_sun = {frac:.2e}, L_2nd = {L2:.2e} W"
               f"   -> {verdict} (over by {frac/bound_frac:.1e})")
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
margin = bound_frac/(theta**2*nu_frac)
print(f"  BASELINE (L_2nd = epsilon^2 L_nu, bound L_2nd < {bound_frac} L_sun):")
print(f"  the worry (epsilon~1) sits at {nu_frac/bound_frac:.2f} of the bound, inside it only by a")
print(f"  factor the unknown O(1) coupling could close; the physics (epsilon~theta_ch)")
print(f"  clears it with margin {bound_frac}/({theta**2:.2e} x {nu_frac:.4f}) = {margin:.1e},")
print(f"  nearly {int(np.round(np.log10(margin)))} orders. The Sun alone does not decide the")
print("  order-unity case; whether burning is chiral-suppressed does, and the")
print("  framework says it is.")
print("="*72)
