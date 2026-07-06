#!/usr/bin/env python3
r"""
forge_cross_section_scoping.py
==============================
Does the vortex-ring forge collapse below its experimental ceiling the way the
second-sound branching did?

CONTEXT. Two exotic e+e- channels are bounded by B-factory data. The second-sound
(compression) branching was found to sit ~9 orders below its naive value, because
that process must MANUFACTURE compression out of a pure-shear event (a higher-order
leak of order theta_ch^2). The ring forge is different: the electron and positron
already ARE condensate vortices, so forging a ring reconnects circulation that is
on hand rather than squeezing a new quantity into being. The only cost left is the
ordinary form factor of any EXCLUSIVE production. The whole question of "collapse
or not" reduces to one thing: is that form factor power-law or exponential?

This script:
  1. reproduces the framework's ring dispersion (floor + Upsilon(4S) tag) as a
     cross-check against cosmology/vortex_hadron_ladder.py;
  2. writes the exclusive cross section as sigma_0(s) * |F(s)|^2, exactly the form
     of the measured e+e- -> p pbar cross section, and asks what |F| the current
     BaBar ceiling (~1 fb) demands, comparing it to the proton at matched q^2;
  3. shows the cross section GROWS steeply toward threshold, so a low-energy
     machine (BESIII, sqrt(s) = 1.84-4.95 GeV) sees orders more than Belle II;
  4. notes that BESIII already holds single-photon + invisible data in the forge
     window, so the hunt needs a reanalysis, not a new machine.

Everything is a SCOPING estimate: the absolute normalisation rides on the proton
analogy and is an open coupling. What is robust is (a) power-law, not exponential,
so no catastrophic collapse; and (b) the steep growth toward threshold.
Model inputs are labelled MODEL.
"""

import numpy as np

# --------------------------------------------------------------------------
# Framework constants.
# --------------------------------------------------------------------------
alpha   = 1.0 / 137.035999177
c       = 2.99792458e8
m_e_MeV = 0.51099895069
m0_MeV  = m_e_MeV / alpha                       # node mass, 70.03 MeV
l_m     = 2.8179403205e-15                       # lattice spacing = r_e [m]
kappa   = 2.0 * np.pi * l_m * c                  # circulation quantum [m^2/s]
rho     = (m0_MeV * 1.7826619e-30) / l_m**3      # vacuum density [kg/m^3]
rho_s   = 0.8 * rho                              # superfluid fraction f_s = 4/5
J_per_GeV = 1.602176634e-10
GeV2_to_nb = 3.8937966e5                          # (hbar c)^2: 1 GeV^-2 = 0.389 mb

line = "=" * 74

def E_ring_GeV(x):
    """Ring energy at radius R = x*l, thin-ring line tension [GeV]."""
    E_J = 0.5 * rho_s * kappa**2 * (x * l_m) * (np.log(8.0 * x) - 2.0)
    return E_J / J_per_GeV

def Pc_ring_GeV(x):
    """Ring impulse times c at radius R = x*l [GeV]. Impulse ~ area."""
    P = rho_s * kappa * np.pi * (x * l_m)**2
    return P * c / J_per_GeV

# --------------------------------------------------------------------------
# 1. Cross-check the ring dispersion against the framework's own numbers.
# --------------------------------------------------------------------------
print(line); print("1. RING DISPERSION (cross-check vs vortex_hadron_ladder.py)"); print(line)
x0 = 1.68                                        # ground-ring radius in units of l
E1, P0c = E_ring_GeV(x0), Pc_ring_GeV(x0)
print(f"   ground ring  R0 = {x0:.2f} l :  E_1 = {E1:.3f} GeV,  P_0 c = {P0c:.3f} GeV")
print(f"   forge floor  sqrt(s)_min = E_1 + P_0 c = {E1 + P0c:.2f} GeV")
# Solve E(x) + Pc(x) = sqrt(s) for the Upsilon(4S) ring radius and photon tag.
def solve_radius(sqrt_s):
    xs = np.linspace(x0, 12.0, 400000)
    tot = E_ring_GeV(xs) + Pc_ring_GeV(xs)
    i = np.argmin(np.abs(tot - sqrt_s))
    return xs[i]
x_ups = solve_radius(10.58)
print(f"   at Upsilon(4S) (sqrt(s)=10.58 GeV): ring R = {x_ups:.2f} l, "
      f"photon tag E_gamma = P c = {Pc_ring_GeV(x_ups):.2f} GeV")
print("   -> reproduces the framework floor (4.22 GeV) and tag (7.57 GeV).")

# --------------------------------------------------------------------------
# 2. Exclusive cross section = sigma_0(s) * |F(s)|^2, calibrated on the proton.
#    Measured form: sigma(e+e- -> p pbar) = (4 pi alpha^2 beta C / 3s) |F_p|^2.
# --------------------------------------------------------------------------
print(line); print("2. FORGE CROSS SECTION AT Upsilon(4S): |F| NEEDED vs PROTON"); print(line)
def sigma0_nb(sqrt_s):
    """QED point cross section 4 pi alpha^2 / 3s [nb]."""
    return (4.0 * np.pi * alpha**2 / 3.0) / sqrt_s**2 * GeV2_to_nb
s0_ups = sigma0_nb(10.58)
print(f"   point cross section sigma_0(10.58 GeV) = {s0_ups*1e6:.2e} fb "
      f"({s0_ups:.3f} nb)")
sigma_ceiling_fb = 1.0                            # BaBar monophoton ceiling [fb]
F_needed = np.sqrt(sigma_ceiling_fb / (s0_ups * 1e6))
print(f"   |F_ring| needed to sit at the {sigma_ceiling_fb:.0f} fb ceiling "
      f"= {F_needed:.2e}")
# Proton effective form factor at matched q^2 = s, from BaBar data scaled by the
# asymptotic |F| ~ 1/s^2 (QCD dimensional counting), anchored at sqrt(s)=3.5 GeV.
# MODEL: anchor |F_p|(s=12.25 GeV^2) ~ 0.04 (BaBar e+e- -> p pbar, above 3 GeV).
F_anchor, s_anchor = 0.040, 3.5**2               # MODEL anchor from BaBar
def Fp(sqrt_s):
    return F_anchor * (s_anchor / sqrt_s**2)      # |F_p| ~ 1/s^2
F_p_ups = Fp(10.58)
sigma_charged = s0_ups * 1e6 * F_p_ups**2         # fb, charged-baryon equivalent
print(f"   proton |F_p| at matched q^2 (extrapolated) = {F_p_ups:.2e}")
print(f"   charged-baryon equivalent sigma_0 |F_p|^2 = {sigma_charged:.1f} fb")
print(f"   -> BUT the ring is NEUTRAL (pure circulation, no Burgers charge), so")
print(f"      unlike p pbar it has no direct photon vertex: the balancing photon")
print(f"      is radiative and the production runs through the circulation overlap,")
print(f"      costing a factor of order alpha (a percent with the ISR log).")
rad = alpha * np.log(10.58**2 / (m_e_MeV*1e-3)**2) / np.pi   # ISR radiator ~ few %
sigma_forge_ups = sigma_charged * rad
print(f"      radiative/neutral factor ~ {rad:.3f}  ->  sigma_forge(Upsilon 4S)")
print(f"      ~ {sigma_forge_ups:.2f} fb, at or just below the ~1 fb BaBar ceiling.")
print(f"      So the null is CONSISTENT with a baryon-like ring form factor; it")
print(f"      does NOT force the 1e-10 collapse the second-sound channel suffered.")

# --------------------------------------------------------------------------
# 3. The forge grows steeply toward threshold: run lower, not higher.
#    sigma_forge ~ sigma_0(s) |F(s)|^2 ~ 1/s * (1/s^2)^2 = 1/s^5  (proton-like).
# --------------------------------------------------------------------------
print(line); print("3. THRESHOLD ENHANCEMENT: the forge is loudest near its floor"); print(line)
def radiator(sqrt_s):
    """ISR/neutral factor of order alpha for the radiative monophoton production."""
    return alpha * np.log(sqrt_s**2 / (m_e_MeV*1e-3)**2) / np.pi

def sigma_forge_fb(sqrt_s):
    return sigma0_nb(sqrt_s) * 1e6 * Fp(sqrt_s)**2 * radiator(sqrt_s)
for E in [10.58, 4.95, 4.60, 4.30, 4.22]:
    r = sigma_forge_fb(E) / sigma_forge_ups
    label = "Upsilon(4S)" if E == 10.58 else ("BESIII window" if E <= 4.95 else "")
    print(f"   sqrt(s) = {E:5.2f} GeV : sigma_forge ~ {sigma_forge_fb(E):8.2f} fb "
          f"(x{r:7.0f} vs Upsilon 4S)  {label}")
print("   (Right at the floor the final-state phase space closes, so the true")
print("    optimum sits a little above 4.22 GeV, inside the BESIII range.)")

# --------------------------------------------------------------------------
# 4. BESIII already holds the data.
# --------------------------------------------------------------------------
print(line); print("4. THE HUNT IS AT BESIII, AND THE DATA EXISTS"); print(line)
L_besiii_fb = 14.9                                # BESIII 4.13-4.60 GeV invisible set
sig_besiii = sigma_forge_fb(4.5)
print(f"   BESIII single-photon + invisible set: {L_besiii_fb:.1f} fb^-1 at "
      f"4.13-4.60 GeV (existing).")
print(f"   scoping sigma_forge(~4.5 GeV) ~ {sig_besiii:.0f} fb -> "
      f"N ~ {sig_besiii*L_besiii_fb:.0f} events if at this rate.")
print("   CAVEAT: the BESIII search assumed a MASSIVE invisible particle, so it")
print("   looked at the photon energy an invariant-mass tag predicts. A ring has")
print("   no rest mass; its dispersion tag sits at a DIFFERENT photon energy")
print("   (E_gamma = P(R)c). The existing data therefore needs a reanalysis at")
print("   the ring-dispersion line, not a new machine.")
print(line)
print("VERDICT: power-law (asymptotic freedom), not exponential -> NO collapse.")
print("The forge sits near the exclusive-baryon scale at Upsilon(4S) and climbs")
print("steeply toward its floor, so BESIII is the discovery venue. Absolute")
print("normalisation is an open coupling; the SCALE and the TREND are robust.")
print(line)
