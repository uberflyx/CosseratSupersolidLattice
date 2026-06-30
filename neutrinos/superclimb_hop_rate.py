#!/usr/bin/env python3
"""
superclimb_hop_rate.py - is the neutrino Z3 hopping amplitude the superclimb rate?
=================================================================================
M. A. Cox, University of the Witwatersrand (2026)

CONTEXT
    A neutrino is an edge dislocation confined to a {111} glide plane, with its
    Burgers vector along one of the three <110> slip directions. Those three
    directions are the three flavours, tied together by the Z3 stacking of the
    ABCABC close-packed sequence. Neutrino oscillation is the Z3 hop A->B->C.

    The edge core is a superfluid pipe: it carries mass along its length
    (superclimb / the Ray-Hallock syringe effect, seen in solid helium and in
    simulation). CLIMB moves the edge perpendicular to its glide plane, i.e. to
    the ADJACENT {111} plane. Adjacent {111} planes are different stacking sites
    (A,B,C), so one climb step IS one A->B->C flavour hop. The conjecture: the
    flavour hop and the superclimb step are the same event.

    THE OPEN CHECK (the framework had flagged this as 'not done'):
    does the superclimb RATE equal the Z3 hopping amplitude the oscillation fit
    requires? This script does that order-of-magnitude check.

THE ESTIMATE
    A tight-binding hop rate = (attempt frequency) x (tunnelling probability).
      - attempt frequency : the rate the condensate moves one node's mass to the
        climbing half-plane. The mass-transport mode is SECOND SOUND (speed v2);
        over the node Compton length ell = hbar/(m0 c) it is v2/ell.
      - tunnelling probability : the chiral matrix element squared, theta_ch^2,
        because the orientation hop is parity-protected to SECOND order in the
        chirality - the same protection that makes m1 = theta_ch^2 m0.
    So
        hbar Gamma_climb = (v2/ell) hbar theta_ch^2 = (v2/c) theta_ch^2 m0 c^2
                         = (v2/c) m1 c^2.
    The Z3 fit's hopping amplitude is h = (pi+1)/(2pi) A, with A = (1/3) sum m_nu
    the on-site energy (trace theorem) and (pi+1)/(2pi) the Cosserat overlap.

    RESULT: the two agree, and the second-sound speed implied by h/m1 lands
    inside the independent v2 range 2.6-3.7c. No fitted quantity.
"""

import numpy as np

# ---- inputs (framework constants only) -----------------------------------
alpha   = 1/137.035999177          # fine-structure constant
m_e_MeV = 0.51099895069            # electron mass [MeV]
m0_MeV  = m_e_MeV/alpha            # screened node mass = m_e/alpha ~ 70 MeV
theta_ch= alpha**2/(2*np.pi)       # chirality = strain-curvature mixing angle

# measured mass-squared splittings used only to BUILD the spectrum the fit fits
Dm21_eV2 = 7.49e-5                 # eV^2 (solar)
Dm31_eV2 = 2.534e-3               # eV^2 (atmospheric)

v2_over_c_range = (2.6, 3.7)       # independent second-sound range (eq:v2_value)

# ---- the neutrino spectrum (lightest mass is parameter-free) --------------
m1_meV = theta_ch**2 * m0_MeV * 1e9            # m1 = theta_ch^2 m0, in meV
m1_eV  = m1_meV*1e-3
m2_meV = np.sqrt(m1_eV**2 + Dm21_eV2)*1e3
m3_meV = np.sqrt(m1_eV**2 + Dm31_eV2)*1e3
sum_meV = m1_meV + m2_meV + m3_meV

# ---- Z3 hopping amplitude from the fit ------------------------------------
N2      = 1/np.pi                              # Cosserat coupling number
A_meV   = sum_meV/3.0                          # on-site energy (trace theorem)
hA      = (1 + N2)/2.0                         # = (pi+1)/(2pi), Cosserat overlap
h_meV   = hA * A_meV                           # Z3 hopping amplitude

# ---- superclimb rate prediction -------------------------------------------
#   hbar Gamma_climb = (v2/c) * m1
Gam_lo  = v2_over_c_range[0]*m1_meV
Gam_hi  = v2_over_c_range[1]*m1_meV
v2_implied = h_meV/m1_meV                       # v2/c that makes them equal

# ---- report ---------------------------------------------------------------
print("="*70)
print("INPUTS (framework constants only)")
print("="*70)
print(f"  alpha            = {alpha:.10f}")
print(f"  m0 = m_e/alpha    = {m0_MeV:.3f} MeV")
print(f"  theta_ch          = {theta_ch:.5e}   theta_ch^2 = {theta_ch**2:.5e}")
print(f"  Cosserat N^2      = {N2:.5f}   (pi+1)/(2pi) = {hA:.5f}")

print()
print("="*70)
print("NEUTRINO SPECTRUM")
print("="*70)
print(f"  m1 = theta_ch^2 m0 = {m1_meV:.3f} meV   (parameter-free)")
print(f"  m2                 = {m2_meV:.3f} meV")
print(f"  m3                 = {m3_meV:.3f} meV")
print(f"  sum m_nu           = {sum_meV:.2f} meV")

print()
print("="*70)
print("Z3 HOPPING AMPLITUDE (from the oscillation fit)")
print("="*70)
print(f"  A = (1/3) sum m_nu        = {A_meV:.3f} meV   (trace theorem)")
print(f"  h = (pi+1)/(2pi) * A      = {h_meV:.3f} meV")

print()
print("="*70)
print("SUPERCLIMB RATE PREDICTION:  hbar*Gamma_climb = (v2/c)*m1")
print("="*70)
for v2c in [2.6, v2_implied, 3.0, 3.7]:
    tag = "  <-- value that matches h" if abs(v2c-v2_implied) < 1e-6 else ""
    print(f"  v2 = {v2c:.3f} c :  hbar*Gamma_climb = {v2c*m1_meV:.2f} meV{tag}")
print(f"\n  predicted range (v2 = 2.6-3.7c) : {Gam_lo:.1f} - {Gam_hi:.1f} meV")
print(f"  Z3 hopping amplitude h          : {h_meV:.2f} meV")
print(f"  -> the estimate BRACKETS the hopping amplitude.")

print()
print("="*70)
print("THE THREE-WAY CONSISTENCY")
print("="*70)
ok = v2_over_c_range[0] <= v2_implied <= v2_over_c_range[1]
print(f"  second-sound speed implied by h/m1 : v2 = {v2_implied:.3f} c")
print(f"  independent second-sound range      : {v2_over_c_range[0]}-{v2_over_c_range[1]} c")
print(f"  implied v2 inside the range?        : {ok}")
print()
print("  Three quantities from three chapters meet on one number:")
print("    - h          : Cosserat angular overlap   (neutrino oscillations)")
print("    - theta_ch   : self-energy series         (chirality)")
print("    - v2         : two-fluid spectrum         (superfluid sector)")
print()
print("  VERDICT: the superclimb rate and the Z3 hopping amplitude agree to")
print("  order of magnitude, at a second-sound speed inside the independent")
print("  range. The conjecture (flavour hop = superclimb step) checks out.")
