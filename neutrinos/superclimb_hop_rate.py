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

    RESULT: with v2 now pinned at sqrt(40/3) c = 3.65c (condensate frozen,
    crystal carries a longitudinal wave against C11 = 8/3 mu), the predicted
    hop rate is (v2/c) m1 ~ 18 meV against the fitted h ~ 14.4 meV: agreement
    to ~25%, inside the order-unity ambiguity of the attempt frequency and the
    tunnelling prefactor. Equivalently h/m1 implies v2 ~ 2.9c vs the two-fluid
    3.65c. No fitted quantity enters the prediction.
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

v2_over_c = np.sqrt(40.0/3.0)      # pinned second sound = sqrt(C11/(f_n rho)), eq:v2_value ~ 3.65c

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
#   hbar Gamma_climb = (v2/c) * m1, with v2 pinned by the two-fluid spectrum
Gam_pred   = v2_over_c*m1_meV                    # predicted hop rate [meV]
v2_implied = h_meV/m1_meV                        # v2/c the fitted h would imply
rel_gap    = abs(Gam_pred - h_meV)/h_meV         # fractional disagreement

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
print(f"  v2 (pinned, two-fluid)          : {v2_over_c:.3f} c")
print(f"  predicted hop rate (v2/c)*m1    : {Gam_pred:.2f} meV")
print(f"  Z3 hopping amplitude h          : {h_meV:.2f} meV")
print(f"  fractional disagreement         : {100*rel_gap:.0f}%")
print(f"  -> the two agree to ~{100*rel_gap:.0f}%, inside the order-unity")
print(f"     ambiguity of the attempt frequency and tunnelling prefactor.")

print()
print("="*70)
print("THE THREE-WAY CONSISTENCY")
print("="*70)
print(f"  second-sound speed implied by h/m1 : v2 = {v2_implied:.3f} c")
print(f"  second-sound speed, two-fluid      : v2 = {v2_over_c:.3f} c")
print(f"  ratio                              : {v2_over_c/v2_implied:.2f}  (~{100*(v2_over_c/v2_implied-1):.0f}% apart)")
print()
print("  Three quantities from three chapters land on one speed to ~25%:")
print("    - h          : Cosserat angular overlap   (neutrino oscillations)")
print("    - theta_ch   : self-energy series         (chirality)")
print("    - v2         : two-fluid spectrum         (superfluid sector)")
print()
print("  VERDICT: the superclimb rate and the Z3 hopping amplitude agree to")
print("  ~25%, well inside the estimate's order-unity ambiguity. The")
print("  conjecture (flavour hop = superclimb step) checks out. A fuller")
print("  attempt-frequency calculation would sharpen the 25% and turn this")
print("  into a genuine test of the C11 assignment for second sound.")
