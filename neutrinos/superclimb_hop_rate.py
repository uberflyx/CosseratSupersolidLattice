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
    The hopping amplitude h is NOT set here. It is DERIVED in the neutrino
    chapter from the Cosserat overlap, h = (pi+1)/(2pi) A, with A = (1/3) sum m_nu
    the on-site energy (trace theorem) and (pi+1)/(2pi) = cos60 + N^2/2 the
    strain-plus-curvature overlap; it matches oscillation data to 0.04%. This
    script only asks whether the superclimb MECHANISM reproduces that h.

    RESULT: with v2 pinned at sqrt(40/3) c = 3.65c (condensate frozen, crystal
    carries a longitudinal wave against C11 = 8/3 mu), the mechanism predicts
    (v2/c) m1 ~ 18 meV against the overlap's h ~ 14.4 meV: agreement to ~25%,
    the precision of a semiclassical attempt frequency. 'v2 = 2.9c' is not an
    independent speed -- it is just h/m1, and h is pinned by the overlap. The
    strain channel of the overlap IS the Peach-Koehler climb force, so the two
    are one coupling. No fitted quantity enters.
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
print("HOW TO READ THIS")
print("="*70)
print(f"  bare hop h0 = h/theta_ch^2          : {h_meV*1e-9/theta_ch**2:.0f} MeV")
print(f"  as a multiple of m0                 : {h_meV*1e-9/theta_ch**2/m0_MeV:.2f} m0  (= h/m1)")
print(f"  two-fluid second sound              : {v2_over_c:.2f} m0")
print(f"  ratio                               : {v2_over_c/v2_implied:.2f}  (~{100*(v2_over_c/v2_implied-1):.0f}% apart)")
print()
print("  h is NOT set by v2. It is derived from the Cosserat overlap")
print("  h/A = (pi+1)/(2pi) (strain channel cos60 + rolling channel N^2/2),")
print("  matching oscillation data to 0.04%. So 'v2 = 2.9c' is not an")
print("  independent measurement -- it is just h/m1, and h is pinned by")
print("  the overlap. The strain channel IS the Peach-Koehler climb force,")
print("  so the overlap and the climb rate are one coupling, not two.")
print()
print("  VERDICT: the superclimb MECHANISM reproduces the overlap hop")
print("  amplitude to ~%d%%, the precision of a semiclassical attempt" % round(100*rel_gap))
print("  frequency. The exact amplitude is the overlap integral itself.")
print("  The bare hop landing at the second-sound scale confirms the hop")
print("  is a longitudinal lattice compression -- i.e. a climb step.")
