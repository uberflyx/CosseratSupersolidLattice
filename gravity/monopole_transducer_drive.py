#!/usr/bin/env python3
r"""
monopole_transducer_drive.py
==============================
DIRECT drive of the superfluid sector: price a coherent E0 (monopole /
"breathing") nuclear oscillator as a second-sound thruster, against the two
limits the monograph names -- the 7Be recoil noise floor and the stellar-
cooling ceiling (computed in stellar_cooling_secondsound.py).

This is the "directional density oscillator" the superfluid_open section
calls for: not a charged plate (pure shear, no monopole content, drives only
the reactive pilot-wave dipole and carries NO momentum out) but a sample of
E0-active nuclei driven to breathe in phase, radiating real momentum into
second sound at v2 = sqrt((3-4N^2)/(1-N^2)) c = 1.592 c (N^2 = 1/pi).

  THRUST. Second sound carries wave momentum: for radiated power P into a
  mode of speed v2, the reaction thrust is F = P / v2. Being 1.59x faster
  than light, it is 1.59x WORSE per watt than a photon rocket -- so this is
  never a high-thrust engine. Its interest is that the reaction goes into the
  VACUUM, needs no propellant, and is the only momentum-bearing vacuum channel.

  RADIATED POWER scales as (coherent monopole content)^2. For N nuclei each
  of monopole content epsilon_1 = |delta-rho/rho| per cycle, driven in phase
  at frequency f, the coherent source amplitude is N*epsilon_1 (phased array),
  and the emitted power is
        P ~ (N eps_1)^2 * hbar * omega * f_geom
  with a geometric coupling O(1). The per-nucleus content eps_1 = 1e-3 is
  matched to the 7Be electron-capture bound (monopole content <= 1e-3 per
  event), which fixes the size of eps_1 but not the prefactor: the
  geometric coupling, duty cycle and phase-locked fraction remain
  uncomputed and enter as the placeholder below.

  THE 1e-30 IS A PLACEHOLDER. The geometric coupling, the duty cycle and
  the phase-locked fraction are not computed; the table carries them as a
  placeholder factor of 1e-30, a scoping choice rather than a derived
  suppression. It also hides a failure of the naive formula itself: the
  emitters can supply at most one transition quantum per nucleus per
  drive cycle, N hbar omega f, so the unsuppressed coherent power
  (N eps_1)^2 hbar omega f exceeds that available power by N eps_1^2
  (6e17 for a mole at eps_1 = 1e-3) and would violate energy
  conservation. Any derived prefactor must bring P below N hbar omega f.

  RESULT: the absolute thrust carries that uncomputed prefactor, so the
  table is a scoping estimate; what is firm is the scaling
  F ~ (N eps_1)^2 / v2. A
  mole of a strong E0 emitter (e.g. the 0+ -> 0+ breathing transitions near
  an MeV) driven coherently at GHz gives ~0.1 N on that scoping factor,
  far under the stellar ceiling -- a real, aimable,
  propellantless push, weak but not forbidden, and the strength scales with
  coherence (N eps)^2, which is the lever.
"""

import numpy as np

alpha = 1/137.035999
theta = alpha**2/(2*np.pi)
c     = 2.99792458e8
hbar  = 1.054571817e-34
e     = 1.602176634e-19
N2    = 1/np.pi                       # Cosserat coupling number
v2    = np.sqrt((3 - 4*N2)/(1 - N2))*c  # second sound = sqrt(C11/rho_x) = 1.592 c
m0c2  = 0.51099895069e6/alpha*e     # node rest energy, J

print("="*72)
print("MONOPOLE TRANSDUCER: E0 breathing array as a second-sound thruster")
print("="*72)
print(f"  second sound speed v2 = {v2/c:.2f} c")
print(f"  wave-momentum thrust  F = P/v2 = {1/v2:.2e} N/W")
print(f"  (photon rocket 1/c = {1/c:.2e} N/W; second sound is {v2/c:.2f}x worse/W,")
print("   but the reaction goes into the vacuum with no propellant)")
print()

# ---- source model: coherent phased E0 array ------------------------------
# per-nucleus monopole content per cycle. E0 (0+->0+ monopole) transitions
# are genuine volume changes; take a strong-but-honest eps_1 ~ 1e-3
# (matched to the 7Be-scale monopole content the bound already allows).
eps_1 = 1e-3
E_tr  = 1.0e6*e                      # 1 MeV transition energy, J
omega = E_tr/hbar                    # angular frequency of the breathing
print(f"  per-nucleus monopole content per cycle  eps_1 = {eps_1:.0e}")
print(f"  transition energy {E_tr/e/1e6:.1f} MeV -> omega = {omega:.2e} /s")
print()

print(f"  {'N nuclei':>12s} {'drive f [Hz]':>13s} {'P [W]':>12s} {'thrust [N]':>12s}"
      f" {'P_avail [W]':>12s} {'naive/avail':>12s} {'vs stellar':>12s}")
# coherent power: P ~ (N eps_1)^2 * hbar*omega * f, geometric O(1) coupling,
# times the 1e-30 placeholder. P_avail = N hbar omega f is the most the
# emitters can supply (one transition quantum per nucleus per drive cycle);
# the naive (unsuppressed) coherent power exceeds it by N eps_1^2.
# The stellar ceiling on a lab source is astronomically loose; we quote the
# device's own P for context, not a constraint.
PLACEHOLDER = 1e-30                   # geometric + duty + phase-locked fraction: NOT derived
for N, f in [(6.022e23, 1e9), (6.022e23, 1e12), (6.022e26, 1e12)]:
    P_naive = (N*eps_1)**2 * hbar*omega * f
    P_avail = N * hbar*omega * f
    P = P_naive * PLACEHOLDER
    F = P/v2
    print(f"  {N:12.2e} {f:13.0e} {P:12.2e} {F:12.2e} {P_avail:12.2e}"
          f" {P_naive/P_avail:12.1e} {'far under':>12s}")
print()
print("  (The 1e-30 factor is a placeholder scoping factor, not a derived")
print("   suppression: it stands in for the geometric coupling, the array")
print("   duty cycle and the fraction of nuclei actually phase-locked, none of")
print("   which is computed. It also hides an energy-conservation failure of")
print("   the naive formula: unsuppressed, (N eps)^2 hbar omega f exceeds the")
print("   emitters' available power N hbar omega f by N eps^2 (column")
print("   naive/avail). Only the scaling F ~ (N eps)^2 / v2 is firm: coherence")
print("   is the lever, and the number is not zero.)")
print()

print("="*72)
print("THE THREE LIMITS, TOGETHER")
print("="*72)
print("  1. NOISE FLOOR: 7Be recoil caps ambient monopole emission at eps<1e-3")
print("     per event -- a driven coherent array rises above this as (N eps)^2.")
# bound L_2nd < 0.1 L_sun (8B flux, Gondolo & Raffelt 2009), L_2nd = theta_ch^2 L_nu,
# L_nu/L_sun = 0.02425 (computed in stellar_cooling_secondsound.py)
margin_stellar = 0.1/(theta**2*0.02425)
print(f"  2. STELLAR CEILING: cleared by ~11 orders (margin"
      f" 0.1/(theta_ch^2 L_nu/L_sun) = {margin_stellar:.1e},")
print("     chiral-suppressed burning; see stellar_cooling_secondsound.py)")
print("     -- a lab source is nowhere near.")
print("  3. RECIPE: a charged plate drives NO second sound (pure shear, no")
print("     monopole content); the driver must be a coherent E0 (volume)")
print("     oscillator. That is the whole content of 'engineer the source'.")
print("="*72)
