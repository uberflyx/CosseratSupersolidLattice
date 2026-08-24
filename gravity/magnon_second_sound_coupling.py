#!/usr/bin/env python3
r"""
magnon_second_sound_coupling.py
================================
The transducer number: how strongly does a precessing electron spin, which
in this framework is a moving vacuum vortex, radiate the condensate's
compression wave (second sound, c_s = sqrt(3) c)?

The clock-engineering appendix leaves the fluid-crystal cross-coupling as
its single named open quantity. This script pins it, with an explicit and
wide error bracket, and shows it appears in the laboratory as an anomalous
contribution to the Gilbert damping of a magnon condensate (YIG). Existing
YIG linewidth data therefore already BOUNDS it.

Two independent estimates bracket the coupling because they make opposite
assumptions about how much real condensate flow a spin precession modulates:

  BRACKET A (magnetic-dipole analog, physical lower estimate). The precessing
    moment radiates EM (the shear wave) at the textbook magnetic-dipole
    Larmor rate P_EM = mu0 w^4 m^2 sin^2/(6 pi c^3). The SAME transverse
    source drives compression only weakly, so second-sound emission is
    bounded above by the shear rate times the compression phase-space
    factor (c/c_s)^3. This uses the ACTUAL magnetic moment, so it is the
    conservative floor.

  BRACKET B (vortex-Larmor analog, upper estimate). Treat the electron as a
    vortex of circulation kappa = h/m0 whose core swings by a core radius
    (d0 ~ ell) as the spin precesses, and apply the sound-Larmor formula
    P/L ~ rho_s kappa^2 a^2/(6 pi c_s^3), a = w^2 d0. The circulation carries
    a large flow, so this over-estimates unless the whole core really moves.

The truth sits between them, and closing the ~12-order gap requires deriving
the magnon's actual vortex structure (how much condensate flow one quantum
of precession modulates). That derivation is the open work; this script
frames it and proves the channel is not already excluded.

Result: both brackets put the second-sound damping far below the measured
YIG Gilbert damping (~1e-5), so the framework is NOT falsified, and YIG is
not yet a sensitive test. The propulsion leverage lives entirely in the
coherent (superradiant, x N) emission of a phased condensate, which is why a
magnon CONDENSATE, not incoherent spins, is the right substrate.
"""

import numpy as np

# ---- framework constants -------------------------------------------------
alpha = 1/137.035999
c     = 2.99792458e8
hbar  = 1.054571817e-34
h     = 2*np.pi*hbar
e     = 1.602176634e-19
ell   = 2.8179403262e-15
mu0   = 4e-7*np.pi
m0    = (0.51099895069e6*e/alpha)/c**2      # node mass [kg]
rho_s = 0.9999*m0/ell**3                    # superfluid density [kg/m^3]
kappa = h/m0                                # circulation quantum = 2 pi c ell [m^2/s]
cs    = np.sqrt(3.0)*c                      # second sound speed [m/s]

# ---- magnon condensate operating point (Demokritov et al., YIG BEC) ------
f_mag = 2.9e9                               # magnon BEC frequency [Hz]
w     = 2*np.pi*f_mag                       # [rad/s]
gamma = 1.76085963e11                       # electron gyromagnetic ratio [rad/s/T]
S, g  = 2.5, 2.0                            # Fe3+ spin and g-factor
m_mom = g*9.2740100783e-24*S               # Fe3+ magnetic moment [J/T]
n_Fe  = 1.05e28                             # Fe spins per m^3 in YIG

def gilbert_from_power(P, sin2=1.0):
    """Anomalous Gilbert damping equivalent to a per-spin radiated power P."""
    return P*gamma/(w**2*m_mom*sin2)

print("="*72)
print("Framework primitives and operating point")
print("="*72)
print(f"  kappa = h/m0 = {kappa:.3e} m^2/s   c_s = sqrt(3) c = {cs:.3e} m/s")
print(f"  rho_s = {rho_s:.3e} kg/m^3   m(Fe3+) = {m_mom:.3e} J/T")
print(f"  magnon BEC: f = {f_mag/1e9:.1f} GHz, omega = {w:.3e} rad/s")
print()

print("="*72)
print("BRACKET A - magnetic-dipole coupling (physical lower estimate)")
print("="*72)
P_EM   = mu0*w**4*m_mom**2/(6*np.pi*c**3)
P_2s_A = P_EM*(c/cs)**3
print(f"  P_EM (shear/EM, exact dipole Larmor) = {P_EM:.3e} W/spin")
print(f"    -> radiative Gilbert alpha_EM      = {gilbert_from_power(P_EM):.2e}")
print(f"  P_2s <= P_EM (c/c_s)^3               = {P_2s_A:.3e} W/spin")
print(f"    -> alpha_2s <=                     = {gilbert_from_power(P_2s_A):.2e}")
print()

print("="*72)
print("BRACKET B - vortex-Larmor, core excursion d0 = ell (upper estimate)")
print("="*72)
d0   = ell
a    = w**2*d0
P_1B = rho_s*kappa**2*a**2*ell/(6*np.pi*cs**3)   # P/L * (line length ~ell)
print(f"  acceleration a = w^2 ell = {a:.3e} m/s^2")
print(f"  P_2s per spin  = rho_s kappa^2 a^2 ell/(6 pi c_s^3) = {P_1B:.3e} W/spin")
print(f"    -> alpha_2s  = {gilbert_from_power(P_1B):.2e}")
print()

print("="*72)
print("Comparison with measured YIG, and the coherence lever")
print("="*72)
for aM,lab in [(1e-5,"best bulk single-crystal YIG"),(1e-4,"good YIG film")]:
    print(f"  measured Gilbert alpha ~ {aM:.0e}   ({lab})")
print(f"  bracket:  {gilbert_from_power(P_2s_A):.0e}  <=  alpha_2s  <=  {gilbert_from_power(P_1B):.0e}")
print(f"  -> below measured by 7 to 19 orders: NOT falsified; not yet a sensitive test.")
print()
lam  = 2*np.pi*cs/w
Vcoh = lam**3
N    = n_Fe*Vcoh
print(f"  second-sound wavelength lam = {lam:.3f} m ; coherence volume ~ {Vcoh:.2e} m^3")
print(f"  N spins phased within lam^3 = {N:.2e}  ->  superradiant emission ~ N x single")
print(f"  incoherent efficiency floor eps = alpha_2s/alpha_meas:")
print(f"     bracket A {gilbert_from_power(P_2s_A)/1e-5:.0e} ... bracket B {gilbert_from_power(P_1B)/1e-5:.0e}")
print(f"  coherent (x N, drive/reaction limited) lifts eps toward O(1).")
print()

print("="*72)
print("FALSIFIABLE SIGNATURE")
print("="*72)
print("  The second-sound channel adds anomalous Gilbert damping with a")
print("  distinctive fingerprint the standard spin-orbit + phonon channels")
print("  do not share:")
print("    - scales as omega^4 (radiative), not omega (intrinsic),")
print("    - scales as 1/c_s^3 with the compression speed sqrt(3) c,")
print("    - GROWS with coherent sample volume (superradiant), unlike a")
print("      local intrinsic damping that is size-independent.")
print("  A magnon-BEC linewidth measured vs drive coherence and sample size,")
print("  isolating the omega^4 super-radiant residual, tests the coupling")
print("  directly and would either detect it or push the bound below 1e-12.")
print("="*72)
