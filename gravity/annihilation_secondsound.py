#!/usr/bin/env python3
r"""
annihilation_secondsound.py
=============================
Second sound from the two natural full-strength monopole sources.

The stellar-cooling gate (stellar_cooling_secondsound.py) established that
second sound couples to a process through its MONOPOLE CONTENT, the net
fractional change in defect (node) count, and that ordinary stellar burning
is chiral-suppressed (eps ~ theta_ch). The next question the framework can
now pose: the LOUDEST clean monopole source is annihilation, which removes
defects outright (delta-rho/rho -> -1). Two epochs run annihilation at full
per-event strength -- cosmological e+e- annihilation (t ~ 1 s, T ~ m_e) and
the supernova core. This script asks what second-sound energy each sheds and
confronts it with the tightest existing bound on each.

KEY PHYSICS (and the reason the answer is not catastrophic):
  A process radiates second sound through the RATE of monopole change, not
  its accumulated value. A single annihilation event has eps_event ~ 1, but
  what a bulk medium radiates coherently is set by how fast the net node
  count changes compared with the mode frequency it can drive. The emitted
  power per unit volume scales as
        P/V ~ K_sf * (d/dt [delta-rho/rho])^2 / (rho c v2)     (dimensional)
  i.e. it is quadratic in the monopole RATE. Slow (equilibrium) annihilation
  is quiet; fast (out-of-equilibrium, shock) annihilation is loud.

EPISODE 1 - cosmological e+e- annihilation.
  Happens in near-thermal equilibrium over a Hubble time, so the monopole
  rate is Hubble-slow: d/dt(delta-rho/rho) ~ H. The natural comparison is the
  radiation energy density, and the bound is BBN / N_eff, which measures any
  exotic energy sink at T ~ MeV to ~few percent (Delta N_eff < 0.2). We show
  the equilibrium rate puts the second-sound drain far below that.

EPISODE 2 - supernova core.
  Neutronization + annihilation in the collapsing core is fast and out of
  equilibrium: the monopole rate is set by the collapse time, not Hubble.
  The bound is the observed neutrino burst of SN 1987A, which carried
  essentially the full binding energy ~3e46 J. Any second-sound channel must
  not have stolen a large fraction of that. We estimate the second-sound
  luminosity and compare.

Both are scoping estimates (order-unity geometric prefactors), but each lands
firmly on the safe side, and each makes a falsifiable statement: a small,
computable second-sound contribution rides along with both episodes.
"""

import numpy as np

# ---- constants -----------------------------------------------------------
alpha  = 1/137.035999
theta  = alpha**2/(2*np.pi)
c      = 2.99792458e8
hbar   = 1.054571817e-34
kB     = 1.380649e-23
e      = 1.602176634e-19
G      = 6.67430e-11
ell    = 2.8179403262e-15
v2     = 3.65*c
m0     = 0.51099895069e6/alpha*e/c**2      # node mass, kg
rho_L  = 0.8*m0/ell**3                      # superfluid density of the lattice
K_sf   = c**4/(G*ell**2)

print("="*72)
print("EPISODE 1: cosmological e+e- annihilation (t ~ 1 s, T ~ m_e c^2)")
print("="*72)
# standard-cosmology anchors at T = m_e
T_MeV   = 0.511
T_K     = T_MeV*1e6*e/kB
t_epoch = 1.0                               # s, ~ age at T ~ m_e
H       = 1/(2*t_epoch)                      # radiation era H ~ 1/(2t)
# radiation energy density at this T (g* ~ 10.75)
gstar   = 10.75
a_rad   = np.pi**2/30 * (kB**4)/(hbar**3 * c**3)   # radiation constant (natural)
rho_rad = gstar * a_rad * T_K**4 / c**2      # kg/m^3 equivalent
print(f"  T = {T_MeV} MeV = {T_K:.2e} K,  age ~ {t_epoch} s,  H ~ {H:.2e} /s")
print(f"  radiation density rho_rad ~ {rho_rad:.2e} kg/m^3")
print()
# monopole rate: near-equilibrium annihilation tracks the expansion, so the
# fractional node-count change per unit time is ~ H (the slowest possible).
mono_rate = H                                # d/dt (delta-rho/rho) ~ H
# second-sound power density ~ K_sf (mono_rate/omega_ref)^2 ... but the honest
# comparison is the DIMENSIONLESS drain per Hubble time relative to radiation.
# Emitted 2nd-sound energy density over a Hubble time, as a fraction of rho_rad,
# scales as (v2/c) * (rho_L-coupling) * (mono_rate * t)^2 with the monopole
# coupling entering as the fraction of the plasma that is annihilating baryons/leptons.
# The leptons annihilating are a fraction ~1 of the plasma, but the RATE is Hubble-slow:
frac_drain = (mono_rate*t_epoch)**2 * (rho_rad/rho_L) * (c/v2)
print(f"  monopole rate ~ H (equilibrium): (H t)^2 = {(mono_rate*t_epoch)**2:.2f}")
print(f"  2nd-sound energy drain / radiation ~ {frac_drain:.2e}")
print(f"  BBN / N_eff bound on exotic sink at this T: ~few percent (dNeff<0.2)")
verdict1 = "PASS" if frac_drain < 0.03 else "FAIL"
print(f"  -> {verdict1}: equilibrium annihilation is Hubble-slow, drain << 1%")
print("  Prediction: a tiny, computable 2nd-sound trickle rides e+e-")
print("  annihilation; it perturbs N_eff at a level far below current reach.")
print()

print("="*72)
print("EPISODE 2: supernova core (fast, out-of-equilibrium annihilation)")
print("="*72)
# SN 1987A anchors
E_bind   = 3.0e46                            # J, gravitational binding released
t_coll   = 1.0e-2                            # s, collapse/neutronization timescale
R_core   = 1.0e4                             # m, ~10 km proto-neutron star
E_nu_obs = 3.0e46                            # J, neutrino burst carried ~all of it
print(f"  binding energy released ~ {E_bind:.1e} J over ~{t_coll} s")
print(f"  SN 1987A neutrino burst carried ~{E_nu_obs:.1e} J (essentially all)")
print()
# Here the monopole rate is fast: a fraction f_ann of core baryons change
# defect character on the collapse time. Second-sound luminosity scales as the
# neutrino luminosity times (monopole content)^2 of the rearrangement, exactly
# as in the stellar case, BUT the SN rearrangement is neutronization
# (p -> n), a beta-type flip -> monopole content ~ theta_ch again (NOT 1),
# because it is the SAME screw/edge flip, chiral-suppressed.
eps_SN   = theta                             # neutronization is a beta flip
L_ratio  = eps_SN**2
E_2nd    = L_ratio * E_bind
print(f"  neutronization p->n is a beta flip -> monopole content ~ theta_ch = {theta:.1e}")
print(f"  2nd-sound energy / binding ~ theta_ch^2 = {L_ratio:.1e}")
print(f"  E_2nd ~ {E_2nd:.1e} J  vs neutrino burst {E_nu_obs:.1e} J")
verdict2 = "PASS" if E_2nd < 0.1*E_nu_obs else "FAIL"
print(f"  -> {verdict2}: chiral suppression again keeps it {E_nu_obs/E_2nd:.0e}x below")
print("     the neutrino burst, safe against SN 1987A energetics.")
print()
print("  The one loud exception the framework allows is TRUE annihilation")
print("  (matter meeting antimatter, eps~1), which does not occur in a SN")
print("  core. Where it does occur -- the early-universe baryon/antibaryon")
print("  annihilation before freeze-out -- the medium is in equilibrium, so")
print("  the RATE is Hubble-slow (Episode 1). Full per-event strength and")
print("  fast bulk rate never coincide in nature, which is why the channel")
print("  is quiet without being forbidden.")
print()
print("="*72)
print("SUMMARY")
print("="*72)
print("  Both natural annihilation episodes shed second sound, and both are")
print("  safe: cosmological e+e- because equilibrium makes it Hubble-slow,")
print("  the SN core because neutronization is a chiral-suppressed beta flip.")
print("  The general rule: full monopole STRENGTH (eps~1) and fast bulk RATE")
print("  never coincide in nature. A built device beats nature precisely by")
print("  making them coincide -- coherent, fast, full-strength E0 breathing.")
print("="*72)
