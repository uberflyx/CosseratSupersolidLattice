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
defects and empties the cells they held. But the per-event strength depends
on WHAT annihilates. Baryonic annihilation removes nodes outright
(delta-rho/rho -> -1, full strength, GeV scale). Light-lepton annihilation
(e+e- -> gamma gamma) merges two screw dislocations into transverse shear
waves (photons); both the incoming defects and the outgoing photons are pure
shear (div u = 0 at linear order), so the monopole content is a small
higher-order leak, NOT of order one. This script treats three arenas: the
LABORATORY (per-event e+e- branching vs the positronium invisible-decay
searches), and two cosmological/astrophysical bulk episodes.

KEY PHYSICS (and the reason the answer is not catastrophic):
  A process radiates second sound through the RATE of monopole change, not
  its accumulated value. A single BARYONIC annihilation event has
  eps_event ~ 1; a single LIGHT-LEPTON event has eps_event << 1, capped by
  energy conservation (2 m_e c^2 cannot pay for a full-cell compression that
  costs ~m_0 c^2) and pinned by the positronium searches to ~1.6e-5 per cell
  (LABORATORY block below). For the bulk episodes, what a medium radiates
  coherently is set by how fast the NET node count changes compared with the
  mode frequency it can drive. The emitted power per unit volume scales as
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
N2     = 1/np.pi                              # Cosserat coupling number
C11_mu = (3 - 4*N2)/(1 - N2)                  # C11/mubar on the D4 shell = 2.533
v2     = np.sqrt(C11_mu)*c                    # second sound = sqrt(C11/rho_x) = 1.592 c
m0     = 0.51099895069e6/alpha*e/c**2      # node mass, kg
F_S    = 0.99932                           # the Leggett ceiling adopted as the working superfluid fraction
rho_L  = F_S*m0/ell**3                      # superfluid density of the lattice
K_sf   = c**4/(G*ell**2)

print("="*72)
print("LABORATORY: per-event e+e- second-sound branching vs positronium bounds")
print("="*72)
# Emission law E_2nd = (1/2) C11 (drho/rho)^2 V over one cell V = ell^3, with
# the longitudinal modulus C11 = 2.533 mubar. Using mubar ell^3 = m0 c^2 and
# m0 = m_e/alpha, the per-event branching into second sound (energy over the
# 2 m_e c^2 released) is
#       Br_2nd = E_2nd / (2 m_e c^2) = (C11/mubar) / (4 alpha) (drho/rho)^2.
coeff = C11_mu/(4.0*alpha)                    # Br_2nd = coeff * (drho/rho)^2
print(f"  Br_2nd = (C11/mubar)/(4 alpha) (drho/rho)^2 = {coeff:.1f} * (drho/rho)^2")
# (1) Energy-conservation ceiling: Br_2nd <= 1 forbids the full-strength value.
drho_ceiling = np.sqrt(1.0/coeff)
print(f"  energy-conservation ceiling: Br<=1 => drho/rho <= {drho_ceiling:.3f}")
print(f"    -> a single e+e- pair CANNOT reach drho/rho ~ 1 (full strength);")
print(f"       full-cell compression costs (1/2) C11 ell^3 = {C11_mu/2:.2f} m0 c^2 = "
      f"{C11_mu/2*0.51099895069/alpha:.1f} MeV >> 2 m_e c^2 = 1.022 MeV.")
# (2) Badertscher 2007 (Phys. Rev. D 75, 032004) invisible-decay limits.
Br_oPs    = 4.2e-7                             # Br(o-Ps -> invisible), 90% CL
Br_direct = 2.1e-8                             # Br(e+e- -> invisible), 90% CL
drho_oPs    = np.sqrt(Br_oPs/coeff)
drho_direct = np.sqrt(Br_direct/coeff)
print(f"  Badertscher o-Ps    bound Br<{Br_oPs:.1e} => drho/rho < {drho_oPs:.2e}")
print(f"  Badertscher direct  bound Br<{Br_direct:.1e} => drho/rho < {drho_direct:.2e}")
print(f"    -> tightest lab limit on the annihilation monopole content:")
print(f"       drho/rho <~ {drho_direct:.1e} per cell (cf. 7Be beta bound ~1e-3).")
# (3) If the leak is the chiral amplitude theta_ch, predict the branching.
Br_chiral = coeff*theta**2
print(f"  IF the leak is chiral (drho/rho ~ theta_ch = {theta:.2e}):")
print(f"       Br_2nd = {Br_chiral:.2e}, a factor {Br_direct/Br_chiral:.1f} below")
print(f"       the current direct bound -- a next-generation target (Vigo 2018).")
print("  Note: the annihilation monopole content is not yet derived from first")
print("  principles; theta_ch is a plausible identification, flagged not asserted.")
print()

# standard-cosmology anchors at k_B T = m_e c^2
T_MeV   = 0.511
T_K     = T_MeV*1e6*e/kB
# radiation energy density at this T, g* = 10.75 (photons, e+e- pairs, three
# neutrino species): u = g* (pi^2/30) (k_B T)^4/(hbar c)^3
gstar   = 10.75
a_rad   = np.pi**2/30 * (kB**4)/(hbar**3 * c**3)   # pi^2/30 k_B^4/(hbar c)^3, per dof
rho_rad = gstar * a_rad * T_K**4 / c**2      # kg/m^3 equivalent
# Friedmann, radiation era: H = sqrt(8 pi G rho_rad/3), a ~ t^(1/2) so
# H = 1/(2t) and H t = 1/2 at every epoch. Equivalently
# t = 2.42/sqrt(g*) (T/MeV)^-2 s.
H       = np.sqrt(8*np.pi*G*rho_rad/3)
t_epoch = 1/(2*H)                            # s, age at k_B T = m_e c^2
t_std   = 2.42/np.sqrt(gstar) * T_MeV**-2    # s, the standard closed form
print("="*72)
print(f"EPISODE 1: cosmological e+e- annihilation (k_B T = m_e c^2,"
      f" t = {t_epoch:.1f} s)")
print("="*72)
print(f"  T = {T_MeV} MeV = {T_K:.2e} K,  g* = {gstar}")
print(f"  age t = 1/(2H) = {t_epoch:.2f} s  (2.42/sqrt(g*) (T/MeV)^-2 s = {t_std:.2f} s),"
      f"  H = {H:.2f} /s,  H t = {H*t_epoch:.1f}")
print(f"  radiation density rho_rad = {rho_rad:.2e} kg/m^3"
      f"  vs rho_s = {rho_L:.2e} kg/m^3")
print()
# monopole rate: near-equilibrium annihilation tracks the expansion, so the
# fractional node-count change per unit time is ~ H (the slowest possible).
mono_rate = H                                # d/dt (delta-rho/rho) ~ H
# second-sound power density ~ K_sf (mono_rate/omega_ref)^2 ... but the honest
# comparison is the DIMENSIONLESS drain per Hubble time relative to radiation.
# Emitted 2nd-sound energy density over a Hubble time, as a fraction of rho_rad,
# scales as (H t)^2 (rho_rad/rho_s) (c/v2), every order-one prefactor set to
# one, with rho_s = f_s m0/ell^3 the condensate density (rho_L here).
# The leptons annihilating are a fraction ~1 of the plasma, but the RATE is Hubble-slow:
frac_drain = (mono_rate*t_epoch)**2 * (rho_rad/rho_L) * (c/v2)
print(f"  monopole rate ~ H (equilibrium): (H t)^2 = {(mono_rate*t_epoch)**2:.2f}")
print(f"  (H t)^2 (rho_rad/rho_s) (c/v2) = {(mono_rate*t_epoch)**2:.2f} x"
      f" {rho_rad/rho_L:.1e} x {c/v2:.2f}")
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
