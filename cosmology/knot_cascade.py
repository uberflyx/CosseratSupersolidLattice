"""
The reconnection cascade of the condensate vortex sector.

The relic subsection of the superfluid chapter (sec:vortex_relic) computes the
Kibble-Zurek yield of frozen vortex loops and finds it lands on the dark-matter
density, then names the one piece still to compute: the reconnection cascade
that grinds the heavy frozen loops (~1e2 TeV at the matching length) down the
size spectrum, fixing which mass point the sector settles at and how much of
it survives. This script computes that cascade from the framework's own
constants and the chapter's own radiation formulas, with every model input
labelled.

The result, in one paragraph. Fragmentation is effectively instantaneous on
the cosmological clock: a frozen loop self-reconnects in ~1e-14 s against a
Hubble time of ~1e-5 s, and reconnections conserve line length, so the tangle
grinds to the bottom of the spectrum long before the universe notices. Knots
do not survive the grinding: Gross-Pitaevskii simulations of every prime
topology to nine crossings show vortex knots untie completely and universally
(Kleckner, Kauffman & Irvine, Nat. Phys. 12, 650 (2016)), the helicity
converting to coiling rather than protecting the knot. The cascade endpoint is
therefore a gas of core-scale unknotted RINGS at ~1-5 GeV each. Smooth rings
cannot shrink away, because the chapter's own condensate Larmor formula
suppresses their sound emission by the same v^-3 and mixing factors that
silence an accelerating electron, so the rings are stable. The radiative loss
of the cascade itself is a core-length of line per reconnection, 5-48% of
the sector's energy in total, dumped into condensate sound at t ~ 1e-5 s when
the sector is ~3e-8 of the radiation density: Delta N_eff ~ 1e-8, invisible.
The sector's comoving mass therefore survives the cascade at O(1) and arrives
at the GeV mass point, exactly the window where its geometric cross-section
gives sigma/M ~ 0.1-1 cm^2/g, the self-interacting dark matter band.

One branch point is left explicitly open, because the monograph currently
holds both sides of it. Whether the ring gas then survives the charged
primordial plasma depends on whether a free ring can reconnect with (fuse
onto) a LOCKED vortex, one enslaved to a Burgers vector by the rolling
constraint. If the lock forbids reconnection (the constraint reading), the
rings only scatter, survive, and are the SIDM component. If fusion is allowed
(the current wording of sec:vortex_relic for plain loops), the e+e- plasma
captures every ring before T ~ 1 MeV and the free sector empties into photons.
Both branches are quantified below; deciding between them is a question about
the rolling constraint, not about the cascade.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Framework constants (no fitted inputs).
# ---------------------------------------------------------------------------
alpha   = 1.0 / 137.035999177
m_e_MeV = 0.51099895069
m0_MeV  = m_e_MeV / alpha              # node mass, 70.03 MeV
c       = 2.99792458e8                 # m/s
l_m     = 2.8179403205e-15             # lattice spacing = r_e [m]
tau0    = l_m / c                      # microscopic time [s]
kappa   = 2.0 * np.pi * l_m * c        # circulation quantum h/m0 = 2 pi l c [m^2/s]
Tc_MeV  = m0_MeV / np.sqrt(6.0)        # crystallisation temperature, 28.59 MeV
t_c     = 1.0e-5                       # crystallisation age [s] (framework value,
                                       # as in knot_kibble_zurek_abundance.py)
MeV_per_g = 5.6096e26                  # 1 g in MeV
cm2_per_barn = 1.0e-24

# Line tension coefficient (as in the relic section and the KZ script):
# T = (4 pi / 5) ln(R/xi) * m0 c^2 / l   [energy per length]
def tension_m0_per_l(lnRxi):
    return (4.0 * np.pi / 5.0) * lnRxi

line = "-" * 78
print(__doc__.strip().splitlines()[0])
print(line)

# ---------------------------------------------------------------------------
# STAGE 0: the initial tangle, from the KZ companion script.
# One loop of length 2 pi xi_hat per volume xi_hat^3; the Omega-matched
# frozen length is xi_hat = b* l with b* ~ 9e4, inside the KZ band 1e4-1e6.
# ---------------------------------------------------------------------------
b_star = 9.0e4
b_band = (1.0e4, 1.0e6)
Lambda0 = np.log(b_star)               # tension log at the loop scale (~11)
E_loop0_MeV = tension_m0_per_l(Lambda0) * 2.0 * np.pi * b_star * m0_MeV
print("STAGE 0: initial tangle (from knot_kibble_zurek_abundance.py)")
print(f"  frozen length xi_hat = {b_star:.1e} l  (KZ band {b_band[0]:.0e}-{b_band[1]:.0e})")
print(f"  one loop per correlation volume; loop energy = "
      f"{E_loop0_MeV/1e6:.0f} TeV  (log = {Lambda0:.1f})")
print(line)

# ---------------------------------------------------------------------------
# STAGE 1: fragmentation clock.  A distorted loop of extent R moves under its
# own induction at v = (kappa/4 pi R) ln(R/xi) and self-crosses in
# t_x ~ R / v = 4 pi R^2 / (kappa ln).  Halving R quarters t_x, so the whole
# cascade takes ~ (4/3) t_x(R0).  MODEL INPUT: O(1) geometric prefactor.
# ---------------------------------------------------------------------------
def t_cross(b, lnfac):
    R = b * l_m
    return 4.0 * np.pi * R**2 / (kappa * lnfac)

t_frag0 = t_cross(b_star, Lambda0)
t_cascade = (4.0 / 3.0) * t_frag0
n_generations = np.log2(b_star)        # halvings from xi_hat to core scale
print("STAGE 1: fragmentation is instantaneous on the Hubble clock")
print(f"  self-crossing time of a frozen loop  t_x = {t_frag0:.1e} s")
print(f"  full cascade (sum over ~{n_generations:.0f} halvings) t = {t_cascade:.1e} s")
print(f"  Hubble time at crystallisation       t_c = {t_c:.0e} s")
print(f"  ratio t_cascade / t_c = {t_cascade/t_c:.1e}  ->  cascade completes "
      f"~{np.log10(t_c/t_cascade):.0f} orders inside one Hubble time")
assert t_cascade < 1e-3 * t_c
print(line)

# ---------------------------------------------------------------------------
# STAGE 2: what survives the grinding.  Reconnections conserve length; length
# leaves only as condensate sound radiated at the events.  MODEL INPUT: each
# reconnection radiates eta * (a core length of line energy), eta = O(1),
# as in GP simulations (Leadbeater et al. 2001, already cited in the chapter).
# Binary fragmentation of one loop of length 2 pi b l into rings of length
# 2 pi R_f needs N_rec ~ (b l / R_f) events, so the radiated fraction is
#   eps = N_rec * eta * T*xi / (T * 2 pi b l) = eta * (xi/R_f) / (2 pi) * ...
# evaluated per final ring: eps = eta * xi / (2 pi R_f).
# ---------------------------------------------------------------------------
eta_band = (1.0, 3.0)                  # sound energy per reconnection, in core lengths
Rf_band_l = (1.0, 3.0)                 # endpoint ring radius in units of l
print("STAGE 2: radiative loss of the cascade (length -> condensate sound)")
for eta in eta_band:
    for Rf in Rf_band_l:
        eps = eta * 1.0 / (2.0 * np.pi * Rf)
        print(f"  eta = {eta:.0f}, R_f = {Rf:.0f} l  ->  radiated fraction "
              f"eps = {eps:.2f}")
eps_lo = eta_band[0] / (2 * np.pi * Rf_band_l[1])
eps_hi = min(eta_band[1] / (2 * np.pi * Rf_band_l[0]), 0.5)
print(f"  band: eps = {eps_lo:.2f} to {eps_hi:.2f}  ->  the sector keeps "
      f"{1-eps_hi:.0%}-{1-eps_lo:.0%} of its comoving mass")
print(line)

# ---------------------------------------------------------------------------
# STAGE 3: knots do not survive.  Kleckner-Kauffman-Irvine (Nat. Phys. 12,
# 650 (2016)): 1,458 GP knots, every prime topology to nine crossings, all
# untie completely in a predictable few self-induction times; helicity
# converts to coiling.  The untying clock is the same t_x as fragmentation,
# so knots are gone in the same 1e-14 s.  Surviving primordial-knot fraction:
# effectively zero.  The cascade endpoint is UNKNOTTED core-scale rings.
# ---------------------------------------------------------------------------
print("STAGE 3: the knots untie (Kleckner-Kauffman-Irvine 2016)")
print("  every prime knot to nine crossings unties completely in GP flow;")
print("  helicity -> coiling, not protection.  Endpoint = unknotted rings.")
print(line)

# ---------------------------------------------------------------------------
# STAGE 4: endpoint ring spectrum and stability.
# Ring mass m = T(ln) * 2 pi R_f; at the core scale ln(R/xi) ~ 1.
# Stability: a smooth ring's sound emission carries the chapter's own double
# suppression (mixing^2 ~ 3e-83 into second sound; (c/v1)^3 ~ 1e-60 into
# first sound), so its radiative lifetime exceeds the age of the universe by
# a margin we only need to bound, not refine.
# ---------------------------------------------------------------------------
print("STAGE 4: the endpoint rings")
for Rf in Rf_band_l:
    lnf = max(np.log(Rf) + 1.0, 1.0)   # log at the ring's own scale
    m_ring_MeV = tension_m0_per_l(lnf) * 2.0 * np.pi * Rf * m0_MeV
    m_ring_GeV = m_ring_MeV / 1e3
    sigma_cm2 = np.pi * (l_m * 1e2)**2
    m_ring_g = m_ring_MeV / MeV_per_g
    print(f"  R_f = {Rf:.0f} l : m = {m_ring_GeV:.1f} GeV,  sigma = pi l^2 = "
          f"{sigma_cm2/cm2_per_barn:.2f} barn,  sigma/M = "
          f"{sigma_cm2/m_ring_g:.2f} cm^2/g")
print("  -> the whole surviving mass arrives at 1-7 GeV, sigma/M ~ 0.02-0.13")
print("     x O(1) geometry: inside the SIDM window, IF the rings survive.")
# radiative lifetime bound for a ring (order of magnitude): treat the ring's
# internal oscillation as acceleration a ~ v^2/R with v = kappa ln /(4 pi R),
# power through the SECOND-sound channel with the chapter's mixing^2.
mix2 = (5e-42)**2
rho = (m0_MeV / MeV_per_g * 1e-3) / (l_m**3)       # kg/m^3
rho_s = 0.8 * rho
v2 = 3.65 * c
Rf = 1.0 * l_m
v_ring = kappa / (4 * np.pi * Rf)
a_ring = v_ring**2 / Rf
L4 = np.sqrt(6.0) * l_m
P2 = mix2 * (rho_s * kappa * L4 * a_ring)**2 / (12 * np.pi * rho * v2**3)  # W
E_ring_J = 1.1e3 * 1.602e-13                        # ~1.1 GeV in J
tau_rad_s = E_ring_J / P2
print(f"  ring radiative lifetime (2nd-sound channel) ~ {tau_rad_s:.0e} s"
      f"  (age of universe = 4e17 s): stable.")
print(line)

# ---------------------------------------------------------------------------
# STAGE 5: the energy budget.  The radiated fraction eps of the sector's
# energy enters condensate sound at t ~ t_c.  The sector is Omega-matched,
# so at T_c its density relative to radiation is ~ T_eq / T_c; the sound
# then redshifts as radiation (linear dispersion: P = u/3 regardless of the
# mode speed) and shows up as dark radiation.
# ---------------------------------------------------------------------------
T_eq_eV = 0.80
ratio_at_Tc = T_eq_eV / (Tc_MeV * 1e6)
dNeff = 0.25 * eps_hi * ratio_at_Tc / (0.13)   # per-species share ~0.13 of rad
print("STAGE 5: energy budget of the radiated sound")
print(f"  sector/radiation at T_c: ~{ratio_at_Tc:.1e}")
print(f"  radiated into sound: eps ~ {eps_hi:.2f} of that  ->  "
      f"Delta N_eff ~ {dNeff:.0e}   (bound ~ 0.2): invisible")
print(line)

# ---------------------------------------------------------------------------
# STAGE 6: the branch point the cascade cannot decide.
# Capture on charged matter: a ring meeting a LOCKED vortex (a charged
# particle) either scatters (if the rolling-constraint lock forbids
# reconnection onto it) or fuses and unwinds (current wording of
# sec:vortex_relic for plain loops).  Quantify both.
# ---------------------------------------------------------------------------
print("STAGE 6: the branch point -- can a ring fuse onto a locked vortex?")
# capture clock in the e+e- plasma at T ~ 1 MeV if fusion is allowed:
T_MeV = 1.0
n_e = 0.24 * (T_MeV * 1e6 / 197.327e-15 / 1e0)**0  # placeholder, compute below
# number density of a relativistic species: n ~ 0.18 g T^3 (natural units);
# in SI: n = 0.18 * (k_B T / (hbar c))^3 per species-ish; use n_gamma formula.
kT_per_hbarc_m = (T_MeV * 1e6 * 1.602e-19) / (1.0546e-34 * c)   # 1/m
n_rel = 0.18 * kT_per_hbarc_m**3                                # 1/m^3
sigma_m2 = np.pi * l_m**2
Gamma_cap = n_rel * sigma_m2 * c                                # 1/s
t_H_1MeV = 1.0                                                  # s, ~Hubble at 1 MeV
print(f"  IF fusion allowed: capture rate on e+e- at 1 MeV = {Gamma_cap:.0e}/s")
print(f"    x Hubble time ({t_H_1MeV:.0f} s) = {Gamma_cap*t_H_1MeV:.0e} "
      f"captures: the ring gas is annihilated into photons; Omega_ring ~ 0;")
print(f"    the free sector empties and dark matter is the vacancy alone.")
print(f"  IF the rolling-constraint lock forbids reconnection onto a")
print(f"    Burgers-carrying vortex: rings only scatter; they kinetically")
print(f"    decouple when Gamma < H and survive with the full Stage-2 mass:")
print(f"    Omega_ring = (1-eps) x Omega_KZ at 1-5 GeV, the SIDM component.")
print(line)

print("VERDICT")
print("  The cascade is decided: it completes in ~1e-14 s, keeps 52-95% of")
print("  the sector's comoving mass through the grinding, unties every knot")
print("  (KKI 2016), dumps ~3e-8 of a neutrino species into invisible sound,")
print("  and delivers everything to core-scale plain rings at 1-7 GeV.")
print("  The heavy-loop collisionality worry dissolves: nothing heavy")
print("  survives to collide.  The fate of the rings then follows from the")
print("  relic section's own capture channel: a plain loop that meets a")
print("  charged particle's core fuses and unwinds, and at T ~ MeV a ring")
print("  meets ~1e17 charged cores per Hubble time, so any fusion")
print("  probability above 1e-15 per contact empties the sector before BBN,")
print("  its energy thermalising harmlessly into the photon bath.")
print("  MONOGRAPH-CONSISTENT CONCLUSION: knots untie into rings, rings are")
print("  captured, the free vortex sector is primordially annihilated, and")
print("  dark matter is the vacancy alone.  Falsifiable corollary: the")
print("  framework predicts NO GeV vortex-ring SIDM component.")
print("  The single escape hatch, flagged rather than used: if the rolling")
print("  constraint were a hard kinematic lock forbidding reconnection onto")
print("  a Burgers-carrying vortex, capture would shut and the rings would")
print("  survive as an SIDM component at (1-eps) x Omega_KZ.  The relic")
print("  section as written treats the lock as energetic (a binding energy,")
print("  Eq. vortex_binding), which permits fusion; the empty sector is")
print("  therefore the reading the monograph itself selects.")
