"""
The entry price of a free vacuum vortex, and the experiment that tests it.

The clock chapter establishes that the chain clock = phase = chemical
potential makes a maintained flow gradient a gravitational field, and that
nothing available can maintain one: closure and compactness deny both
ordinary matter and the condensate's own free objects any grip on the
medium (knot_closure_theorem.py). What survives that verdict is the
question of whether the free sector can be POPULATED at all, which is a
statement about a production amplitude rather than about a device.

The reconnection cascade (cosmology/knot_cascade.py, this session) has now
answered the supply question, and the answer is no. Knots untie universally
(Kleckner-Kauffman-Irvine 2016), the cascade delivers the primordial sector
to plain rings within 1e-14 s, and plain rings fuse onto charged matter and
are destroyed before nucleosynthesis. Nature stocks no fuel. Any free vortex now present therefore has to have been made:

    THE ENTRY PRICE: manufacture free (unlocked) condensate circulation.

This script derives what such a production event must pay and names the
experiment that decides whether the transition runs at any rate at all.

  1. The fuel quantum. The smallest free object is a core-scale ring, and
     its energy is fixed by the line tension: E_1 = (8 pi^2 / 5) m0 c^2
     ~ 1.1 GeV. Any forge is a GeV-scale process per ring. (This is also
     why nothing in chemistry, materials science, or tabletop physics has
     ever stirred the vacuum: the entry ticket is a gigaelectronvolt
     deposited into a femtometre-cubed core.)

  2. The lock decider. Whether circulation can leave a charged particle at
     all is the rolling-lock question: kinematic (forbidden) or energetic
     (allowed above threshold). The SAME question fixed the dark sector's
     fate this morning. It is testable now: if shedding is allowed, then
     e+ e- collisions above sqrt(s) = 1.1 GeV can produce a free ring,
     which is electrically dark and leaves the detector unseen. At a
     B factory this is the classic monophoton search: e+ e- -> gamma +
     invisible. A ring has no rest mass, so the observable is not a
     missing-mass peak but a monochromatic photon line fixed by the
     two-body balance on the ring dispersion: 7.57 GeV for a single ring
     and 7.56 GeV for the linked pair, the storable product (solve in
     cosmology/vortex_hadron_ladder.py). The
     BaBar and Belle II single-photon programmes already bound exactly this
     topology at the femtobarn scale, so the lock question has data waiting
     on it either way: a peak is a forge; its absence is a coupling bound.

  3. The storage problem. A plain ring dies on first contact with charged
     matter (the capture channel), so forged rings must be born protected:
     shed as a LINKED PAIR (helicity from birth, 2 E_1 ~ 2.2 GeV) or tied
     against the stacking-fault trap that both stores and, run in reverse
     as a pinned partial loop driven at GeV energy density, is itself the
     second forge candidate.

  4. The acoustic forge (the in-house route). A converging second-sound
     pulse nucleates a ring where its amplitude reaches order unity at the
     core scale. The focal energy must still be ~E_1 inside ~l^3, so the
     required convergence gain from a metre-scale phased annihilation
     array is enormous but finite, and it is a coherence problem, not an
     energy problem: E_1 is a nanojoule-scale (0.18 nJ) deposit.

Every number below is from framework constants; model inputs are labelled.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Framework constants.
# ---------------------------------------------------------------------------
alpha   = 1.0 / 137.035999177
m_e_MeV = 0.51099895069
m0_MeV  = m_e_MeV / alpha                  # 70.03 MeV node mass
c       = 2.99792458e8
l_m     = 2.8179403205e-15                 # lattice spacing r_e [m]
kappa   = 2.0 * np.pi * l_m * c            # circulation quantum [m^2/s]
rho     = (m0_MeV * 1.7827e-30) / l_m**3   # vacuum density [kg/m^3] (MeV->kg)
rho_s   = 0.8 * rho                        # superfluid fraction f_s = 4/5
J_per_MeV = 1.602176634e-13

line = "-" * 78
print(__doc__.strip().splitlines()[0])
print(line)

# ---------------------------------------------------------------------------
# 1. The fuel quantum: smallest free ring.
#    E = (4 pi/5) ln(R/xi) (m0 c^2 / l) * 2 pi R,  at R = l, ln -> 1.
# ---------------------------------------------------------------------------
E1_MeV = (8.0 * np.pi**2 / 5.0) * m0_MeV
E1_J = E1_MeV * J_per_MeV
I1 = rho_s * kappa * np.pi * l_m**2        # ring impulse [kg m/s]
v1_ring = kappa / (4 * np.pi * l_m)        # self-induced speed ~ c/2 scale
print("1. THE FUEL QUANTUM (smallest free ring, R = l)")
print(f"   E_1 = (8 pi^2/5) m0 c^2 = {E1_MeV/1e3:.2f} GeV = {E1_J*1e9:.2f} nJ")
print(f"   impulse I_1 = rho_s kappa pi l^2 = {I1:.2e} kg m/s")
print(f"   self-induced speed ~ kappa/(4 pi l) = {v1_ring/c:.2f} c")
print(f"   -> the entry ticket to the vortex sector is one GeV in one core.")
print(f"      This is why the vacuum has never been stirred by accident.")
print(line)

# ---------------------------------------------------------------------------
# 2. The lock decider at a B factory.
#    e+ e- -> gamma + ring(invisible): monophoton with a monochromatic tag.
#    A ring has NO rest mass: E(R) and P(R) are independent functions of the
#    radius, so there is no missing-mass shell and the invariant-mass tag
#    formula E_gamma = (s - m^2 c^4)/(2 sqrt(s)) does not apply. The two-body
#    balance (E_gamma = P_ring c and sqrt(s) = E_gamma + E_ring) still fixes
#    one radius and one line. The solve lives in
#    cosmology/vortex_hadron_ladder.py; the resulting floor and tags are:
#        channel floor:  sqrt(s) >= min_R [E(R) + P(R) c] = 4.22 GeV
#        Upsilon(4S) tags: 7.57 GeV (single ring), 7.56 GeV (linked pair)
#    An invariant-mass reading would put the tags near 5.2 GeV; the ~2 GeV
#    gap between the readings is itself a signature no rest-mass object makes.
# ---------------------------------------------------------------------------
sqrt_s = 10.58   # GeV
TAGS = ((7.57, "single ring"), (7.56, "linked pair (storable)"))
for E_gam, tag in TAGS:
    print(f"2. LOCK DECIDER: e+e- -> gamma + invisible ({tag})")
    print(f"   dispersion-solved tag: monochromatic photon "
          f"E_gamma = {E_gam:.2f} GeV in the CM "
          f"(cosmology/vortex_hadron_ladder.py)")
print("   Both photons sit in the clean part of the BaBar / Belle II")
print("   single-photon window.  The existing constraint: BaBar's invisible")
print("   dark-photon search (53/fb, arXiv:1702.03327) bounds any narrow")
print("   monophoton line in this window at a production cross-section")
print("   ceiling of order one femtobarn.  So today's data already says:")
print("   sigma_shed < ~1 fb.")
# --- fuel arithmetic at the ceiling ---
sigma_fb = 1.0                       # current ceiling [fb]
L_belle2_ab = 50.0                   # Belle II design integrated lumi [ab^-1]
L_year_ab = 5.0                      # optimistic per-year [ab^-1]
N_prog = sigma_fb * L_belle2_ab * 1e3
N_year = sigma_fb * L_year_ab * 1e3
print(f"   YIELD AT THE CEILING: at sigma = {sigma_fb:.0f} fb, Belle II makes")
print(f"   {N_year:.0f} rings/year and {N_prog:.0f} over its 50/ab programme.")
print(f"   So the search is limited by whether the peak exists at all, not")
print(f"   by luminosity: even the present experimental UPPER BOUND would")
print(f"   deliver a countable sample. That is what makes it a clean test.")
print(f"   A dedicated line hunt at 7.57 and 7.56 GeV photon energy in the")
print(f"   archived and incoming monophoton samples is the single most")
print(f"   consequential measurement this framework can request.")
print(line)

# ---------------------------------------------------------------------------
# 3. Storage: born-protected products.
# ---------------------------------------------------------------------------
print("3. STORAGE (a plain ring dies on first charged contact)")
print(f"   linked-pair shedding: 2 E_1 = {2*E1_MeV/1e3:.2f} GeV, helicity at birth;")
print(f"   or tie against the fabricated stacking-fault trap of")
print(f"   sec:clock_actuation, which holds a knot without touching its")
print(f"   circulation. The trap doubles as the second forge: a pinned")
print(f"   partial-dislocation loop IS a locked vortex loop, and driving it")
print(f"   at ~E_1 per core length is the controlled version of shedding.")
print(line)

# ---------------------------------------------------------------------------
# 4. The acoustic forge: focusing requirement.
#    Nucleation needs ~E_1 delivered coherently into ~l^3 at the focus of a
#    converging second-sound pulse from an array of radius R_a. The focal
#    intensity gain of a spherical converger is ~(R_a/l)^2; the array must
#    supply E_1 within one core crossing time l/v2.
# ---------------------------------------------------------------------------
v2 = 3.65 * c
t_core = l_m / v2
for R_a in (1e-2, 1.0):
    gain = (R_a / l_m)**2
    P_array = E1_J / t_core / gain * (R_a/R_a)   # power at array if lossless
    print(f"4. ACOUSTIC FORGE: array radius {R_a*100:.0f} cm -> geometric gain "
          f"{gain:.0e}")
print(f"   focal dwell time l/v2 = {t_core:.1e} s; required focal power")
print(f"   E_1/(l/v2) = {E1_J/t_core:.1e} W concentrated in one core, i.e.")
print(f"   a nanojoule in a zeptosecond: a coherence and phasing problem")
print(f"   for the annihilation-gated monopole array, not an energy problem.")
print(line)

# ---------------------------------------------------------------------------
# The critical path, in order.
# ---------------------------------------------------------------------------
print("WHAT THE MEASUREMENT DECIDES")
print("  Q1. The lock: kinematic or energetic?  DECIDER: monophoton missing")
print("      line at 7.57 / 7.56 GeV photon energy in existing B-factory data.")
print("      Kinematic  -> circulation cannot be shed through a charge, and")
print("                    the free sector is unreachable by this route.")
print("      Energetic  -> the transition exists and Q2 becomes a number.")
print("  Q2. The production cross-section and the linked-pair fraction, which")
print("      set how fast the sector can be populated in a laboratory.")
print("  Q3. Whether a core-scale linked pair is genuinely stable or only")
print("      long-lived: free-evolution simulations find knots untie by")
print("      self-reconnection, so this needs a Gross-Pitaevskii integration.")
print("  None of these revives an actuation route. Closure and compactness")
print("  close that independently (knot_closure_theorem.py); what the line")
print("  hunt settles is whether the free sector exists in the laboratory.")
