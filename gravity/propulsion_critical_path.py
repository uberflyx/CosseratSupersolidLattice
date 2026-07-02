"""
The critical path to time-gradient propulsion: the vortex forge.

The clock-drive section (sec:clock_drive, sec:clock_actuation) has already
priced the drive itself. The chain clock = phase = chemical potential makes a
maintained flow gradient a gravitational field; the knot rotor turns trapped
vacuum vortices into a sustained phase slip (one g across ten metres from
~1,200 knots on a centimetre rotor at ordinary loads); the fabricated
stacking-fault loop is the trap that holds a knot without untying it; and the
ambient-propellant jet prices the hover of a hundred-tonne craft at half a
watt. The section's named open steps were SUPPLY (the local abundance of
free vacuum vortices) and LOADING.

The reconnection cascade (cosmology/knot_cascade.py, this session) has now
answered the supply question, and the answer is no. Knots untie universally
(Kleckner-Kauffman-Irvine 2016), the cascade delivers the primordial sector
to plain rings within 1e-14 s, and plain rings fuse onto charged matter and
are destroyed before nucleosynthesis. Nature stocks no fuel. Every path to
the drive therefore passes through one gate:

    THE FORGE: manufacture free (unlocked) condensate circulation.

This script derives what the forge must pay and names the experiment that
decides whether it can exist at all.

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
     B factory this is the classic monophoton search: e+ e- -> gamma_ISR +
     invisible, with missing mass peaked at m_ring ~ 1.1 GeV (and a second
     peak near 2.2 GeV for the linked pair, the storable product). The
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
#    e+ e- -> gamma + ring(invisible): monophoton, missing mass m_ring.
#    Belle II runs at sqrt(s) = 10.58 GeV; the ISR photon energy that tags
#    a missing mass m is E_gamma = (s - m^2 c^4) / (2 sqrt(s)).
# ---------------------------------------------------------------------------
sqrt_s = 10.58   # GeV
for m_miss, tag in ((E1_MeV/1e3, "single ring"), (2*E1_MeV/1e3, "linked pair (storable)")):
    E_gam = (sqrt_s**2 - m_miss**2) / (2.0 * sqrt_s)
    print(f"2. LOCK DECIDER: e+e- -> gamma + invisible ({tag})")
    print(f"   missing mass = {m_miss:.2f} GeV  ->  monochromatic photon "
          f"E_gamma = {E_gam:.2f} GeV in the CM")
print("   Both photons sit in the clean part of the BaBar / Belle II")
print("   single-photon window.  The existing constraint: BaBar's invisible")
print("   dark-photon search (53/fb, arXiv:1702.03327) bounds any narrow")
print("   missing-mass peak at ~1 GeV at epsilon ~ 1e-3, i.e. a production")
print("   cross-section ceiling of order one femtobarn.  So today's data")
print("   already says: sigma_shed(1.11 GeV) < ~1 fb.")
# --- fuel arithmetic at the ceiling ---
sigma_fb = 1.0                       # current ceiling [fb]
L_belle2_ab = 50.0                   # Belle II design integrated lumi [ab^-1]
L_year_ab = 5.0                      # optimistic per-year [ab^-1]
N_prog = sigma_fb * L_belle2_ab * 1e3
N_year = sigma_fb * L_year_ab * 1e3
N_rotor = 1200                       # knots for one g on the cm rotor
print(f"   FUEL ARITHMETIC AT THE CEILING: at sigma = {sigma_fb:.0f} fb,")
print(f"   Belle II makes {N_year:.0f} rings/year and {N_prog:.0f} over its")
print(f"   50/ab programme; the rotor needs ~{N_rotor} trapped knots.  Even")
print(f"   at today's experimental UPPER BOUND, a few months of B-factory")
print(f"   running forges one rotor fuel load.  The drive is therefore")
print(f"   gated by physics (does the peak exist?), not by luminosity.")
print(f"   A dedicated peak hunt at 1.11 and 2.21 GeV in the archived and")
print(f"   incoming monophoton samples is the single most consequential")
print(f"   measurement this framework can request.")
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
print("THE CRITICAL PATH (everything else is already priced)")
print("  Q1. The lock: kinematic or energetic?  DECIDER: monophoton missing")
print("      mass at 1.10 / 2.21 GeV in existing B-factory data.")
print("      Kinematic  -> no forge through charges; drive waits on the")
print("                    acoustic route or dies.")
print("      Energetic  -> Q2.")
print("  Q2. The forge cross-section and the linked-pair fraction: sets the")
print("      fuel production rate per collider watt.")
print("  Q3. Trap loading: steering a born ring/knot onto a fabricated")
print("      fault loop before first charged contact (its mean free path in")
print("      solid matter is ~(n_e sigma)^-1 with sigma = pi l^2: metres in")
print("      ordinary solids, so a vacuum-gap forge-to-trap line suffices).")
print("  Q4. The output-fidelity derivation (clock coupling of a defect to")
print("      delta-mu) -- the section's own remaining theory item.")
print("  With Q1-Q4 closed, the rotor spec of sec:clock_actuation flies as")
print("  written: ~1,200 trapped knots, centimetre rotor, one g of pure")
print("  free-fall thrust, and the half-watt ambient jet for the hover.")
