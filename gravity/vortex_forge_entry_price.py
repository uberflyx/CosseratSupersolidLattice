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
     its energy is fixed by the line tension over one lattice spacing:
     E_1 = 2 pi^2 f_s m0 c^2 = 1.38 GeV at the adopted superfluid fraction
     f_s = 0.99932, a fifth of a nanojoule. Any forge is a GeV-scale
     process per ring. (This is also why nothing in chemistry, materials
     science, or tabletop physics has ever stirred the vacuum: the entry
     ticket is a gigaelectronvolt deposited into a femtometre-cubed core.)

  2. The lock decider. Whether circulation can leave a charged particle at
     all is the rolling-lock question: kinematic (forbidden) or energetic
     (allowed above threshold). The SAME question fixed the dark sector's
     fate. It is testable now: if shedding is allowed, then e+ e-
     collisions above the ring's dispersion floor,
     sqrt(s) = E_1 [1 + (R0/ell)^2] = 5.27 GeV, can produce a free ring,
     which is electrically dark and leaves the detector unseen. At a
     B factory this is the classic monophoton search: e+ e- -> gamma +
     invisible. A ring has no rest mass, so the observable is not a
     missing-mass peak but a monochromatic photon line fixed by the
     two-body balance on the ring dispersion: 7.57 GeV for a single ring
     at the Upsilon(4S) (solve in cosmology/vortex_hadron_ladder.py). The
     linked pair, the storable product, cannot be made there: its floor
     2 E_1 + E_12 + 2 P_0 c = 10.73 GeV lies above the resonance. The
     BaBar and Belle II single-photon programmes bound monophoton peaks at
     the femtobarn scale, but a photon recoiling against anything with rest
     mass cannot exceed sqrt(s)/2 = 5.29 GeV and the ring's line lies above
     that endpoint, so whether those bounds reach the ring is unchecked. A
     peak is a forge; its absence is a coupling bound.

  3. The storage problem. A plain ring dies on first contact with charged
     matter (the capture channel), so forged rings must be born protected:
     shed as a LINKED PAIR (helicity from birth, production floor 10.73 GeV,
     above the Upsilon(4S)) or tied against the stacking-fault trap that
     both stores and, run in reverse as a pinned partial loop driven at GeV
     energy density, is itself the second forge candidate. At the
     Upsilon(4S) only the single ring, the messenger, is within reach.

  4. The acoustic forge (the in-house route). A converging second-sound
     pulse nucleates a ring where its amplitude reaches order unity at the
     core scale. The focal energy must still be ~E_1 inside ~l^3, so the
     required convergence gain from a metre-scale phased annihilation
     array is enormous but finite, and it is a coherence problem, not an
     energy problem: E_1 is a nanojoule-scale (0.22 nJ) deposit.

Every number below is from framework constants; model inputs are labelled.
"""

import numpy as np
from scipy.optimize import brentq
from scipy.special import ellipk, ellipe

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
F_S     = 0.99932  # the Leggett ceiling adopted as the working superfluid fraction
rho_s   = F_S * rho                        # superfluid density [kg/m^3]
J_per_MeV = 1.602176634e-13

line = "-" * 78
print(__doc__.strip().splitlines()[0])
print(line)

# ---------------------------------------------------------------------------
# 1. The fuel quantum: smallest free ring.
#    Line tension rho_s kappa^2/(4 pi) = (pi f_s) m0 c^2 / l at ln -> 1,
#    times one lattice spacing of circumference 2 pi l: E_1 = 2 pi^2 f_s m0 c^2.
#    The ground ring carries that energy at R0, where x [ln(8x) - 2] = 1
#    with x = R0/l (core xi = l), so R0 carries no f_s.
# ---------------------------------------------------------------------------
E1_MeV = 2.0 * np.pi**2 * F_S * m0_MeV
E1_J = E1_MeV * J_per_MeV
x0 = brentq(lambda x: x * (np.log(8.0 * x) - 2.0) - 1.0, 1.0, 3.0)
I0 = rho_s * kappa * np.pi * (x0 * l_m)**2  # ground-ring impulse [kg m/s]
P0c_MeV = E1_MeV * x0**2                   # the same impulse x c, = E_1 x0^2
v1_ring = kappa / (4 * np.pi * l_m)        # self-induced speed ~ c/2 scale
print("1. THE FUEL QUANTUM (smallest free ring)")
print(f"   f_s = {F_S}")
print(f"   E_1 = 2 pi^2 f_s m0 c^2 = {E1_MeV/1e3:.3f} GeV = {E1_J*1e9:.2f} nJ")
print(f"   ground radius R0 = {x0:.2f} l; impulse P_0 = rho_s kappa pi R0^2 = "
      f"{I0:.2e} kg m/s (P_0 c = E_1 x0^2 = {P0c_MeV/1e3:.2f} GeV)")
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
#    cosmology/vortex_hadron_ladder.py and is repeated here in closed form:
#        channel floor:  sqrt(s) >= min_R [E(R) + P(R) c] = E_1 (1 + x0^2)
#        Upsilon(4S) tag: E_gamma = P(R) c at E(R) + P(R) c = sqrt(s)
#        linked pair floor: 2 E_1 + E_12 + 2 P_0 c (above the Upsilon(4S))
#    An invariant-mass reading would put the tag near 5.2 GeV; the ~2 GeV
#    gap between the readings is itself a signature no rest-mass object makes.
# ---------------------------------------------------------------------------
sqrt_s = 10.58   # GeV
E_GeV = lambda x: E1_MeV / 1e3 * x * (np.log(8.0 * x) - 2.0)   # ring energy
Pc_GeV = lambda x: E1_MeV / 1e3 * x**2                         # ring impulse x c
floor = E_GeV(x0) + Pc_GeV(x0)
x_ups = brentq(lambda x: E_GeV(x) + Pc_GeV(x) - sqrt_s, x0, 12.0)
E_gam = Pc_GeV(x_ups)
# coaxial pair energy at spacing 3 R0 (Neumann form), rho_s kappa^2 l = 2 E_1
m_k = 4 * x0**2 / ((2 * x0)**2 + (3 * x0)**2)
k_k = np.sqrt(m_k)
E12 = 2 * E1_MeV / 1e3 * x0 * ((2 / k_k - k_k) * ellipk(m_k) - (2 / k_k) * ellipe(m_k))
pair_floor = 2 * E_GeV(x0) + E12 + 2 * Pc_GeV(x0)
print("2. LOCK DECIDER: e+e- -> gamma + invisible (single ring)")
print(f"   channel floor sqrt(s) >= E_1 (1 + x0^2) = {floor:.2f} GeV")
print(f"   dispersion-solved tag at the Upsilon(4S): monochromatic photon "
      f"E_gamma = {E_gam:.2f} GeV in the CM")
print(f"   linked pair: floor 2 E_1 + E_12 + 2 P_0 c = {pair_floor:.2f} GeV "
      f"(E_12 = {E12*1e3:.0f} MeV) > {sqrt_s} GeV -> not made there")
print("   The existing constraint: BaBar's invisible dark-photon search")
print("   (53/fb, arXiv:1702.03327) bounds monophoton peaks at a production")
print("   cross-section of order one femtobarn, but a photon recoiling against")
print(f"   anything with rest mass cannot exceed sqrt(s)/2 = {sqrt_s/2:.2f} GeV, and the")
print("   ring's line lies above that endpoint, so whether the bound reaches")
print("   the ring is not yet checked. If it does, sigma_shed < ~1 fb.")
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
print(f"   A dedicated line hunt at {E_gam:.2f} GeV photon energy in the")
print(f"   archived and incoming monophoton samples is the single most")
print(f"   consequential measurement this framework can request.")
print(line)

# ---------------------------------------------------------------------------
# 3. Storage: born-protected products.
# ---------------------------------------------------------------------------
print("3. STORAGE (a plain ring dies on first charged contact)")
print(f"   linked-pair shedding: floor {pair_floor:.2f} GeV, helicity at birth,")
print(f"   out of reach at the Upsilon(4S), where only the messenger ring forms;")
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
# Second sound rides the crystal's longitudinal modulus over its transverse
# one, C11/C44 = (3 - 4 N^2)/(1 - N^2) at the vacuum's N^2 = 1/pi, with the
# inertia cancelling: v2 = 1.592 c, independent of f_s.
N2 = 1.0 / np.pi
v2 = np.sqrt((3.0 - 4.0 * N2) / (1.0 - N2)) * c
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
print(f"      line at {E_gam:.2f} GeV photon energy in existing B-factory data.")
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
