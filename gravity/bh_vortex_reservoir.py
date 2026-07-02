"""
The black-hole vortex reservoir: inventory, dynamics, and escape flux.

The black-hole chapter establishes a shear-melted, superfluid interior in
which "defects persist topologically as zero-energy vortex lines but carry
no mass" (sec:interior_geometry). Read against the superfluid chapter's two
topological sectors, that sentence has a sharp meaning: the horizon is a
lock-breaker. An infalling charge's dislocation half (Burgers vector, mass,
coupling to light) is stripped adiabatically and paid into the exterior
strain field; its circulation half is a winding of the condensate phase,
and the condensate survives the melt, so the winding persists inside as a
FREE vortex line. The free vortex sector, primordially emptied outside
horizons (cosmology/knot_cascade.py), survives inside them.

This script prices the reservoir in four stages, all from framework
constants, model inputs labelled.

  1. INVENTORY: circulation quanta delivered per accreted solar mass.
     Neutral matter delivers +/- in equal measure, so the delivered tangle
     is net-neutral and can annihilate.

  2. TANGLE DYNAMICS: the interior cascade clock scales as R^2/kappa, which
     at horizon radius is ~10^13 s. A quiescent (Schwarzschild) interior
     therefore relaxes to empty on ~Myr timescales, consistent with the
     chapter's mass-outside bookkeeping; an accreting hole sustains a
     steady-state tangle set by feed rate against Vinen decay.

  3. THE KERR ARRAY: a superfluid can only co-rotate by quantised vortices
     (Feynman rule n_v = 2 Omega / kappa), so a spinning interior holds a
     permanent, spin-protected net array. Its line count and its share of
     the hole's angular momentum and rotational energy are computed.

  4. ESCAPE FLUX: outside the horizon, accretion collisions above the ring
     threshold E_1 = 1.11 GeV shed free rings at the same sigma_shed the
     B-factory monophoton search bounds. Per-proton escape probability is
     ~ sigma_shed / sigma_capture; the resulting dark-ring luminosity of
     an Eddington source and its flux at Earth are given as functions of
     sigma_shed at the current experimental ceiling.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Framework constants.
# ---------------------------------------------------------------------------
alpha  = 1.0 / 137.035999177
m_e    = 9.1093837015e-31            # kg
m0     = m_e / alpha                 # node mass [kg]
c      = 2.99792458e8
G      = 6.6743e-11
l_m    = 2.8179403205e-15            # lattice spacing r_e [m]
kappa  = 2.0 * np.pi * l_m * c       # circulation quantum [m^2/s]
rho_s  = 0.8 * (m0 / l_m**3)         # superfluid density [kg/m^3]
m_p    = 1.67262192e-27              # proton mass [kg]
Msun   = 1.989e30                    # kg
J_per_MeV = 1.602176634e-13
E1_J   = (8*np.pi**2/5) * (m0*c**2)  # smallest-ring energy [J] ~ 1.11 GeV
sigma_cap = np.pi * l_m**2           # geometric capture cross-section [m^2]

line = "-" * 78
print(__doc__.strip().splitlines()[0]); print(line)

# ---------------------------------------------------------------------------
# 1. INVENTORY per accreted solar mass (solar composition: ~0.87 electrons
#    and protons per nucleon; MODEL INPUT, composition-dependent O(1)).
# ---------------------------------------------------------------------------
f_charge = 0.87
N_pairs_per_Msun = f_charge * Msun / m_p
print("1. INVENTORY: circulation delivered by accretion")
print(f"   charge pairs per accreted solar mass: {N_pairs_per_Msun:.1e}")
print(f"   each pair leaves one + and one - circulation quantum in the melt:")
print(f"   a 10 Msun hole built from ~10 Msun of matter has swallowed")
print(f"   ~{10*N_pairs_per_Msun:.0e} quanta of each sign. Net: ~zero.")
print(f"   The delivered tangle is vortex-antivortex, so it CAN decay.")
print(line)

# ---------------------------------------------------------------------------
# 2. TANGLE DYNAMICS in the interior.
#    Cascade/self-induction clock: t(R) ~ 4 pi R^2 / (kappa ln), as in
#    knot_cascade.py, now at horizon scale. Vinen steady state under feed:
#    dL/dt = S - chi kappa L^2, chi ~ 0.1 (MODEL INPUT, superfluid
#    turbulence value), L = line length per volume.
# ---------------------------------------------------------------------------
M_bh   = 10 * Msun
r_s    = 2*G*M_bh/c**2
lnfac  = 10.0                         # MODEL INPUT: tension/induction log
t_casc = 4*np.pi*r_s**2/(kappa*lnfac)
print("2. TANGLE DYNAMICS (10 Msun, r_s = %.1f km)" % (r_s/1e3))
print(f"   horizon-scale cascade clock: t ~ 4 pi r_s^2/(kappa ln) = "
      f"{t_casc:.1e} s ~ {t_casc/3.15e13:.0f} Myr")
print(f"   -> a quiescent interior relaxes to empty in ~Myr: Schwarzschild")
print(f"      holes hold no standing tangle, as the chapter's mass-outside")
print(f"      bookkeeping requires. The relaxation energy leaves in first")
print(f"      sound and is absorbed at the boundary shell.")
# steady state under Eddington feed
L_Edd  = 1.26e31 * (M_bh/Msun)        # W
Mdot   = L_Edd / (0.1 * c**2)         # kg/s (10% efficiency; MODEL INPUT)
S_feed = (f_charge * Mdot/m_p) * (2*np.pi*l_m*3) / ((4/3)*np.pi*r_s**3)
# each pair injects ~ one small loop of length ~ 2 pi * few l per quantum
chi    = 0.1                          # MODEL INPUT (Vinen)
L_star = np.sqrt(S_feed/(chi*kappa))
E_line = (4*np.pi/5)*lnfac*(m0*c**2)/l_m      # J per metre of line
E_tangle = E_line * L_star * (4/3)*np.pi*r_s**3
print(f"   Eddington-fed steady state: L* = {L_star:.1e} m/m^3,")
print(f"   stored tangle energy ~ {E_tangle:.1e} J = "
      f"{E_tangle/(M_bh*c**2):.1e} of the hole's mass-energy: a trace.")
print(line)

# ---------------------------------------------------------------------------
# 3. THE KERR ARRAY (spin-protected, permanent).
#    Feynman rule: n_v = 2 Omega / kappa. Near-extremal: Omega_H = c/(2 r_+),
#    r_+ = G M / c^2. J_Kerr(extremal) = G M^2 / c.
# ---------------------------------------------------------------------------
r_plus  = G*M_bh/c**2
Omega_H = c/(2*r_plus)
n_v     = 2*Omega_H/kappa
N_lines = n_v * np.pi * r_plus**2
# condensate share of J: rigid-body-mimicking array, I_s = (2/5) M_s r^2
M_s = rho_s * (4/3)*np.pi*r_plus**3
J_array = (2/5) * M_s * r_plus**2 * Omega_H
J_kerr  = G*M_bh**2/c
E_array = 0.5 * (2/5) * M_s * (Omega_H*r_plus)**2
E_spin  = (1 - 1/np.sqrt(2)) * M_bh * c**2
print("3. THE KERR ARRAY (near-extremal 10 Msun)")
print(f"   Omega_H = {Omega_H:.1e} rad/s  ->  Feynman density n_v = "
      f"{n_v:.1e} lines/m^2")
print(f"   total quantised vortex lines threading the interior: "
      f"N = {N_lines:.1e}")
print(f"   array share of angular momentum: J_array/J_Kerr = "
      f"{J_array/J_kerr:.1e}")
print(f"   array rotational energy / extractable spin energy = "
      f"{E_array/E_spin:.1e}")
print(f"   -> a ~1e19-line quantised bundle, holding ~1e-4 of the spin:")
print(f"      real, permanent, spin-protected, and a minority shareholder.")
print(line)

# ---------------------------------------------------------------------------
# 4. ESCAPE FLUX from the exterior forge.
#    Per accreted proton, P(shed a ring that escapes) ~ sigma_shed/sigma_cap
#    (MODEL INPUT: order-of-magnitude branching against capture).
# ---------------------------------------------------------------------------
sigma_shed_fb = 1.0                   # current B-factory ceiling [fb]
sigma_shed = sigma_shed_fb * 1e-43    # m^2
P_escape = sigma_shed / sigma_cap
N_dot_rings = P_escape * (Mdot/m_p)
d = 3.086e19                          # 1 kpc in m
flux = N_dot_rings/(4*np.pi*d**2)
print("4. ESCAPE FLUX (Eddington source, sigma_shed at the BaBar ceiling)")
print(f"   branching per proton: sigma_shed/sigma_cap = {P_escape:.1e}")
print(f"   dark-ring luminosity: {N_dot_rings:.1e} rings/s")
print(f"   flux at 1 kpc: {flux:.2e} rings/m^2/s "
      f"({flux*3.15e7:.1e} per m^2 per year)")
print(f"   Each ring scatters on nuclei with sigma = 0.25 barn, so any")
print(f"   detected flux is conspicuous; the luminosity scales linearly")
print(f"   with sigma_shed, tying every microquasar's dark-ring output to")
print(f"   the same number Belle II measures.")
print(line)

print("VERDICT")
print("  Outside horizons the free vortex sector is empty (the primordial")
print("  cascade + capture). Inside, the melt breaks the rolling lock for")
print("  free: Schwarzschild interiors relax to empty in ~Myr; Kerr")
print("  interiors permanently hold a ~1e19-line quantised Feynman array")
print("  carrying ~1e-4 of the spin; accreting holes add a live tangle of")
print("  ~1e-12 of the mass-energy. The universe's free vortices exist,")
print("  and they are all inside or beside black holes -- fuel that only")
print("  a forge, not a harvest, can practically replace.")
