"""
Primordial uncrystallised pockets: seed masses, growth budget, and the
overmassive-ratio trend.

This script reproduces, from first principles and a small set of lattice
constants, the quantitative claims about black holes left behind by an
incomplete vacuum-crystallisation transition:

  1. the minimum stable pocket mass M_min for several readings of the
     pre-crystallisation fluid density;
  2. the pocket vs Schwarzschild radius across the observed mass range;
  3. the Eddington-limited growth of a floor-mass seed to the redshift of
     the most distant luminous quasars, using the correct mass e-folding
     (Salpeter) time rather than the Eddington luminosity timescale;
  4. the dark-accretion lever in the lightless fluid phase and why the
     fluid's compressional stiffness bounds it;
  5. the overmassive black-hole-to-stellar-mass ratio trend with redshift.

All inputs are physical constants or quantities derived elsewhere in the
framework. No observational quantity is fitted.
"""

import numpy as np

# --------------------------------------------------------------------------
# Physical constants (CODATA 2018/2022; SI units unless noted)
# --------------------------------------------------------------------------
c = 2.99792458e8          # speed of light [m/s]
G = 6.67430e-11           # Newtonian constant [m^3 kg^-1 s^-2]
hbar = 1.054571817e-34    # reduced Planck constant [J s]
m_e = 9.1093837015e-31    # electron mass [kg]
alpha = 7.2973525693e-3   # fine-structure constant [-]
sigma_T = 6.6524587321e-29  # Thomson cross-section [m^2]
m_p = 1.67262192369e-27   # proton mass [kg]

# Derived lattice quantities
m_0 = m_e / alpha                 # fundamental node mass m_e/alpha [kg]
ell = 2.8179403262e-15            # lattice spacing = classical electron radius [m]
ell_P = np.sqrt(hbar * G / c**3)  # Planck length [m]

# Astronomical units
M_sun = 1.98892e30        # solar mass [kg]
Myr = 3.15576e13          # 1 Myr [s]


def schwarzschild_radius(M):
    """Schwarzschild radius 2GM/c^2 [m]."""
    return 2.0 * G * M / c**2


def pocket_radius(M, rho_fluid):
    """Uniform-sphere radius (3M / 4 pi rho)^(1/3) [m]."""
    return (3.0 * M / (4.0 * np.pi * rho_fluid)) ** (1.0 / 3.0)


def M_min_pocket(rho_fluid):
    """Minimum stable pocket mass from R_pocket = R_S (analytic), [kg]."""
    return np.sqrt(3.0 * c**6 / (32.0 * np.pi * G**3 * rho_fluid))


# --------------------------------------------------------------------------
# 1. Minimum stable pocket mass for three density readings
# --------------------------------------------------------------------------
rho_crystal = 4.0 * m_0 / (np.sqrt(2.0) * ell) ** 3   # FCC: a = sqrt(2) * ell
rho_fluid_D4 = 3.0 * rho_crystal                       # D4 projection (3 stacking layers)
# Formal "strong-coupling" density chosen so that M_min coincides with the
# one-spacing horizon (shown for completeness; not the physical fluid density):
rho_strong = 3.0 * c**2 / (8.0 * np.pi * G * ell**2)

print("=" * 70)
print("1. MINIMUM STABLE POCKET MASS")
print("=" * 70)
print(f"  rho_crystal (FCC)          = {rho_crystal:.3e} kg/m^3")
print(f"  rho_fluid (D4, 3*crystal)  = {rho_fluid_D4:.3e} kg/m^3")
for label, rho in [("D4 fluid (3 rho_cryst)", rho_fluid_D4),
                   ("crystal density", rho_crystal),
                   ("formal strong-coupling", rho_strong)]:
    Mmin = M_min_pocket(rho)
    Rs = schwarzschild_radius(Mmin)
    print(f"  {label:24s}: M_min = {Mmin/M_sun:10.3e} Msun "
          f"({Mmin:9.3e} kg),  R_S = {Rs/1e3:9.3e} km")

M_star = c**2 * ell / (2.0 * G)   # evaporation floor: r_s = ell
print(f"\n  evaporation floor M_star    = {M_star:.3e} kg "
      f"(r_s = {schwarzschild_radius(M_star)/ell:.3f} * ell)")

# --------------------------------------------------------------------------
# 2. Pocket vs Schwarzschild radius across the observed mass range
# --------------------------------------------------------------------------
print("\n" + "=" * 70)
print("2. POCKET RADIUS vs SCHWARZSCHILD RADIUS (rho_fluid = 3 rho_crystal)")
print("=" * 70)
print(f"  {'M [Msun]':>10} {'R_S [m]':>12} {'R_pocket [m]':>14} {'R_pocket/R_S':>14}")
for M_solar in [M_min_pocket(rho_fluid_D4) / M_sun, 1e2, 1e6, 1e9, 1e10]:
    M = M_solar * M_sun
    Rs = schwarzschild_radius(M)
    Rp = pocket_radius(M, rho_fluid_D4)
    print(f"  {M_solar:10.2e} {Rs:12.3e} {Rp:14.3e} {Rp/Rs:14.3e}")

# --------------------------------------------------------------------------
# 3. Eddington-limited growth: the CORRECT e-folding time
# --------------------------------------------------------------------------
# The Eddington luminosity is L_Edd = 4 pi G M m_p c / sigma_T. The black-hole
# MASS e-folds on the Salpeter time, NOT on t_Edd = M c^2 / L_Edd:
#     t_Salpeter = (eps / (1 - eps)) * sigma_T * c / (4 pi G m_p)
# where eps is the radiative efficiency. Confusing t_Edd (~450 Myr) with the
# mass e-folding time understates the number of e-foldings by ~1/eps.
print("\n" + "=" * 70)
print("3. EDDINGTON-LIMITED GROWTH OF A FLOOR-MASS SEED")
print("=" * 70)
t_Edd = sigma_T * c / (4.0 * np.pi * G * m_p)   # Eddington timescale M c^2 / L_Edd
print(f"  Eddington timescale t_Edd (= M c^2 / L_Edd) = {t_Edd/Myr:.1f} Myr")

age_z7 = 750.0 * Myr      # approximate cosmic age at z ~ 7
target = 1e9 * M_sun      # observed luminous-quasar mass at z ~ 7

for eps in [0.057, 0.1, 0.2]:
    t_S = (eps / (1.0 - eps)) * t_Edd    # mass e-folding (Salpeter) time
    n_efold = age_z7 / t_S
    print(f"\n  radiative efficiency eps = {eps:.3f}: "
          f"Salpeter time t_S = {t_S/Myr:5.1f} Myr, "
          f"{n_efold:5.1f} e-folds in 750 Myr")
    for M_seed_solar in [28.0, 48.0, 300.0]:
        M_final = M_seed_solar * M_sun * np.exp(n_efold)
        print(f"      seed {M_seed_solar:6.1f} Msun -> "
              f"{M_final/M_sun:9.2e} Msun by z=7")
    # time to reach 1e9 from the floor seed
    M_seed = 28.0 * M_sun
    t_reach = t_S * np.log(target / M_seed)
    print(f"      time for 28 Msun -> 1e9 Msun: {t_reach/Myr:.0f} Myr "
          f"(available: 750 Myr)")

# --------------------------------------------------------------------------
# 4. Dark-accretion lever in the lightless fluid phase
# --------------------------------------------------------------------------
# In the fluid phase mu = 0, so there are no photons and no Eddington limit.
# Accretion is instead set by the Bondi rate, with the relevant signal speed
# being the longitudinal (compression) speed v_p, which is enormous in the
# stiff superfluid. The Bondi radius r_B = G M / v_p^2 is therefore tiny.
print("\n" + "=" * 70)
print("4. DARK ACCRETION: BONDI SUPPRESSION IN THE STIFF FLUID")
print("=" * 70)
vp_over_c = 1e20                     # longitudinal speed ~ sqrt(K/mu) in fluid
v_p = vp_over_c * c
M_demo = M_min_pocket(rho_fluid_D4)  # a floor-mass seed
r_B = G * M_demo / v_p**2            # Bondi radius with compression speed
print(f"  fluid compression speed v_p     = {vp_over_c:.0e} c")
print(f"  floor seed mass                 = {M_demo/M_sun:.1f} Msun")
print(f"  Bondi radius r_B = G M / v_p^2  = {r_B:.3e} m")
print(f"  Schwarzschild radius R_S        = {schwarzschild_radius(M_demo):.3e} m")
print(f"  r_B / R_S                       = {r_B/schwarzschild_radius(M_demo):.3e}")
print("  -> the lightless phase removes the Eddington cap, but the stiffness")
print("     keeps the accretable region small, so dark accretion is not fast.")

# --------------------------------------------------------------------------
# 5. Overmassive ratio trend with redshift
# --------------------------------------------------------------------------
print("\n" + "=" * 70)
print("5. OVERMASSIVE RATIO M_BH / M_star vs REDSHIFT")
print("=" * 70)
ratio_local = 1e-3        # ~0.1% in the local Universe
print("  mean trend, (1+z)^{5/2}:")
for z in [0, 2, 4, 5, 7, 10]:
    ratio = ratio_local * (1.0 + z) ** 2.5    # (1+z)^{5/2} evolution
    print(f"    z = {z:4.1f}:  M_BH/M_star = {ratio:.3e}  "
          f"(overmassive factor {ratio/ratio_local:6.1f})")

# Individual measured systems plotted in the figure (M_BH and M_star in Msun).
# Sources: UHZ1 (Bogdan 2024; Natarajan 2024), GN-z11 (Maiolino 2024),
# A2744-QSO1 (Furtak 2024), the dormant BH GN-1001830 (Juodzbalis 2024),
# J0313-1806 (Wang 2021; ratio uses the host dynamical mass).
print("\n  measured high-z systems:")
systems = [
    # name, z, M_BH, M_star (or dynamical), note
    ("UHZ1",          10.1, 4.0e7,  4.0e7,  "M_BH ~ host stellar mass"),
    ("GN-z11",        10.6, 1.6e6,  8.0e8,  "near the local relation"),
    ("A2744-QSO1",     7.05, 4.0e7,  4.0e8,  "ratio >= 3%, possibly ~1"),
    ("GN-1001830",     6.68, 4.0e8,  1.0e9,  "dormant, ratio ~ 0.4"),
    ("J0313-1806",     7.64, 1.6e9,  4.0e10, "quasar, M_dyn proxy"),
]
print(f"  {'object':12} {'z':>5} {'M_BH[Msun]':>11} {'M_*[Msun]':>11} {'ratio':>8}  note")
for name, z, mbh, mstar, note in systems:
    print(f"  {name:12} {z:5.2f} {mbh:11.2e} {mstar:11.2e} {mbh/mstar:8.3f}  {note}")
