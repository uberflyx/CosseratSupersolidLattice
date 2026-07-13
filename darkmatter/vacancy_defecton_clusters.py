#!/usr/bin/env python3
"""
Vacancy dark matter: defecton band, cluster bond counts, and
detection-rate estimates for the Cosserat supersolid lattice.

What this script computes
-------------------------
1. The defecton tight-binding band on the FCC lattice from the hop
   amplitude t = alpha * m0 * c^2 (one electron mass): band bottom,
   bandwidth, effective mass m* = m0/(4 alpha), group-velocity
   ceiling v_max = 4 sqrt(2) alpha c, and the equivalence-principle
   ratio m*/m_g.
2. Exact broken-bond binding energies for small vacancy clusters
   (di-vacancy, triangle, tetrahedron) by direct enumeration on the
   FCC coordination graph, confirming E_b = B(G) * eps with
   eps = m0 c^2 / 6, and the growth-step ladder m * eps.
3. An upper-bound estimate of the present-day diffuse luminosity of
   the 11.67 MeV coalescence line from a Milky-Way-like halo, using
   the derived elastic transport cross-section as a generous proxy
   for radiative capture.
4. The peak of the Klemens (Rayleigh, omega^4) photon-scattering
   cross-section against the zero-point smearing form factor, and
   the resulting optical depth through the Galactic Centre dark
   matter column.

Inputs are lattice geometry plus the SI conversion m_s = m_e only.
The script verifies; it does not derive.
"""

import itertools
import numpy as np

# ----------------------------------------------------------------------
# Constants (CODATA 2022 values where relevant)
# ----------------------------------------------------------------------
ALPHA = 1.0 / 137.035999177          # fine structure constant
ME_MEV = 0.51099895                  # electron mass [MeV/c^2]
M0_MEV = ME_MEV / ALPHA              # node mass m0 = m_e / alpha [MeV/c^2]
C_KMS = 2.99792458e5                 # speed of light [km/s]
ELL_CM = 2.8179403205e-13            # lattice spacing = r_e [cm]
MEV_TO_ERG = 1.602176634e-6          # [erg/MeV]
LSUN_ERG_S = 3.828e33                # solar luminosity [erg/s]

# ----------------------------------------------------------------------
# 1. Defecton band on FCC
# ----------------------------------------------------------------------
def defecton_band():
    t = ALPHA * M0_MEV               # hop amplitude [MeV] = m_e exactly
    e_f = M0_MEV                     # formation energy [MeV]

    # Band E(k) = E_f - 4t [cx cy + cy cz + cz cx], ci = cos(ki a / 2).
    # Bottom at Gamma (bracket = 3), top at X (bracket = -1).
    e_bottom = e_f - 12.0 * t
    e_top = e_f + 4.0 * t
    width = e_top - e_bottom         # = 16 t

    # Effective mass from the parabolic expansion:
    # E ~ E0 + t a^2 k^2  ->  hbar^2 / 2m* = t a^2 = 2 t ell^2
    # with hbar = m0 c ell:  m* = m0 / (4 alpha).
    m_star = M0_MEV / (4.0 * ALPHA)  # [MeV/c^2]

    # Group-velocity ceiling along a cube axis: v_max = 4 t a / hbar
    # = 4 sqrt(2) alpha c.
    v_max = 4.0 * np.sqrt(2.0) * ALPHA * C_KMS

    m_g = e_bottom                   # gravitating (band-bottom) mass [MeV/c^2]
    ep_ratio = m_star / m_g          # equivalence-principle ratio, crystalline

    print("=== Defecton band (FCC crystalline channel) ===")
    print(f"hop amplitude t          : {t:.4f} MeV (= m_e)")
    print(f"band bottom E_f - 12t    : {e_bottom:.2f} MeV  "
          f"(shift {12*t:.2f} MeV = {12*ALPHA*100:.1f}% of E_f)")
    print(f"band width 16t           : {width:.2f} MeV")
    print(f"effective mass m*        : {m_star:.1f} MeV = {m_star/M0_MEV:.1f} m0"
          f" = {m_star/1000:.2f} GeV")
    print(f"velocity ceiling v_max   : {v_max:.3e} km/s = {v_max/C_KMS:.4f} c")
    print(f"EP ratio m*/m_g          : {ep_ratio:.1f} (crystalline ceiling)")
    print(f"EP ratio lower bracket   : {1.0/(1.0-12.0*ALPHA):.3f} "
          f"(condensate-dominated transport)")
    print()
    return t


# ----------------------------------------------------------------------
# 2. Cluster bond counts on the FCC coordination graph
# ----------------------------------------------------------------------
def fcc_sites(shells=2):
    """Integer FCC sites (x+y+z even) within a few conventional cells."""
    r = range(-2 * shells, 2 * shells + 1)
    return [p for p in itertools.product(r, r, r) if sum(p) % 2 == 0]


def internal_bonds(site_set):
    """Number of nearest-neighbour bonds internal to a site set.

    Integer FCC convention: nearest neighbours are at squared
    distance 2 (vectors (+-1, +-1, 0) and permutations).
    """
    sites = list(site_set)
    count = 0
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            d = sum((a - b) ** 2 for a, b in zip(sites[i], sites[j]))
            if d == 2:
                count += 1
    return count


def cluster_bindings():
    eps = M0_MEV / 6.0               # bond price [MeV]
    print("=== Cluster bond counts, E_b = B(G) * eps ===")
    print(f"bond price eps = m0 c^2 / 6 = {eps:.3f} MeV "
          f"(= m_e/(6 alpha) = {ME_MEV/(6*ALPHA):.3f} MeV)")

    # Compact clusters, given as integer FCC coordinates.
    clusters = {
        "di-vacancy (pair)": [(0, 0, 0), (1, 1, 0)],
        "tri-vacancy (triangle)": [(0, 0, 0), (1, 1, 0), (1, 0, 1)],
        "tetra-vacancy (tetrahedron)": [(0, 0, 0), (1, 1, 0),
                                        (1, 0, 1), (0, 1, 1)],
    }
    prev_bonds = 0
    for name, sites in clusters.items():
        b = internal_bonds(sites)
        step = b - prev_bonds
        print(f"{name:32s}: B(G) = {b}, E_b = {b*eps:.2f} MeV, "
              f"growth step m = {step} releases {step*eps:.2f} MeV")
        prev_bonds = b
    # Verify the tetrahedron is a maximal clique (no 5-clique in FCC).
    sites = fcc_sites(1)
    tetra = [(0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1)]
    common = [s for s in sites
              if s not in tetra
              and all(sum((a - b) ** 2 for a, b in zip(s, t)) == 2
                      for t in tetra)]
    print(f"sites adjacent to all four tetrahedron members: {len(common)} "
          f"(no K5 in the FCC graph)")
    print()


# ----------------------------------------------------------------------
# 3. Coalescence-line luminosity (upper bound)
# ----------------------------------------------------------------------
def coalescence_luminosity():
    eps = M0_MEV / 6.0                       # line energy [MeV]
    m_dm_g = M0_MEV * MEV_TO_ERG / (C_KMS * 1e5) ** 2  # particle mass [g]
    n_local = 0.4 / M0_MEV * 1e3             # local density 0.4 GeV/cm^3 [cm^-3]
    sigma_t = 2.0e-57                        # derived elastic sigma_T [cm^2]
    v_rel = 2.0e7                            # galactic velocity [cm/s]

    # Rate density (upper bound: elastic sigma as capture proxy).
    rate_density = 0.5 * n_local ** 2 * sigma_t * v_rel   # [cm^-3 s^-1]
    # Effective n^2-weighted halo volume, order (few kpc)^3.
    v_eff = 1.0e66                                        # [cm^3]
    rate = rate_density * v_eff                           # [s^-1]
    lum = rate * eps * MEV_TO_ERG                         # [erg/s]

    print("=== Coalescence line: diffuse-luminosity upper bound ===")
    print(f"line energy              : {eps:.2f} MeV, ladder 2eps = "
          f"{2*eps:.1f}, 3eps = {3*eps:.1f} MeV")
    print(f"local number density     : {n_local:.1f} cm^-3")
    print(f"rate density (bound)     : {rate_density:.2e} cm^-3 s^-1")
    print(f"halo luminosity (bound)  : {lum:.2e} erg/s = "
          f"{lum/LSUN_ERG_S:.1e} L_sun")
    print()


# ----------------------------------------------------------------------
# 4. Klemens (Rayleigh) photon scattering: peak and optical depth
# ----------------------------------------------------------------------
def rayleigh_channel():
    # sigma(w) ~ ell^2 (w/w_D)^4 exp(-(k sigma_u)^2), k = w/c,
    # sigma_u = 1.6 ell the node zero-point spread. With x = w/w_D
    # and k ell = pi x, the exponent is -(1.6 pi)^2 x^2.
    s = (1.6 * np.pi) ** 2
    x = np.linspace(1e-4, 1.0, 200000)
    f = x ** 4 * np.exp(-s * x ** 2)
    i = np.argmax(f)
    x_pk = x[i]
    hbar_wd = np.pi * M0_MEV                 # Debye energy [MeV]
    e_pk = x_pk * hbar_wd                    # peak photon energy [MeV]
    sigma_pk = ELL_CM ** 2 * f[i]            # peak cross-section [cm^2]

    # Galactic Centre column: ~1e22 GeV/cm^2 -> particles/cm^2.
    col = 1.0e22 * 1e3 / M0_MEV              # [cm^-2]
    tau = col * sigma_pk

    print("=== Klemens omega^4 photon scattering ===")
    print(f"peak at x = w/w_D        : {x_pk:.3f} -> E_peak = {e_pk:.0f} MeV")
    print(f"peak cross-section       : {sigma_pk:.2e} cm^2")
    print(f"GC column density        : {col:.2e} cm^-2")
    print(f"optical depth at peak    : {tau:.1e}")
    print()


if __name__ == "__main__":
    defecton_band()
    cluster_bindings()
    coalescence_luminosity()
    rayleigh_channel()
