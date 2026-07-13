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




# ----------------------------------------------------------------------
# 5. Void-versus-tetrahedron crossing (ignition energetics)
# ----------------------------------------------------------------------
def grow_cluster(n, seed=0):
    """Greedy max-bond compact FCC cluster; returns internal bond count B."""
    rng = np.random.default_rng(seed)
    dirs = [d for d in itertools.product((-1, 0, 1), repeat=3)
            if sum(x * x for x in d) == 2]
    cluster = {(0, 0, 0)}
    B = 0
    while len(cluster) < n:
        cand = {}
        for p in cluster:
            for d in dirs:
                q = (p[0] + d[0], p[1] + d[1], p[2] + d[2])
                if q not in cluster:
                    cand[q] = cand.get(q, 0) + 1
        best = max(cand.values())
        picks = [q for q, v in cand.items() if v == best]
        cen = np.mean([list(p) for p in cluster], axis=0)
        picks.sort(key=lambda q: (np.sum((np.array(q) - cen) ** 2),
                                  rng.random()))
        cluster.add(picks[0])
        B += best
    return B


def void_energy(n):
    """Broken-bond formation energy [MeV] of the best compact n-void found."""
    B = max(grow_cluster(n, s) for s in range(6))
    return (6 * n - B) * (M0_MEV / 6.0)


def sft_energy(k, r0=0.5, c_core=1.0, gamma_isf=0.0):
    """Stacking-fault tetrahedron energy [MeV] for platelet edge k sites.

    n = k(k+1)/2 vacancies, edge length L = (k-1) ell.  Faces are
    colour-closed (self-screening), so the area term carries only the
    residual fault energy gamma_isf [m0/ell^2], near zero.  Edges are
    six stair-rods, b = ell/3, nu = 1/2, giving a line prefactor
    sqrt(2) mu ell^2 b^2-scaled = sqrt(2) m0 / (18 pi) per ell.
    """
    L = k - 1
    if L < 1:
        return None
    pre = np.sqrt(2.0) * M0_MEV / (18.0 * np.pi)
    line = 6.0 * L * pre * (np.log(L / r0) + c_core)
    fault = np.sqrt(3.0) * L * L * gamma_isf * M0_MEV
    return line + fault


def ignition_crossing():
    print("=== Void vs stacking-fault tetrahedron (ignition) ===")
    print(f"{'k':>3} {'n':>4} {'E_void':>8} {'E_SFT':>8} {'burst':>8} "
          f"{'burst frac':>10}")
    for k in range(2, 13):
        n = k * (k + 1) // 2
        ev = void_energy(n)
        es = sft_energy(k)
        print(f"{k:>3} {n:>4} {ev:8.1f} {es:8.1f} {ev - es:8.1f} "
              f"{(ev - es) / ev:10.1%}")
    # Crossing sweep: robustness against core cutoff, core constant,
    # and residual fault energy.
    ncs = set()
    for r0 in (0.25, 0.5, 1.0):
        for c in (0.0, 1.0, 2.0):
            for g in (0.0, 0.01, 0.05):
                for k in range(2, 13):
                    if sft_energy(k, r0, c, g) < void_energy(k * (k + 1) // 2):
                        ncs.add(k * (k + 1) // 2)
                        break
    print(f"thermodynamic crossing n_c^eq across sweep: {sorted(ncs)}")


# ----------------------------------------------------------------------
# 6. Trigger cross-section and population depletion
# ----------------------------------------------------------------------
def trigger_and_depletion():
    """Athermal trigger of void collapse by an impinging matter defect.

    Regime A (fault-cancelling path): threshold is line-tension-limited,
    tau_c ~ mu b_p/(2 pi R); a screw of Burgers vector ell delivers that
    stress out to r_c = sqrt(3) R, so sigma = 3 pi R^2.
    Regime B (open-fault path): tau_c ~ gamma_USF/b_p, core-on impact
    only, sigma = pi R^2.  Both are geometric.
    """
    v0 = ELL_CM ** 3 / np.sqrt(2.0)          # FCC volume per site [cm^3]
    barn = 1e-24                              # [cm^2]
    print("=== Trigger cross-section (athermal, both regimes) ===")
    for n in (6, 21, 55, 78):
        R = (3.0 * n * v0 / (4.0 * np.pi)) ** (1.0 / 3.0)
        sA = 3.0 * np.pi * R * R
        sB = np.pi * R * R
        print(f"n={n:3d}: R = {R/ELL_CM:4.2f} ell, sigma = "
              f"{sB/barn:4.2f}-{sA/barn:4.2f} barn")

    # Depletion timescales, sigma ~ 1 barn, v_rel ~ 250 km/s
    sigma, v = 1e-24, 2.5e7
    print("\nDepletion time 1/(n_H sigma v), sigma = 1 barn:")
    for name, nH in (("halo gas", 1e-4), ("warm ISM", 0.5),
                     ("atomic gas", 1.0), ("molecular cloud", 1e3),
                     ("CMZ", 1e4)):
        t = 1.0 / (nH * sigma * v)
        print(f"  {name:16s} n_H={nH:8.1e} cm^-3 -> {t/3.15e16:10.3g} Gyr")

    # Relic vacancy site fraction: clustering cannot happen by chance
    n_vac = 0.4 / M0_MEV * 1e3                # local DM number density [cm^-3]
    c_v = n_vac * v0
    print(f"\nvacancy site fraction c_v = {c_v:.1e} "
          f"(chance adjacency ~ 12 c_v = {12*c_v:.0e}: never)")

    # Cap on the surviving cluster mass fraction f_c from the observed
    # GC GeV excess (~1e37 erg/s ~ 6e39 GeV/s), r < 1.5 kpc, n_H ~ 10.
    L_obs = 6e39                              # [GeV/s]
    V = 4.0 / 3.0 * np.pi * (1.5 * 3.086e21) ** 3
    M_cl, E_burst = 0.7, 0.5                  # [GeV], n ~ 21 cluster
    for name, rho in (("halo-tracing (rho_DM ~ 1 GeV/cm^3)", 1.0),
                      ("comoving-uniform (mean rho_DM)", 1.4e-6)):
        rate_den = (rho / M_cl) * 10.0 * sigma * v      # per f_c [cm^-3 s^-1]
        f_cap = L_obs / (rate_den * E_burst * V)
        print(f"  f_c cap, {name}: {f_cap:.0e}")


if __name__ == "__main__":
    defecton_band()
    cluster_bindings()
    coalescence_luminosity()
    rayleigh_channel()
    ignition_crossing()
    trigger_and_depletion()
