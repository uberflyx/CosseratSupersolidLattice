#!/usr/bin/env python3
"""
Vacancy dark matter: defecton band, cluster bond counts, and
detection-rate estimates for the Cosserat supersolid lattice.

The lattice is D_4, whose kissing shell has Z = 24 neighbours; the FCC
shell of 12 is its spatial slice.  Every bond count here is taken on the
D_4 contact graph, and the bond quantum is eps = m0 c^2 / 12, fixed by
requiring that a single vacancy (24 bonds severed, 12 recovered on
re-deposition at a remote perfect site) cost exactly the node mass m0.

What this script computes
-------------------------
1. The defecton tight-binding band on the 24-neighbour shell from the hop
   amplitude t = alpha * m0 * c^2 (one electron mass): band bottom
   E_f - 24t, band top E_f + 8t, and full width 32t.  The band bottom is
   the particle's rest energy.  Its INERTIA is not read off this band: an
   equal-time hop to one of the twelve crossing bonds is a step in
   imaginary time, not a step in space, so the worldline register settles
   the dispersion instead (see defecton_worldline.py).
2. Broken-bond binding energies for small vacancy clusters (di-vacancy,
   triangle, tetrahedron) by direct enumeration on the D_4 coordination
   graph, confirming E_b = B(G) * eps and the growth-step ladder m * eps.
   Cross-checks d4_conversion.py, which enumerates B_max(n) independently.
3. An upper-bound estimate of the present-day diffuse luminosity of the
   5.835 MeV coalescence line from a Milky-Way-like halo, using the
   derived elastic transport cross-section as a generous proxy for
   radiative capture.
4. The peak of the Klemens (Rayleigh, omega^4) photon-scattering
   cross-section against the zero-point smearing form factor, and the
   resulting optical depth through the Galactic Centre dark matter column.
5. The void-versus-tetrahedron energetics behind the ignition channel,
   with the void ledger given on both graphs.  The tetrahedron side is
   still a slice construction, so the burst column is the slice pair and
   is labelled as such; see the note in ignition_crossing().

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
# 1. Defecton band on the D_4 kissing shell
# ----------------------------------------------------------------------
def defecton_band():
    """Rest energy from the band bottom; inertia deliberately not taken here.

    The band Hamiltonian offers the vacancy 24 equal-time destinations, but
    twelve of them are crossing bonds with a leg along the compact axis.  In
    the Matsubara reading that axis is imaginary time, so an equal-time hop to
    a different time slice is not an available move, and pricing those twelve
    amplitudes as spatial mobility makes the inertia far too heavy.  The band
    bottom is still the rest energy, because it is a static shift.  The
    dispersion belongs to the worldline register (defecton_worldline.py).
    """
    t = ALPHA * M0_MEV               # hop amplitude [MeV] = m_e exactly
    e_f = M0_MEV                     # formation energy [MeV]

    # Structure factor S(k) = sum over the 24 neighbours of cos(k.delta) runs
    # from S = 24 at Gamma down to S = -8, so the band spans [E_f - 24t,
    # E_f + 8t] with full width 32t, all of it reachable at k_4 = 0.
    e_bottom = e_f - 24.0 * t
    e_top = e_f + 8.0 * t
    width = e_top - e_bottom         # = 32 t

    print("=== Defecton band on the D_4 shell (rest energy only) ===")
    print(f"hop amplitude t          : {t:.4f} MeV (= m_e)")
    print(f"band bottom E_f - 24t    : {e_bottom:.2f} MeV  "
          f"(shift {24*t:.2f} MeV = {24*ALPHA*100:.1f}% of E_f)")
    print(f"band top    E_f +  8t    : {e_top:.2f} MeV")
    print(f"band width 32t           : {width:.2f} MeV")
    print("inertia                  : not from this band; the worldline "
          "register gives m*/m_g = 1 - (m_g/m0)^4/720")
    print("                           (see defecton_worldline.py)")
    print()
    return t


# ----------------------------------------------------------------------
# 2. Cluster bond counts on the D_4 coordination graph
# ----------------------------------------------------------------------
def d4_sites(reach=2):
    """Integer D_4 sites (coordinate sum even) within a small box."""
    r = range(-reach, reach + 1)
    return [p for p in itertools.product(r, repeat=4) if sum(p) % 2 == 0]


def internal_bonds(site_set):
    """Number of nearest-neighbour bonds internal to a site set.

    Integer D_4 convention: nearest neighbours sit at squared distance 2,
    the 24 vectors (+-1, +-1, 0, 0) and their permutations.  The 12 with a
    zero fourth component are the FCC slice shell; the other 12 are the
    crossing bonds, each with one leg along the compact axis.
    """
    sites = list(site_set)
    count = 0
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            d = sum((a - b) ** 2 for a, b in zip(sites[i], sites[j]))
            if d == 2:
                count += 1
    return count


def bond_quantum():
    """eps from the single-vacancy ledger: 24 severed, 12 recovered."""
    return M0_MEV / 12.0


def cluster_bindings():
    eps = bond_quantum()             # bond price [MeV]
    print("=== Cluster bond counts, E_b = B(G) * eps ===")
    print(f"bond price eps = m0 c^2 / 12 = {eps:.3f} MeV "
          f"(= m_e/(12 alpha) = {ME_MEV/(12*ALPHA):.3f} MeV)")
    print(f"single-vacancy check: 12 eps = {12*eps:.3f} MeV "
          f"= m0 = {M0_MEV:.3f} MeV")

    # Compact clusters, given as integer D_4 coordinates.  All four lie in
    # the slice (fourth coordinate zero); the tetrahedron is the largest
    # clique of the D_4 graph, as it is of the slice's.
    clusters = {
        "di-vacancy (pair)": [(0, 0, 0, 0), (1, 1, 0, 0)],
        "tri-vacancy (triangle)": [(0, 0, 0, 0), (1, 1, 0, 0), (1, 0, 1, 0)],
        "tetra-vacancy (tetrahedron)": [(0, 0, 0, 0), (1, 1, 0, 0),
                                        (1, 0, 1, 0), (0, 1, 1, 0)],
    }
    prev_bonds = 0
    for name, sites in clusters.items():
        b = internal_bonds(sites)
        step = b - prev_bonds
        print(f"{name:32s}: B(G) = {b}, E_b = {b*eps:.2f} MeV, "
              f"growth step m = {step} releases {step*eps:.2f} MeV")
        prev_bonds = b
    # Verify the tetrahedron is a maximal clique (no K5 in the D_4 graph).
    tetra = [(0, 0, 0, 0), (1, 1, 0, 0), (1, 0, 1, 0), (0, 1, 1, 0)]
    common = [s for s in d4_sites(2)
              if s not in tetra
              and all(sum((a - b) ** 2 for a, b in zip(s, t)) == 2
                      for t in tetra)]
    print(f"sites adjacent to all four tetrahedron members: {len(common)} "
          f"(no K5 in the D_4 graph)")
    print(f"tetrahedron binding 6 eps = {6*eps:.2f} MeV = m0/2; bound "
          f"tetra-vacancy weighs 4m0 - m0/2 = {3.5*M0_MEV:.1f} MeV")
    print("B_max(5) = 9 and B_max(6) = 13 are enumerated exhaustively in "
          "d4_conversion.py")
    print()


# ----------------------------------------------------------------------
# 3. Coalescence-line luminosity (upper bound)
# ----------------------------------------------------------------------
def coalescence_luminosity():
    eps = bond_quantum()                     # line energy [MeV]
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
SHELL = {3: [d for d in itertools.product((-1, 0, 1), repeat=3)
             if sum(x * x for x in d) == 2],
         4: [d for d in itertools.product((-1, 0, 1), repeat=4)
             if sum(x * x for x in d) == 2 and sum(d) % 2 == 0]}


def grow_cluster(n, dim, seed=0):
    """Greedy max-bond compact cluster; returns the site set.

    At each step the site adjacent to the most current members is added,
    ties broken toward the centroid.  Greedy alone leaves a bond or two on
    the table at larger n, so anneal_cluster() polishes the result.
    """
    rng = np.random.default_rng(seed)
    dirs = SHELL[dim]
    cluster = {tuple([0] * dim)}
    while len(cluster) < n:
        cand = {}
        for p in cluster:
            for d in dirs:
                q = tuple(p[i] + d[i] for i in range(dim))
                if q not in cluster:
                    cand[q] = cand.get(q, 0) + 1
        best = max(cand.values())
        picks = [q for q, v in cand.items() if v == best]
        cen = np.mean([list(p) for p in cluster], axis=0)
        picks.sort(key=lambda q: (np.sum((np.array(q) - cen) ** 2),
                                  rng.random()))
        cluster.add(picks[0])
    return cluster


def _degree(q, cluster, dirs, dim):
    return sum(1 for d in dirs
               if tuple(q[i] + d[i] for i in range(dim)) in cluster)


def anneal_cluster(cluster, dim, iters, rng):
    """Move the least connected site to the best free site, keeping gains."""
    dirs = SHELL[dim]
    B = sum(_degree(p, cluster, dirs, dim) for p in cluster) // 2
    best = B
    for _ in range(iters):
        degs = {p: _degree(p, cluster, dirs, dim) for p in cluster}
        lo = min(degs.values())
        movers = [p for p, v in degs.items() if v <= lo + 1]
        p = movers[rng.integers(len(movers))]
        trial = set(cluster)
        trial.discard(p)
        cand = {}
        for site in trial:
            for d in dirs:
                q = tuple(site[i] + d[i] for i in range(dim))
                if q not in trial:
                    cand[q] = cand.get(q, 0) + 1
        gain = max(cand.values())
        picks = [q for q, v in cand.items() if v == gain]
        trial.add(picks[rng.integers(len(picks))])
        Bt = B - degs[p] + gain
        if Bt >= B:                      # accept sideways moves to escape
            cluster, B = trial, Bt       # shallow local optima
            best = max(best, B)
    return best


def best_bonds(n, dim, seeds=6, iters=150):
    """Best internal bond count found for a compact n-site cluster."""
    return max(anneal_cluster(grow_cluster(n, dim, s), dim, iters,
                              np.random.default_rng(1000 + s))
               for s in range(seeds))


def void_energy(n, dim=4):
    """Broken-bond formation energy [MeV] of the best compact n-void found.

    Removing n nodes severs Z*n - B distinct bonds (each of the B internal
    bonds joins two removed nodes and would otherwise be counted from both
    ends); re-depositing them at remote perfect sites recovers (Z/2)*n.  The
    net cost is (Z/2 * n - B) * eps, which at n = 1 returns (Z/2)*eps = m0
    on either graph.
    """
    half_Z = {3: 6, 4: 12}[dim]
    eps = M0_MEV / half_Z
    return (half_Z * n - best_bonds(n, dim)) * eps


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
    """Void versus tetrahedron, with the void ledger given on both graphs.

    OPEN, and the reason is worth stating plainly.  The void is a set of n
    missing nodes, and n is a D_4 count.  The tetrahedron below is still a
    slice construction: n = k(k+1)/2 is a triangular platelet on ONE
    stacking layer, and its six stair-rod edges are priced per unit length
    in three dimensions.  A platelet wrapping all three layers of the
    compact ring is a different number of missing nodes, so until that
    extent is settled the two curves are not indexed by the same n and
    their difference is not the burst.  The burst column below is therefore
    the slice pair throughout, which is what the monograph quotes; the D_4
    void column is printed beside it to show the size of the pending shift.
    """
    print("=== Void vs stacking-fault tetrahedron (ignition) ===")
    print("burst column = slice pair (as quoted); E_void(D4) = pending recount")
    print(f"{'k':>3} {'n':>4} {'Evoid(sl)':>10} {'Evoid(D4)':>10} "
          f"{'ratio':>6} {'E_SFT':>8} {'burst(sl)':>10} {'frac':>7}")
    for k in range(2, 13):
        n = k * (k + 1) // 2
        ev3 = void_energy(n, dim=3)
        ev4 = void_energy(n, dim=4)
        es = sft_energy(k)
        print(f"{k:>3} {n:>4} {ev3:10.1f} {ev4:10.1f} {ev4/ev3:6.3f} "
              f"{es:8.1f} {ev3 - es:10.1f} {(ev3 - es) / ev3:7.1%}")
    # Crossing sweep: robustness against core cutoff, core constant, and
    # residual fault energy.  Run on the slice pair for consistency with the
    # burst column; the D_4 void is dearer, so the crossing only moves down.
    ncs = set()
    for r0 in (0.25, 0.5, 1.0):
        for c in (0.0, 1.0, 2.0):
            for g in (0.0, 0.01, 0.05):
                for k in range(2, 13):
                    n = k * (k + 1) // 2
                    if sft_energy(k, r0, c, g) < void_energy(n, dim=3):
                        ncs.add(n)
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
    # Slice volume per site.  The trigger is a spatial cross-section, so the
    # void's spatial footprint is what matters; if the collapsing cluster
    # turns out to span three stacking layers, R_void is set by its slice
    # projection rather than by n, which is the same open question the
    # ignition table carries.
    v0 = ELL_CM ** 3 / np.sqrt(2.0)          # slice volume per site [cm^3]
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
    # Cluster mass and burst for n ~ 21, on the slice pair (the cluster mass
    # is just the void energy, n*m0 - B*eps).  Both shift once the platelet
    # extent is settled; the cap is an order-of-magnitude statement and
    # survives the shift.
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
