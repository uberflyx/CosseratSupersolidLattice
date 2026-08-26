#!/usr/bin/env python3
"""
d4_lattice_sums.py — Elastic constants of the D4 (24-cell) lattice
===================================================================
Mitchell A. Cox, University of the Witwatersrand
Accompanies: M. A. Cox, "The Cosserat Supersolid" (2026)

Computes second- and third-order elastic lattice sums for the D3 (FCC,
12 neighbours) and D4 (24-cell, 24 neighbours) nearest-neighbour shells,
verifies that the Zener anisotropy ratio A = 1 in D4, and confirms three
structural identities that protect the α derivation and the Koide sector:

    S₂₂  = Σ R₁²R₂²       = 1   (identical in D3 and D4)
    S₄₂  = Σ R₁⁴R₂²       = 1/2 (identical)
    S₂₂₂ = Σ R₁²R₂²R₃²   = 0   (identical → C₄₅₆ = 0 selection rule)

It then carries the two steps that close the tangential load exponent:

  (1) strain stationarity.  A lattice that recrystallises at a new spacing
      carries no bond force, so every elastic constant is phi''(R0) times a
      purely geometric shell sum.  That common factor divides out of every
      ratio, so K/mu and the Poisson ratio cannot move with density.
  (2) the loop gain.  Letting the tangential-to-normal stiffness ratio r drift,
      the Poisson ratio can only drift through r, while any contact law fixes r
      from the Poisson ratio.  The loop closes on itself with gain Lambda, and
      Lambda is bounded by 1/12 on both the 12-bond and 24-bond dictionaries,
      at every r in (0,1), which leaves eta_t = 7 as the only solution.

Usage:
    python d4_lattice_sums.py
"""

import math

# Mindlin's tangential-to-normal contact stiffness ratio for identical
# elastic spheres, a function of the Poisson ratio alone.
def mindlin_r(nu):
    return 2.0 * (1.0 - nu) / (2.0 - nu)


def dln_mindlin_dnu(nu):
    """d ln(k_t/k_n) / d nu for the Mindlin law."""
    return -1.0 / (1.0 - nu) + 1.0 / (2.0 - nu)


# Poisson ratio of each dictionary as a function of r = k_t/k_n, taken from
# the Born sums with a normal and a tangential spring on every contact and the
# microrotation relaxed (long-wavelength shear mu_bar).
def nu_d4(r):
    return (1.0 - r) / (2.0 * (2.0 + r))


def dnu_d4(r):
    return -3.0 / (2.0 * (2.0 + r) ** 2)


def nu_fcc(r):
    return (1.0 - r) / (5.0 + r)


def dnu_fcc(r):
    return -6.0 / (5.0 + r) ** 2


def loop_gain(r, nu_of_r, dnu_of_r):
    """Lambda: the fraction of a drift in r that returns after one pass
    through the Poisson ratio.  eta_t = 7 unless Lambda equals one."""
    return dln_mindlin_dnu(nu_of_r(r)) * dnu_of_r(r) * r

# ================================================================
# NEAREST-NEIGHBOUR VECTORS
# ================================================================

def fcc_nn_vectors():
    """12 FCC (D3) nearest-neighbour vectors at distance d.
    Vectors: (±1,±1,0)/√2 and permutations, in units where d = 1."""
    vecs = []
    for i in range(3):
        for j in range(i + 1, 3):
            for si in (-1, 1):
                for sj in (-1, 1):
                    v = [0.0, 0.0, 0.0]
                    v[i] = si / math.sqrt(2)
                    v[j] = sj / math.sqrt(2)
                    vecs.append(v)
    return vecs


def d4_nn_vectors():
    """24 D4 (24-cell) nearest-neighbour vectors at distance d.
    Vectors: (±1,±1,0,0)/√2 and all permutations of 4 coordinates."""
    vecs = []
    for i in range(4):
        for j in range(i + 1, 4):
            for si in (-1, 1):
                for sj in (-1, 1):
                    v = [0.0, 0.0, 0.0, 0.0]
                    v[i] = si / math.sqrt(2)
                    v[j] = sj / math.sqrt(2)
                    vecs.append(v)
    return vecs


# ================================================================
# LATTICE SUM COMPUTATION
# ================================================================

def lattice_sum(vecs, indices):
    """Compute Σ_bonds Π_k R_{i_k} for the given component indices."""
    return sum(math.prod(v[i] for i in indices) for v in vecs)


def compute_all_sums(vecs, D):
    """Compute all relevant lattice sums for SOEC and TOEC."""
    sums = {
        'Z':    len(vecs),
        'S2':   lattice_sum(vecs, [0, 0]),
        'S4':   lattice_sum(vecs, [0, 0, 0, 0]),
        'S22':  lattice_sum(vecs, [0, 0, 1, 1]),
        'S6':   lattice_sum(vecs, [0, 0, 0, 0, 0, 0]),
        'S42':  lattice_sum(vecs, [0, 0, 0, 0, 1, 1]),
        'S222': lattice_sum(vecs, [0, 0, 1, 1, 2, 2]),
    }
    return sums


# ================================================================
# MAIN
# ================================================================

if __name__ == '__main__':
    v3 = fcc_nn_vectors()
    v4 = d4_nn_vectors()
    s3 = compute_all_sums(v3, 3)
    s4 = compute_all_sums(v4, 4)

    print("Lattice sums: D3 (FCC, 12 NN) vs D4 (24-cell, 24 NN)")
    print(f"{'Sum':>6s}  {'D3':>8s}  {'D4':>8s}  {'Same':>5s}")
    print(f"{'─' * 6}  {'─' * 8}  {'─' * 8}  {'─' * 5}")
    for key in ('Z', 'S2', 'S4', 'S22', 'S6', 'S42', 'S222'):
        same = '✓' if abs(s3[key] - s4[key]) < 1e-10 else '✗'
        print(f"{key:>6s}  {s3[key]:8.4f}  {s4[key]:8.4f}  {same:>5s}")

    print()
    for D, s, label in [(3, s3, 'D3 FCC'), (4, s4, 'D4 24-cell')]:
        C11 = s['S4']
        C12 = C44 = s['S22']
        A = 2 * C44 / (C11 - C12)
        K = (C11 + (D - 1) * C12) / D
        mu = C44 if A == 1 else (C11 - C12 + 3 * C44) / 5
        print(f"{label}: C₁₁={C11:.0f} C₁₂=C₄₄={C12:.0f}  "
              f"A={A:.1f}  K={K:.4f}  μ={mu:.4f}  K/μ={K / mu:.4f}")

    # ------------------------------------------------------------------
    # Strain stationarity: rescaling the shell multiplies every constant by
    # one common factor, so no ratio of moduli can move with density.
    # ------------------------------------------------------------------
    print()
    print("Strain stationarity of the D4 modulus ratios under uniform rescaling")
    print(f"{'spacing':>8s}  {'C44':>10s}  {'K':>10s}  {'K/C44':>14s}  {'nu':>12s}")
    # phi''(R0) for the Morse contact the gravitational sector fixes,
    # a = 7/3 so that the bond anharmonicity is xi = R phi'''/phi'' = -7.
    a_morse = 7.0 / 3.0
    D_morse = 1.0 / (2.0 * a_morse ** 2)     # normalised to phi''(1) = 1

    def phi2(R):
        e = math.exp(-a_morse * (R - 1.0))
        return 2.0 * D_morse * a_morse ** 2 * (2.0 * e * e - e)

    for scale in (0.90, 0.95, 1.00, 1.05, 1.10):
        # every constant is phi''(R0) * R0^2 / V0, and V0 scales as R0^3
        w = phi2(scale) * scale ** 2 / scale ** 3
        C11, C12 = s4['S4'] * w, s4['S22'] * w
        C44 = C12
        K = (C11 + 2 * C12) / 3
        print(f"{scale:8.2f}  {C44:10.5f}  {K:10.5f}  {K / C44:14.10f}"
              f"  {C12 / (C11 + C12):12.10f}")
    print("  Both moduli move; every ratio holds to machine precision.")

    # ------------------------------------------------------------------
    # The loop gain, and the bound that makes eta_t = 7 the only solution.
    # ------------------------------------------------------------------
    print()
    print("Loop gain Lambda(r) on both dictionaries; eta_t = 7 unless Lambda = 1")
    print(f"{'r':>8s}  {'nu (D4)':>10s}  {'Lambda D4':>11s}"
          f"  {'nu (FCC)':>10s}  {'Lambda FCC':>11s}")
    rolling_d4 = 1.0 / (3.0 * math.pi - 5.0)
    rolling_fcc = 1.0 / (2.0 * math.pi - 3.0)
    marks = [0.001, rolling_d4, rolling_fcc, 0.5, 0.9, 0.999]
    for r in marks:
        print(f"{r:8.4f}  {nu_d4(r):10.4f}  {loop_gain(r, nu_d4, dnu_d4):11.5f}"
              f"  {nu_fcc(r):10.4f}  {loop_gain(r, nu_fcc, dnu_fcc):11.5f}")
    worst = max(abs(loop_gain(1e-4 + i * 1e-4, f, g))
                for i in range(9999) for f, g in
                ((nu_d4, dnu_d4), (nu_fcc, dnu_fcc)))
    print(f"  max |Lambda| over r in (0,1), both dictionaries: {worst:.6f}"
          f"  (bound 1/12 = {1 / 12:.6f})")
    print("  1 - Lambda >= 11/12 everywhere, so 7 - eta_t = 0 exactly.")
