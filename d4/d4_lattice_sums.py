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

Usage:
    python d4_lattice_sums.py
"""

import math

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
