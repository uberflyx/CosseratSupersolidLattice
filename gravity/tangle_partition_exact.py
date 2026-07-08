"""
Exact canonical partition function of the interior tangle, from first
principles, with brute-force verification and its two expansions.

Physics context (monograph, black-holes chapter; paper "Black holes as
superfluid droplets"): the tangle of n line defects carries one mode per
unordered pair, whose occupation number is that pair's winding class,
truncated to k values.  The Hamiltonian is harmonic in the total winding,
H = eps * sum_{i<j} n_ij with eps(r) = eps_0 * psi(r) the elastic winding
gap, and the relabelling of indistinguishable lines is a gauged S_N.
Projecting the thermal trace onto gauge invariants and evaluating it
cycle by cycle gives the closed form

    Z_N(x; k) = sum_{lambda |- N} (1/z_lambda)
                 prod_{cycles c} (1 - x^{k l_c}) / (1 - x^{l_c}),

with x = exp(-beta * eps) the Boltzmann weight per winding quantum and
l_c the pair-cycle lengths.  This is the k-truncated Molien-Weyl series
of O'Connor-Ramgoolam (JHEP 04 (2024) 080; JHEP 07 (2024) 152) for the
pair representation, now read as a physical thermal object.

Checks and outputs:
  1. BRUTE FORCE: for small N, k the closed form equals the explicit
     gauge-projected trace over all k^(N(N-1)/2) occupation patterns.
  2. LIMITS: Z(x -> 1) = Omega(N, k) (the microstate count);
     Z(x -> 0) = 1 (frozen record); low-temperature coefficients are the
     stabilised multigraph numbers 1, 1, 3, 8, 23, 66, ... for N >= 2g.
  3. SKIN BAND: with x(r) = exp(-2 ln N (xi/xi_*)^3) from the Tolman law
     and the winding gap, the per-pair entropy rises from 0 to ln k
     between the outer edge xi_* (first thermal winding) and the inner
     edge xi_* (2 ln N)^(-1/3) (full mixing): the earlier robustness
     bracket is the band structure of the exact Z.
"""

from fractions import Fraction
from itertools import product
from math import log, exp, lgamma

from tangle_matrix_model_matching import (partitions, sym_factor,
                                          pair_cycle_lengths, omega)


# ----------------------------------------------------------------------
# Closed form
# ----------------------------------------------------------------------

def Z_closed(N, k, x):
    """Exact gauge-projected canonical partition function (float x)."""
    total = 0.0
    for p in partitions(N):
        term = 1.0
        for l in pair_cycle_lengths(p):
            if abs(1.0 - x) < 1e-14:
                term *= k
            else:
                term *= (1.0 - x ** (k * l)) / (1.0 - x ** l)
        total += term / sym_factor(p)
    return total


# ----------------------------------------------------------------------
# Brute force: explicit trace over occupation patterns
# ----------------------------------------------------------------------

def Z_brute(N, k, x):
    """Tr[P exp(-beta H)] evaluated state by state: sum x^(total winding)
    over one representative per S_N orbit, or equivalently the projected
    trace sum_states x^W * |Orb|^{-1}... implemented as the plain orbit
    sum: (1/N!) * sum over ALL patterns of x^W * (number of sigma fixing
    the pattern) is the same as summing x^W over orbit representatives.
    We use the direct projected trace: (1/N!) sum_sigma Tr[U(sigma) x^W],
    with the sigma-trace evaluated pattern by pattern."""
    from itertools import permutations
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    idx = {pq: a for a, pq in enumerate(pairs)}
    d = len(pairs)
    total = 0.0
    nfact = 1
    for i in range(2, N + 1):
        nfact *= i
    for sigma in permutations(range(N)):
        # induced permutation on pair slots
        perm = [idx[tuple(sorted((sigma[i], sigma[j])))] for (i, j) in pairs]
        # trace of U(sigma) x^W: patterns constant on cycles of perm
        # decompose perm into cycles
        seen, cycles = [False] * d, []
        for a in range(d):
            if not seen[a]:
                l, b = 0, a
                while not seen[b]:
                    seen[b] = True
                    b = perm[b]
                    l += 1
                cycles.append(l)
        term = 1.0
        for l in cycles:
            term *= sum(x ** (m * l) for m in range(k))
        total += term
    return total / nfact


# ----------------------------------------------------------------------
# Thermodynamic profile through the skin (identity-dominant regime)
# ----------------------------------------------------------------------

def per_pair_entropy(x, k):
    """Entropy per pair mode, s/k_B, of one truncated oscillator at
    Boltzmann weight x (identity-dominant, valid for macroscopic N where
    the ln N! rebate per pair is negligible)."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return log(k)
    zs = (1 - x ** k) / (1 - x)
    # mean occupation of one mode
    n = x / (1 - x) - k * x ** k / (1 - x ** k)
    return log(zs) - n * log(x)


if __name__ == "__main__":
    print("Brute force vs closed form (max |rel. err.| over x grid):")
    for N, k in ((4, 2), (4, 3), (5, 2), (5, 3)):
        errs = []
        for x in (0.0, 0.1, 0.35, 0.7, 0.95, 1.0):
            zb, zc = Z_brute(N, k, x), Z_closed(N, k, x)
            errs.append(abs(zb - zc) / zc)
        status = "PASS" if max(errs) < 1e-12 else "FAIL"
        print(f"  N={N} k={k}: {max(errs):.2e}  {status}")
        assert max(errs) < 1e-12

    print("\nLimits: Z(1) = Omega, Z(0) = 1:")
    for N, k in ((5, 2), (6, 2), (6, 3)):
        z1, z0 = Z_closed(N, k, 1.0), Z_closed(N, k, 0.0)
        print(f"  N={N} k={k}: Z(1)={z1:.6g} Omega={omega(N,k)}  Z(0)={z0:.3g}")
        assert abs(z1 - omega(N, k)) < 1e-6 * omega(N, k) and abs(z0 - 1) < 1e-12

    print("\nSkin band (stellar hole, 2 ln N = 178, k = 2): per-pair entropy")
    lnN2 = 178.0
    print(f"  {'xi/xi_*':>8} {'x(r)':>12} {'s/ln k':>8}")
    for r in (2.0, 1.0, 0.5, 0.178, 0.05):
        x = exp(-lnN2 * r ** 3)
        s = per_pair_entropy(x, 2) / log(2)
        print(f"  {r:>8.3f} {x:>12.3e} {s:>8.4f}")
    print("  outer edge xi_*: first thermal winding; inner edge "
          f"{(1/lnN2)**(1/3):.3f} xi_*: full mixing")
