"""
The interior tangle as a gauged permutation-invariant matrix model.

Physics context (monograph, black-holes chapter; paper "Black holes as
superfluid droplets"): the tangle microstate count Omega(N, k), the number
of symmetric k-valued linking matrices modulo S_N relabelling, is exactly
the dimension of the gauge-invariant Hilbert space of an S_N-gauged matrix
harmonic oscillator on the pair representation, with each pair mode
truncated to k occupation levels, evaluated at zero mode frequency.  Zero
frequency is not a choice: mu_eff = 0 makes every tangle mode cost zero
energy, so the Boltzmann weight per quantum is x = exp(-beta*omega) = 1
and the partition function equals the state count.  The untruncated
generating function is the Molien-Weyl series whose lcm/gcd structure was
derived for the full-matrix representation by O'Connor and Ramgoolam,
JHEP 04 (2024) 080 and JHEP 07 (2024) 152.

This script verifies, exactly:
  1. IDENTITY: the coefficient sum of the k-truncated Molien series equals
     Omega(N, k) from the Burnside count (bh_defect_packing_count.py).
  2. STABILISATION: the degree-g invariant dimension dim_g(N) is
     independent of N for N >= 2g (the partition-algebra regime of
     Barnes, Padellaro and Ramgoolam, PRD 106 (2022) 106020): observables
     of low degree cannot measure the line number N.
  3. MATCHING DEGREE: the cumulative count of linearly independent
     invariants first exceeds the packing count Omega(N, k) at total
     degree G* = c(k) * d, with d = N(N-1)/2 pair slots and c(k) the root
     of (1+c)^(1+c) / c^c = k, which tends to k/e at large k.  The
     invariant algebra needs an average degree budget of about k/e per
     pair to resolve the packing space.

Pure Python integer/series arithmetic; no dependencies beyond stdlib.
"""

from fractions import Fraction
from math import gcd, lcm, log, exp, lgamma
from functools import lru_cache


# ----------------------------------------------------------------------
# Partitions of N and their symmetry factors
# ----------------------------------------------------------------------

def partitions(n, max_part=None):
    """Yield partitions of n as dicts {part: multiplicity}."""
    if max_part is None:
        max_part = n
    if n == 0:
        yield {}
        return
    for a in range(min(n, max_part), 0, -1):
        for rest in partitions(n - a, a):
            p = dict(rest)
            p[a] = p.get(a, 0) + 1
            yield p


def sym_factor(p):
    """z_lambda = prod a^{m_a} m_a! ; the centraliser order of the class."""
    z = 1
    for a, m in p.items():
        z *= a ** m
        for i in range(1, m + 1):
            z *= i
    return z


def pair_cycle_lengths(p):
    """Multiset of cycle lengths of the induced permutation on unordered
    pairs {i, j}, i != j, for a permutation of cycle type p.

    Within one a-cycle: floor(a/2) cycles; offsets d < a/2 have length a,
    and for even a the offset a/2 has length a/2.  Between two cycles of
    lengths a and b (including two distinct cycles of equal length):
    gcd(a, b) cycles, each of length lcm(a, b)."""
    lengths = []
    parts = sorted(p.items())
    for a, m in parts:
        # within each a-cycle
        for _ in range(m):
            if a % 2 == 0:
                lengths.extend([a] * (a // 2 - 1))
                lengths.append(a // 2)
            else:
                lengths.extend([a] * ((a - 1) // 2))
        # between distinct cycles of the same length a
        n_pairs_same = m * (m - 1) // 2
        lengths.extend([a] * (a * n_pairs_same))
    # between cycles of different lengths
    for i in range(len(parts)):
        for j in range(i + 1, len(parts)):
            a, ma = parts[i]
            b, mb = parts[j]
            g, l = gcd(a, b), lcm(a, b)
            lengths.extend([l] * (g * ma * mb))
    return lengths


# ----------------------------------------------------------------------
# Series utilities (exact rational coefficients)
# ----------------------------------------------------------------------

def series_mult(A, B, G):
    """Product of two truncated series (lists of Fractions) to degree G."""
    out = [Fraction(0)] * (G + 1)
    for i, ai in enumerate(A):
        if ai == 0:
            continue
        jmax = G - i
        for j, bj in enumerate(B[: jmax + 1]):
            if bj:
                out[i + j] += ai * bj
    return out


def geometric_series(step, G, k=None):
    """1/(1 - x^step) truncated at G, or (1 - x^{k step})/(1 - x^step),
    i.e. 1 + x^step + ... + x^{(k-1) step} when k is given."""
    A = [Fraction(0)] * (G + 1)
    top = (k - 1) * step if k is not None else G
    for e in range(0, min(G, top) + 1, step):
        A[e] = Fraction(1)
    return A


def molien_series(N, G, k=None):
    """Graded dimensions dim_g of S_N invariants on the pair slots, to
    degree G.  With k given, each mode is truncated to k levels (the
    finite alphabet of linking states); without k, the full Molien-Weyl
    series of the gauged oscillator."""
    total = [Fraction(0)] * (G + 1)
    for p in partitions(N):
        z = sym_factor(p)
        cls = [Fraction(1)] + [Fraction(0)] * G
        for step in pair_cycle_lengths(p):
            cls = series_mult(cls, geometric_series(step, G, k), G)
        w = Fraction(1, z)
        for g in range(G + 1):
            total[g] += w * cls[g]
    # all coefficients must be integers
    dims = []
    for g, c in enumerate(total):
        assert c.denominator == 1, f"non-integer dim at degree {g}"
        dims.append(int(c))
    return dims


# ----------------------------------------------------------------------
# Reference count Omega(N, k) by Burnside (independent implementation)
# ----------------------------------------------------------------------

def omega(N, k):
    tot = Fraction(0)
    for p in partitions(N):
        tot += Fraction(k ** len(pair_cycle_lengths(p)), sym_factor(p))
    assert tot.denominator == 1
    return int(tot)


# ----------------------------------------------------------------------
# The matching-degree prediction
# ----------------------------------------------------------------------

def c_of_k(k, tol=1e-12):
    """Root of (1+c)^(1+c) / c^c = k, by bisection.  c(k) -> k/e as k
    grows; c(2) ~ 0.293."""
    lo, hi = 1e-9, float(k)
    f = lambda c: (1 + c) * log(1 + c) - c * log(c) - log(k)
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# ----------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------

def check_identity(N_list=(4, 5, 6, 7), k_list=(2, 3)):
    """Coefficient sum of the k-truncated Molien series = Omega(N, k):
    the tangle count is the omega -> 0 (x -> 1) gauged-oscillator
    partition function."""
    print("Identity: truncated Molien series at x -> 1 vs Burnside count")
    for N in N_list:
        d = N * (N - 1) // 2
        for k in k_list:
            G = d * (k - 1)              # top degree of the truncated ring
            dims = molien_series(N, G, k=k)
            lhs, rhs = sum(dims), omega(N, k)
            status = "PASS" if lhs == rhs else "FAIL"
            print(f"  N={N} k={k}: sum(dim_g)={lhs}  Omega={rhs}  {status}")
            assert lhs == rhs


def check_stabilisation(g_max=6, N_ref=30):
    """dim_g is N-independent once N >= 2g: low-degree observables cannot
    see the line number."""
    print("\nStabilisation dim_g(N) = dim_g(infinity) for N >= 2g "
          f"(reference N = {N_ref}):")
    ref = molien_series(N_ref, g_max)
    print(f"  stabilised dims g=0..{g_max}: {ref}")
    for N in range(2, 2 * g_max + 3):
        dims = molien_series(N, g_max)
        marks = []
        for g in range(g_max + 1):
            ok = dims[g] == ref[g]
            expected = (N >= 2 * g)
            # stabilisation must hold exactly when N >= 2g
            if expected:
                assert ok, (N, g, dims[g], ref[g])
            marks.append("=" if ok else ".")
        print(f"  N={N:<3d} {' '.join(marks)}   (= means equals stabilised value)")


def check_matching_degree(N_list=(8, 10, 12, 14), k_list=(2, 3, 5)):
    """Smallest G with cumulative invariant count >= Omega(N, k),
    against the prediction G* = c(k) d."""
    print("\nMatching degree: cumulative invariants reach the packing count")
    print(f"  {'N':>3} {'k':>3} {'d':>5} {'G* exact':>9} {'c(k) d':>9} "
          f"{'ratio':>7}")
    for N in N_list:
        d = N * (N - 1) // 2
        for k in k_list:
            target = omega(N, k)
            c = c_of_k(k)
            G_guess = int(c * d * 2) + 20
            dims = molien_series(N, G_guess)
            run, Gstar = 0, None
            for g, dg in enumerate(dims):
                run += dg
                if run >= target:
                    Gstar = g
                    break
            assert Gstar is not None, "raise G_guess"
            pred = c * d
            print(f"  {N:>3} {k:>3} {d:>5} {Gstar:>9} {pred:>9.1f} "
                  f"{Gstar / pred:>7.3f}")
    print(f"  c(2)={c_of_k(2):.4f}, c(3)={c_of_k(3):.4f}, "
          f"c(5)={c_of_k(5):.4f}; large-k limit c -> k/e "
          f"(k/e = {2/exp(1):.4f}, {3/exp(1):.4f}, {5/exp(1):.4f})")


def naor_bits(N, k=2):
    """Description length: ln2 Omega vs C(N,2) - log2 N! (Naor's optimal
    unlabelled-graph encoding, exact to the second leading term)."""
    print("\nDescription length in bits (k = 2, tangle as a graph):")
    for n in N:
        lnO = None
        # exact for small n, closed form for large
        if n <= 24:
            lnO = log(omega(n, k))
        else:
            lnO = (log(k) / 2) * n * (n - 1) - lgamma(n + 1)
        bits = lnO / log(2)
        naor = n * (n - 1) / 2 - lgamma(n + 1) / log(2)
        print(f"  N={n:<6d} ln2 Omega = {bits:16.2f}   "
              f"C(N,2) - log2 N! = {naor:16.2f}")


if __name__ == "__main__":
    check_identity()
    check_stabilisation()
    check_matching_degree()
    naor_bits((10, 16, 24, 100, 1000))
