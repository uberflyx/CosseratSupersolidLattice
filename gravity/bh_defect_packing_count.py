"""
Exact count of the black-hole interior tangle degeneracy.

Physics context (monograph, black-holes chapter; paper "Black holes as
superfluid droplets"): the shear-free interior keeps its defects as
zero-energy vortex lines, and the pairwise record of a tangle of n lines
is the linking matrix, a symmetric integer matrix modulo relabelling of
the (indistinguishable) lines.  This script enumerates those microstates
exactly and checks the closed form

    ln Omega = (ln k / 2) n(n-1) - ln n!  + O(n^2 k^-(n-2)),

with k the number of linking states per pair, against three anchors:
  1. the unlabelled-graph sequence (k = 2 must reproduce Polya counting);
  2. the exact finite-size formula above, to machine precision;
  3. the invariant-theory counts (Krull dimension, Hironaka secondary
     invariants, Procesi's matrix-tuple dimension), which supply the same
     n^2 exponent from an independent direction.
It closes with the Bekenstein-Hawking mapping n = sqrt(8 pi / ln k) M/M_P.

Only pairwise line topology produces a quadratic exponent: point-defect
positions give n ln(N/n) and internal states give n ln q, both linear.
"""


from math import gcd, comb, lgamma, log, pi, sqrt
from sympy.utilities.iterables import partitions

# ----------------------------------------------------------------------
# Burnside orbit count over cycle types of S_n
# ----------------------------------------------------------------------

def z_lambda(lam):
    """Size of the centraliser z_lambda = prod a^{m_a} m_a! for cycle type lam
    (dict {part_length: multiplicity}).  n!/z_lambda permutations share this type."""
    z = 1
    for a, m in lam.items():
        z *= a ** m
        for i in range(1, m + 1):
            z *= i
    return z


def pair_cycles(lam):
    """Number of cycles of the induced permutation on unordered pairs {i,j}, i<j.
    Standard graph-counting decomposition (Harary & Palmer):
      - pairs inside one cycle of length a: floor(a/2) cycles;
      - pairs across two distinct cycles of lengths a, b: gcd(a, b) cycles."""
    parts = list(lam.items())
    total = 0
    for a, m in parts:
        total += m * (a // 2)          # within each cycle of length a
        total += comb(m, 2) * a        # between two cycles of equal length a
    for i in range(len(parts)):
        for j in range(i + 1, len(parts)):
            a, ma = parts[i]
            b, mb = parts[j]
            total += ma * mb * gcd(a, b)
    return total


def diag_cycles(lam):
    """Cycles of the action on diagonal entries = number of cycles of the permutation."""
    return sum(lam.values())


def count_orbits_offdiag(n, states):
    """Exact number of symmetric off-diagonal matrices (n x n, zero diagonal)
    with entries from a set of `states` values, modulo S_n relabelling.
    states=2 reproduces the number of simple graphs on n unlabelled vertices."""
    total = 0
    nfact = 1
    for i in range(2, n + 1):
        nfact *= i
    for lam in partitions(n):
        total += (nfact // z_lambda(lam)) * states ** pair_cycles(lam)
    assert total % nfact == 0
    return total // nfact


def count_orbits_full(n, states):
    """As above but including the diagonal (self-linking / framing), so the
    label is the full symmetric matrix with n(n+1)/2 entries."""
    total = 0
    nfact = 1
    for i in range(2, n + 1):
        nfact *= i
    for lam in partitions(n):
        total += (nfact // z_lambda(lam)) * states ** (pair_cycles(lam) + diag_cycles(lam))
    assert total % nfact == 0
    return total // nfact


# ----------------------------------------------------------------------
# Checks and asymptotics
# ----------------------------------------------------------------------

def check_a000088():
    """Unlabelled simple graphs on n vertices, OEIS A000088."""
    known = [1, 1, 2, 4, 11, 34, 156, 1044, 12346, 274668, 12005168]
    got = [count_orbits_offdiag(n, 2) for n in range(len(known))]
    ok = got == known
    print("A000088 check (unlabelled graphs = binary linking patterns):", "PASS" if ok else "FAIL")
    print("  n     :", list(range(len(known))))
    print("  Omega :", got)
    return ok


def asymptotic_table(states_list=(2, 3, 5), n_max=40):
    """ln Omega / n^2 should converge to ln(states)/2: the n! relabelling and
    the diagonal are subleading (n ln n and n respectively)."""
    print("\nConvergence of ln(Omega)/n^2 toward ln(states)/2 (off-diagonal count):")
    header = "  n    " + "".join(f"k={s:<12d}" for s in states_list)
    print(header)
    for n in list(range(4, n_max + 1, 4)):
        row = f"  {n:<4d}"
        for s in states_list:
            om = count_orbits_offdiag(n, s)
            # om is an exact (huge) integer; take its log via string length guard
            row += f" {log_of_int(om) / n**2:<12.5f}"
        print(row)
    print("  limit" + "".join(f" {log(s)/2:<12.5f}" for s in states_list))


def log_of_int(m):
    """Natural log of a (possibly enormous) exact integer."""
    if m < 10**300:
        return log(m)
    s = str(m)
    return log(float("0." + s[:15])) + len(s) * log(10.0)


def point_defect_comparison(n_max=40):
    """Point defects: n defects on N = (2n)^3 sites (droplet comfortably larger
    than the defect count).  ln C(N, n) grows like n ln N, i.e. linearly in n
    up to logs; per n^2 it dies."""
    print("\nPoint-defect packing for contrast, ln C(N, n) / n^2 with N = (2n)^3:")
    for n in range(4, n_max + 1, 8):
        N = (2 * n) ** 3
        lnC = lgamma(N + 1) - lgamma(n + 1) - lgamma(N - n + 1)
        print(f"  n = {n:<4d}  ln C/n = {lnC/n:8.3f}   ln C/n^2 = {lnC/n**2:8.5f}")


def finite_size_formula_check(n_list=(10, 20, 30, 40), k_list=(2, 3, 5)):
    """The two-term formula ln Omega = (ln k/2) n(n-1) - ln n! is exact to
    machine precision once the transposition class k^{-(n-2)} dies."""
    print("\nExact vs closed form, ln Omega - [(ln k/2)n(n-1) - ln n!]:")
    for n in n_list:
        row = f"  n = {n:<4d}"
        for k in k_list:
            ex = log_of_int(count_orbits_offdiag(n, k))
            pred = (log(k) / 2) * n * (n - 1) - lgamma(n + 1)
            row += f"  k={k}: {ex - pred:11.3e}"
        print(row)


def invariant_theory_checks(n_list=(4, 8, 16, 32, 64)):
    """Question 2: independent invariants and secondary invariants.

    (a) Krull dimension of C[Sym_n(C)]^{S_n} = n(n+1)/2  (the action permutes
        coordinates, is finite, so the invariant ring has full dimension).
        This is the choice-free count of independent invariant labels.
    (b) Hironaka decomposition with canonical primaries e_1..e_d in the
        d = n(n+1)/2 matrix entries (degrees 1..d): number of secondary
        invariants  m_sec = d! / n!.
    (c) Procesi/Razmyslov: m-tuples of n x n matrices mod simultaneous GL_n
        conjugation form a variety of dimension (m-1) n^2 + 1  (m >= 2).
    """
    print("\nInvariant-theory counts against the microstate exponent:")
    print("  n     d=n(n+1)/2   ln(m_sec)=ln(d!/n!)   ln(m_sec)/n^2    ln Omega(k=1)/n^2  Procesi dim (m=2)")
    for n in n_list:
        d = n * (n + 1) // 2
        ln_msec = lgamma(d + 1) - lgamma(n + 1)
        om = count_orbits_full(n, 2)
        print(f"  {n:<5d} {d:<12d} {ln_msec:<21.2f} {ln_msec/n**2:<16.4f} {log_of_int(om)/n**2:<18.5f} {n*n + 1}")
    print("  ln(m_sec)/n^2 grows like (1/2) ln(n^2/2) - i.e. the SAME n^2 skeleton,")
    print("  with the per-pair state count promoted from a constant to ~ d/e.")
    print("  Setting Omega = m_sec:  d ln(k) = ln(d!)  =>  k = d/e  (Stirling), so the")
    print("  secondary count equals the packing count when each pair explores ~ n^2/(2e) linking values.")


def black_hole_mapping():
    """Map S = (ln 2 / 2) n^2 (binary linking) onto S_BH = 4 pi (M/M_P)^2 k_B."""
    print("\nBlack-hole mapping (Schwarzschild):")
    coeff = log(2.0) / 2.0
    n_per_planck = sqrt(4.0 * pi / coeff)
    print(f"  S_BH/k_B = 4 pi (M/M_P)^2;  S_tangle/k_B = (ln2/2) n^2")
    print(f"  Equal counts  =>  n = sqrt(8 pi / ln 2) (M/M_P) = {n_per_planck:.3f} (M/M_P)")
    M_P = 2.176434e-8       # kg, CODATA Planck mass
    M_sun = 1.98841e30      # kg
    for msun in (3.0, 30.0, 150.0):
        M = msun * M_sun
        n = n_per_planck * M / M_P
        print(f"  M = {msun:6.1f} M_sun:  n = {n:.3e} lines,  S_BH/k_B = {4*pi*(M/M_P)**2:.3e}")
    print("  For comparison, the Onsager-Feynman vortex count of a near-extremal")
    print("  spinning droplet is J/hbar ~ (M/M_P)^2, the same order as S_BH itself.")


if __name__ == "__main__":
    check_a000088()
    asymptotic_table()
    finite_size_formula_check()
    point_defect_comparison()
    invariant_theory_checks()
    black_hole_mapping()
