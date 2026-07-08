"""
The tangle partition function as a holonomy-field integral, and its
large-N saddle analysis.

Physics context (monograph, black-holes chapter; paper "Black holes as
superfluid droplets"): the exact gauged partition function

    Z_N(x;k) = (1/N!) sum_sigma prod_{pair cycles c} (1-x^{k l_c})/(1-x^{l_c})

can be rewritten, by the plethystic logarithm, as an exponential of the
S_N holonomy traces:

    Z_N(x;k) = (1/N!) sum_sigma exp[ sum_{u>=1} (1/u) chi2(sigma^u) (x^u - x^{ku}) ],

where chi2(tau) = number of unordered pairs fixed by tau
              = C(F1(tau),2) + F2(tau),
F1 = fixed points of tau, F2 = number of 2-cycles of tau, and for
tau = sigma^u these follow from the cycle type {m_l} of sigma:

    F1(sigma^u) = sum_{l | u} l m_l,
    F2(sigma^u) = [F1(sigma^{2u}) - F1(sigma^u)] / 2.

Since chi2 ~ t_u^2/2 with t_u = Tr sigma^u, the exponent is quadratic in
the holonomy traces: the S_N counterpart of the |Tr U^n|^2 effective
action of the Hagedorn/deconfinement analysis in large-N gauge theory
(Sundborg; Aharony-Marsano-Minwalla-Papadodimas-Van Raamsdonk), with the
permutation eigenvalues (roots of unity clustered by cycle) playing the
role of the Polyakov-loop eigenvalue density.

Checks and outputs:
  1. IDENTITY: the holonomy-field form equals the pair-cycle form for
     random permutations, and Z_N summed both ways agrees, to machine
     precision.
  2. SADDLES: the identity condensate gives ln Z ~ (N^2/2) ln f(x)
     (the N^2, black-hole saddle); the long-cycle phase gives Z ~ O(1)
     (the trivial saddle).  The exchange of dominance defines
     x_c ~ 2 ln N / ((k-1) N), the S_N deconfinement point, located
     numerically and compared with the asymptote.
  3. GEOMETRY: with x(r) = N^{-2(xi/xi_*)^3}, the deconfinement point
     sits at xi_c = 2^{-1/3} xi_* ~ 0.79 xi_*: inside the Hagedorn skin,
     between first winding (xi_*) and full mixing (0.18 xi_*).
"""

from math import log, exp, lgamma, gcd
from itertools import permutations
import random

from tangle_matrix_model_matching import (partitions, sym_factor,
                                          pair_cycle_lengths)
from tangle_partition_exact import Z_closed


# ----------------------------------------------------------------------
# chi2 from the cycle type, and the two forms of the sigma-term
# ----------------------------------------------------------------------

def F1(m, u):
    """Fixed points of sigma^u for cycle type m = {l: m_l}."""
    return sum(l * ml for l, ml in m.items() if u % l == 0)


def chi2(m, u):
    """Number of unordered pairs (i<j) fixed (as a pair) by sigma^u."""
    f1, f2 = F1(m, u), (F1(m, 2 * u) - F1(m, u)) // 2
    return f1 * (f1 - 1) // 2 + f2


def term_pair_cycles(p, k, x):
    """prod over pair cycles of (1-x^{k l})/(1-x^l) for partition p."""
    t = 1.0
    for l in pair_cycle_lengths(p):
        t *= k if abs(1 - x) < 1e-14 else (1 - x ** (k * l)) / (1 - x ** l)
    return t


def term_holonomy(p, k, x, umax=400):
    """exp[sum_u (1/u) chi2(sigma^u) (x^u - x^{ku})] for partition p,
    where p is a cycle type given as a {length: multiplicity} dict."""
    m = dict(p)
    if x <= 0:
        return 1.0
    s, u = 0.0, 1
    while u <= umax:
        w = (x ** u - x ** (k * u)) / u
        if abs(w) < 1e-18 and u > 8:
            break
        s += chi2(m, u) * w
        u += 1
    return exp(s)


def Z_holonomy(N, k, x):
    """Z from the holonomy-field (plethystic) form."""
    return sum(term_holonomy(p, k, x) / sym_factor(p) for p in partitions(N))


# ----------------------------------------------------------------------
# Saddle bookkeeping: identity condensate vs long-cycle phase
# ----------------------------------------------------------------------

def lnZ_identity_saddle(N, k, x):
    """ln of the identity-class contribution: (N(N-1)/2) ln f(x) - ln N!."""
    d = N * (N - 1) // 2
    f = k if abs(1 - x) < 1e-14 else (1 - x ** k) / (1 - x)
    return d * log(f) - (lgamma(N + 1))


def x_deconfine(N, k, tol=1e-12):
    """Numerical deconfinement point: where the identity-class term
    crosses O(1), i.e. lnZ_identity = 0 (trivial phase has Z ~ 1)."""
    lo, hi = 1e-12, 1.0 - 1e-12
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if lnZ_identity_saddle(N, k, mid) < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    print("1. Holonomy form vs pair-cycle form, per class and summed:")
    for N, k in ((5, 2), (6, 2), (6, 3), (7, 2)):
        worst = 0.0
        for x in (0.05, 0.3, 0.6, 0.9):
            for p in partitions(N):
                a, b = term_pair_cycles(p, k, x), term_holonomy(p, k, x)
                worst = max(worst, abs(a - b) / max(a, 1e-300))
            za, zb = Z_closed(N, k, x), Z_holonomy(N, k, x)
            worst = max(worst, abs(za - zb) / za)
        status = "PASS" if worst < 1e-9 else "FAIL"
        print(f"   N={N} k={k}: worst rel. err. {worst:.2e}  {status}")
        assert worst < 1e-9

    print("\n2. Deconfinement point x_c vs asymptote 2 ln N / ((k-1) N), k=2:")
    print(f"   {'N':>8} {'x_c (numeric)':>14} {'2lnN/N':>10} {'ratio':>7}")
    for N in (10, 30, 100, 300, 1000, 10000):
        xc = x_deconfine(N, 2)
        asym = 2 * log(N) / N
        print(f"   {N:>8} {xc:>14.6g} {asym:>10.4g} {xc/asym:>7.4f}")

    print("\n3. Geometric location: x(xi) = N^(-2(xi/xi_*)^3) crosses x_c at")
    print("   (xi_c/xi_*)^3 = [ln N - ln(2 ln N)] / (2 ln N)  ->  1/2 at large N")
    for lnN in (45, 89.2, 200):     # ln N for PBH-ish, solar, SMBH scales
        cube = (lnN - log(2 * lnN)) / (2 * lnN)
        print(f"   ln N = {lnN:>6}: xi_c/xi_* = {cube ** (1/3):.4f}"
              f"   (2^(-1/3) = {2 ** (-1/3):.4f})")
