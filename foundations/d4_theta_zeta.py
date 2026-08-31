#!/usr/bin/env python3
"""
d4_theta_zeta.py
================
Two exact statements about the D_4 lattice's shells, checked by direct
enumeration.

1. SHELL POPULATIONS ARE DIVISOR SUMS.  Writing D_4 as the integer vectors
   of even coordinate sum, every squared length is an even integer and the
   number of lattice points at squared length 2m is

       r(2m) = 24 * sigma_odd(m),

   with sigma_odd(m) the sum of the odd divisors of m.  This is Jacobi's
   four-square theorem restricted to the even branch: a vector's norm and
   its coordinate sum have the same parity, so D_4 is exactly the set of
   integer vectors of even norm.  The kissing number 24 is the m = 1 case.

2. THE EPSTEIN ZETA FUNCTION FACTORISES.  Summing an inverse power of the
   squared length over the lattice then has a closed form in the Riemann
   zeta function,

       zeta_{D4}(s) = 24 * 2^(-s) * (1 - 2^(1-s)) * zeta(s) * zeta(s-1),

   convergent for Re s > 2.  The pole of zeta(s-1) at s = 2 is the physical
   divergence of an inverse-fourth-power sum in four dimensions.  Analytic
   continuation assigns finite values below it, which is the exact form of
   the regularisation a Casimir-type lattice sum needs.

   The three-dimensional slice has no counterpart: counting FCC shells means
   counting representations as a sum of three squares, which run on class
   numbers rather than divisor sums, with no elementary closed form.  The
   script prints the first FCC shells alongside for the contrast.

No calculation in the framework currently consumes the closed form; it is
recorded as an exact tool for the first convergent lattice sum that needs it.

Run:  python3 d4_theta_zeta.py
"""

import itertools
from collections import Counter

from mpmath import mp, zeta

mp.dps = 30

N_SHELLS = 12
ENUMERATION_REACH = 9      # integer coordinate range for the direct count


def sigma_odd(m):
    """Sum of the odd divisors of m."""
    return sum(d for d in range(1, m + 1) if m % d == 0 and d % 2 == 1)


def enumerate_d4_shells(reach=ENUMERATION_REACH):
    """Direct count of D_4 points by squared length, out to the reach."""
    pts = (v for v in itertools.product(range(-reach, reach + 1), repeat=4)
           if sum(v) % 2 == 0)
    return Counter(sum(x * x for x in v) for v in pts)


def enumerate_fcc_shells(reach=ENUMERATION_REACH):
    """Direct count of FCC (= D_3) points by squared length, for contrast."""
    pts = (v for v in itertools.product(range(-reach, reach + 1), repeat=3)
           if sum(v) % 2 == 0)
    return Counter(sum(x * x for x in v) for v in pts)


def epstein_closed_form(s):
    """zeta_{D4}(s) from the factorisation."""
    return float(24 * 2 ** (-s) * (1 - 2 ** (1 - s)) * zeta(s) * zeta(s - 1))


def epstein_direct(s, shells):
    """Truncated direct sum over enumerated shells, zero excluded."""
    return sum(count / (norm ** s) for norm, count in shells.items() if norm > 0)


def main():
    d4 = enumerate_d4_shells()
    print("1. D_4 SHELL POPULATIONS AGAINST 24 * sigma_odd(m)")
    print(f"   {'m':>3s}{'norm':>7s}{'enumerated':>13s}{'24 sigma_odd':>15s}")
    mismatches = 0
    for m in range(1, N_SHELLS + 1):
        predicted = 24 * sigma_odd(m)
        found = d4[2 * m]
        mismatches += (found != predicted)
        print(f"   {m:3d}{2 * m:7d}{found:13d}{predicted:15d}"
              f"{'' if found == predicted else '   MISMATCH'}")

    print("\n   for contrast, the three-dimensional slice has no such formula")
    fcc = enumerate_fcc_shells()
    counts = [fcc[2 * m] for m in range(1, 9)]
    print(f"   FCC shell populations: {counts}")
    print("   (these are three-square counts, governed by class numbers)")

    print("\n2. EPSTEIN ZETA: DIRECT SUM AGAINST THE CLOSED FORM")
    print("   convergence is slow near the pole at s = 2 and fast above it")
    print(f"   {'s':>5s}{'direct (truncated)':>22s}{'closed form':>18s}{'ratio':>12s}")
    for s in (2.5, 3.0, 4.0, 6.0, 8.0):
        direct = epstein_direct(s, d4)
        closed = epstein_closed_form(s)
        print(f"   {s:5.1f}{direct:22.12f}{closed:18.12f}{direct / closed:12.8f}")

    print(f"\n   truncation reaches squared length {max(d4)}, so the residual "
          f"tail\n   explains the shortfall at small s and vanishes at large s")
    print(f"\n{'SHELL COUNTS EXACT' if mismatches == 0 else f'{mismatches} MISMATCHES'}")
    return mismatches


if __name__ == "__main__":
    raise SystemExit(main())
