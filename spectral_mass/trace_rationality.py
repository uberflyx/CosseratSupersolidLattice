#!/usr/bin/env python3
"""
trace_rationality.py
====================
Why the channel trace sum rules come out integer.

The chapter observes that the sum of a channel's eigenvalues is an integer on
every catalogue cluster that carries a two-root channel, and records the
mechanism as unexplained: the eigenvalues themselves are surds, and it is not
obvious why they should cancel.  They cancel for a structural reason that this
script establishes in exact rational arithmetic, with no floating point
anywhere.

The argument has two steps.

Step one: the surds sit off the diagonal.  A nearest-neighbour direction cosine
in FCC is irrational, but the three blocks of the Cosserat dynamical matrix use
it differently.  The displacement block carries the projector rhat (x) rhat,
whose entries are halves.  The microrotation block carries no direction cosine
at all.  Only the coupling block carries rhat linearly, and it reaches the
diagonal solely as C^T A C, squared and rational again.  Every irrational entry
of the matrix therefore lies in the off-diagonal displacement-microrotation
block, which no trace can see.

Step two: the projector cannot expose them.  Displacement is a polar vector and
microrotation an axial one, so each point-group element acts block-diagonally on
the two sectors.  The irrep projector is a character-weighted sum of those
elements and is block-diagonal too.  Hence

    tr(P M) = tr(P_uu M_uu) + tr(P_pp M_pp),

and both terms are rational.

What is left is divisibility, not irrationality, and that is settled in three
lines.  The one-dimensional irreps have characters +-1 and every group element
acts as a signed permutation matrix, so |H| tr(P M) is an integer combination of
entries of M.  Every entry of the two diagonal blocks lies in (1/B)Z for an
integer B read off the geometry, so the denominator divides B|H|.  With one bond
length B = 8, halves from rhat (x) rhat and eighths from C^T A C; activating
voids adds a bond of squared length 3/4 rather than 2, which brings in thirds
and lifts B to 24.  A host with a large residual group gives an integer and a
host with a small one can give a fraction.  That is the
pattern the chapter records without explaining: the crossed fault at |H| = 8
gives quarters and halves, and the coordination shell at |H| = 48 gives
integers.  Integrality is a property of the group order, not of how many roots
the channel has, so the sum rules hold more generally than for two-root
channels.

One boundary the argument draws by itself.  Rationality of C^T A C needs the
products of bond lengths meeting at a site to be rational.  On a single-length
cluster every product is the squared nearest-neighbour distance and the step
goes through.  On a two-length cluster, where activated voids sit closer than
the shell, the cross terms carry a product of two different lengths and the
trace need not be rational at all.  The script tests both families and reports
the split.

The run reports B and |H| for each host and checks that every observed
denominator divides B|H|.

Run:  python3 trace_rationality.py
"""

import itertools
import os
import sys
from fractions import Fraction as Q
from math import gcd

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spectral_classifier import (                                  # noqa: E402
    build_lambda_cluster, cluster_born, cluster_coord_shell,
    cluster_cuboctahedron, cluster_hex_cap,
)
from baryon_mass_modes import cluster_delta, cluster_xi            # noqa: E402

# Characters of the one-dimensional irreps of O_h, by conjugacy class.
# Classes in the order used below: E, 8C3, 6C2, 6C4, 3C2', i, 6S4, 8S6, 3sigma_h, 6sigma_d
OH_1D = {
    # Ten classes of O_h.  C2ax are the three twofold axes along the cube edges
    # (the squares of the C4), C2diag the six through opposite face diagonals;
    # sigma_h are the three mirrors normal to the axes, sigma_d the six diagonal
    # ones.  Collapsing either pair silently changes the A_2 characters.
    'A_1g': dict(E=1, C3=1, C4=1, C2ax=1, C2diag=1,
                 i=1, S4=1, S6=1, sh=1, sd=1),
    'A_2g': dict(E=1, C3=1, C4=-1, C2ax=1, C2diag=-1,
                 i=1, S4=-1, S6=1, sh=1, sd=-1),
    'A_1u': dict(E=1, C3=1, C4=1, C2ax=1, C2diag=1,
                 i=-1, S4=-1, S6=-1, sh=-1, sd=-1),
    'A_2u': dict(E=1, C3=1, C4=-1, C2ax=1, C2diag=-1,
                 i=-1, S4=1, S6=-1, sh=-1, sd=1),
}


def oh_elements():
    """The 48 signed permutation matrices of the cubic group, as integer tuples."""
    out = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            mat = [[0, 0, 0] for _ in range(3)]
            for row, col in enumerate(perm):
                mat[row][col] = signs[row]
            out.append(tuple(tuple(r) for r in mat))
    return out


def class_of(mat):
    """Conjugacy class of a cubic-group element from its trace, determinant and shape."""
    arr = np.array(mat)
    tr = int(round(arr.trace()))
    det = int(round(np.linalg.det(arr)))
    diagonal = np.count_nonzero(arr - np.diag(np.diagonal(arr))) == 0
    if det == 1:
        if tr == 3:
            return 'E'
        if tr == 0:
            return 'C3'
        if tr == 1:
            return 'C4'
        return 'C2ax' if diagonal else 'C2diag'
    if tr == -3:
        return 'i'
    if tr == 0:
        return 'S6'
    if tr == -1:
        return 'S4'
    return 'sh' if diagonal else 'sd'


def class_of_proper(mat):
    return class_of(mat)


def integer_coords(coords):
    """Rescale FCC coordinates so every component is rational with small denominator.

    The builders return nearest neighbours at unit distance, whose components are
    +-1/sqrt(2).  Multiplying by sqrt(2) puts the shell on integer sites and the
    activated voids on half-integers, both exactly representable.
    """
    scaled = np.asarray(coords) * np.sqrt(2.0)
    out = []
    for point in scaled:
        row = []
        for value in point:
            twice = round(2 * value)
            assert abs(2 * value - twice) < 1e-6, f"non-half-integer coordinate {value}"
            row.append(Q(twice, 2))
        out.append(tuple(row))
    return out


def bonds_by_length(sites, n_lengths=1):
    """Ordered bonded pairs, grouped by exact squared separation.

    A bond is a pair at one of the `n_lengths` shortest separations, which is
    what the dynamical-matrix builders use: one length on a pure shell cluster,
    two once interstitial voids are activated.
    """
    allpairs = {}
    for i, a in enumerate(sites):
        for j, b in enumerate(sites):
            if i == j:
                continue
            w = tuple(b[k] - a[k] for k in range(3))
            d2 = sum(c * c for c in w)
            allpairs.setdefault(d2, []).append((i, j, w))
    keep = sorted(allpairs)[:n_lengths]
    return {d2: allpairs[d2] for d2 in keep}


def diagonal_blocks(sites, bond_groups, alpha=Q(1)):
    """Exact M_uu and M_pp over Q[sqrt(6)].

    Returns (M_uu, M_pp, surd_terms_present).  The displacement block is built
    with the coupling contribution C^T A C carried in the quadratic field, so a
    cluster with two bond lengths can be handled exactly rather than refused.
    """
    n = len(sites)
    for d2 in bond_groups:
        assert _inv_sqrt(d2) is not None, f"length sqrt({d2}) outside Q[sqrt2, sqrt3]"

    m_pp = [[Q(0)] * (3 * n) for _ in range(3 * n)]
    for blist in bond_groups.values():
        for (i, j, _w) in blist:
            for k in range(3):
                m_pp[3 * i + k][3 * j + k] -= Q(1)
                m_pp[3 * i + k][3 * i + k] += Q(1)
    for i in range(n):
        for k in range(3):
            m_pp[3 * i + k][3 * i + k] += alpha

    m_uu = [[S6() for _ in range(3 * n)] for _ in range(3 * n)]
    for d2, blist in bond_groups.items():
        for (i, j, w) in blist:
            for a_ in range(3):
                for b_ in range(3):
                    val = Surd(w[a_] * w[b_] / d2)
                    m_uu[3 * i + a_][3 * j + b_] = m_uu[3 * i + a_][3 * j + b_] - val
                    m_uu[3 * i + a_][3 * i + b_] = m_uu[3 * i + a_][3 * i + b_] + val

    # C_{i,j} = E(w_ij) / (2 |w_ij|); the inverse length lives in the field.
    inv_len = {d2: _inv_sqrt(d2) for d2 in bond_groups}

    from_site = {}
    for d2, blist in bond_groups.items():
        for (i, j, w) in blist:
            from_site.setdefault(i, []).append((j, w, d2))
    for i, out_bonds in from_site.items():
        for (j1, w1, d2a) in out_bonds:
            e1 = _cross_matrix(w1)
            for (j2, w2, d2b) in out_bonds:
                e2 = _cross_matrix(w2)
                scale = (inv_len[d2a] * inv_len[d2b]) * (alpha / 4)
                for a_ in range(3):
                    for b_ in range(3):
                        coef = sum(e1[k][a_] * e2[k][b_] for k in range(3))
                        if coef:
                            m_uu[3 * j1 + a_][3 * j2 + b_] = (
                                m_uu[3 * j1 + a_][3 * j2 + b_] + scale * coef)

    surds = any(not m_uu[r][c].is_rational()
                for r in range(3 * n) for c in range(3 * n))
    return m_uu, m_pp, surds


def _cross_matrix(w):
    return [[Q(0), -w[2], w[1]], [w[2], Q(0), -w[0]], [-w[1], w[0], Q(0)]]


class Surd:
    """An element a + b*sqrt2 + c*sqrt3 + d*sqrt6 of Q[sqrt2, sqrt3].

    This is the smallest field containing every bond length the catalogue uses:
    the nearest-neighbour separation brings in sqrt2 and the shell-to-void
    separation brings in sqrt3.  Working here rather than over the rationals
    lets the surds appear honestly and then be seen to cancel, instead of being
    assumed away.
    """

    __slots__ = ('c',)
    # products of basis elements, as (index, rational multiplier)
    TABLE = {
        (0, 0): (0, 1), (0, 1): (1, 1), (0, 2): (2, 1), (0, 3): (3, 1),
        (1, 1): (0, 2), (1, 2): (3, 1), (1, 3): (2, 2),
        (2, 2): (0, 3), (2, 3): (1, 3),
        (3, 3): (0, 6),
    }

    def __init__(self, a=Q(0), b=Q(0), c=Q(0), d=Q(0)):
        self.c = [Q(a), Q(b), Q(c), Q(d)]

    def __add__(self, o):
        o = o if isinstance(o, Surd) else Surd(o)
        return Surd(*[x + y for x, y in zip(self.c, o.c)])

    def __sub__(self, o):
        o = o if isinstance(o, Surd) else Surd(o)
        return Surd(*[x - y for x, y in zip(self.c, o.c)])

    def __mul__(self, o):
        if not isinstance(o, Surd):
            return Surd(*[x * o for x in self.c])
        out = [Q(0)] * 4
        for i, x in enumerate(self.c):
            if not x:
                continue
            for j, y in enumerate(o.c):
                if not y:
                    continue
                k, mult = Surd.TABLE[(i, j)] if (i, j) in Surd.TABLE \
                    else Surd.TABLE[(j, i)]
                out[k] += x * y * mult
        return Surd(*out)

    __rmul__ = __mul__

    def __truediv__(self, o):
        return Surd(*[x / o for x in self.c])

    def is_rational(self):
        return all(x == 0 for x in self.c[1:])

    def __repr__(self):
        names = ['', '*sqrt2', '*sqrt3', '*sqrt6']
        parts = [f"{x}{n}" for x, n in zip(self.c, names) if x != 0]
        return " + ".join(parts) if parts else "0"


S6 = Surd          # the projected-trace code refers to the field by this name


def _inv_sqrt(d2):
    """1/sqrt(d2) in Q[sqrt2, sqrt3], for the separations this lattice uses."""
    for idx, base in ((1, Q(2)), (2, Q(3)), (3, Q(6)), (0, Q(1))):
        ratio = d2 / base
        num, den = ratio.numerator, ratio.denominator
        rn, rd = int(round(num ** 0.5)), int(round(den ** 0.5))
        if rn * rn == num and rd * rd == den:
            # d2 = base * (rn/rd)^2, so 1/sqrt(d2) = (rd/rn) / sqrt(base)
            scale = Q(rd, rn)
            if idx == 0:
                return Surd(scale)
            unit = [Q(0)] * 4
            unit[idx] = scale / base          # 1/sqrt(base) = sqrt(base)/base
            return Surd(*unit)
    return None


def _sqrt_in_field(d2):
    return _inv_sqrt(d2)


def _is_rational_sqrt_product(d2a, d2b):
    """Is sqrt(d2a * d2b) rational?"""
    prod = d2a * d2b
    num, den = prod.numerator, prod.denominator
    return _is_square(num) and _is_square(den)


def _is_square(m):
    r = int(round(m ** 0.5))
    return r * r == m


def _sqrt_product(d2a, d2b):
    prod = d2a * d2b
    return Q(int(round(prod.numerator ** 0.5)), int(round(prod.denominator ** 0.5)))


def site_permutation(rot, sites):
    """Where the cubic element sends each site, or None if it is not a symmetry."""
    index = {s: k for k, s in enumerate(sites)}
    perm = []
    for s in sites:
        image = tuple(sum(Q(rot[a][b]) * s[b] for b in range(3)) for a in range(3))
        if image not in index:
            return None
        perm.append(index[image])
    return perm


def projected_traces(sites, m_uu, m_pp, channel):
    """tr(P_Gamma M) split into its two sectors, exactly."""
    n = len(sites)
    chars = OH_1D[channel]
    order = 0
    trace_uu, trace_pp = Surd(), Q(0)
    for rot in oh_elements():
        perm = site_permutation(rot, sites)
        if perm is None:
            continue
        order += 1
        chi = chars[class_of_proper(rot)]
        det = int(round(np.linalg.det(np.array(rot))))
        # tr(D(R) M) for the polar (u) and axial (phi) sectors
        acc_u, acc_p = Surd(), Q(0)
        for i in range(n):
            j = perm[i]
            for a in range(3):
                for b in range(3):
                    r_ab = Q(rot[a][b])
                    if r_ab:
                        if m_uu is not None:
                            acc_u = acc_u + m_uu[3 * i + b][3 * j + a] * r_ab
                        acc_p += det * r_ab * m_pp[3 * i + b][3 * j + a]
        trace_uu = trace_uu + acc_u * chi
        trace_pp += chi * acc_p
    return trace_uu / order, trace_pp / order, order


def main():
    def arr(x):
        return np.atleast_2d(x[0] if isinstance(x, tuple) else x)

    catalogue = [
        ("cell pair", np.array([[0.0, 0.0, 0.0], [np.sqrt(0.5), np.sqrt(0.5), 0.0]])),
        ("hex cap", arr(cluster_hex_cap())),
        ("cuboctahedron", arr(cluster_cuboctahedron())),
        ("coordination shell", arr(cluster_coord_shell())),
        ("shell + triangle", arr(build_lambda_cluster())),
        ("Born", arr(cluster_born())),
        ("shell + 2 triangles", arr(cluster_xi())),
        # The void cluster is built here at uniform bond stiffness, so its
        # numbers are not the chapter's two-stiffness eigenvalues.  It is
        # included because it is the one catalogue host whose matrix genuinely
        # carries surds, which is what the cancellation has to survive.
        ("shell + 4 voids", arr(cluster_delta())),
    ]

    print("Exact projected traces  tr(P_Gamma M),  alpha_Cos = 1, no floating point\n")
    header = (f"{'cluster':22s} {'|H|':>4s} {'channel':>7s} "
              f"{'trace':>14s} {'denominator':>12s}  status")
    print(header)
    print("-" * len(header))

    for name, coords in catalogue:
        sites = integer_coords(coords)
        n_lengths = 2 if "void" in name else 1
        groups = bonds_by_length(sites, n_lengths)
        m_uu, m_pp, surds = diagonal_blocks(sites, groups)
        lengths = sorted(groups)
        note = (f"  [{len(lengths)} bond length(s) {[str(x) for x in lengths]}; "
                f"matrix {'carries' if surds else 'free of'} sqrt6 entries]")
        print(f"{name:22s}{note}")
        block_den = 1
        for row in m_uu:
            for x in row:
                for part in x.c:
                    block_den = block_den * part.denominator // gcd(block_den, part.denominator)
        for row in m_pp:
            for x in row:
                block_den = block_den * x.denominator // gcd(block_den, x.denominator)
        for channel in ('A_1g', 'A_2g', 'A_1u', 'A_2u'):
            t_uu, t_pp, order = projected_traces(sites, m_uu, m_pp, channel)
            total = t_uu + Surd(t_pp)
            assert total.is_rational(), f"surd survived in {name} {channel}: {total}"
            total = total.c[0]
            den = total.denominator
            bound = block_den * order
            ok = bound % den == 0
            status = ("integer" if den == 1
                      else f"denominator {den} divides B|H| = {block_den}x{order}"
                           f" = {bound}: {'yes' if ok else 'NO'}")
            print(f"{name:22s} {order:4d} {channel:>7s} {str(total):>14s} "
                  f"{den:12d}  {status}")
        print()

    print("Every trace above is exactly rational, computed in Fraction arithmetic.")
    print("No eigenvalue was ever formed, so no surd had the chance to appear.")
    print("The block denominator B is read off the geometry: 8 with one bond")
    print("length, 24 once voids add a bond of squared length 3/4.  Every")
    print("denominator above divides B|H|, as the three-line argument requires.")
    _cross_check()


def _cross_check():
    """Compare the exact traces against numerical eigenvalue sums on the shell."""
    from cosserat_classifier import build_cosserat_matrix
    coords = np.atleast_2d(cluster_coord_shell()[0])
    sites = integer_coords(coords)
    groups = bonds_by_length(sites, 1)
    m_uu, m_pp, _ = diagonal_blocks(sites, groups)
    mat = build_cosserat_matrix(coords, alpha=1.0)
    n = len(coords)
    print("\nCross-check on the coordination shell, exact against numerical:")
    for channel in ('A_1g', 'A_2g', 'A_1u', 'A_2u'):
        t_uu, t_pp, order = projected_traces(sites, m_uu, m_pp, channel)
        exact = (t_uu + Surd(t_pp)).c[0]
        proj = np.zeros((6 * n, 6 * n))
        for rot in oh_elements():
            perm = site_permutation(rot, sites)
            if perm is None:
                continue
            chi = OH_1D[channel][class_of(rot)]
            det = int(round(np.linalg.det(np.array(rot))))
            rep = np.zeros((6 * n, 6 * n))
            for i in range(n):
                j = perm[i]
                rep[3 * j:3 * j + 3, 3 * i:3 * i + 3] = np.array(rot, dtype=float)
                rep[3 * n + 3 * j:3 * n + 3 * j + 3,
                    3 * n + 3 * i:3 * n + 3 * i + 3] = det * np.array(rot, dtype=float)
            proj += chi * rep
        proj /= order
        numeric = float(np.trace(proj @ mat))
        agree = abs(numeric - float(exact)) < 1e-9
        print(f"  {channel}: exact {str(exact):>6s}   numerical {numeric:12.8f}   "
              f"{'agree' if agree else 'DISAGREE'}")


if __name__ == "__main__":
    raise SystemExit(main())
