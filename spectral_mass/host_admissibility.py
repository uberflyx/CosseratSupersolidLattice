#!/usr/bin/env python3
"""
host_admissibility.py
=====================
Which FCC clusters can host a defect.

The spectral mass formula charges each disturbed node the Peierls-Nabarro
barrier height m_0, so a node only belongs to a host if the cluster actually
holds it at the top of that barrier.  This script turns that requirement into a
criterion on the cluster's geometry and enumerates the clusters that satisfy it.

The criterion
-------------
Under the lattice's contact law a bond carries a central spring k_n and a
tangential spring k_t = r_t k_n, with r_t = 1/(3 pi - 5) = 0.22600 the rolling
point of the D4 dictionary.  A node p bonded to cluster nodes q therefore sits
in the on-site stiffness tensor

    K_p = sum_b [ rhat_b rhat_b^T  +  r_t (I - rhat_b rhat_b^T) ] ,

where rhat_b runs over the unit vectors of p's bonds into the cluster.  A fully
coordinated FCC node has twelve bonds and K = (Z1/3 + 8 r_t) I = 5.808 I,
isotropic, whose central part Z1/3 = 4 is the mass formula's reference
eigenvalue.  If p's bond directions fail to span R^3 then K_p has an eigenvalue
of order r_t alone: along that direction nothing but the weak tangential spring
opposes the node, it slides back into its own valley at almost no cost, and
charging it a full m_0 is wrong.

    ADMISSIBILITY.  Every node of an added orbit must have bond directions
    spanning R^3.

Two further conditions come from the physics rather than from this tensor.
Orbit completeness: the disturbed set is a level set of the defect's own
disregistry, which is invariant under the defect's residual point group H, so
nodes fall in or out an H-orbit at a time.  Radial closure: for a compact
winding the disregistry falls with distance, so no orbit may be occupied while a
closer one stands empty.  Radial closure does NOT apply to a Shockley partial,
whose disregistry is a shear across a {111} habit plane rather than a radial
decay; the strange sector is flagged in the output for that reason.

Units: a/2 throughout, so an FCC site is an integer triple of even coordinate
sum, the nearest-neighbour separation squared is ell^2 = 2, and a tetrahedral
interstice sits at half-odd coordinates with void-corner separation squared 3/4.

Author: M. Cox, with Claude (Anthropic).  License: MIT.
"""
import numpy as np
from itertools import product, combinations
from fractions import Fraction as F
from collections import defaultdict

R_T = 1.0 / (3.0 * np.pi - 5.0)   # tangential/central contact ratio, D4 rolling point
ELL2 = 2.0                        # nearest-neighbour separation squared, (a/2)^2
VOID2 = 0.75                      # void-to-corner separation squared, (a/2)^2
Z1_OVER_3 = 4.0                   # the mass formula's on-site reference eigenvalue
FULL_ONSITE = Z1_OVER_3 + 8.0 * R_T
ORIGIN = (F(0), F(0), F(0))


# --------------------------------------------------------------- geometry
def r2(p):
    """Squared distance from the origin, in (a/2)^2."""
    return float(sum(F(x) ** 2 for x in p))


def lattice_sites(rmax):
    """FCC sites (integer triples of even coordinate sum) inside rmax, origin excluded."""
    n = int(np.ceil(rmax)) + 1
    return [tuple(F(c) for c in cc)
            for cc in product(range(-n, n + 1), repeat=3)
            if sum(cc) % 2 == 0 and 0 < np.dot(cc, cc) <= rmax ** 2 + 1e-9]


def tetrahedral_sites():
    """The eight tetrahedral interstices of the origin's own coordination shell."""
    return [tuple(F(s[i], 2) for i in range(3)) for s in product((1, -1), repeat=3)]


def bond_directions(p, cluster):
    """Unit vectors of p's bonds into `cluster`, at either allowed bond length."""
    out = []
    for q in cluster:
        if q == p:
            continue
        d = np.array([float(p[i] - q[i]) for i in range(3)])
        d2 = float(d @ d)
        if abs(d2 - ELL2) < 1e-9 or abs(d2 - VOID2) < 1e-9:
            out.append(d / np.sqrt(d2))
    return out


def onsite_stiffness(p, cluster):
    """(smallest eigenvalue of K_p, bond count, span of the bond directions)."""
    dirs = bond_directions(p, cluster)
    if not dirs:
        return 0.0, 0, 0
    I3 = np.eye(3)
    K = sum(np.outer(u, u) + R_T * (I3 - np.outer(u, u)) for u in dirs)
    span = int(np.linalg.matrix_rank(np.array(dirs), tol=1e-8))
    return float(np.linalg.eigvalsh(K)[0]), len(dirs), span


def admissible(orbit, host):
    """True when every node of `orbit` is held in all three directions."""
    assembled = list(host) + list(orbit)
    return all(onsite_stiffness(p, assembled)[2] == 3 for p in orbit)


# ------------------------------------------------------------ point group
def oh_group():
    """The 48 signed permutation matrices, i.e. the full octahedral group."""
    out = []
    for perm in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]:
        for sg in product((1, -1), repeat=3):
            M = np.zeros((3, 3), dtype=int)
            for i, pi in enumerate(perm):
                M[i, pi] = sg[i]
            out.append(tuple(map(tuple, M)))
    return out


OH = oh_group()


def act(M, p):
    return tuple(M[i][0] * p[0] + M[i][1] * p[1] + M[i][2] * p[2] for i in range(3))


def stabiliser_order(cluster):
    S = set(cluster)
    return sum(1 for M in OH if set(act(M, p) for p in S) == S)


# ------------------------------------------------------------------ report
def shells_by_radius(rmax):
    out = defaultdict(list)
    for p in lattice_sites(rmax):
        out[r2(p)].append(p)
    return out


def main():
    shells = shells_by_radius(4.5)
    shell13 = [ORIGIN] + shells[2.0]
    born19 = shell13 + shells[4.0]

    print("=" * 78)
    print("On-site stiffness of an added orbit")
    print(f"r_t = {R_T:.5f};  a fully coordinated node has mu = {FULL_ONSITE:.3f}")
    print("=" * 78)
    print(f"{'host':<22}{'orbit r^2':>10}{'size':>6}{'bonds':>7}{'span':>6}"
          f"{'mu_min':>9}  verdict")
    for name, host in (("coordination shell", shell13), ("Born cluster", born19)):
        rmax_host = max(r2(p) for p in host)
        for k in sorted(shells):
            if k <= rmax_host + 1e-9 or k > 10.5:
                continue
            orb = shells[k]
            cand = list(host) + orb
            vals = [onsite_stiffness(p, cand) for p in orb]
            mu = min(v[0] for v in vals)
            nb = min(v[1] for v in vals)
            sp = min(v[2] for v in vals)
            print(f"{name:<22}{k:>10g}{len(orb):>6}{nb:>7}{sp:>6}{mu:>9.4f}  "
                  f"{'admissible' if sp == 3 else 'rejected: bonds do not span'}")

    print("\n" + "=" * 78)
    print("The cubic ladder: complete O_h orbits that pass both conditions")
    print("=" * 78)
    V, ladder = [ORIGIN], [1]
    for k in sorted(shells):
        orb = shells[k]
        if not admissible(orb, V):
            print(f"   r^2 = {k:g} ({len(orb):2d} sites): rejected, "
                  f"span {min(onsite_stiffness(p, list(V) + orb)[2] for p in orb)}")
            continue
        V = V + orb
        ladder.append(len(V))
        print(f"   r^2 = {k:g} ({len(orb):2d} sites): N -> {len(V)}, "
              f"|H| = {stabiliser_order(V)}")
        if len(V) > 60:
            break
    print(f"\n   ladder: {ladder}")
    print(f"   the N = 13 + 6k reading would allow: {[13 + 6 * j for j in range(8)]}")

    print("\n" + "=" * 78)
    print("Void activation: the tetrahedral orbit on the shell and on Born")
    print("=" * 78)
    voids = tetrahedral_sites()
    pauli = [p for p in voids if float(p[0] * p[1] * p[2]) > 0]
    for tag, host in (("shell", shell13), ("Born", born19)):
        for label, sel in (("one tetrahedron", pauli), ("both tetrahedra", voids)):
            cand = list(host) + sel
            vals = [onsite_stiffness(p, cand) for p in sel]
            print(f"   {tag} + {label:<16} N = {len(cand):<3} "
                  f"mu_min = {min(v[0] for v in vals):.4f}  "
                  f"|H| = {stabiliser_order(cand)}")

    print("\n" + "=" * 78)
    print("Why the strange extension arrives as a triangle")
    print("=" * 78)
    facet_tri = sorted(shells[6.0], key=lambda p: -float(sum(p)))[:3]
    for m in (1, 2, 3):
        sub = facet_tri[:m]
        vals = [onsite_stiffness(p, shell13 + sub) for p in sub]
        sp = min(v[2] for v in vals)
        print(f"   {m} node(s) of the facet triangle: bonds "
              f"{min(v[1] for v in vals)}, span {sp}, "
              f"mu_min {min(v[0] for v in vals):.4f}  "
              f"{'admissible' if sp == 3 else 'rejected'}")
    print("   A lone node fails the span test and a pair is not a complete orbit")
    print("   of the facet's three-fold axis, so the triangle is the smallest")
    print("   addition that is both held and orbit-complete.")

    print("\n" + "=" * 78)
    print("Soft-mode check: the rejected orbit seen in the spectrum")
    print("=" * 78)
    for tag, coords in (("coordination shell", shell13),
                        ("Born (shell + r^2=4)", born19),
                        ("shell + r^2=8 orbit", shell13 + shells[8.0])):
        X = np.array([[float(c) for c in p] for p in coords])
        n = len(X)
        P = np.zeros((3 * n, 3 * n))
        I3 = np.eye(3)
        for i in range(n):
            for j in range(i + 1, n):
                d = X[i] - X[j]
                L2 = float(d @ d)
                if abs(L2 - ELL2) > 1e-9:
                    continue
                u = d / np.sqrt(L2)
                K = np.outer(u, u) + R_T * (I3 - np.outer(u, u))
                for a, b, s in ((i, i, 1), (j, j, 1), (i, j, -1), (j, i, -1)):
                    P[3 * a:3 * a + 3, 3 * b:3 * b + 3] += s * K
        w = np.sort(np.linalg.eigvalsh(P))
        internal = w[int(np.sum(w < 1e-9)):]
        print(f"   {tag:<24} N = {n:<3} lowest internal lambda = "
              f"{internal[0]:.4f}, modes below 0.5: {int(np.sum(internal < 0.5))}")


if __name__ == "__main__":
    main()
