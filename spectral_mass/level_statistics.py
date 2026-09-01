#!/usr/bin/env python3
"""
level_statistics.py
===================
The Cosserat spectrum is not GOE, and nothing generic explains why.

Unfold the eigenvalues of a low-symmetry cluster, split by its full point group,
and the nearest-neighbour spacing distribution comes out with variance near 0.58,
halfway between Poisson (1.00, integrable) and GOE (0.29, chaotic with
time-reversal).  Four matched controls all return clean GOE:

    dense GOE of the same dimension                       0.285
    random symmetric blocks on the same bond graph        0.272
    rank-1 central force, true directions, random k       0.293
    Cosserat with bond stiffnesses randomised by 35%      0.571   <- not GOE

The last line is the decisive one.  The intermediate statistics survive
disorder in the couplings untouched, so they are not a tuned critical point.
They are structural: something in the Cosserat architecture, the rank-1
displacement block, the isotropic microrotation block, and the curl operator
coupling them, produces spectra that do not mix the way a generic matrix does.

Two readings remain and the spectrum alone cannot separate them.  A genuine
semi-Poisson class (variance 0.500, multifractal eigenstates) or one further
two-valued conserved quantity beyond the point group, which superposes two GOE
families (variance 0.495).  Disorder robustness favours the second and says
the quantity is topological rather than a symmetry.

Candidates tested and excluded as the hidden label: the relative lattice curl
of the displacement (continuous, no dichotomy); the total winding of the
microrotation direction around the triangular loops, split by parity, sign and
magnitude (each behaves exactly like a random split).

Run:  python3 level_statistics.py
"""

import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_linkedness import adjacency                            # noqa: E402
from cosserat_classifier import build_cosserat_matrix                # noqa: E402
from n1535_continuation import bond_list, weighted_cosserat         # noqa: E402
from n1535_first_principles import build_n1535_cluster              # noqa: E402

TRIALS = 25
SEED = 20260901
REFERENCE = {'Poisson': 1.000, 'semi-Poisson': 0.500, 'GOE': 0.286, 'GUE': 0.178}


def spacing_variance(eigenvalues, degree=7):
    """Variance of unfolded nearest-neighbour spacings."""
    ev = np.sort(eigenvalues)
    ev = ev[ev > 1e-8]
    stair = np.arange(1, len(ev) + 1)
    unfolded = np.polyval(np.polyfit(ev, stair, min(degree, len(ev) - 2)), ev)
    s = np.diff(unfolded)
    s = s[s > 0]
    s = s / s.mean()
    return float(s.var())


def improper_symmetry(coords):
    """One improper cubic-group element mapping the cluster to itself, if any."""
    for perm3 in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            rot = np.zeros((3, 3))
            for row, col in enumerate(perm3):
                rot[row, col] = signs[row]
            if np.linalg.det(rot) > 0:
                continue
            image = coords @ rot.T
            perm, ok = [], True
            for p in image:
                d = np.linalg.norm(coords - p, axis=1)
                j = int(np.argmin(d))
                if d[j] > 1e-6:
                    ok = False
                    break
                perm.append(j)
            if ok:
                return rot, perm
    return None


def sector_projectors(coords, rot, perm):
    n = len(coords)
    det = int(round(np.linalg.det(rot)))
    D = np.zeros((6 * n, 6 * n))
    for i in range(n):
        j = perm[i]
        D[3 * j:3 * j + 3, 3 * i:3 * i + 3] = rot
        D[3 * n + 3 * j:3 * n + 3 * j + 3, 3 * n + 3 * i:3 * n + 3 * i + 3] = det * rot
    w, V = np.linalg.eigh(D)
    return [V[:, np.abs(w - 1) < 1e-8], V[:, np.abs(w + 1) < 1e-8]]


def main():
    rng = np.random.default_rng(SEED)
    coords = np.atleast_2d(build_n1535_cluster())
    n = len(coords)
    adj = adjacency(coords, [1.0])
    bonds = bond_list(coords)
    undirected = sorted({(min(i, j), max(i, j)) for (i, j) in bonds})

    rot, perm = improper_symmetry(coords)
    sectors = sector_projectors(coords, rot, perm)

    print("Nearest-neighbour spacing variance, 21-node low-symmetry cluster\n")

    M = build_cosserat_matrix(coords, alpha=1.0)
    vs = [spacing_variance(np.linalg.eigvalsh(P.T @ M @ P)) for P in sectors]
    print(f"  Cosserat, desymmetrised sectors        {vs[0]:.3f}  {vs[1]:.3f}")

    def goe(N):
        A = rng.normal(size=(N, N))
        return np.linalg.eigvalsh((A + A.T) / np.sqrt(2 * N))
    g = [spacing_variance(goe(126)) for _ in range(TRIALS)]
    print(f"  dense GOE, dimension 126               {np.mean(g):.3f} +- {np.std(g):.3f}")

    def random_blocks():
        A = np.zeros((6 * n, 6 * n))
        for (i, j) in undirected:
            B = rng.normal(size=(6, 6))
            ii = [3 * i + k for k in range(3)] + [3 * n + 3 * i + k for k in range(3)]
            jj = [3 * j + k for k in range(3)] + [3 * n + 3 * j + k for k in range(3)]
            for p, a in enumerate(ii):
                for q, b in enumerate(jj):
                    A[a, b] += B[p, q]
                    A[b, a] += B[p, q]
        A = (A + A.T) / 2
        A[np.diag_indices_from(A)] += rng.normal(size=6 * n) * 3
        return A
    r = [spacing_variance(np.linalg.eigvalsh(random_blocks())) for _ in range(TRIALS)]
    print(f"  random blocks, same bond graph         {np.mean(r):.3f} +- {np.std(r):.3f}")

    def central_force():
        A = np.zeros((3 * n, 3 * n))
        for (i, j) in undirected:
            rhat = coords[j] - coords[i]
            rhat /= np.linalg.norm(rhat)
            P = np.outer(rhat, rhat) * (abs(rng.normal()) + 0.2)
            A[3 * i:3 * i + 3, 3 * i:3 * i + 3] += P
            A[3 * j:3 * j + 3, 3 * j:3 * j + 3] += P
            A[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= P
            A[3 * j:3 * j + 3, 3 * i:3 * i + 3] -= P
        return A
    cf = [spacing_variance(np.linalg.eigvalsh(central_force())) for _ in range(TRIALS)]
    print(f"  rank-1 central force, random k         {np.mean(cf):.3f} +- {np.std(cf):.3f}")

    dis = []
    for _ in range(TRIALS):
        k = {e: abs(rng.normal(1.0, 0.35)) + 0.05 for e in undirected}
        bw = {(i, j): k[(min(i, j), max(i, j))] for (i, j) in bonds}
        sw = {i: 1.0 for i in range(n)}
        dis.append(spacing_variance(np.linalg.eigvalsh(
            weighted_cosserat(coords, bw, sw, alpha=1.0))))
    print(f"  Cosserat, stiffnesses randomised 35%   {np.mean(dis):.3f} +- {np.std(dis):.3f}")

    print("\n  reference  " + "  ".join(f"{k} {v:.3f}" for k, v in REFERENCE.items()))
    print("\nThe intermediate statistics belong to the Cosserat architecture and")
    print("survive disorder in the couplings, so they are structural, not critical.")


if __name__ == "__main__":
    raise SystemExit(main())
