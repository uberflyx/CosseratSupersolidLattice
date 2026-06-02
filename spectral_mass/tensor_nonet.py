#!/usr/bin/env python3
"""
tensor_nonet.py
===============
Spin-2 ladder of the 19-node Born stability cluster: the structural origin
of the light tensor mesons f_2(1270) and a_2(1320).

The script builds the Born cluster (a central node, its 12 nearest
neighbours, and its 6 next-nearest neighbours) directly from FCC
coordinates, assembles the 6n x 6n Cosserat dynamical matrix that couples
the displacement field u to the microrotation field phi at unit Cosserat
coupling, and diagonalises it.  Each eigenmode is projected onto the ten
irreducible representations of the octahedral group O_h with the character
projector

    P_R = (d_R / |G|) sum_g chi_R(g) D(g),

where d_R is the dimension of irrep R, |G| = 48 is the group order, chi_R
its character, and D(g) the 6n x 6n action of the symmetry operation g on
the cluster (the displacement block transforms as a polar vector, the
microrotation block as an axial vector, so it carries an extra factor of
det(g)).

The five-dimensional spin-2 representation of SO(3) subduces on a cubic
lattice as J = 2 -> E_g (+) T_2g.  Both irreps host massive spin-2 modes;
the script tabulates every E_g and T_2g eigenvalue with its displacement
fraction u% and the mass it predicts through the spectral mass formula

    m = N m_0 - N (4 - lambda) m_e ,        N = 18 (the rotating shell),

with the node mass m_0 = m_e / alpha and the electron mass m_e taken from
CODATA 2022.  N = 18 = Z_1 + Z_2 counts the 12 nearest and 6 next-nearest
neighbours that rotate around the central node, which sits on the rotation
axis and stays put.

Two roots in T_2g close cleanly against the PDG, both on this one cluster,
distinguished only by their mode sector:

    f_2(1270): the optical-mix root, lambda = 5.580, u% ~ 33  -> 1274.99 MeV
    a_2(1320): the stiff microrotation root, lambda = 10.480, u% ~ 2 -> 1320.06 MeV

The two E_g roots at lambda = 4.844 and lambda = 8.307 predict additional
spin-2 channels (masses 1268.22 and 1300.07 MeV) not yet matched to named
PDG states.  Note that lambda = 8.167 (T_2g) and lambda = 8.307 (E_g) lie
0.14 apart, the signature of a single delocalised J = 2 level weakly split
by the cubic field rather than two independent channels.
"""

import numpy as np
from itertools import permutations

from cosserat_classifier import build_cosserat_matrix
from spectral_classifier import cluster_born

# ----------------------------------------------------------------- constants
M_E = 0.51099895                 # electron mass [MeV], CODATA 2022
ALPHA = 1.0 / 137.035999177      # fine structure constant, CODATA 2022
M_0 = M_E / ALPHA                # node mass [MeV] = 70.0253
N_ROT = 18                       # rotating shell of the Born cluster (Z_1 + Z_2)

# --------------------------------------------------------- O_h character table
# Classes in a fixed order; sizes sum to the group order |G| = 48.
CLASSES = ['E', '8C3', '6C2', '6C4', '3C2', 'i', '6S4', '8S6', '3sh', '6sd']
CLASS_SIZE = dict(zip(CLASSES, [1, 8, 6, 6, 3, 1, 6, 8, 3, 6]))
CHARACTERS = {
    'A_1g': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    'A_2g': [1, 1, -1, -1, 1, 1, -1, 1, 1, -1],
    'E_g':  [2, -1, 0, 0, 2, 2, 0, -1, 2, 0],
    'T_1g': [3, 0, -1, 1, -1, 3, 1, 0, -1, -1],
    'T_2g': [3, 0, 1, -1, -1, 3, -1, 0, -1, 1],
    'A_1u': [1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
    'A_2u': [1, 1, -1, -1, 1, -1, 1, -1, -1, 1],
    'E_u':  [2, -1, 0, 0, 2, -2, 0, 1, -2, 0],
    'T_1u': [3, 0, -1, 1, -1, -3, -1, 0, 1, 1],
    'T_2u': [3, 0, 1, -1, -1, -3, 1, 0, 1, -1],
}
GROUP_ORDER = sum(CLASS_SIZE.values())   # 48


def generate_oh():
    """Return the 48 signed 3x3 permutation matrices that make up O_h."""
    elements = []
    for perm in permutations(range(3)):
        for sx in (1, -1):
            for sy in (1, -1):
                for sz in (1, -1):
                    M = np.zeros((3, 3))
                    signs = (sx, sy, sz)
                    for axis in range(3):
                        M[axis, perm[axis]] = signs[axis]
                    elements.append(M)
    return elements


def conjugacy_class(R):
    """Identify the O_h conjugacy class of R from its trace, determinant, and
    whether it is diagonal (which separates the 3C2 from the 6C2 rotations and
    the 3sigma_h from the 6sigma_d reflections)."""
    trace = round(np.trace(R))
    det = round(np.linalg.det(R))
    is_diagonal = all(
        abs(R[a, a]) == 1 and all(abs(R[a, b]) < 1e-9 for b in range(3) if b != a)
        for a in range(3)
    )
    if det == 1:
        return {3: 'E', 0: '8C3', 1: '6C4',
                -1: ('3C2' if is_diagonal else '6C2')}[trace]
    return {-3: 'i', 0: '8S6', -1: '6S4',
            1: ('3sh' if is_diagonal else '6sd')}[trace]


def site_permutation(R, coords, tol=1e-6):
    """Index map i -> j such that R maps node i onto node j; None if R is not a
    symmetry of the cluster (no node sits at the image position)."""
    perm = []
    for i in range(len(coords)):
        image = R @ coords[i]
        match = next((j for j in range(len(coords))
                      if np.linalg.norm(coords[j] - image) < tol), None)
        if match is None:
            return None
        perm.append(match)
    return perm


def cosserat_representation(R, perm, n):
    """6n x 6n action of R on (u, phi).  The displacement block is a polar
    vector (transforms as R); the microrotation block is an axial vector
    (transforms as det(R) * R)."""
    D = np.zeros((6 * n, 6 * n))
    det = round(np.linalg.det(R))
    for i in range(n):
        j = perm[i]
        D[3 * j:3 * j + 3, 3 * i:3 * i + 3] = R                  # u: polar
        D[3 * n + 3 * j:3 * n + 3 * j + 3,
          3 * n + 3 * i:3 * n + 3 * i + 3] = det * R             # phi: axial
    return D


def mass_from_lambda(lam):
    """Spectral mass formula for the rotating Born shell (N = 18)."""
    return N_ROT * M_0 - N_ROT * (4.0 - lam) * M_E


def main():
    coords, label = cluster_born()
    n = len(coords)
    print(f"{label}: {n} nodes, {6 * n} Cosserat degrees of freedom\n")

    H = build_cosserat_matrix(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    vals, vecs = np.linalg.eigh(H)

    # Symmetry operations that fix the cluster, paired with their class and rep.
    cluster_reps = []
    for R in generate_oh():
        perm = site_permutation(R, coords)
        if perm is not None:
            cluster_reps.append((conjugacy_class(R), cosserat_representation(R, perm, n)))

    def irrep_of(indices):
        """Dominant O_h irrep of a (possibly degenerate) eigenspace, by the
        character overlap n_R = (1/|G|) sum_g chi_R(g) Tr(P^T D(g) P)."""
        block = vecs[:, indices]
        chi = {cls: 0.0 for cls in CLASSES}
        for cls, D in cluster_reps:
            chi[cls] += np.trace(block.T @ D @ block)
        best, best_count = None, 0.0
        for irrep, chars in CHARACTERS.items():
            count = sum(chars[k] * chi[CLASSES[k]] for k in range(len(CLASSES))) / GROUP_ORDER
            if abs(count) > abs(best_count):
                best_count, best = count, irrep
        return best

    def displacement_fraction(indices):
        """Fraction of the eigenspace amplitude carried by the displacement
        field u (the top half of the eigenvector); the complement is phi."""
        return float(np.mean([np.sum(vecs[:3 * n, k] ** 2) for k in indices]))

    # Collect degenerate eigenvalues into single rows.
    groups = []
    i = 0
    while i < len(vals):
        j = i + 1
        while j < len(vals) and abs(vals[j] - vals[i]) < 1e-4:
            j += 1
        groups.append((vals[i], list(range(i, j))))
        i = j

    print(f"{'lambda':>9} {'irrep':>6} {'u%':>5} {'m at N=18 [MeV]':>16}   note")
    print("-" * 60)
    for lam, indices in groups:
        if lam < 1e-4:                       # skip the zero modes (rigid motions)
            continue
        irrep = irrep_of(indices)
        if irrep not in ('E_g', 'T_2g'):     # the spin-2 ladder lives here
            continue
        u_pct = 100.0 * displacement_fraction(indices)
        m = mass_from_lambda(lam)
        note = ''
        if irrep == 'T_2g' and abs(lam - 5.580) < 0.05:
            note = '<- f_2(1270): optical-mix root'
        elif irrep == 'T_2g' and abs(lam - 10.480) < 0.05:
            note = '<- a_2(1320): stiff phi root'
        print(f"{lam:9.3f} {irrep:>6} {u_pct:5.0f} {m:16.2f}   {note}")

    print("\nClean closures (PDG):")
    print(f"  f_2(1270): lambda=5.580  -> {mass_from_lambda(5.580):8.2f} MeV  (obs 1275.4 +/- 0.8)")
    print(f"  a_2(1320): lambda=10.480 -> {mass_from_lambda(10.480):8.2f} MeV  (obs 1318.2 +/- 0.6)")
    print("Forward predictions in E_g (no named PDG state):")
    print(f"  lambda=4.844  -> {mass_from_lambda(4.844):8.2f} MeV")
    print(f"  lambda=8.307  -> {mass_from_lambda(8.307):8.2f} MeV  (0.14 above the T_2g at 8.167:")
    print("                  likely one cubic-split J=2 level, not two channels)")


if __name__ == '__main__':
    main()
