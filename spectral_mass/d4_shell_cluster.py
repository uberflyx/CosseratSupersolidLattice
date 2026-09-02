#!/usr/bin/env python3
"""
d4_shell_cluster.py
===================
The coordination shell as a D_4 cluster: four displacement components and six
rotation planes per node, all twenty-four D_4 nearest-neighbour bonds, the
clamped-frame curl on every plane, single-scale convention. Two treatments of
the twelve inter-layer (compact-direction) bonds that leave the cluster:

  severed : the defect is out of registry with every neighbour, so the bonds
            are removed. The A_2u channel returns the slice pair {4.697, 8.303}
            exactly, plus a root at 6 from the mixed rotations phi_i4 that no
            state reads. The four-dimensional cluster is the slice cluster with
            extra channels, not shifted ones.
  clamped : the bonds are kept with the exterior held at zero. The pair moves
            to {8.94, 27.06}; the nucleon would weigh about 1050 MeV. Excluded.

Depends on spectral_classifier, proton_first_principles and baryon_mass_modes.
"""
import itertools
import numpy as np
from spectral_classifier import cluster_coord_shell
from proton_first_principles import build_oh_elements, find_perm, class_of
from baryon_mass_modes import OH_CHAR

PAIRS = [(a, b) for a in range(4) for b in range(a + 1, 4)]
NN4 = []
for i, j in itertools.combinations(range(4), 2):
    for si in (1, -1):
        for sj in (1, -1):
            v = np.zeros(4); v[i] = si; v[j] = sj; NN4.append(v / np.sqrt(2))
NN4 = np.array(NN4)


def build(coords, mode):
    n = len(coords)
    du, dp = 4 * n, 6 * n
    M = np.zeros((du + dp, du + dp))
    Uu, Pp = M[:du, :du], M[du:, du:]
    C = np.zeros((dp, du))
    key = {tuple(np.round(c * 4).astype(int)): k for k, c in enumerate(coords)}
    for i in range(n):
        for r in NN4:
            j = key.get(tuple(np.round((coords[i] + r) * 4).astype(int)))
            o = np.outer(r, r)
            if j is not None:
                Uu[4*i:4*i+4, 4*i:4*i+4] += o; Uu[4*i:4*i+4, 4*j:4*j+4] -= o
                Pp[6*i:6*i+6, 6*i:6*i+6] += np.eye(6); Pp[6*i:6*i+6, 6*j:6*j+6] -= np.eye(6)
                for p, (a, b) in enumerate(PAIRS):
                    C[6*i+p, 4*j+b] += r[a] / 2; C[6*i+p, 4*j+a] -= r[b] / 2
            elif mode == 'clamped':
                Uu[4*i:4*i+4, 4*i:4*i+4] += o; Pp[6*i:6*i+6, 6*i:6*i+6] += np.eye(6)
    M[:du, :du] += C.T @ C; M[du:, du:] += np.eye(dp); M[:du, du:] = -C.T; M[du:, :du] = -C
    return M


def channel_roots(M, coords3, chi):
    n = len(coords3); du = 4 * n
    w, V = np.linalg.eigh(M)
    P = np.zeros_like(M)
    for R in build_oh_elements():
        perm = find_perm(R, coords3)
        R4 = np.eye(4); R4[:3, :3] = R
        D = np.zeros_like(M)
        for i in range(n):
            j = perm[i]
            D[4*j:4*j+4, 4*i:4*i+4] = R4
            T = np.zeros((6, 6))
            for p, (a, b) in enumerate(PAIRS):
                for q, (c, d) in enumerate(PAIRS):
                    T[p, q] = R4[a, c] * R4[b, d] - R4[a, d] * R4[b, c]
            D[du+6*j:du+6*j+6, du+6*i:du+6*i+6] = T
        P += chi[class_of(R)] * D
    P /= 48
    return sorted(round(w[k], 4) for k in range(len(w)) if np.linalg.norm(P @ V[:, k]) > 0.5)


if __name__ == "__main__":
    shell3, _ = cluster_coord_shell()
    coords = np.hstack([shell3, np.zeros((len(shell3), 1))])
    for mode in ("severed", "clamped"):
        M = build(coords, mode)
        print(f"{mode:8s}: A_2u roots {channel_roots(M, shell3, OH_CHAR['A_2u'])}, "
              f"A_2g roots {channel_roots(M, shell3, OH_CHAR['A_2g'])}")
