#!/usr/bin/env python3
"""
clamped_frame_curl.py
=====================
Why the Cosserat coupling's discrete curl carries no self-term, and what
happens if one is added.

The coupling energy is (alpha/2) |phi_i - omega_i|^2, with omega_i the local
rotation of the displacement field.  The exact lattice curl at site i is

    omega_i = (1/2l) sum_{j in NN(i)} rhat_ij x (u_j - u_i).

In the infinite FCC lattice every site has twelve nearest neighbours whose
bond directions sum to zero, so the self-term -(sum rhat) x u_i vanishes
identically.  Restricting to a finite defect cluster with the environment
held at u = 0 therefore gives

    omega_i = (1/2l) sum_{j in cluster NN(i)} rhat_ij x u_j,

with no self-term: this is the CLAMPED-FRAME convention used throughout the
framework.  Physically, the defect's nodes measure their microrotation
against the ambient crystal's orientation frame.  The defect severs its
translational bonds to the environment (that cost is the N m_0 assembly
term), but the orientational reference of the surrounding crystal survives.

The alternative FREE-CLUSTER convention keeps the self-term restricted to
in-cluster neighbours, so a rigid translation of the cluster costs nothing.
This script computes both on the coordination shell:

  clamped frame:  A2u = {(13-sqrt13)/2, (13+sqrt13)/2} = {4.697, 8.303}
                  A2g = {(7-sqrt29)/2, (7+sqrt29)/2}   = {0.807, 6.193}
                  -> proton isoscalar 938.912 MeV (obs 938.919, -7 keV)
  free cluster:   A2u = {4.185, 8.066}, A2g = {3.622, 8.628}
                  -> proton 937.34 MeV (misses by 1.6 MeV) and no soft
                     baryon-winding mode near zero.

The data therefore select the clamped frame at the 6 keV level: the defect
keeps the crystal's rotational frame.
"""
import numpy as np
from itertools import permutations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cosserat_classifier import (build_uu_block, build_phiphi_block, ELL)
from spectral_classifier import fcc_nn_vectors

M_E = 0.51099895069
INV_ALPHA = 137.035999177
M_0 = M_E * INV_ALPHA

def curl_operator(coords, include_self, ell=ELL, tol=1e-6):
    n = len(coords)
    C = np.zeros((3*n, 3*n))
    eps = np.zeros((3, 3, 3))
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1
    eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            r = coords[j] - coords[i]
            d = np.linalg.norm(r)
            if abs(d - ell) > tol:
                continue
            rh = r / d
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        if eps[a, b, c]:
                            C[3*i+a, 3*j+c] += eps[a, b, c]*rh[b]/(2*ell)
                            if include_self:
                                C[3*i+a, 3*i+c] -= eps[a, b, c]*rh[b]/(2*ell)
    return C

def cosserat(coords, include_self, alpha=1.0):
    n = len(coords)
    M = np.zeros((6*n, 6*n))
    M[:3*n, :3*n] = build_uu_block(coords)
    M[3*n:, 3*n:] = build_phiphi_block(coords)
    C = curl_operator(coords, include_self)
    M[:3*n, :3*n] += alpha * C.T @ C
    M[3*n:, 3*n:] += alpha * np.eye(3*n)
    M[:3*n, 3*n:] -= alpha * C.T
    M[3*n:, :3*n] -= alpha * C
    return M

def irrep_roots(coords, M, chi):
    from integer_blocks import oh_elements, classify, cosserat_rep
    n = len(coords)
    P = np.zeros((6*n, 6*n))
    G = [(R, cosserat_rep(R, coords)) for R in oh_elements()]
    G = [(R, D) for R, D in G if D is not None]
    for R, D in G:
        P += chi[classify(R)] * D
    P /= len(G)
    w, v = np.linalg.eigh(P)
    B = v[:, w > 0.5]
    return np.round(np.linalg.eigvalsh(B.T @ M @ B), 4)

def main():
    # (1) bond directions of a bulk FCC site sum to zero
    nn = fcc_nn_vectors()
    print("Sum of the twelve FCC bond directions:", np.round(nn.sum(axis=0), 12))

    shell = np.vstack([[[0., 0., 0.]], fcc_nn_vectors()])
    CLASSES = ['E', '8C3', '6C4', '6C2', '3C2', 'i', '8S6', '6S4', '6sd', '3sh']
    A2U = dict(zip(CLASSES, [1, 1, -1, -1, 1, -1, -1, 1, 1, -1]))
    A2G = dict(zip(CLASSES, [1, 1, -1, -1, 1, 1, 1, -1, -1, 1]))

    for label, self_term in [("clamped frame (no self-term)", False),
                             ("free cluster (with self-term)", True)]:
        M = cosserat(shell, self_term)
        a2u = irrep_roots(shell, M, A2U)
        a2g = irrep_roots(shell, M, A2G)
        mN = 13*M_0 - 13*(4 - a2u[-1])*M_E
        # rigid-translation energy
        t = np.zeros(6*13)
        for i in range(13):
            t[3*i] = 1
        t /= np.linalg.norm(t)
        print(f"\n{label}:")
        print(f"  A2u roots {a2u}, A2g roots {a2g}")
        print(f"  proton isoscalar -> {mN:.3f} MeV  (obs 938.919)")
        print(f"  rigid-translation energy <t|M|t> = {t @ M @ t:.4f}")

if __name__ == "__main__":
    main()
