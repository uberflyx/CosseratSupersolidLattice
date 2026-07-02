#!/usr/bin/env python3
"""
integer_blocks.py
=================
The Cosserat singlet eigenvalues are quadratic surds, not numerical outputs.

For any O_h irrep that appears exactly once in the displacement sector and
once in the microrotation sector of a cluster, the coupled Cosserat matrix
restricted to that irrep is a 2x2 symmetric matrix.  In the decoupled
(alpha = 0) eigenbasis the entries are integers fixed by bond counting:

  coordination shell A2u:  [[5, 1], [1, 8]]  ->  lambda = (13 +- sqrt(13))/2
      5 = universal displacement value 4 + one curl unit <C^T C> = 1
      8 = microrotation integer 7 (= cuboctahedron 6 + one centre bond) + 1
      off-diagonal exactly 1 (the curl maps the u-mode to a unit-norm
      phi-pattern)

  coordination shell A2g:  [[1, 1], [1, 6]]  ->  lambda = (7 -+ sqrt(29))/2
      the baryon winding is the lower root (7 - sqrt29)/2 = 0.80742

  Born cluster A1u: both modes are pure microrotation, one on each shell
      (cuboctahedron and NNN octahedron), mixing with trace 9, det 16
      ->  (9 +- sqrt(17))/2, lifted by the coupling to (11 +- sqrt(17))/2.
      The dilaton is the upper root (11 + sqrt17)/2 = 7.56155.

Closed-form masses (m = N m_0 - N (4 - lambda) m_e):
  proton (isoscalar):  m_N = 13 m_0 + (13/2)(5 + sqrt13) m_e = 938.9116 MeV
  dilaton eta(1295):   m   = 18 m_0 + 9 (3 + sqrt17) m_e     = 1293.2136 MeV
"""
import numpy as np
from itertools import permutations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cosserat_classifier import build_cosserat_matrix
from spectral_classifier import fcc_nn_vectors

M_E = 0.51099895069
INV_ALPHA = 137.035999177
M_0 = M_E * INV_ALPHA

CLASSES = ['E', '8C3', '6C4', '6C2', '3C2', 'i', '8S6', '6S4', '6sd', '3sh']
CH = {'A2u': dict(zip(CLASSES, [1, 1, -1, -1, 1, -1, -1, 1, 1, -1])),
      'A2g': dict(zip(CLASSES, [1, 1, -1, -1, 1, 1, 1, -1, -1, 1])),
      'A1u': dict(zip(CLASSES, [1, 1, 1, 1, 1, -1, -1, -1, -1, -1]))}

def oh_elements():
    out = []
    for p in permutations(range(3)):
        for s in np.ndindex(2, 2, 2):
            R = np.zeros((3, 3))
            for row, (col, sg) in enumerate(zip(p, s)):
                R[row, col] = 1 - 2*sg
            out.append(R)
    return out

def classify(R):
    tr, det = round(np.trace(R)), round(np.linalg.det(R))
    diag = np.allclose(R, np.diag(np.diag(R)))
    if det == 1:
        return {3: 'E', 0: '8C3', 1: '6C4'}.get(tr, '3C2' if diag else '6C2')
    return {-3: 'i', 0: '8S6', -1: '6S4'}.get(tr, '3sh' if diag else '6sd')

def cosserat_rep(R, coords):
    n = len(coords)
    perm = []
    for c in coords:
        d = np.linalg.norm(coords - (R @ c), axis=1)
        k = int(np.argmin(d))
        if d[k] > 1e-6:
            return None
        perm.append(k)
    det = np.linalg.det(R)
    M = np.zeros((6*n, 6*n))
    for i, j in enumerate(perm):
        M[3*j:3*j+3, 3*i:3*i+3] = R
        M[3*n+3*j:3*n+3*j+3, 3*n+3*i:3*n+3*i+3] = det * R
    return M

def block(coords, irrep):
    """Coupled Cosserat matrix restricted to the irrep, expressed in the
    decoupled (alpha = 0) eigenbasis."""
    M0 = build_cosserat_matrix(coords, 1.0, 1.0, alpha=0.0)
    M1 = build_cosserat_matrix(coords, 1.0, 1.0, alpha=1.0)
    n = len(coords)
    P = np.zeros((6*n, 6*n))
    G = [(R, cosserat_rep(R, coords)) for R in oh_elements()]
    G = [(R, D) for R, D in G if D is not None]
    for R, D in G:
        P += CH[irrep][classify(R)] * D
    P /= len(G)
    w, v = np.linalg.eigh(P)
    B = v[:, w > 0.5]
    w0, v0 = np.linalg.eigh(B.T @ M0 @ B)
    U = B @ v0
    return w0, U.T @ M1 @ U

def report(name, coords, irrep):
    w0, blk = block(coords, irrep)
    s, det = np.trace(blk), np.linalg.det(blk)
    D = s*s - 4*det
    print(f"{name} {irrep}:")
    print(f"  decoupled eigenvalues: {np.round(w0, 4)}")
    print(f"  coupled block:\n{np.round(blk, 6)}")
    print(f"  trace = {s:.4f}, det = {det:.4f}  ->  "
          f"lambda = ({s:.0f} +- sqrt({D:.0f}))/2 = "
          f"{(s+np.sqrt(D))/2:.6f}, {(s-np.sqrt(D))/2:.6f}\n")

def main():
    shell = np.vstack([[[0., 0., 0.]], fcc_nn_vectors()])
    a = np.sqrt(2.0)
    born = np.vstack([shell, a*np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0],
                                         [0, -1, 0], [0, 0, 1], [0, 0, -1]],
                                        float)])
    report("coordination shell", shell, 'A2u')
    report("coordination shell", shell, 'A2g')
    report("Born cluster", born, 'A1u')

    mN = 13*M_0 + (13/2)*(5 + np.sqrt(13))*M_E
    md = 18*M_0 + 9*(3 + np.sqrt(17))*M_E
    print(f"closed-form proton isoscalar: {mN:.4f} MeV   (PDG 938.9187)")
    print(f"closed-form dilaton eta(1295): {md:.4f} MeV  (PDG 1294)")

if __name__ == "__main__":
    main()
