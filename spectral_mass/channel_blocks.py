#!/usr/bin/env python3
"""
channel_blocks.py
=================
The one-dimensional channels of an O_h host as exact 2x2 blocks.

For a one-dimensional irrep that appears once in the displacement sector and
once in the microrotation sector, O_h-equivariance of the discrete curl C forces
C u = c phi with u the displacement eigenvector (Phi_uu u = d u) and phi the
unique twist pattern of the irrep, a graph-Laplacian eigenvector (T phi = mu phi).
The coupled Cosserat matrix then closes on the pair,

    M = [[d + alpha c^2, alpha c], [alpha c, mu + alpha]],

and its two roots are the channel's coupled eigenvalues exactly. This script
extracts d, mu and c for the A_2g and A_2u channels of the coordination shell
and of the Born cluster, checks that C u is a Laplacian eigenvector, and
compares the block's roots with the full coupled spectrum. On every channel
c^2 = 1 to machine precision, so the blocks are the integer matrices
[[1,1],[1,6]], [[5,1],[1,8]] (shell) and [[3,1],[1,8]], [[5,1],[1,10]] (Born),
with roots (7 +/- sqrt29)/2, (13 +/- sqrt13)/2, (11 +/- sqrt29)/2 and
(15 +/- sqrt29)/2.
"""
import numpy as np
from delta_first_principles import build_cosserat_matrix_two_d
from spectral_classifier import cluster_coord_shell
from proton_first_principles import build_oh_elements, find_perm, build_rep, class_of
from baryon_mass_modes import OH_CHAR

OH = build_oh_elements()
A1G = {'E':1,'8C_3':1,'6C_4':1,'6C_2':1,'3C_2':1,'i':1,'8S_6':1,'6S_4':1,'3σ_h':1,'6σ_d':1}
A1U = {k: (v if k in ('E','8C_3','6C_4','6C_2','3C_2') else -v) for k, v in A1G.items()}
CHANNELS = {'A_1g': A1G, 'A_2g': OH_CHAR['A_2g'], 'A_1u': A1U, 'A_2u': OH_CHAR['A_2u']}


def projector(chi, coords):
    n = len(coords)
    P = np.zeros((6 * n, 6 * n))
    for R in OH:
        P += chi[class_of(R)] * build_rep(R, find_perm(R, coords), n)
    return P / len(OH)


def analyse(coords, label, alpha=1.0):
    n = len(coords)
    M = build_cosserat_matrix_two_d(coords, alpha=alpha)
    M0 = build_cosserat_matrix_two_d(coords, alpha=0.0)
    Tuu, Tpp = M0[:3 * n, :3 * n], M0[3 * n:, 3 * n:]
    C = -M[3 * n:, :3 * n] / alpha
    w, V = np.linalg.eigh(M)
    print(f"--- {label}")
    for name, chi in CHANNELS.items():
        P = projector(chi, coords)
        Pu = P[:3 * n, :3 * n]
        k = int(round(np.trace(Pu)))
        if k == 0:
            print(f"{name}: no displacement mode")
            continue
        B = np.linalg.svd(Pu)[0][:, :k]
        dvals, dvecs = np.linalg.eigh(B.T @ Tuu @ B)
        U = B @ dvecs
        coupled = sorted(set(round(w[i], 4) for i in range(len(w)) if np.linalg.norm(P @ V[:, i]) > 0.5))
        for j in range(k):
            u = U[:, j]
            cu = C @ u
            c = np.linalg.norm(cu)
            if c < 1e-9:
                print(f"{name}: d = {dvals[j]:.4f}, curl-free (single-sector channel); coupled roots {coupled}")
                continue
            phi = cu / c
            mu = phi @ Tpp @ phi
            resid = np.linalg.norm(Tpp @ phi - mu * phi)
            block = np.array([[dvals[j] + alpha * c * c, alpha * c], [alpha * c, mu + alpha]])
            roots = np.linalg.eigvalsh(block)
            print(f"{name}: d = {dvals[j]:.4f}  mu = {mu:.4f}  c^2 = {c*c:.6f}  |T C u - mu C u| = {resid:.1e}"
                  f"  block roots {np.round(roots, 4)}  coupled roots {coupled}")


if __name__ == "__main__":
    shell, _ = cluster_coord_shell()
    analyse(shell, "coordination shell (13 nodes)")
    s2 = np.sqrt(2.0)
    born = np.vstack([shell, [[s2,0,0],[-s2,0,0],[0,s2,0],[0,-s2,0],[0,0,s2],[0,0,-s2]]])
    analyse(born, "Born cluster (19 nodes)")
