#!/usr/bin/env python3
"""
sigma_trace.py
==============
Close the Sigma(1660) (J^P = 1/2+) by adiabatic continuation of the Sigma
ground mass mode onto the second-shell cluster.

The positive-parity excited Sigma is a radial excitation: the same A_2
phi-dominant mode that gives the Sigma ground (lambda = 4.624 on the shell
plus four voids, see sigma_ground.py), now read on the cluster extended by
the six second-neighbour nodes (23 nodes total).  The selection rule is the
adiabatic continuation used for the Delta(1600): identify the parent mode
with the extension decoupled, then follow it by eigenvector overlap as the
extension bonds turn on.  Nothing in the trace consults the measured mass.

The parent rises monotonically from lambda = 4.624 to lambda = 7.545 at full
coupling, giving m = 1652.2 MeV against PDG 1660, a residual of -0.47%.

Adiabatic method.  The turn-on parameter tau scales every bond that touches
one of the six new nodes; bonds internal to the 17-node ground cluster are
held fixed.  At tau = 0 the new nodes are decoupled and the parent is the
bare Sigma ground mode; at tau = 1 the matrix is the full two-length
Cosserat operator on 23 nodes.  This is the physical adiabatic turn-on,
distinct from a linear interpolation between the tau = 0 and tau = 1
matrices.  A sanity check confirms matrix_tau(tau = 1) reproduces
build_cosserat_matrix_two_d on the same cluster.

Matrix convention: unit stiffness on every nearest-neighbour (d = 1) and
void-corner (d = sqrt6/4) bond, with length entering only through the
discrete curl, weighted 1/(2d) per bond.

Run:  python3 sigma_trace.py
"""

import numpy as np
from delta_first_principles import (void_positions_Td, build_cosserat_matrix_two_d,
                                     cluster_delta)
from spectral_classifier import fcc_nn_vectors, ELL, A_LAT
from hadron_spectral_mass import mass_from_lambda

VOID_D = np.sqrt(6.0) / 4.0                          # void-to-corner bond length
LEVI_CIVITA = {(0, 1, 2): 1.0, (1, 2, 0): 1.0, (2, 0, 1): 1.0,
               (0, 2, 1): -1.0, (2, 1, 0): -1.0, (1, 0, 2): -1.0}

# Six second-neighbour (NNN) nodes at distance A_LAT along the cube axes.
NNN = np.array([[s * A_LAT if axis == k else 0.0 for axis in range(3)]
                for k in range(3) for s in (+1, -1)])


def matrix_tau(coords, new_idx, tau, alpha=1.0, tol=1e-6):
    """Two-length Cosserat matrix with every bond touching a node in
    new_idx scaled by tau.  At tau = 1 this equals
    build_cosserat_matrix_two_d(coords)."""
    n = len(coords)
    P = np.zeros((6 * n, 6 * n))
    I3 = np.eye(3)
    new = set(new_idx)

    bonds = []
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(coords[i] - coords[j])
            if abs(d - ELL) < tol or abs(d - VOID_D) < tol:
                scale = tau if (i in new or j in new) else 1.0
                bonds.append((i, j, d, scale))

    # Displacement (Puu) and microrotation (Ppp) Laplacian blocks.
    Puu = P[:3 * n, :3 * n]
    Ppp = P[3 * n:, 3 * n:]
    for i, j, d, s in bonds:
        rhat = (coords[j] - coords[i]) / d
        outer = s * np.outer(rhat, rhat)
        Puu[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= outer
        Puu[3 * j:3 * j + 3, 3 * i:3 * i + 3] -= outer
        Puu[3 * i:3 * i + 3, 3 * i:3 * i + 3] += outer
        Puu[3 * j:3 * j + 3, 3 * j:3 * j + 3] += outer
        Ppp[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= s * I3
        Ppp[3 * j:3 * j + 3, 3 * i:3 * i + 3] -= s * I3
        Ppp[3 * i:3 * i + 3, 3 * i:3 * i + 3] += s * I3
        Ppp[3 * j:3 * j + 3, 3 * j:3 * j + 3] += s * I3

    # Discrete curl operator, weighted 1/(2d) per bond.
    C = np.zeros((3 * n, 3 * n))
    for i, j, d, s in bonds:
        rhat = (coords[j] - coords[i]) / d
        for (k, m, sign) in [(i, j, +1), (j, i, -1)]:
            rh = sign * rhat
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        eps = LEVI_CIVITA.get((a, b, c), 0.0)
                        if eps:
                            C[3 * k + a, 3 * m + c] += s * eps * rh[b] / (2.0 * d)

    P[:3 * n, :3 * n] += alpha * (C.T @ C)
    P[3 * n:, 3 * n:] += alpha * np.eye(3 * n)
    P[:3 * n, 3 * n:] = -alpha * C.T
    P[3 * n:, :3 * n] = -alpha * C
    return P


def main():
    centre = np.array([[0.0, 0.0, 0.0]])
    shell = fcc_nn_vectors()
    voids = void_positions_Td()
    ground = np.vstack([centre, shell, voids])               # 17 nodes
    second_shell = np.vstack([centre, shell, voids, NNN])     # 23 nodes
    n_ground = len(ground)
    n = len(second_shell)
    new_idx = list(range(n_ground, n))                        # the 6 NNN nodes

    # Sanity check: matrix_tau at tau = 1 matches the canonical builder.
    diff = np.max(np.abs(matrix_tau(second_shell, new_idx, 1.0)
                         - build_cosserat_matrix_two_d(second_shell)))
    assert diff < 1e-9, f"matrix_tau(1) mismatch: {diff}"

    # Identify the A_2 parent (phi-dominant, localised on the ground cluster)
    # with the second shell decoupled.
    P0 = matrix_tau(second_shell, new_idx, 0.0)
    ev0, vec0 = np.linalg.eigh(P0)
    ev0[np.abs(ev0) < 1e-9] = 0.0
    phi = np.sum(vec0[3 * n:, :] ** 2, axis=0)
    ground_wt = (np.sum(vec0[:3 * n_ground, :] ** 2, axis=0)
                 + np.sum(vec0[3 * n:3 * n + 3 * n_ground, :] ** 2, axis=0))
    candidates = [c for c in range(2 * 3 * n)
                  if 4.0 < ev0[c] < 5.2 and phi[c] > 0.9 and ground_wt[c] > 0.9]
    parent = min(candidates, key=lambda c: abs(ev0[c] - 4.624))
    print(f"Sigma A_2 parent (second shell decoupled): lambda = {ev0[parent]:.3f}, "
          f"phi {100 * phi[parent]:.0f}%, ground weight {ground_wt[parent]:.2f}")

    # Follow the parent by eigenvector overlap as the second shell turns on.
    prev = vec0[:, parent].copy()
    current = parent
    trajectory = []
    for tau in np.linspace(0, 1, 61):
        P = matrix_tau(second_shell, new_idx, tau)
        ev, vec = np.linalg.eigh(P)
        ev[np.abs(ev) < 1e-9] = 0.0
        if tau > 0:
            current = int(np.argmax(np.abs(vec.T @ prev)))
        prev = vec[:, current].copy()
        trajectory.append((tau, ev[current]))

    print("\nAdiabatic trajectory (every 15th step):")
    for tau, lam in trajectory[::15] + [trajectory[-1]]:
        print(f"  tau = {tau:.2f}   lambda = {lam:.3f}   m(N=23) = {mass_from_lambda(23, lam):.1f}")

    lam_final = trajectory[-1][1]
    m = mass_from_lambda(23, lam_final)
    print(f"\nSigma(1660): lambda = {lam_final:.3f}, m = {m:.1f} MeV "
          f"(PDG 1660, {100 * (m - 1660) / 1660:+.2f}%)")


if __name__ == "__main__":
    main()
