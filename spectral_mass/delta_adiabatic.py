#!/usr/bin/env python3
"""
delta_adiabatic.py
==================

Adiabatic interpolation from Delta(1232) cluster (N=17, T_d) to the
dual-orbit Delta(1600) cluster (N=21, O_h).

Parameter tau in [0, 1] scales the bonds involving the 4 "new" voids
(negative T_d orbit) by tau.  At tau=0, the new voids are decoupled
(their 24 d.o.f. are zero modes; the 17-atom Delta spectrum is recovered).
At tau=1, the full dual-orbit Cosserat matrix is reached.

Track each eigenvalue continuously by maximum eigenvector overlap between
adjacent tau steps.  Identify which eigenvalue trajectory descends from
the Delta(1232)'s T_1 mass mode at lambda=9.05.

The question: does the parent T_1 mode trace to lambda=7.17 (gerade T_1g
at tau=1), lambda=22.36 (shell-localised T_1g at tau=1), or somewhere
else?  This is the principled answer to which mode is the Delta(1600).
"""

import numpy as np
import scipy.linalg
import sys

from hadron_spectral_mass import mass_from_lambda, M_0, M_E
from spectral_classifier import cluster_coord_shell
from delta_first_principles import build_cosserat_matrix_two_d


def build_dual_orbit_cluster():
    """13 coord shell + 8 tetrahedral voids (both T_d orbits)."""
    shell, _ = cluster_coord_shell()
    s = np.sqrt(2.0) / 4.0
    voids_plus = np.array([
        [+s, +s, +s], [+s, -s, -s], [-s, +s, -s], [-s, -s, +s]
    ])
    voids_minus = np.array([
        [-s, -s, -s], [-s, +s, +s], [+s, -s, +s], [+s, +s, -s]
    ])
    return np.vstack([shell, voids_plus, voids_minus])


def build_interpolating_matrix(coords, tau, n_old=17):
    """Build Cosserat dynamical matrix with bonds to new atoms (indices
    >= n_old) scaled by tau.

    At tau=0: bonds within atoms 0..n_old-1 active; atoms n_old..N-1
              decoupled (zero modes).
    At tau=1: full dual-orbit Cosserat matrix.

    We accomplish this by placing the new atoms at their dual-orbit
    positions but scaling the spring constants of any bond that touches
    one of them by tau.
    """
    N = len(coords)
    # Get the unscaled Cosserat matrix at tau=1 (full coupling)
    H_full = build_cosserat_matrix_two_d(coords, K_u=1.0, K_phi=1.0, alpha=1.0)

    # Get the Cosserat matrix at tau=0: cluster of only the first n_old
    # atoms, padded with zero rows/columns for the decoupled new atoms.
    H_partial_inner = build_cosserat_matrix_two_d(
        coords[:n_old], K_u=1.0, K_phi=1.0, alpha=1.0
    )
    # Embed in the 6N-dim space.  The Cosserat matrix is structured as
    # 3N u-block + 3N phi-block.  For each block we copy the n_old-block
    # submatrix into the right positions.
    H_partial = np.zeros((6 * N, 6 * N))
    # u block: H_partial_inner has 3n_old rows for u and 3n_old for phi.
    # The layout depends on how build_cosserat_matrix_two_d orders things.
    # Standard convention: [u_0, u_1, ..., u_{n-1}, phi_0, ..., phi_{n-1}],
    # 3 components each.
    # So u-u block: H_inner[:3n_old, :3n_old] -> H_partial u-u positions
    # phi-phi block: H_inner[3n_old:, 3n_old:] -> H_partial phi-phi
    # u-phi block: H_inner[:3n_old, 3n_old:] -> H_partial u-phi cross
    # The new atoms (indices n_old..N-1) get zero entries.
    H_partial[:3*n_old, :3*n_old] = H_partial_inner[:3*n_old, :3*n_old]  # u-u
    H_partial[3*N:3*N+3*n_old, 3*N:3*N+3*n_old] = H_partial_inner[3*n_old:, 3*n_old:]  # phi-phi
    H_partial[:3*n_old, 3*N:3*N+3*n_old] = H_partial_inner[:3*n_old, 3*n_old:]  # u-phi
    H_partial[3*N:3*N+3*n_old, :3*n_old] = H_partial_inner[3*n_old:, :3*n_old]  # phi-u

    # Interpolate: H(tau) = (1-tau) H_partial + tau H_full
    # This isn't quite "scale bonds touching new atoms by tau" but it's
    # a clean smooth interpolation with the correct endpoints.
    return (1 - tau) * H_partial + tau * H_full


def per_atom_amplitude(eigvec, n_atoms):
    """Total amplitude per atom |u_i|^2 + |phi_i|^2."""
    amps = np.zeros(n_atoms)
    for i in range(n_atoms):
        u = np.sum(eigvec[3*i:3*(i+1)]**2)
        phi = np.sum(eigvec[3*n_atoms+3*i:3*n_atoms+3*(i+1)]**2)
        amps[i] = u + phi
    return amps


def main():
    print("=" * 78)
    print("Adiabatic interpolation: Delta(1232) cluster -> dual-orbit cluster")
    print("=" * 78)

    coords = build_dual_orbit_cluster()
    N = len(coords)
    print(f"\nCluster: N = {N} (17 single-orbit + 4 new voids)")
    print(f"  Atoms 0..12: inner shell")
    print(f"  Atoms 13..16: positive T_d orbit voids (active at all tau)")
    print(f"  Atoms 17..20: negative T_d orbit voids (gated by tau)")

    # Reference: spectrum at tau=0 (should reproduce Delta(1232) spectrum
    # plus 24 zero modes for the decoupled atoms)
    print("\n--- tau = 0 (decoupled new voids) ---")
    H0 = build_interpolating_matrix(coords, tau=0.0)
    vals0, vecs0 = scipy.linalg.eigh(H0)
    print(f"Spectrum: {len(vals0)} eigenvalues")
    # The 24 lowest should be zero modes; the rest match Delta(1232) spectrum
    n_zero = np.sum(np.abs(vals0) < 1e-6)
    print(f"Near-zero eigenvalues: {n_zero} (expected 24 from decoupled atoms)")
    # Find Delta T_1 mode: should be near lambda=9.05
    candidates_delta = np.where((vals0 > 8.8) & (vals0 < 9.3))[0]
    print(f"Eigenvalues in [8.8, 9.3] (parent T_1): {len(candidates_delta)}")
    for k in candidates_delta:
        # Check support on inner 17 atoms
        amps = per_atom_amplitude(vecs0[:, k], N)
        inner = amps[:17].sum() / amps.sum()
        print(f"  k={k}: lambda={vals0[k]:.4f}, inner-17 weight={inner*100:.1f}%")

    # Track eigenvalues through tau
    print("\n--- Adiabatic tracking from tau=0 to tau=1 ---")
    n_steps = 41
    taus = np.linspace(0.0, 1.0, n_steps)

    # We'll track ALL eigenvalues but report only those that pass through
    # the parent T_1 mode.  Use eigenvector overlap to match indices.

    H_prev = build_interpolating_matrix(coords, tau=taus[0])
    vals_prev, vecs_prev = scipy.linalg.eigh(H_prev)

    # Initial labels: 0..(6N-1)
    n_dim = 6 * N
    labels = np.arange(n_dim)

    # Store trajectories: trajectory[label] = list of (tau, lambda) pairs
    trajectories = {i: [(taus[0], vals_prev[i])] for i in range(n_dim)}
    # Also store final eigenvector for each label
    eigenvectors = {i: vecs_prev[:, i].copy() for i in range(n_dim)}

    for step_idx in range(1, n_steps):
        tau = taus[step_idx]
        H = build_interpolating_matrix(coords, tau=tau)
        vals, vecs = scipy.linalg.eigh(H)

        # Match each new eigenvector to the previous step's eigenvectors
        # by maximum overlap.  Use Hungarian-style assignment.
        # Overlap matrix: O[i,j] = |<vecs_prev[:,i], vecs[:,j]>|^2
        overlap = np.abs(vecs_prev.T @ vecs) ** 2

        # Greedy matching: for each row (old eigenvector), find best col
        # (new eigenvector); use cost = -overlap so minimum-cost matching
        # = maximum overlap.
        # Use scipy.optimize.linear_sum_assignment for proper assignment.
        from scipy.optimize import linear_sum_assignment
        row_ind, col_ind = linear_sum_assignment(-overlap)

        # Permutation: new index col_ind[i] corresponds to old index row_ind[i]
        # Build permutation: perm[i] = j means new column j is old column i
        perm = np.zeros(n_dim, dtype=int)
        for i, j in zip(row_ind, col_ind):
            perm[i] = j
        # Reorder vals and vecs so that index i is the continuation of
        # previous index i
        vals_ordered = vals[perm]
        vecs_ordered = vecs[:, perm]

        # Update trajectories
        for i in range(n_dim):
            trajectories[i].append((tau, vals_ordered[i]))
            eigenvectors[i] = vecs_ordered[:, i].copy()

        vals_prev = vals_ordered
        vecs_prev = vecs_ordered

    # Identify the trajectories that started near the parent T_1 mode
    print(f"\nFinal spectrum at tau=1: lambda range "
          f"[{vals_prev[0]:.4f}, {vals_prev[-1]:.4f}]")

    print(f"\n--- Trajectories descending from parent T_1 (lambda ~ 9.05) ---")
    parent_indices = []
    for i in range(n_dim):
        traj = trajectories[i]
        initial_lambda = traj[0][1]
        if 8.8 < initial_lambda < 9.3:
            parent_indices.append(i)

    print(f"Found {len(parent_indices)} indices starting at lambda ~ 9.05")
    print(f"\n{'idx':>4} {'tau=0':>8} {'tau=0.25':>10} {'tau=0.5':>10} "
          f"{'tau=0.75':>10} {'tau=1':>10} {'final m':>10}")
    for i in parent_indices:
        traj = trajectories[i]
        # Sample at specific tau values
        lambda_0 = traj[0][1]
        lambda_25 = traj[int(0.25 * (n_steps - 1))][1]
        lambda_50 = traj[int(0.50 * (n_steps - 1))][1]
        lambda_75 = traj[int(0.75 * (n_steps - 1))][1]
        lambda_1 = traj[-1][1]
        m_final = mass_from_lambda(N, lambda_1)
        print(f"  {i:>3} {lambda_0:>8.4f} {lambda_25:>10.4f} "
              f"{lambda_50:>10.4f} {lambda_75:>10.4f} {lambda_1:>10.4f} "
              f"{m_final:>10.2f}")

    # Final eigenvector characterisation
    print(f"\n--- Final character of descendant modes ---")
    print(f"{'idx':>4} {'lambda_final':>12} {'phi%':>5} {'inner17%':>9} "
          f"{'new4%':>6}")
    for i in parent_indices:
        v = eigenvectors[i]
        # phi fraction
        u_tot = np.sum(v[:3*N]**2)
        phi_tot = np.sum(v[3*N:]**2)
        phi_frac = phi_tot / (u_tot + phi_tot)
        # localisation: inner 17 vs new 4
        amps = per_atom_amplitude(v, N)
        inner17 = amps[:17].sum() / amps.sum()
        new4 = amps[17:].sum() / amps.sum()
        print(f"  {i:>3} {trajectories[i][-1][1]:>12.4f} "
              f"{phi_frac*100:>4.1f} "
              f"{inner17*100:>8.1f}% {new4*100:>5.1f}%")

    # Summary
    print(f"\n--- Summary: where does the parent's T_1 mode go? ---")
    sorted_indices = sorted(parent_indices,
                            key=lambda i: trajectories[i][-1][1])
    for i in sorted_indices:
        lam = trajectories[i][-1][1]
        m = mass_from_lambda(N, lam)
        print(f"  lambda_final = {lam:7.4f}, m = {m:.2f} MeV")


if __name__ == '__main__':
    main()
