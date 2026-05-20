#!/usr/bin/env python3
"""
backward_trace.py
=================

For the N(1535) on proton+bilayer (N=21):
The chapter's claim is m = 1513.95 MeV at lambda = 8.046, characterised
as a shell-concentrated stiff phi-dominant singlet.

My forward adiabatic check showed that psi_p (lambda = 8.303 on the bare
13-shell) traces to lambda = 11.55 (m = 1551 MeV), bilayer-localised
73%.  This is NOT the chapter's N(1535).

So which bare-shell mode IS the parent of lambda = 8.046?  Trace
backwards from tau=1 to tau=0 and report the starting eigenvalue.

We also do the same for Roper (proton + second-shell, N=19) to verify
whether psi_p truly descends to lambda=10.193 (chapter) or 10.48 (my
forward tracking).
"""

import numpy as np
import scipy.linalg
from scipy.optimize import linear_sum_assignment
import sys

from hadron_spectral_mass import mass_from_lambda
from spectral_classifier import cluster_coord_shell
from cosserat_classifier import build_cosserat_matrix
from n1535_first_principles import build_n1535_cluster


def add_second_shell(cluster):
    s = np.sqrt(2.0)
    extra = np.array([
        [+s, 0, 0], [-s, 0, 0],
        [0, +s, 0], [0, -s, 0],
        [0, 0, +s], [0, 0, -s],
    ])
    return np.vstack([cluster, extra])


def trace_all(coords, n_inner, matrix_builder, n_steps=81):
    """Run adiabatic interpolation; return (initial_vals, final_vals, perm)
    where perm[i] = j means initial index i ended up at final index j."""
    N = len(coords)
    n_dim = 6 * N

    H_full = matrix_builder(coords)
    H_inner = matrix_builder(coords[:n_inner])
    H_partial = np.zeros((6 * N, 6 * N))
    H_partial[:3*n_inner, :3*n_inner] = H_inner[:3*n_inner, :3*n_inner]
    H_partial[3*N:3*N+3*n_inner, 3*N:3*N+3*n_inner] = H_inner[3*n_inner:, 3*n_inner:]
    H_partial[:3*n_inner, 3*N:3*N+3*n_inner] = H_inner[:3*n_inner, 3*n_inner:]
    H_partial[3*N:3*N+3*n_inner, :3*n_inner] = H_inner[3*n_inner:, :3*n_inner]

    def H_of_tau(tau):
        return (1 - tau) * H_partial + tau * H_full

    taus = np.linspace(0.0, 1.0, n_steps)
    H_prev = H_of_tau(taus[0])
    vals_prev, vecs_prev = scipy.linalg.eigh(H_prev)
    initial_vals = vals_prev.copy()
    initial_vecs = vecs_prev.copy()

    cumulative_perm = np.arange(n_dim)
    for step_idx in range(1, n_steps):
        H = H_of_tau(taus[step_idx])
        vals, vecs = scipy.linalg.eigh(H)
        overlap = np.abs(vecs_prev.T @ vecs) ** 2
        row_ind, col_ind = linear_sum_assignment(-overlap)
        perm = np.zeros(n_dim, dtype=int)
        for i, j in zip(row_ind, col_ind):
            perm[i] = j
        cumulative_perm = perm[cumulative_perm]
        vals_prev = vals[cumulative_perm]
        vecs_prev = vecs[:, cumulative_perm]

    final_vals = vals_prev
    final_vecs = vecs_prev

    return initial_vals, final_vals, initial_vecs, final_vecs


def describe(eigvec, N, n_inner):
    """Phi fraction, inner localisation."""
    u = np.sum(eigvec[:3*N]**2)
    phi = np.sum(eigvec[3*N:]**2)
    phi_frac = phi / (u + phi)
    amps = np.zeros(N)
    for i in range(N):
        amps[i] = (np.sum(eigvec[3*i:3*(i+1)]**2) +
                   np.sum(eigvec[3*N+3*i:3*N+3*(i+1)]**2))
    inner = amps[:n_inner].sum() / amps.sum()
    return phi_frac, inner


def main():
    print("=" * 78)
    print("Backward-trace audit of N(1535) and Roper")
    print("=" * 78)

    # ------------------------------------------------------------
    # Case 1: N(1535) on proton+bilayer (N=21)
    # ------------------------------------------------------------
    print("\n" + "=" * 78)
    print("N(1535) on proton+bilayer (N=21)")
    print("=" * 78)

    coords = build_n1535_cluster()
    N = len(coords)
    initial_vals, final_vals, initial_vecs, final_vecs = trace_all(
        coords, n_inner=13,
        matrix_builder=lambda c: build_cosserat_matrix(c, K_u=1.0, K_phi=1.0, alpha=1.0),
        n_steps=81
    )

    print(f"\nFinal spectrum near chapter's N(1535) claim (lambda = 8.046):")
    print(f"{'idx':>4} {'init lambda':>12} {'final lambda':>13} "
          f"{'m_final':>9} {'phi%':>5} {'inner%':>7}")
    for i in range(len(final_vals)):
        if 7.9 < final_vals[i] < 8.5:
            phi, inner = describe(final_vecs[:, i], N, 13)
            m = mass_from_lambda(N, final_vals[i])
            tag = ""
            if abs(final_vals[i] - 8.046) < 0.01:
                tag = " <-- chapter's N(1535)"
            print(f"  {i:>3} {initial_vals[i]:>11.4f} {final_vals[i]:>12.4f} "
                  f"{m:>9.2f} {phi*100:>4.1f}  {inner*100:>5.1f}% {tag}")

    print(f"\nWhere does the bare proton's psi_p (lambda=8.303) descend to?")
    # Find initial index closest to 8.303 with phi > 90% and shell concentration
    for i in range(len(initial_vals)):
        if abs(initial_vals[i] - 8.303) < 0.01:
            phi_init, inner_init = describe(initial_vecs[:, i], N, 13)
            phi_fin, inner_fin = describe(final_vecs[:, i], N, 13)
            m = mass_from_lambda(N, final_vals[i])
            print(f"  i={i}: lambda(tau=0) = {initial_vals[i]:.4f} -> "
                  f"lambda(tau=1) = {final_vals[i]:.4f}, m={m:.2f}")
            print(f"     phi: {phi_init*100:.0f}% -> {phi_fin*100:.0f}%, "
                  f"inner: {inner_init*100:.0f}% -> {inner_fin*100:.0f}%")

    print(f"\nWhere do the bare T_1g stiff phi-dom triplet "
          f"(lambda~8.07) descend to?")
    for i in range(len(initial_vals)):
        if 8.05 < initial_vals[i] < 8.10:
            phi_init, inner_init = describe(initial_vecs[:, i], N, 13)
            phi_fin, inner_fin = describe(final_vecs[:, i], N, 13)
            m = mass_from_lambda(N, final_vals[i])
            print(f"  i={i}: lambda(tau=0) = {initial_vals[i]:.4f} -> "
                  f"lambda(tau=1) = {final_vals[i]:.4f}, m={m:.2f}, "
                  f"final phi={phi_fin*100:.0f}%, inner={inner_fin*100:.0f}%")

    # ------------------------------------------------------------
    # Case 2: Roper on proton + second-shell (N=19)
    # ------------------------------------------------------------
    print("\n" + "=" * 78)
    print("N(1440) Roper on proton+second-shell (N=19)")
    print("=" * 78)

    shell, _ = cluster_coord_shell()
    coords = add_second_shell(shell)
    N = len(coords)
    initial_vals, final_vals, initial_vecs, final_vecs = trace_all(
        coords, n_inner=13,
        matrix_builder=lambda c: build_cosserat_matrix(c, K_u=1.0, K_phi=1.0, alpha=1.0),
        n_steps=81
    )

    print(f"\nFinal spectrum near chapter's Roper claim (lambda = 10.193):")
    print(f"{'idx':>4} {'init lambda':>12} {'final lambda':>13} "
          f"{'m_final':>9} {'phi%':>5} {'inner%':>7}")
    for i in range(len(final_vals)):
        if 10.0 < final_vals[i] < 10.6:
            phi, inner = describe(final_vecs[:, i], N, 13)
            m = mass_from_lambda(N, final_vals[i])
            tag = ""
            if abs(final_vals[i] - 10.193) < 0.01:
                tag = " <-- chapter's Roper"
            print(f"  {i:>3} {initial_vals[i]:>11.4f} {final_vals[i]:>12.4f} "
                  f"{m:>9.2f} {phi*100:>4.1f}  {inner*100:>5.1f}% {tag}")

    print(f"\nWhere does the bare proton's psi_p (lambda=8.303) descend to?")
    for i in range(len(initial_vals)):
        if abs(initial_vals[i] - 8.303) < 0.01:
            phi_init, inner_init = describe(initial_vecs[:, i], N, 13)
            phi_fin, inner_fin = describe(final_vecs[:, i], N, 13)
            m = mass_from_lambda(N, final_vals[i])
            print(f"  i={i}: lambda(tau=0) = {initial_vals[i]:.4f} -> "
                  f"lambda(tau=1) = {final_vals[i]:.4f}, m={m:.2f}")
            print(f"     phi: {phi_init*100:.0f}% -> {phi_fin*100:.0f}%, "
                  f"inner: {inner_init*100:.0f}% -> {inner_fin*100:.0f}%")


if __name__ == '__main__':
    main()
