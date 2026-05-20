#!/usr/bin/env python3
"""
parity_flip_rule.py
===================

Test the hypothesis that the mass-mode irrep has parity OPPOSITE to the
hadron's spinor parity.  For each baryon, trace candidate bare-shell
irreps forward through the relevant cluster perturbation and verify
which lands at the right mass.

Tests:
  - Proton (1/2+): A_2u stiff phi-dom at lambda=8.303. No perturbation.
  - N(1535) (1/2-): A_2g stiff phi-dom at lambda=6.193 traced through
    bilayer to lambda~8.046 (chapter claim, m=1514).
  - Delta(1232) (3/2+): T_1u (8.054), T_2u (8.036), T_1g (8.071), or
    T_2g (8.188) traced through single-orbit voids to lambda~9.05.
  - N(1520) (3/2-): which gerade triplet (T_1g lower at 6.159, T_1g
    upper at 8.071, T_2g lower at 6.155, T_2g upper at 8.188, or
    T_1g very stiff at 14.120) traces to lambda~8.32 on bilayer?
"""

import numpy as np
import scipy.linalg
from scipy.optimize import linear_sum_assignment
import sys

from hadron_spectral_mass import mass_from_lambda
from spectral_classifier import cluster_coord_shell
from cosserat_classifier import build_cosserat_matrix
from n1535_first_principles import build_n1535_cluster
from delta_first_principles import build_cosserat_matrix_two_d
from delta1600_dual_orbit import (
    build_oh_elements, build_cosserat_oh_action,
    project_onto_irrep, O_H_IRREPS, CLASS_ORDERS, G_ORDER
)


def build_single_orbit_delta_cluster():
    """13 inner shell + 4 voids of one T_d orbit (positive orbit)."""
    shell, _ = cluster_coord_shell()
    s = np.sqrt(2.0) / 4.0
    voids = np.array([
        [+s, +s, +s], [+s, -s, -s], [-s, +s, -s], [-s, -s, +s]
    ])
    return np.vstack([shell, voids])


def trace_all(coords, n_inner, matrix_builder, n_steps=81):
    """Run adiabatic interpolation; return all eigenvalue trajectories."""
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

    return initial_vals, vals_prev, initial_vecs, vecs_prev


def describe(eigvec, N, n_inner):
    u = np.sum(eigvec[:3*N]**2)
    phi = np.sum(eigvec[3*N:]**2)
    phi_frac = phi / (u + phi)
    amps = np.zeros(N)
    for i in range(N):
        amps[i] = (np.sum(eigvec[3*i:3*(i+1)]**2) +
                   np.sum(eigvec[3*N+3*i:3*N+3*(i+1)]**2))
    inner = amps[:n_inner].sum() / amps.sum()
    return phi_frac, inner


def assign_irrep(eigvec, N, projectors, irrep_names):
    """Return the irrep with highest projection content."""
    max_name = None
    max_content = 0
    for name in irrep_names:
        Pv = projectors[name] @ eigvec
        content = np.dot(eigvec, Pv) / np.dot(eigvec, eigvec)
        if content > max_content:
            max_content = content
            max_name = name
    return max_name, max_content


def main():
    print("=" * 78)
    print("Parity-flip rule test for baryon mass modes")
    print("=" * 78)

    # ============================================================
    # Test 1: N(1535) - confirm A_2g parent traces to lambda=8.046
    # ============================================================
    print("\n" + "-" * 78)
    print("TEST 1: N(1535) (J^P = 1/2-) - parity-flip rule predicts A_2g parent")
    print("-" * 78)

    coords = build_n1535_cluster()
    N = len(coords)
    print(f"Cluster: N = {N} (proton + bilayer)")
    print(f"Bare-shell A_2g stiff phi-dom: lambda = 6.193")
    print(f"Bare-shell A_2u stiff phi-dom (psi_p): lambda = 8.303")

    initial_vals, final_vals, initial_vecs, final_vecs = trace_all(
        coords, n_inner=13,
        matrix_builder=lambda c: build_cosserat_matrix(c, 1.0, 1.0, 1.0),
        n_steps=81
    )

    # Identify the A_2g singlet at initial lambda = 6.193
    print(f"\nForward trace results:")
    for i in range(len(initial_vals)):
        if 6.18 < initial_vals[i] < 6.20:
            final_phi, final_inner = describe(final_vecs[:, i], N, 13)
            m = mass_from_lambda(N, final_vals[i])
            tag = ""
            if abs(final_vals[i] - 8.046) < 0.05:
                tag = "  <-- MATCHES chapter N(1535)"
            print(f"  A_2g (lambda_init = {initial_vals[i]:.4f}) -> "
                  f"lambda_final = {final_vals[i]:.4f}, m = {m:.2f} MeV "
                  f"(phi={final_phi*100:.0f}%, "
                  f"shell={final_inner*100:.0f}%){tag}")
        if 8.30 < initial_vals[i] < 8.31:
            final_phi, final_inner = describe(final_vecs[:, i], N, 13)
            m = mass_from_lambda(N, final_vals[i])
            print(f"  A_2u/psi_p (lambda_init = {initial_vals[i]:.4f}) -> "
                  f"lambda_final = {final_vals[i]:.4f}, m = {m:.2f} MeV "
                  f"(phi={final_phi*100:.0f}%, shell={final_inner*100:.0f}%)")

    print(f"\nPDG N(1535) mass: 1510 +/- 20 MeV (range 1490-1530)")

    # ============================================================
    # Test 2: Delta(1232) - which T-type irrep is the parent?
    # ============================================================
    print("\n" + "-" * 78)
    print("TEST 2: Delta(1232) (J^P = 3/2+) - which T-type irrep is the parent?")
    print("-" * 78)

    coords = build_single_orbit_delta_cluster()
    N = len(coords)
    print(f"Cluster: N = {N} (proton + 4 T_d voids of one orbit)")
    print(f"Bare-shell T-type stiff phi-dom candidates:")
    print(f"  T_2u (lambda=8.036, parity u): predicted by 'parity-flip' for J=3/2+")
    print(f"  T_1u (lambda=8.054, parity u): predicted by 'parity-flip' for J=3/2+")
    print(f"  T_1g (lambda=8.071, parity g): chapter's claim")
    print(f"  T_2g (lambda=8.188, parity g)")

    initial_vals, final_vals, initial_vecs, final_vecs = trace_all(
        coords, n_inner=13,
        matrix_builder=lambda c: build_cosserat_matrix_two_d(c, 1.0, 1.0, 1.0),
        n_steps=81
    )

    print(f"\nForward trace results for triplets near the chapter's Delta T_1 "
          f"at lambda=9.05:")
    print(f"{'init lambda':>12} {'final lambda':>13} {'m_final':>10} "
          f"{'phi%':>5}")
    for i in range(len(initial_vals)):
        if 7.95 < initial_vals[i] < 8.30:
            phi, inner = describe(final_vecs[:, i], N, 13)
            m = mass_from_lambda(N, final_vals[i])
            tag = ""
            if abs(final_vals[i] - 9.052) < 0.05:
                tag = "  <-- MATCHES Delta T_1 at 9.05"
            print(f"  {initial_vals[i]:>11.4f} {final_vals[i]:>12.4f} "
                  f"{m:>9.2f} {phi*100:>4.1f}{tag}")

    # ============================================================
    # Test 3: N(1520) - which gerade triplet is the parent?
    # ============================================================
    print("\n" + "-" * 78)
    print("TEST 3: N(1520) (J^P = 3/2-) - which gerade triplet is the parent?")
    print("-" * 78)

    coords = build_n1535_cluster()
    N = len(coords)
    print(f"Cluster: N = {N} (proton + bilayer)")
    print(f"Bare-shell gerade triplet candidates:")
    print(f"  T_2g (lambda=6.155, parity g): lower")
    print(f"  T_1g (lambda=6.159, parity g): lower")
    print(f"  T_1g (lambda=8.071, parity g): upper, same as Delta parent")
    print(f"  T_2g (lambda=8.188, parity g): upper")

    initial_vals, final_vals, initial_vecs, final_vecs = trace_all(
        coords, n_inner=13,
        matrix_builder=lambda c: build_cosserat_matrix(c, 1.0, 1.0, 1.0),
        n_steps=81
    )

    print(f"\nForward trace results for triplets near N(1520) at lambda~8.32:")
    print(f"{'init lambda':>12} {'final lambda':>13} {'m_final':>10} "
          f"{'phi%':>5} {'inner%':>7}")
    for i in range(len(initial_vals)):
        if 6.10 < initial_vals[i] < 6.25:
            phi, inner = describe(final_vecs[:, i], N, 13)
            m = mass_from_lambda(N, final_vals[i])
            tag = ""
            if abs(final_vals[i] - 8.32) < 0.10:
                tag = "  <-- near N(1520)"
            print(f"  {initial_vals[i]:>11.4f} {final_vals[i]:>12.4f} "
                  f"{m:>9.2f} {phi*100:>4.1f}  {inner*100:>5.1f}%{tag}")

    # ============================================================
    # Summary
    # ============================================================
    print("\n" + "=" * 78)
    print("Summary of parity-flip rule test")
    print("=" * 78)

    print(f"\nParity-flip rule (proposed):")
    print(f"  Mass-mode irrep has parity OPPOSITE to spinor parity.")
    print(f"  J=1/2+: A_2u stiff phi-dom (proton, lambda=8.303)")
    print(f"  J=1/2-: A_2g stiff phi-dom (N(1535), lambda=6.193 -> 8.046)")
    print(f"  J=3/2+: T_1u or T_2u stiff phi-dom (Delta family)")
    print(f"  J=3/2-: T_1g or T_2g stiff phi-dom (N(1520) family)")


if __name__ == '__main__':
    main()
