#!/usr/bin/env python3
"""
retrospective_adiabatic.py
==========================

Retrospective adiabatic-continuation check for the earlier closures:
  - N(1535): proton + bilayer (N=21), chapter claim lambda = 8.046
  - N(1440) Roper: proton + second-shell (N=19), chapter claim lambda = 10.193

For each, identify the proton's psi_p mode at lambda = 8.303 in the
tau=0 spectrum (extension decoupled), then trace it through to tau=1.
The final lambda is the principled prediction; it should match the
chapter's claim if the original picking rule happened to identify the
adiabatic descendant.
"""

import sys

import numpy as np
import scipy.linalg
from scipy.optimize import linear_sum_assignment

from hadron_spectral_mass import mass_from_lambda
from spectral_classifier import cluster_coord_shell
from cosserat_classifier import build_cosserat_matrix
from delta_first_principles import build_cosserat_matrix_two_d
from n1535_first_principles import build_n1535_cluster


def add_second_shell(cluster):
    """Append six octahedral second-shell atoms at +/-sqrt(2) along axes."""
    s = np.sqrt(2.0)
    extra = np.array([
        [+s, 0, 0], [-s, 0, 0],
        [0, +s, 0], [0, -s, 0],
        [0, 0, +s], [0, 0, -s],
    ])
    return np.vstack([cluster, extra])


def adiabatic_track(coords, n_inner, matrix_builder, n_steps=41,
                    parent_lambda=8.303, parent_label="psi_p"):
    """Trace eigenvalues from the n_inner-atom cluster to the full cluster.

    Parameters
    ----------
    coords : (N, 3) array of all atom positions; first n_inner are the
        parent cluster atoms, rest are the extension atoms.
    n_inner : number of parent-cluster atoms.
    matrix_builder : function (coords, ...) -> H_Cosserat.
    n_steps : number of tau values in [0, 1].
    parent_lambda : eigenvalue to track from tau=0.

    Returns
    -------
    final_lambda : eigenvalue of the parent's descendant at tau=1.
    trajectory : list of (tau, lambda) pairs.
    final_eigvec : the descendant eigenvector at tau=1.
    """
    N = len(coords)
    n_dim = 6 * N

    # H_full at tau=1
    H_full = matrix_builder(coords)
    # H_partial at tau=0: only the inner cluster has bonds; extension atoms
    # are decoupled (zero modes).  Embed the inner Cosserat matrix into
    # the 6N space with zeros for extension atoms.
    H_inner = matrix_builder(coords[:n_inner])
    H_partial = np.zeros((6 * N, 6 * N))
    # H_inner layout: [u_0..u_{n-1}, phi_0..phi_{n-1}], each 3D.
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

    # Find parent eigenvalue index in initial spectrum
    parent_idx = int(np.argmin(np.abs(initial_vals - parent_lambda)))
    if abs(initial_vals[parent_idx] - parent_lambda) > 0.05:
        print(f"  WARNING: closest initial eigenvalue to {parent_lambda} is "
              f"{initial_vals[parent_idx]}")

    final_lambda = vals_prev[parent_idx]
    final_eigvec = vecs_prev[:, parent_idx]

    # Build trajectory by re-doing the interpolation and tracking only the
    # parent's index; simpler: just report initial and final
    trajectory = [(0.0, initial_vals[parent_idx]), (1.0, final_lambda)]

    return final_lambda, trajectory, final_eigvec, initial_vals, vals_prev


def describe_eigvec(eigvec, n_atoms, n_inner):
    """Compute phi fraction and inner/outer localisation."""
    u_tot = np.sum(eigvec[:3*n_atoms]**2)
    phi_tot = np.sum(eigvec[3*n_atoms:]**2)
    phi_frac = phi_tot / (u_tot + phi_tot)
    amps = np.zeros(n_atoms)
    for i in range(n_atoms):
        amps[i] = (np.sum(eigvec[3*i:3*(i+1)]**2) +
                   np.sum(eigvec[3*n_atoms+3*i:3*n_atoms+3*(i+1)]**2))
    inner = amps[:n_inner].sum() / amps.sum()
    outer = amps[n_inner:].sum() / amps.sum()
    return phi_frac, inner, outer


def main():
    print("=" * 78)
    print("Retrospective adiabatic check for earlier closures")
    print("=" * 78)

    # ===========================================================
    # Check 1: N(1535) on proton + bilayer
    # Chapter claim: lambda = 8.046, m = 1513.95 MeV
    # ===========================================================
    print("\n" + "-" * 78)
    print("Check 1: N(1535) = proton + bilayer (N=21)")
    print("-" * 78)
    print("Chapter claim: lambda = 8.046, m = 1513.95 MeV")

    coords = build_n1535_cluster()
    N = len(coords)
    print(f"Cluster: {N} atoms (13 shell + 8 bilayer)")

    # The bilayer uses NN bonds only (single bond length)
    final_lambda, _, eigvec, initial_vals, final_vals = adiabatic_track(
        coords, n_inner=13,
        matrix_builder=lambda c: build_cosserat_matrix(
            c, K_u=1.0, K_phi=1.0, alpha=1.0),
        n_steps=41, parent_lambda=8.303
    )

    m_pred = mass_from_lambda(N, final_lambda)
    phi_frac, inner_frac, outer_frac = describe_eigvec(eigvec, N, 13)
    print(f"  Adiabatic continuation of psi_p (lambda=8.303 at tau=0):")
    print(f"    -> final lambda = {final_lambda:.4f}")
    print(f"    -> m = {m_pred:.2f} MeV")
    print(f"    -> phi content = {phi_frac*100:.1f}%, "
          f"shell = {inner_frac*100:.1f}%, bilayer = {outer_frac*100:.1f}%")
    print(f"  Chapter's lambda = 8.046")
    print(f"  Match: {'YES' if abs(final_lambda - 8.046) < 0.01 else 'NO'} "
          f"(difference {final_lambda - 8.046:+.4f})")

    # ===========================================================
    # Check 2: Roper N(1440) on proton + second-shell
    # Chapter claim: lambda = 10.193, m = 1390.60 MeV
    # ===========================================================
    print("\n" + "-" * 78)
    print("Check 2: N(1440) Roper = proton + second-shell (N=19)")
    print("-" * 78)
    print("Chapter claim: lambda = 10.193, m = 1390.60 MeV")

    shell, _ = cluster_coord_shell()
    coords = add_second_shell(shell)
    N = len(coords)
    print(f"Cluster: {N} atoms (13 shell + 6 second-shell)")

    # Second-shell uses two bond lengths (1 and sqrt(2))
    final_lambda, _, eigvec, initial_vals, final_vals = adiabatic_track(
        coords, n_inner=13,
        matrix_builder=lambda c: build_cosserat_matrix(
            c, K_u=1.0, K_phi=1.0, alpha=1.0),
        n_steps=41, parent_lambda=8.303
    )

    m_pred = mass_from_lambda(N, final_lambda)
    phi_frac, inner_frac, outer_frac = describe_eigvec(eigvec, N, 13)
    print(f"  Adiabatic continuation of psi_p (lambda=8.303 at tau=0):")
    print(f"    -> final lambda = {final_lambda:.4f}")
    print(f"    -> m = {m_pred:.2f} MeV")
    print(f"    -> phi content = {phi_frac*100:.1f}%, "
          f"shell = {inner_frac*100:.1f}%, second-shell = {outer_frac*100:.1f}%")
    print(f"  Chapter's lambda = 10.193")
    print(f"  Match: {'YES' if abs(final_lambda - 10.193) < 0.01 else 'NO'} "
          f"(difference {final_lambda - 10.193:+.4f})")

    # ===========================================================
    # Check 3: N(1520) on proton + bilayer
    # Chapter claim: lambda = 8.325, m = 1516.94 MeV
    # The N(1520) is bilayer-LOCALISED; its parent isn't psi_p
    # (which is shell-localised).  What's its parent?
    # ===========================================================
    print("\n" + "-" * 78)
    print("Check 3: N(1520) = proton + bilayer, J=3/2-")
    print("-" * 78)
    print("Chapter claim: lambda = 8.325 (bilayer-localised triplet)")
    print("Parent identification needed: bilayer-localised modes have no")
    print("clean parent in the tau=0 spectrum (the decoupled bilayer atoms")
    print("contribute only zero modes).  This needs separate analysis.")

    # Try: identify the bilayer-localised modes in the full spectrum and
    # report their adiabatic origin (where they come from at tau=0)
    coords = build_n1535_cluster()
    N = len(coords)
    H_full = build_cosserat_matrix(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    vals_full, vecs_full = scipy.linalg.eigh(H_full)
    print(f"\nFinal spectrum: searching for bilayer-localised modes "
          f"near lambda = 8.325")
    for k in range(len(vals_full)):
        if 8.2 < vals_full[k] < 8.5:
            phi_frac, inner_frac, outer_frac = describe_eigvec(
                vecs_full[:, k], N, 13)
            if outer_frac > 0.5:  # bilayer-dominated
                print(f"  k={k}: lambda={vals_full[k]:.4f}, "
                      f"phi={phi_frac*100:.1f}%, "
                      f"shell={inner_frac*100:.1f}%, "
                      f"bilayer={outer_frac*100:.1f}%")

    # ===========================================================
    # Check 4: Lambda(1600)
    # Chapter claim: lambda = 10.380, m = 1612.28 on N=22
    # Cluster: Lambda ground state + 6 second-shell
    # Parent: Lambda's psi_p-analog on the N=16 cluster
    # ===========================================================
    print("\n" + "-" * 78)
    print("Check 4: Lambda(1600) = Lambda + second-shell (N=22)")
    print("-" * 78)
    print("Chapter claim: lambda = 10.380, m = 1612.28 MeV")
    print("Need: Lambda ground state's mass mode at lambda = 3.295 (Class 1)")
    print("BUT Lambda(1600) was previously picked at high lambda via")
    print("max-projection-onto-psi_p, which is NOT Lambda's own mass mode.")
    print("If Lambda(1600) is the second-shell radial excitation of Lambda,")
    print("its adiabatic parent should be Lambda's lambda=3.295 mode.")

    from spectral_classifier import build_lambda_cluster
    try:
        lambda_ground = build_lambda_cluster()
        n_inner = len(lambda_ground)
        coords = add_second_shell(lambda_ground)
        N = len(coords)
        print(f"\nCluster: {N} atoms ({n_inner} Lambda ground + 6 second-shell)")

        # First: identify Lambda's mass mode at lambda=3.295 in the inner
        # spectrum
        H_inner = build_cosserat_matrix(lambda_ground, 1.0, 1.0, 1.0)
        vals_inner, vecs_inner = scipy.linalg.eigh(H_inner)
        # Look for lambda ~ 3.295
        candidates = np.where((vals_inner > 3.2) & (vals_inner < 3.4))[0]
        print(f"\n  Lambda ground state spectrum near lambda = 3.295:")
        for k in candidates:
            print(f"    k={k}: lambda={vals_inner[k]:.4f}")

        # Now adiabatic continuation of lambda = 3.295
        final_lambda, _, eigvec, _, _ = adiabatic_track(
            coords, n_inner=n_inner,
            matrix_builder=lambda c: build_cosserat_matrix(
                c, K_u=1.0, K_phi=1.0, alpha=1.0),
            n_steps=41, parent_lambda=3.295
        )
        m_pred = mass_from_lambda(N, final_lambda)
        print(f"\n  Adiabatic continuation of Lambda's lambda = 3.295:")
        print(f"    -> final lambda = {final_lambda:.4f}")
        print(f"    -> m = {m_pred:.2f} MeV")
        print(f"  vs PDG Lambda(1600) = 1600 +/- 70 MeV (range 1560-1700)")

        # Also try psi_p (lambda=8.303) as parent
        final_lambda_p, _, eigvec_p, _, _ = adiabatic_track(
            coords, n_inner=n_inner,
            matrix_builder=lambda c: build_cosserat_matrix(
                c, K_u=1.0, K_phi=1.0, alpha=1.0),
            n_steps=41, parent_lambda=8.303
        )
        m_pred_p = mass_from_lambda(N, final_lambda_p)
        print(f"\n  Alternative: adiabatic continuation of psi_p (lambda=8.303):")
        print(f"    -> final lambda = {final_lambda_p:.4f}")
        print(f"    -> m = {m_pred_p:.2f} MeV")
        print(f"  Chapter's earlier pick: lambda = 10.380, m = 1612.28 MeV")
    except Exception as e:
        print(f"  Lambda cluster build failed: {e}")


if __name__ == '__main__':
    main()
