#!/usr/bin/env python3
"""
roper_continuation.py
=====================
Do the two selection rules for extension clusters agree?

The chapter names two ways of picking the mass mode on a cluster carrying an
extension beyond the inner coordination shell.  One is adiabatic continuation:
turn the extension on continuously and follow the parent cluster's mass mode
across.  The other is the primary-cluster heuristic, the lowest stiff
microrotation-dominant shell-concentrated root, which is what fixes the proton,
the Sigma and the Delta on their own hosts.

On the N(1535) bilayer the two happen to agree.  On the Roper they do not, and
this script is the demonstration.

The Roper's host is the 19-node Born cluster: the 13-node coordination shell
plus the six second-shell sites on the <100> axes.  Continuing the proton's
A_2u mass mode at lambda = 8.303 across the turn-on gives lambda = 10.1926 with
a worst step overlap of 1.0000, so the tracking never comes near a crossing.
That is the value the mass closure uses.  The heuristic applied to the 19-node
cluster instead returns lambda = 8.1673, some twenty MeV lower in mass.

The conclusion is that continuation is the rule and the heuristic is not an
independent route to the same answer.  On a primary cluster the two coincide,
because there is no extension to continue across and the parent is the cluster
itself.  On an extension cluster the heuristic can return a mode that never
descends from the parent at all.

Run:  python3 roper_continuation.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cosserat_classifier import build_cosserat_matrix            # noqa: E402
from spectral_classifier import cluster_born, cluster_coord_shell  # noqa: E402

# --- CODATA 2022 -----------------------------------------------------------
M_E = 0.51099895069           # electron mass [MeV]
ALPHA_FS = 7.2973525643e-3    # fine-structure constant [-]
M_0 = M_E / ALPHA_FS          # node assembly quantum [MeV]

ELL = 1.0
N_NODES = 19
STEPS = 600
PARENT_LAMBDA = 8.3028        # the proton's A_2u stiff mass mode on the shell


def mass(lam, n=N_NODES):
    """Master spectral formula."""
    return n * (M_0 + (lam - 4.0) * M_E)


def as_coords(builder_output):
    if isinstance(builder_output, tuple):
        builder_output = builder_output[0]
    return np.atleast_2d(builder_output)


def bond_list(coords, tol=1e-6):
    """Ordered pairs at nearest-neighbour separation."""
    return [(i, j)
            for i in range(len(coords)) for j in range(len(coords))
            if i != j and abs(np.linalg.norm(coords[j] - coords[i]) - ELL) < tol]


def weighted_cosserat(coords, bonds, bond_weight, site_weight, alpha=1.0):
    """Cosserat dynamical matrix with per-bond and per-site coupling weights.

    Reduces to build_cosserat_matrix when every weight is one.
    """
    n = len(coords)
    uu = np.zeros((3 * n, 3 * n))
    pp = np.zeros((3 * n, 3 * n))
    curl = np.zeros((3 * n, 3 * n))
    for (i, j) in bonds:
        rhat = (coords[j] - coords[i]) / np.linalg.norm(coords[j] - coords[i])
        w = bond_weight[(i, j)]
        block = np.outer(rhat, rhat) * w
        uu[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= block
        uu[3 * i:3 * i + 3, 3 * i:3 * i + 3] += block
        pp[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= np.eye(3) * w
        pp[3 * i:3 * i + 3, 3 * i:3 * i + 3] += np.eye(3) * w
        cross = np.array([[0.0, -rhat[2], rhat[1]],
                          [rhat[2], 0.0, -rhat[0]],
                          [-rhat[1], rhat[0], 0.0]])
        curl[3 * i:3 * i + 3, 3 * j:3 * j + 3] += cross * w / (2.0 * ELL)
    a3 = np.diag(np.repeat([alpha * site_weight[i] for i in range(n)], 3))
    mat = np.zeros((6 * n, 6 * n))
    mat[:3 * n, :3 * n] = uu + curl.T @ a3 @ curl
    mat[3 * n:, 3 * n:] = pp + a3
    mat[:3 * n, 3 * n:] = -(a3 @ curl).T
    mat[3 * n:, :3 * n] = -(a3 @ curl)
    return mat


def heuristic_candidates(matrix, n, shell_idx,
                         min_phi=0.95, min_shell=0.95, stiff_above=4.0):
    """Stiff, microrotation-dominant, shell-concentrated roots, ascending."""
    vals, vecs = np.linalg.eigh(matrix)
    out = []
    for a in range(len(vals)):
        if vals[a] <= stiff_above:
            continue
        psi = vecs[:, a]
        phi_w = float(np.sum(psi[3 * n:] ** 2))
        shell_w = sum(psi[3 * k + m] ** 2 + psi[3 * n + 3 * k + m] ** 2
                      for k in shell_idx for m in range(3))
        if phi_w >= min_phi and shell_w >= min_shell:
            out.append((vals[a], phi_w, shell_w))
    return sorted(out)


def main():
    parent = as_coords(cluster_coord_shell())
    child = as_coords(cluster_born())
    n = len(child)
    shell_idx = [int(np.argmin(np.linalg.norm(child - p, axis=1))) for p in parent]
    extension = [j for j in range(n) if j not in shell_idx]
    bonds = bond_list(child)

    print(f"Roper host: {len(parent)}-node shell inside the {n}-node Born cluster")
    print(f"extension nodes (second shell, <100> axes): {extension}")

    unit = weighted_cosserat(child, bonds,
                             {b: 1.0 for b in bonds},
                             {i: 1.0 for i in range(n)})
    check = np.abs(unit - build_cosserat_matrix(child, alpha=1.0)).max()
    print(f"weighted builder reproduces the standard one: max abs diff {check:.2e}")

    def bw(tau):
        return {(i, j): (tau if (i in extension or j in extension) else 1.0)
                for (i, j) in bonds}

    def sw(tau):
        return {i: (tau if i in extension else 1.0) for i in range(n)}

    taus = np.linspace(1e-4, 1.0, STEPS)
    vals, vecs = np.linalg.eigh(weighted_cosserat(child, bonds,
                                                  bw(taus[0]), sw(taus[0])))
    seed = int(np.argmin(np.abs(vals - PARENT_LAMBDA)))
    print(f"\nseed at the decoupled end: lambda = {vals[seed]:.4f}")
    psi, worst = vecs[:, seed], 1.0
    for tau in taus[1:]:
        vals, vecs = np.linalg.eigh(weighted_cosserat(child, bonds,
                                                      bw(tau), sw(tau)))
        overlap = (vecs.T @ psi) ** 2
        k = int(np.argmax(overlap))
        worst = min(worst, overlap[k])
        psi, lam_cont = vecs[:, k], vals[k]

    print(f"CONTINUATION   lambda = {lam_cont:8.4f}   m = {mass(lam_cont):8.2f} MeV"
          f"   (worst step overlap {worst:.4f})")

    cands = heuristic_candidates(unit, n, shell_idx)
    lam_heur = cands[0][0]
    print(f"HEURISTIC      lambda = {lam_heur:8.4f}   m = {mass(lam_heur):8.2f} MeV"
          f"   ({len(cands)} qualifying roots)")
    print(f"\ndisagreement: {abs(mass(lam_cont) - mass(lam_heur)):.2f} MeV")
    print("The mass closure uses the continuation value.")

    print("\nlowest qualifying roots under the heuristic:")
    for lam, phi_w, shell_w in cands[:6]:
        print(f"   lambda = {lam:8.4f}  phi = {phi_w:.3f}  shell = {shell_w:.3f}"
              f"  m = {mass(lam):8.2f} MeV")


if __name__ == "__main__":
    raise SystemExit(main())
