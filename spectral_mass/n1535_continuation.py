#!/usr/bin/env python3
"""
n1535_continuation.py
=====================
Which bare-shell mode does each mode of the bilayer-extended cluster descend
from?

The 21-node cluster carries eight phi-dominant, shell-concentrated, stiff
modes packed between lambda = 8.046 and 8.130, and the mass closure reads the
lowest of them.  Their eigenvalues sit close to the bare 13-node shell's own
group near 8.07, which invites the reading that they are that group slightly
split.  This script shows that reading is wrong.  The lineage is settled by
turning the extension on continuously and tracking each eigenvector by maximum
overlap between adjacent steps, and it puts the 8.07 group somewhere else
entirely.

Two independent paths are run, because a continuation result that depends on
the path is not a result:

  path A   linear interpolation of the dynamical matrix between the
           block-diagonal limit (bare shell plus isolated bilayer) and the
           full 21-node matrix;
  path B   scaling every bond that touches a bilayer node, together with the
           Cosserat coupling on those nodes, by the same parameter.

Both agree.  The shell's 8.07 multiplet lands between 10.9 and 11.5, and the
proton's own A_2u mass mode at 8.303 lands at 11.55.  The 8.03 to 8.13 band on
the extended cluster descends instead from the shell's stiff microrotation
group near lambda = 6.2, which is where the A_2g singlet at 6.193 sits, the
mode the parity-flip rule assigns to a spin-half negative-parity state.

One thing the script does not settle.  Within the 6.19 to 6.23 cluster the
step-to-step overlaps fall as low as 0.52, so those modes mix on the way
across and no individual member can be followed to an individual descendant.
Only the group-level lineage is firm.

Run:  python3 n1535_continuation.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cosserat_classifier import build_cosserat_matrix          # noqa: E402
from spectral_classifier import cluster_coord_shell            # noqa: E402
from n1535_first_principles import build_n1535_cluster         # noqa: E402

# --- CODATA 2022 -----------------------------------------------------------
M_E = 0.51099895069           # electron mass [MeV]
ALPHA_FS = 7.2973525643e-3    # fine-structure constant [-]
M_0 = M_E / ALPHA_FS          # node assembly quantum [MeV]

ELL = 1.0
N_NODES = 21
STEPS = 500


def mass(lam, n=N_NODES):
    """Master spectral formula."""
    return n * (M_0 + (lam - 4.0) * M_E)


def clusters():
    """Bare shell, extended cluster, and the index maps between them."""
    out = cluster_coord_shell()
    shell = np.atleast_2d(out[0] if isinstance(out, tuple) else out)
    child = np.atleast_2d(build_n1535_cluster())
    shell_idx = [int(np.argmin(np.linalg.norm(child - p, axis=1))) for p in shell]
    bilayer_idx = [j for j in range(len(child)) if j not in shell_idx]
    return shell, child, shell_idx, bilayer_idx


def bond_list(coords, tol=1e-6):
    """Ordered pairs at nearest-neighbour separation."""
    out = []
    for i in range(len(coords)):
        for j in range(len(coords)):
            if i != j and abs(np.linalg.norm(coords[j] - coords[i]) - ELL) < tol:
                out.append((i, j))
    return out


def weighted_cosserat(coords, bond_weight, site_weight, alpha=1.0):
    """Cosserat dynamical matrix with per-bond and per-site coupling weights.

    Reduces to build_cosserat_matrix when every weight is one.
    """
    n = len(coords)
    uu = np.zeros((3 * n, 3 * n))
    pp = np.zeros((3 * n, 3 * n))
    curl = np.zeros((3 * n, 3 * n))
    for (i, j) in bond_list(coords):
        r = coords[j] - coords[i]
        rhat = r / np.linalg.norm(r)
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


def embed_block(sub, nodes, n_total):
    """Place a 6m x 6m Cosserat block for `nodes` into the 6n x 6n space."""
    m = len(nodes)
    out = np.zeros((6 * n_total, 6 * n_total))
    for a, ja in enumerate(nodes):
        for b, jb in enumerate(nodes):
            out[3 * ja:3 * ja + 3, 3 * jb:3 * jb + 3] = \
                sub[3 * a:3 * a + 3, 3 * b:3 * b + 3]
            out[3 * ja:3 * ja + 3, 3 * n_total + 3 * jb:3 * n_total + 3 * jb + 3] = \
                sub[3 * a:3 * a + 3, 3 * m + 3 * b:3 * m + 3 * b + 3]
            out[3 * n_total + 3 * ja:3 * n_total + 3 * ja + 3, 3 * jb:3 * jb + 3] = \
                sub[3 * m + 3 * a:3 * m + 3 * a + 3, 3 * b:3 * b + 3]
            out[3 * n_total + 3 * ja:3 * n_total + 3 * ja + 3,
                3 * n_total + 3 * jb:3 * n_total + 3 * jb + 3] = \
                sub[3 * m + 3 * a:3 * m + 3 * a + 3, 3 * m + 3 * b:3 * m + 3 * b + 3]
    return out


def track(matrices, seed_vector):
    """Follow one eigenvector along a path by maximum overlap at each step."""
    psi = seed_vector
    worst = 1.0
    lam = None
    for mat in matrices:
        vals, vecs = np.linalg.eigh(mat)
        overlap = (vecs.T @ psi) ** 2
        k = int(np.argmax(overlap))
        worst = min(worst, overlap[k])
        psi, lam = vecs[:, k], vals[k]
    return lam, worst


def path_a(shell, child, shell_idx, bilayer_idx):
    """Linear interpolation between the decoupled limit and the full matrix."""
    n = len(child)
    full = build_cosserat_matrix(child, alpha=1.0)
    decoupled = (embed_block(build_cosserat_matrix(shell, alpha=1.0), shell_idx, n)
                 + embed_block(build_cosserat_matrix(child[bilayer_idx], alpha=1.0),
                               bilayer_idx, n))
    taus = np.linspace(0.0, 1.0, STEPS + 1)
    return decoupled, [(1 - t) * decoupled + t * full for t in taus[1:]]


def path_b(child, bilayer_idx):
    """Scale every bond and coupling that touches a bilayer node."""
    bonds = bond_list(child)
    n = len(child)

    def build(tau):
        bw = {(i, j): (tau if (i in bilayer_idx or j in bilayer_idx) else 1.0)
              for (i, j) in bonds}
        sw = {i: (tau if i in bilayer_idx else 1.0) for i in range(n)}
        return weighted_cosserat(child, bw, sw, alpha=1.0)

    taus = np.linspace(1e-4, 1.0, STEPS)
    return build(taus[0]), [build(t) for t in taus[1:]]


def stiff_phi_shell_seeds(matrix, shell_idx, n, min_shell=0.98, min_phi=0.90):
    """Stiff, microrotation-dominant, shell-concentrated modes at the path start."""
    vals, vecs = np.linalg.eigh(matrix)
    seeds = []
    for a in range(len(vals)):
        if vals[a] <= 4.0:
            continue
        psi = vecs[:, a]
        shell_w = sum(psi[3 * j + m] ** 2 + psi[3 * n + 3 * j + m] ** 2
                      for j in shell_idx for m in range(3))
        phi_w = float(np.sum(psi[3 * n:] ** 2))
        if shell_w > min_shell and phi_w > min_phi:
            seeds.append((vals[a], vecs[:, a]))
    return seeds


def main():
    shell, child, shell_idx, bilayer_idx = clusters()
    n = len(child)

    check = np.abs(
        weighted_cosserat(child,
                          {b: 1.0 for b in bond_list(child)},
                          {i: 1.0 for i in range(n)}, alpha=1.0)
        - build_cosserat_matrix(child, alpha=1.0)).max()
    print(f"weighted builder reproduces the standard one at unit weights: "
          f"max abs difference {check:.2e}")

    for name, (start, steps) in (("A (matrix interpolation)",
                                  path_a(shell, child, shell_idx, bilayer_idx)),
                                 ("B (bond scaling)",
                                  path_b(child, bilayer_idx))):
        print(f"\nPATH {name}")
        seeds = stiff_phi_shell_seeds(start, shell_idx, n)
        print(f"  {len(seeds)} stiff phi-dominant shell seeds at the path start")
        print(f"  {'seed':>8s} -> {'descendant':>11s} {'mass (MeV)':>11s} "
              f"{'worst overlap':>14s}")
        for lam0, vec in seeds:
            lam, worst = track(steps, vec)
            print(f"  {lam0:8.4f} -> {lam:11.4f} {mass(lam):11.2f} {worst:14.4f}")


if __name__ == "__main__":
    raise SystemExit(main())
