"""
N(1535) negative-parity baryon: first-principles spectral derivation.

Cluster: 13-node cuboctahedral coordination shell + 8-node bilayer attached
along one <111> stacking direction.  The bilayer is a hexagonal cap (7 nodes
on the adjacent {111} layer above the shell's top triangle) plus 1 apex node
on the next {111} layer.  The bilayer attachment breaks the shell's O_h
symmetry down to C_s.

Selection rule: lowest stiff (lambda > 4) phi-dominant (> 95% microrotation)
mode on the 21-node cluster.  This is the parity-flipped analog of the
proton's stiff A_{2u} mass mode (lambda = 8.303 on the bare 13-node shell),
slightly shifted by the bilayer perturbation.  The bilayer's role is to flip
the cluster's parity (defining the J^P = 1/2^- assignment); the eigenvalue
shift gives the small additional mass increment beyond N m_0.

Output (with N = 21, m_0 = 70.0253 MeV, m_e = 0.510999 MeV):

  lambda = 8.0458
  m_N(1535) = 21 * m_0 - 21 * (4 - 8.0458) * m_e = 1513.95 MeV

  PDG N(1535) pole position: 1510 +/- 20 MeV
  Residual: +0.26% (within electromagnetic correction band)
"""

import numpy as np
from itertools import product
from scipy.linalg import eigh
import sys
from cosserat_classifier import build_cosserat_matrix


# ---------------------------------------------------------------------------
# 1. Build the 21-node cluster from raw FCC geometry
# ---------------------------------------------------------------------------

def build_n1535_cluster():
    """Construct shell + bilayer cluster from FCC positions.  Returns (21,3) array."""
    a = np.sqrt(2.0)                                # cubic edge for NN = 1
    sublat = [np.array([0, 0, 0]),
              np.array([a/2, a/2, 0]),
              np.array([a/2, 0, a/2]),
              np.array([0, a/2, a/2])]
    fcc = set()
    for ni, nj, nk in product(range(-3, 4), repeat=3):
        for s in sublat:
            p = np.array([ni*a, nj*a, nk*a]) + s
            fcc.add(tuple(np.round(p, 8)))
    fcc = np.array(sorted(fcc, key=lambda p: (p[0]**2 + p[1]**2 + p[2]**2,
                                              p[0], p[1], p[2])))

    n111 = np.array([1, 1, 1]) / np.sqrt(3)
    plane_dz = a / np.sqrt(3)                       # {111} plane spacing
    projs = fcc @ n111

    # Shell: centre + 12 nearest neighbours of origin
    d0 = np.linalg.norm(fcc, axis=1)
    shell_idx = np.concatenate([[np.argmin(d0)],
                                np.where(np.abs(d0 - 1.0) < 1e-6)[0]])
    shell = fcc[shell_idx]

    # Hex cap centre: FCC site on plane 2 (above shell's top triangle) that
    # is closest to the (111) axis through the origin
    p2_mask = np.abs(projs - 2*plane_dz) < 1e-6
    p2_nodes = fcc[p2_mask]
    perp = np.linalg.norm(p2_nodes - np.outer(p2_nodes @ n111, n111), axis=1)
    cap_centre = p2_nodes[np.argmin(perp)]

    # Hex cap ring: 6 FCC NN of cap centre that also live on plane 2
    d_cap = np.linalg.norm(p2_nodes - cap_centre, axis=1)
    ring = p2_nodes[np.argsort(d_cap)[1:7]]         # skip the centre itself

    # Apex: the FCC NN of cap centre on the next {111} plane up (plane 3)
    # nearest to the (111) axis through the origin
    p3_mask = np.abs(projs - 3*plane_dz) < 1e-6
    p3_nodes = fcc[p3_mask]
    perp3 = np.linalg.norm(p3_nodes - np.outer(p3_nodes @ n111, n111), axis=1)
    apex = p3_nodes[np.argmin(perp3)]

    return np.vstack([shell, [cap_centre], ring, [apex]])


# ---------------------------------------------------------------------------
# 2. Solve the Cosserat eigenvalue problem and report the mass closure
# ---------------------------------------------------------------------------

def main():
    cluster = build_n1535_cluster()
    n = len(cluster)

    # Sanity checks
    nn_bonds = sum(1 for i in range(n) for j in range(i+1, n)
                   if abs(np.linalg.norm(cluster[i] - cluster[j]) - 1.0) < 1e-6)
    print(f"Cluster:        {n} nodes")
    print(f"NN bonds:       {nn_bonds}")
    print(f"Node breakdown: 13 shell + 1 cap centre + 6 hex ring + 1 apex = 21")

    # Cosserat matrix at alpha = 1 with NN-only convention
    M = build_cosserat_matrix(cluster, K_u=1.0, K_phi=1.0, alpha=1.0)
    eigvals, eigvecs = eigh(M)

    # Selection: lowest stiff phi-dominant mode
    shell_indices  = list(range(13))
    bilayer_indices = list(range(13, 21))

    print(f"\nCandidate mass-mode search (lambda > 4, phi > 95%):")
    print(f"  {'lambda':>8s}  {'phi%':>5s}  {'shell%':>6s}  {'bilayer%':>8s}  {'m@N=21 (MeV)':>13s}")

    chosen = None
    for k, lam in enumerate(eigvals):
        if lam <= 4.0:
            continue
        mode = eigvecs[:, k] / np.linalg.norm(eigvecs[:, k])
        phi_frac = float(np.sum(mode[3*n:]**2))
        shell_frac = sum(np.sum(mode[3*i:3*i+3]**2) +
                         np.sum(mode[3*n+3*i:3*n+3*i+3]**2)
                         for i in shell_indices)
        if phi_frac < 0.95:
            continue
        m = 21 * 70.0253 - 21 * (4.0 - lam) * 0.510999
        print(f"  {lam:8.4f}  {phi_frac*100:4.1f}%  {shell_frac*100:5.1f}%  "
              f"{(1.0-shell_frac)*100:7.1f}%  {m:13.2f}")
        if chosen is None:
            chosen = (lam, phi_frac, shell_frac, m)

    lam, phi_frac, shell_frac, m = chosen
    pdg_pole = 1510.0                               # PDG N(1535) pole, MeV
    pull_pole = (m - pdg_pole) / pdg_pole * 100
    print(f"\nSelection: lowest stiff phi-dominant mode")
    print(f"  lambda     = {lam:.4f}")
    print(f"  phi-content = {phi_frac*100:.1f}%")
    print(f"  shell-content = {shell_frac*100:.1f}%")
    print(f"  m_N(1535)  = {m:.2f} MeV")
    print(f"  PDG pole   = {pdg_pole:.0f} +/- 20 MeV")
    print(f"  residual   = {pull_pole:+.3f}%  ({m-pdg_pole:+.2f} MeV)")


if __name__ == "__main__":
    main()
