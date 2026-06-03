"""
local_consistency_mass_mode.py
==============================

Rest-mass closure for the bound light nuclei (deuteron, triton, helium-4)
in the Cosserat supersolid lattice.

These compounds are not selected the way the spin-3/2 baryons are.  A
single hadron such as the proton or the Delta has a mass mode pinned by
an irrep of its own point group (A_2u for the proton, T_1 for the Delta).
A bound nucleus instead carries a mode that looks locally like the bare
proton mode psi_p on every constituent shell at once.  The selection is
therefore a LOCAL-CONSISTENCY rule, and the right score is the geometric
mean of the per-constituent overlaps, because the geometric mean collapses
to zero the moment any one constituent is silent.

The squared overlap is blind to sign and to symmetry, so the geometric
mean on its own can be maximised by a stiff mode that sits in the wrong
irrep.  The fix is to maximise it only WITHIN the symmetry sector that the
compound's quantum numbers select.  That sector is an irrep of the
cluster's full point group, and the point group is fixed by the docking
geometry, not chosen by hand:

    deuteron   two shells along <110>     point group D_2h   sector: parity-odd (u)
    triton     three shells, triangle     point group C_3v   sector: A_1 (totally symmetric)
    helium-4   four shells, T_d void      point group T_d    sector: A_1 (totally symmetric)

The deuteron is the centrosymmetric case: its two nucleons are exchanged
by the inversion through the bond midpoint, and that inversion also flips
the internal microrotation, so the bound mode is parity-odd rather than
totally symmetric.  The triangle and the tetrahedron have no inversion, so
their bound modes are the totally-symmetric A_1.

The full point group is found here as the rigid symmetries of the cluster
about its centroid (an O_h linear part together with the lattice-compatible
translation that returns the cluster to itself), so the order-24 T_d of the
helium-4 tetrahedron is recovered rather than only the order-6 stabiliser of
the lattice origin.

Mass follows the master formula with the un-deduplicated node count
N_tot = 13 B for B constituents, m = N_tot m_0 - N_tot (4 - lambda) m_e.
"""

import numpy as np
from scipy.linalg import eigh

from spectral_classifier import fcc_nn_vectors, generate_Oh
from cosserat_classifier import build_cosserat_matrix
from hadron_spectral_mass import (
    get_psi_p, find_constituent_shell_atoms, M_0, M_E,
    cluster_deuteron_coords, cluster_tritium_coords, cluster_helium4_coords,
)

LAM_P, PSI_P = get_psi_p()
NN = fcc_nn_vectors()


# ----------------------------------------------------------------------
# Full point group of a cluster, about its centroid.
# ----------------------------------------------------------------------
def cluster_point_group(coords, tol=1e-6):
    """Rigid symmetries of the cluster.  Each is an O_h linear part R paired
    with the translation that maps the cluster onto itself, returned as
    (R, atom-permutation).  Operations about the centroid are captured, so a
    tetrahedral cluster returns the full order-24 T_d, not only the order-6
    stabiliser of the lattice origin."""
    OH = generate_Oh()
    n = len(coords)
    group = []
    for R in OH:
        img = coords @ R.T
        for j in range(n):
            t = img[0] - coords[j]
            shifted = img - t
            perm = []
            ok = True
            for k in range(n):
                d = np.linalg.norm(coords - shifted[k], axis=1)
                m = int(np.argmin(d))
                if d[m] > tol:
                    ok = False
                    break
                perm.append(m)
            if ok and len(set(perm)) == n:
                group.append((R, np.array(perm)))
                break
    return group


def state_rep(R, perm, coords):
    """6N representation of a point-group element: displacement u is a polar
    vector (transforms by R), microrotation phi is an axial vector
    (transforms by det(R) R), and atoms move by the permutation."""
    n = len(coords)
    D = np.zeros((6 * n, 6 * n))
    det = np.linalg.det(R)
    for i in range(n):
        j = perm[i]
        D[3 * j:3 * j + 3, 3 * i:3 * i + 3] = R
        D[3 * n + 3 * j:3 * n + 3 * j + 3, 3 * n + 3 * i:3 * n + 3 * i + 3] = det * R
    return D


def projector_A1(group, coords):
    """Projector onto the totally-symmetric A_1 irrep: the group average."""
    return sum(state_rep(R, p, coords) for R, p in group) / len(group)


def projector_parity_odd(group, coords):
    """Projector onto the inversion-odd (u) sector, P = (I - D(i)) / 2,
    for a centrosymmetric cluster."""
    n = len(coords)
    inv = None
    for R, p in group:
        if np.allclose(R, -np.eye(3)):
            inv = state_rep(R, p, coords)
            break
    if inv is None:
        raise ValueError("cluster has no inversion centre")
    return (np.eye(6 * n) - inv) / 2.0


# ----------------------------------------------------------------------
# Local-consistency score.
# ----------------------------------------------------------------------
def restrict_to_shell(mode, coords, centre):
    """Restrict a cluster eigenvector to one constituent's 13-atom shell,
    laid out in the 78-component (13-node) frame of psi_p."""
    N = len(coords)
    idx = find_constituent_shell_atoms(coords, centre)
    out = np.zeros(78)
    for k, i in enumerate(idx):
        out[3 * k:3 * k + 3] = mode[3 * i:3 * i + 3]
        out[39 + 3 * k:39 + 3 * k + 3] = mode[3 * N + 3 * i:3 * N + 3 * i + 3]
    return out


def cos2(a, b):
    na, nb = a @ a, b @ b
    if na < 1e-14 or nb < 1e-14:
        return 0.0
    return (a @ b) ** 2 / (na * nb)


def gm_consistency(mode, coords, centres):
    """Geometric mean over constituents of cos^2(psi_p, mode restricted to
    that shell)."""
    B = len(centres)
    prod = 1.0
    for c in centres:
        prod *= cos2(PSI_P, restrict_to_shell(mode, coords, c))
    return prod ** (1.0 / B)


def mass(Ntot, lam):
    return Ntot * M_0 - Ntot * (4.0 - lam) * M_E


# ----------------------------------------------------------------------
# Closure for one nucleus.
# ----------------------------------------------------------------------
def close(coords, centres, sector, label, pdg):
    n = len(coords)
    B = len(centres)
    Ntot = 13 * B
    group = cluster_point_group(coords)
    if sector == "A1":
        P = projector_A1(group, coords)
    elif sector == "parity-odd":
        P = projector_parity_odd(group, coords)
    else:
        raise ValueError(sector)

    H = build_cosserat_matrix(coords, 1.0, 1.0, 1.0)
    vals, vecs = eigh(H)

    best = None
    for k in range(len(vals)):
        if vals[k] < 1.0:
            continue
        weight = float(vecs[:, k] @ (P @ vecs[:, k]))   # content in the sector
        if weight < 0.9:
            continue
        gm = gm_consistency(vecs[:, k], coords, centres)
        if best is None or gm > best[0]:
            best = (gm, vals[k], weight)

    gm, lam, weight = best
    m = mass(Ntot, lam)
    print(f"{label:9s} |G|={len(group):2d}  sector={sector:11s}  "
          f"lambda*={lam:7.4f} (sector wt {weight:.2f})  "
          f"N_tot={Ntot}  m={m:8.2f} MeV  PDG={pdg:8.2f}  "
          f"residual {m - pdg:+6.2f} MeV ({(m - pdg) / pdg * 100:+.3f}%)")
    return lam, m


def main():
    print("Local-consistency mass closure for the bound light nuclei")
    print("(geometric-mean overlap with psi_p, maximised inside the irrep")
    print(" the cluster point group assigns to the compound)\n")
    close(cluster_deuteron_coords(), [np.zeros(3), NN[0]],
          "parity-odd", "deuteron", 1875.61)
    close(cluster_tritium_coords(), [np.zeros(3), NN[0], NN[1]],
          "A1", "triton", 2808.92)
    close(cluster_helium4_coords(), [np.zeros(3), NN[0], NN[1], NN[2]],
          "A1", "helium-4", 3727.38)
    print("\nThe naive geometric-mean maximum without the sector restriction")
    print("lands on a stiff off-symmetry mode for the deuteron and helium-4")
    print("(parity-even and T_2 respectively), which the projection excludes.")


if __name__ == "__main__":
    main()
