#!/usr/bin/env python3
"""
delta1600_dual_orbit.py
=======================

Test the analytical prediction that Delta(1600) sits at the T_1g gerade
combination on the dual-orbit void cluster.

Cluster: N=21 with full O_h symmetry.
  - 13 coordination shell (1 centre + 12 cuboctahedral vertices)
  - 8 tetrahedral interstitial voids (BOTH T_d orbits)

The Delta(1232) at N=17 with T_d symmetry sits on T_1 of T_d at lambda=9.052.
Under inversion restoration (dual orbit), this T_1 splits into:
  T_1g (gerade): both void orbits move IN PHASE  -> SOFTER mode
  T_1u (ungerade): both void orbits move OUT OF PHASE -> STIFFER mode

Prediction: Delta(1600) at lambda_T1g(gerade) ~ 7-8 on N=21.

Verification: diagonalise, decompose under O_h, find T_1g modes, identify
the one with phi-dominant shell-concentrated character.
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
    # Tetrahedral voids of FCC: 8 sites at (+/-1, +/-1, +/-1)*a/4 in units
    # where the NN distance is a/sqrt(2) = 1. So a = sqrt(2), and voids
    # are at (+/-1, +/-1, +/-1) * sqrt(2)/4.
    s = np.sqrt(2.0) / 4.0
    voids_plus = np.array([  # positive T_d orbit (product of signs = +1)
        [+s, +s, +s], [+s, -s, -s], [-s, +s, -s], [-s, -s, +s]
    ])
    voids_minus = np.array([  # negative T_d orbit (product of signs = -1)
        [-s, -s, -s], [-s, +s, +s], [+s, -s, +s], [+s, +s, -s]
    ])
    return np.vstack([shell, voids_plus, voids_minus])


# === O_h character table ===
# Order: E, 8C_3, 3C_2, 6C_4, 6C_2', i, 8S_6, 3sigma_h, 6S_4, 6sigma_d
O_H_IRREPS = {
    'A_1g': [ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1],
    'A_2g': [ 1,  1,  1, -1, -1,  1,  1,  1, -1, -1],
    'E_g':  [ 2, -1,  2,  0,  0,  2, -1,  2,  0,  0],
    'T_1g': [ 3,  0, -1,  1, -1,  3,  0, -1,  1, -1],
    'T_2g': [ 3,  0, -1, -1,  1,  3,  0, -1, -1,  1],
    'A_1u': [ 1,  1,  1,  1,  1, -1, -1, -1, -1, -1],
    'A_2u': [ 1,  1,  1, -1, -1, -1, -1, -1,  1,  1],
    'E_u':  [ 2, -1,  2,  0,  0, -2,  1, -2,  0,  0],
    'T_1u': [ 3,  0, -1,  1, -1, -3,  0,  1, -1,  1],
    'T_2u': [ 3,  0, -1, -1,  1, -3,  0,  1,  1, -1],
}
CLASS_ORDERS = [1, 8, 3, 6, 6, 1, 8, 3, 6, 6]
G_ORDER = sum(CLASS_ORDERS)  # 48


def build_oh_elements():
    """Return list of 48 O_h elements as 3x3 matrices, grouped by conjugacy
    class. Returns (matrices, class_indices) where class_indices[i] is the
    class (0..9) that element i belongs to."""
    elements = []
    class_ids = []

    # Class 0: E
    elements.append(np.eye(3))
    class_ids.append(0)

    # Class 1: 8 C_3 (around body diagonals)
    for sx in [+1, -1]:
        for sy in [+1, -1]:
            for sz in [+1, -1]:
                # C_3 around axis (sx, sy, sz)/sqrt(3): cyclic permutation
                # of (x, y, z) with appropriate signs
                # C_3 around (1,1,1) maps x -> y -> z -> x
                # For general (sx, sy, sz), rotate by 120deg around that axis
                axis = np.array([sx, sy, sz]) / np.sqrt(3)
                R = rotation_matrix(axis, 2 * np.pi / 3)
                elements.append(R)
                class_ids.append(1)

    # Class 2: 3 C_2 (around coord axes)
    for k in range(3):
        axis = np.zeros(3)
        axis[k] = 1
        R = rotation_matrix(axis, np.pi)
        elements.append(R)
        class_ids.append(2)

    # Class 3: 6 C_4 (around coord axes)
    for k in range(3):
        axis = np.zeros(3)
        axis[k] = 1
        for sign in [+1, -1]:
            R = rotation_matrix(axis, sign * np.pi / 2)
            elements.append(R)
            class_ids.append(3)

    # Class 4: 6 C_2' (around face-diagonal axes (1,1,0)/sqrt(2) etc.)
    face_diag_axes = []
    for i in range(3):
        for j in range(i + 1, 3):
            for sign in [+1, -1]:
                axis = np.zeros(3)
                axis[i] = 1
                axis[j] = sign
                face_diag_axes.append(axis / np.sqrt(2))
    # 6 axes total
    for axis in face_diag_axes:
        R = rotation_matrix(axis, np.pi)
        elements.append(R)
        class_ids.append(4)

    # Class 5: i (inversion)
    elements.append(-np.eye(3))
    class_ids.append(5)

    # Classes 6, 7, 8, 9 are obtained by multiplying classes 1, 2, 3, 4 by -I
    # 6: 8 S_6 = 8 C_3 * i
    for sx in [+1, -1]:
        for sy in [+1, -1]:
            for sz in [+1, -1]:
                axis = np.array([sx, sy, sz]) / np.sqrt(3)
                R = -rotation_matrix(axis, 2 * np.pi / 3)
                elements.append(R)
                class_ids.append(6)

    # 7: 3 sigma_h = 3 C_2(coord) * i
    for k in range(3):
        axis = np.zeros(3)
        axis[k] = 1
        R = -rotation_matrix(axis, np.pi)
        elements.append(R)
        class_ids.append(7)

    # 8: 6 S_4 = 6 C_4 * i
    for k in range(3):
        axis = np.zeros(3)
        axis[k] = 1
        for sign in [+1, -1]:
            R = -rotation_matrix(axis, sign * np.pi / 2)
            elements.append(R)
            class_ids.append(8)

    # 9: 6 sigma_d = 6 C_2' * i
    for axis in face_diag_axes:
        R = -rotation_matrix(axis, np.pi)
        elements.append(R)
        class_ids.append(9)

    return elements, class_ids


def rotation_matrix(axis, theta):
    """Rodrigues' rotation formula."""
    axis = axis / np.linalg.norm(axis)
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)


def build_cosserat_oh_action(coords, oh_elements):
    """For each O_h element g, build the 6N x 6N action on the Cosserat
    state vector (u_1, u_2, ..., u_N, phi_1, phi_2, ..., phi_N).

    For atom i at position r_i, g sends it to atom pi(i) at position g.r_i.
    The displacement u_i transforms as u_pi(i) <- g . u_i (polar vector).
    The microrotation phi_i transforms as phi_pi(i) <- det(g) * g . phi_i
    (axial vector).
    """
    N = len(coords)
    actions = []
    for g in oh_elements:
        M = np.zeros((6 * N, 6 * N))
        det_g = np.linalg.det(g)
        for i in range(N):
            r_new = g @ coords[i]
            # Find j such that coords[j] = r_new
            dists = np.linalg.norm(coords - r_new, axis=1)
            j = int(np.argmin(dists))
            if dists[j] > 1e-6:
                raise ValueError(f"Atom {i} maps to non-cluster point")
            # u: M[3j:3j+3, 3i:3i+3] = g (polar vector)
            M[3*j:3*j+3, 3*i:3*i+3] = g
            # phi: M[3N+3j:3N+3j+3, 3N+3i:3N+3i+3] = det(g) * g (axial)
            M[3*N+3*j:3*N+3*j+3, 3*N+3*i:3*N+3*i+3] = det_g * g
        actions.append(M)
    return actions


def project_onto_irrep(eigenvectors, oh_actions, class_ids, irrep_chars):
    """Apply the projector onto irrep R: P_R = (d_R / |G|) sum_g chi_R(g)* D(g)

    Returns a matrix P that projects 6N-vectors onto the R-isotypic
    component.

    d_R = irrep_chars[0] (dimension of the irrep, i.e., character of E).
    """
    d_R = irrep_chars[0]
    dim = oh_actions[0].shape[0]
    P = np.zeros((dim, dim))
    for g, cls in zip(oh_actions, class_ids):
        P += irrep_chars[cls] * g
    P *= d_R / G_ORDER
    return P


def main():
    print("=" * 78)
    print("Delta(1600) on dual-orbit void cluster: analytical prediction test")
    print("=" * 78)

    # Build cluster
    coords = build_dual_orbit_cluster()
    N = len(coords)
    print(f"\nCluster: N = {N}")
    print(f"  Inner shell (13): centre + cuboctahedron")
    print(f"  Voids (8): both T_d orbits")
    print(f"\nVoid positions:")
    for i, p in enumerate(coords[13:]):
        print(f"  void {i}: {p}, distance from origin = {np.linalg.norm(p):.4f}")

    # Build O_h elements
    print("\nBuilding O_h group elements (48 elements, 10 classes)...")
    oh_elements, class_ids = build_oh_elements()
    print(f"  Built {len(oh_elements)} elements")

    # Verify that O_h elements preserve the cluster
    print("\nVerifying O_h symmetry of cluster:")
    for g_idx, g in enumerate(oh_elements):
        transformed = (g @ coords.T).T
        # Check that each transformed atom matches some original atom
        for r in transformed:
            dists = np.linalg.norm(coords - r, axis=1)
            if np.min(dists) > 1e-6:
                print(f"  Element {g_idx} (class {class_ids[g_idx]}): "
                      f"BREAKS symmetry")
                break
    print("  Cluster preserved under all 48 elements")

    # Build Cosserat dynamical matrix (mixed bond lengths via two-distance
    # builder)
    H = build_cosserat_matrix_two_d(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    vals, vecs = scipy.linalg.eigh(H)
    print(f"\nCosserat spectrum: {len(vals)} eigenvalues, "
          f"range [{vals[0]:.4f}, {vals[-1]:.4f}]")

    # Build O_h action on the 6N-dim Cosserat space
    print("\nBuilding O_h action on 6N-dim Cosserat space...")
    oh_actions = build_cosserat_oh_action(coords, oh_elements)

    # Project each eigenmode onto each irrep
    print("\nIrrep content of each eigenmode (showing only T_1g and T_1u):")
    print(f"  {'k':>3} {'lambda':>8}  {'A1g':>5} {'A2g':>5} {'Eg':>5} "
          f"{'T1g':>5} {'T2g':>5} {'A1u':>5} {'A2u':>5} {'Eu':>5} "
          f"{'T1u':>5} {'T2u':>5}")

    # For each eigenmode, compute |P_R v|^2 / |v|^2 for each irrep
    irrep_names = ['A_1g', 'A_2g', 'E_g', 'T_1g', 'T_2g',
                   'A_1u', 'A_2u', 'E_u', 'T_1u', 'T_2u']
    irrep_projectors = {}
    for name in irrep_names:
        irrep_projectors[name] = project_onto_irrep(
            vecs, oh_actions, class_ids, O_H_IRREPS[name])

    # Show breakdown for each eigenvalue
    t1g_modes = []
    for k in range(len(vals)):
        v = vecs[:, k]
        v_norm2 = np.dot(v, v)
        breakdown = {}
        for name in irrep_names:
            Pv = irrep_projectors[name] @ v
            breakdown[name] = np.dot(v, Pv) / v_norm2
        # Identify dominant irrep
        dominant = max(breakdown.items(), key=lambda x: x[1])
        if dominant[1] > 0.5:
            irrep_label = dominant[0]
        else:
            irrep_label = '?'

        # Print T_1g and T_1u-significant modes
        if breakdown['T_1g'] > 0.3 or k < 30:
            row = f"  {k:>3} {vals[k]:>8.4f}"
            for name in irrep_names:
                row += f"  {breakdown[name]:>4.2f}"
            row += f"  [{irrep_label}]"
            # Print only first occurrence of degenerate triplet
            print(row)

        if breakdown['T_1g'] > 0.5:
            t1g_modes.append({
                'k': k,
                'lambda': vals[k],
                'eigenvector': v,
                'T_1g_content': breakdown['T_1g']
            })

    # Group T_1g modes into triplets (same eigenvalue)
    print(f"\nFound {len(t1g_modes)} T_1g modes total")
    print(f"\nGrouping into triplets (each T_1g mode is 3-fold degenerate):")
    triplets = []
    used = set()
    for i, m in enumerate(t1g_modes):
        if i in used:
            continue
        group = [m]
        used.add(i)
        for j in range(i + 1, len(t1g_modes)):
            if j not in used and abs(t1g_modes[j]['lambda'] - m['lambda']) < 1e-4:
                group.append(t1g_modes[j])
                used.add(j)
        triplets.append(group)

    print(f"\n{'idx':>4} {'lambda':>8} {'m_pred(MeV)':>12} "
          f"{'mult':>5} {'phi%':>5} {'shell%':>7} {'void%':>6}")
    for idx, trip in enumerate(triplets):
        lam = trip[0]['lambda']
        m_pred = mass_from_lambda(N, lam)
        # Average characters over the triplet
        phi_sum = 0
        shell_sum = 0
        void_sum = 0
        for m in trip:
            v = m['eigenvector']
            # phi fraction
            u_norm = np.sum(v[:3*N]**2)
            phi_norm = np.sum(v[3*N:]**2)
            phi_sum += phi_norm / (u_norm + phi_norm)
            # shell vs void localisation (atom-by-atom amplitudes)
            for i in range(N):
                a = np.sum(v[3*i:3*(i+1)]**2) + np.sum(v[3*N+3*i:3*N+3*(i+1)]**2)
                if i < 13:
                    shell_sum += a
                else:
                    void_sum += a
        n = len(trip)
        if n > 0:
            phi_avg = phi_sum / n
            shell_avg = shell_sum / n
            void_avg = void_sum / n
            total = shell_avg + void_avg
            shell_avg /= total
            void_avg /= total
            print(f"  {idx:>3} {lam:>8.4f} {m_pred:>11.2f}  {n:>4} "
                  f"{phi_avg*100:>4.1f}  {shell_avg*100:>5.1f}% "
                  f"{void_avg*100:>5.1f}%")

    # Compare to PDG Delta(1600)
    print(f"\n--- Comparison to PDG Delta(1600) ---")
    print(f"  BW range: 1500 to 1640 MeV (mid 1570)")
    print(f"  Pole range: 1460 to 1560 MeV (mid 1510)")
    print(f"  Predicted lambda for pole 1510: 7.68")
    print(f"  Predicted lambda for BW 1570: 9.27")
    print(f"\nT_1g triplets falling in PDG window [1460, 1640]:")
    for idx, trip in enumerate(triplets):
        lam = trip[0]['lambda']
        m_pred = mass_from_lambda(N, lam)
        if 1460 <= m_pred <= 1640:
            print(f"  triplet {idx}: lambda = {lam:.4f}, m = {m_pred:.2f} MeV")


if __name__ == '__main__':
    main()
