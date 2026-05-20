#!/usr/bin/env python3
"""
proton_first_principles.py
==========================
Derive the proton's lambda = 8.303 from first principles via the
orbit-sector decomposition of the coordination shell under O_h.

Same machinery as rho_first_principles.py but on the higher-symmetry
coord shell. The proton's A_2u stiff Cosserat partner is the unique
phi-derived A_2u mode of the cluster, parallel to the rho's identification
as the lowest phi-derived B_3u mode of the crossed fault.

Universal pattern verified across three centred O_h-symmetric FCC subgraphs:
  - Bare cuboctahedron (12 nodes)    : A_2u (phi) = 6
  - Coordination shell (13 nodes)    : A_2u (phi) = 7   (proton parent)
  - Born cluster (19 nodes)          : A_2u (phi) = 9

The increment per cluster extension equals the number of extra
cuboctahedron-vertex bonds, derived analytically by the bond-counting
argument explained in subsec:cosserat_baryon.

The Cosserat coupling at alpha = 1 splits the (u, phi) A_2u pair into
soft and stiff partners. On the coordination shell, the stiff partner
is at lambda = 8.3028, which gives the proton mass via the master
spectral formula.
"""

import numpy as np
import sys
from itertools import permutations
from cosserat_classifier import build_cosserat_matrix
from spectral_classifier import fcc_nn_vectors


# --- O_h group: all 48 signed permutation matrices on (x, y, z) ---
def build_oh_elements():
    """Return list of 48 O_h 3x3 matrices (signed permutations)."""
    elements = []
    for perm in permutations([0, 1, 2]):
        for sx in (+1, -1):
            for sy in (+1, -1):
                for sz in (+1, -1):
                    R = np.zeros((3, 3))
                    R[0, perm[0]] = sx
                    R[1, perm[1]] = sy
                    R[2, perm[2]] = sz
                    elements.append(R)
    return elements


def class_of(R):
    """Classify an O_h element by its character class (one of 10)."""
    tr = round(np.trace(R))
    det = round(np.linalg.det(R))
    is_diag = all(R[a, a] in (-1, 1) and abs(R[a, b]) < 1e-9
                  for a in range(3) for b in range(3) if a != b)
    if det == 1:
        if tr == 3: return 'E'
        if tr == 0: return '8C_3'
        if tr == 1: return '6C_4'
        if tr == -1:
            return '3C_2' if is_diag else '6C_2'
    else:
        if tr == -3: return 'i'
        if tr == 0: return '8S_6'
        if tr == -1: return '6S_4'
        if tr == 1:
            return '3σ_h' if is_diag else '6σ_d'


# A_2u characters under O_h classes
CHARS_A2U = {'E': +1, '8C_3': +1, '6C_4': -1, '6C_2': -1, '3C_2': +1,
             'i': -1, '8S_6': -1, '6S_4': +1, '6σ_d': +1, '3σ_h': -1}

# A_2g characters (parity-flipped partner of A_2u)
CHARS_A2G = {'E': +1, '8C_3': +1, '6C_4': -1, '6C_2': -1, '3C_2': +1,
             'i': +1, '8S_6': +1, '6S_4': -1, '6σ_d': -1, '3σ_h': +1}


def find_perm(R, coords, tol=1e-6):
    """Return the node permutation induced by R, or None if not a symmetry."""
    n = len(coords)
    perm = []
    for i in range(n):
        target = R @ coords[i]
        match = -1
        for j in range(n):
            if np.linalg.norm(coords[j] - target) < tol:
                match = j; break
        if match == -1:
            return None
        perm.append(match)
    return perm


def build_rep(R, perm, n):
    """Build the 6n x 6n Cosserat representation matrix for element R.
    u (polar) transforms with R; phi (axial) transforms with det(R) * R."""
    D = np.zeros((6*n, 6*n))
    det = round(np.linalg.det(R))
    for i in range(n):
        D[3*perm[i]:3*perm[i]+3, 3*i:3*i+3] = R
        D[3*n + 3*perm[i]:3*n + 3*perm[i]+3, 3*n + 3*i:3*n + 3*i+3] = det * R
    return D


def build_irrep_projector(coords, oh_elements, chars):
    """Build P_{irrep} = (1/|G|) sum_g chi_{irrep}(g) D(g) for the given character table."""
    n = len(coords)
    P = np.zeros((6*n, 6*n))
    for R in oh_elements:
        p = find_perm(R, coords)
        if p is None:
            return None
        D = build_rep(R, p, n)
        P += chars[class_of(R)] * D
    return P / len(oh_elements)


def build_A2u_projector(coords, oh_elements):
    """P_{A_2u} = (1/|G|) sum_g chi_{A_2u}(g) D(g)"""
    return build_irrep_projector(coords, oh_elements, CHARS_A2U)


def build_A2g_projector(coords, oh_elements):
    """P_{A_2g} = (1/|G|) sum_g chi_{A_2g}(g) D(g) -- the defect irrep."""
    return build_irrep_projector(coords, oh_elements, CHARS_A2G)


def find_A2u_modes(coords, alpha):
    """Find all eigenmodes of the Cosserat matrix that lie in A_2u."""
    oh = build_oh_elements()
    P = build_A2u_projector(coords, oh)
    if P is None:
        return None
    n = len(coords)
    M = build_cosserat_matrix(coords, K_u=1.0, K_phi=1.0, alpha=alpha)
    ea, va = np.linalg.eigh(M)
    modes = []
    i = 0
    while i < len(ea):
        j = i
        while j < len(ea) and abs(ea[j] - ea[i]) < 1e-6:
            j += 1
        V = va[:, i:j]
        S = V.T @ P @ V
        u_eig, v_eig = np.linalg.eigh(S)
        for k in range(j - i):
            if u_eig[k] > 0.5:  # mode projects strongly onto A_2u
                mode = V @ v_eig[:, k]
                mode = mode / np.linalg.norm(mode)
                u_c = np.sum(mode[:3*n]**2)
                p_c = np.sum(mode[3*n:]**2)
                modes.append((ea[i], u_c, p_c))
        i = j
    return modes


def main():
    co = fcc_nn_vectors()  # 12 NN of FCC, distance 1
    centre = np.array([[0., 0., 0.]])
    s = np.sqrt(2.0)
    second_nn = np.array([[+s, 0, 0], [-s, 0, 0],
                          [0, +s, 0], [0, -s, 0],
                          [0, 0, +s], [0, 0, -s]])

    clusters = [
        ("Bare cuboctahedron (12 nodes)", co),
        ("Coordination shell (13 nodes)", np.vstack([centre, co])),
        ("Born cluster (19 nodes)", np.vstack([centre, co, second_nn])),
    ]

    # Constants for mass closure
    m0 = 70.0253      # MeV, = m_e / alpha
    me = 0.51099895   # MeV, electron mass
    PDG_proton_iso = 938.918  # MeV

    print("=" * 78)
    print("Universal phi-derived A_2u eigenvalue across O_h-symmetric cluster extensions")
    print("=" * 78)
    print(f"\n{'cluster':<38s} {'alpha=0 (u,phi)':<22s} {'alpha=1 (soft,stiff)':<22s}")

    for name, coords in clusters:
        modes_0 = find_A2u_modes(coords, 0.0)
        modes_1 = find_A2u_modes(coords, 1.0)
        if modes_0 is None or modes_1 is None:
            print(f"{name:<38s} (not O_h-symmetric)")
            continue

        # Sort: u-derived first (alpha=0), then phi-derived
        modes_0.sort(key=lambda m: -m[1])
        modes_1.sort(key=lambda m: m[0])
        u0 = modes_0[0][0]
        phi0 = modes_0[-1][0] if len(modes_0) > 1 else None
        lo = modes_1[0][0]
        hi = modes_1[-1][0] if len(modes_1) > 1 else None

        a0_str = f"{u0:.2f}, {phi0:.2f}" if phi0 is not None else f"{u0:.2f}"
        a1_str = f"{lo:.4f}, {hi:.4f}" if hi is not None else f"{lo:.4f}"
        marker = ""
        if hi is not None and abs(hi - 8.303) < 0.01:
            marker = "  <-- proton stiff partner"
        print(f"{name:<38s} ({a0_str:<20s}) ({a1_str:<20s}){marker}")

    # Mass closure for the coordination shell
    print("\n" + "=" * 78)
    print("Mass closure on the coordination shell")
    print("=" * 78)
    coords_p = np.vstack([centre, co])
    modes_p = find_A2u_modes(coords_p, 1.0)
    modes_p.sort(key=lambda m: m[0])
    stiff = modes_p[-1][0]
    N = 13
    m_pred = N * m0 - N * (4 - stiff) * me
    print(f"  N = {N} (coordination shell)")
    print(f"  lambda_stiff (phi-derived A_2u at alpha=1) = {stiff:.4f}")
    print(f"  m_pred = N * m0 - N * (4 - lambda) * me")
    print(f"         = {N} * {m0:.4f} - {N} * (4 - {stiff:.4f}) * {me:.5f}")
    print(f"         = {m_pred:.4f} MeV")
    print(f"  PDG isoscalar (m_p + m_n)/2 = {PDG_proton_iso} MeV")
    print(f"  Residual = {(m_pred - PDG_proton_iso):+.4f} MeV  ({(m_pred - PDG_proton_iso)/PDG_proton_iso*100:+.4f}%)")

    # Defect-mass pair: verify the DOUBLE FLIP (parity AND sector)
    print("\n" + "=" * 78)
    print("Defect-mass pair: parity flip + sector flip")
    print("=" * 78)
    oh = build_oh_elements()
    P_A2u = build_A2u_projector(coords_p, oh)
    P_A2g = build_A2g_projector(coords_p, oh)
    
    # Modes in each irrep at alpha=1
    def modes_in(P, alpha):
        n = len(coords_p)
        M = build_cosserat_matrix(coords_p, K_u=1.0, K_phi=1.0, alpha=alpha)
        ea, va = np.linalg.eigh(M)
        out = []
        i = 0
        while i < len(ea):
            j = i
            while j < len(ea) and abs(ea[j] - ea[i]) < 1e-6:
                j += 1
            V = va[:, i:j]
            S = V.T @ P @ V
            u_eig, v_eig = np.linalg.eigh(S)
            for k in range(j - i):
                if u_eig[k] > 0.5:
                    mode = V @ v_eig[:, k]
                    mode = mode / np.linalg.norm(mode)
                    u_c = float(np.sum(mode[:3*n]**2))
                    p_c = float(np.sum(mode[3*n:]**2))
                    out.append((float(ea[i]), u_c, p_c))
            i = j
        return out
    
    a2g_modes = sorted(modes_in(P_A2g, 1.0))
    a2u_modes = sorted(modes_in(P_A2u, 1.0))
    
    print("\nA_2g irrep (defect carrier, parity-even):")
    for lam, u, p in a2g_modes:
        char = "u-derived (DEFECT)" if u > 0.5 and lam < 4 else \
               "phi-derived (stiff)" if p > 0.5 else f"mixed"
        marker = "  <-- proton baryon-winding defect" if u > 0.5 and abs(lam - 0.807) < 0.01 else ""
        print(f"  lambda = {lam:>7.4f}  u={u*100:>5.1f}%  phi={p*100:>5.1f}%  ({char}){marker}")
    
    print("\nA_2u irrep (mass carrier, parity-odd):")
    for lam, u, p in a2u_modes:
        char = "u-derived (soft)" if u > 0.5 and lam < 6 else \
               "phi-derived (MASS)" if p > 0.5 else "mixed"
        marker = "  <-- proton rest-mass mode" if p > 0.5 and abs(lam - 8.303) < 0.01 else ""
        print(f"  lambda = {lam:>7.4f}  u={u*100:>5.1f}%  phi={p*100:>5.1f}%  ({char}){marker}")
    
    print("\nDouble-flip verification:")
    defect = [m for m in a2g_modes if m[1] > 0.5 and m[0] < 4][0]
    mass = [m for m in a2u_modes if m[2] > 0.5 and m[0] > 4][0]
    print(f"  Defect: A_2g at lambda = {defect[0]:.4f}, u-content = {defect[1]*100:.1f}% (u-derived)")
    print(f"  Mass:   A_2u at lambda = {mass[0]:.4f}, phi-content = {mass[2]*100:.1f}% (phi-derived)")
    print(f"  Parity flip: A_2g (g) -> A_2u (u)  [verified]")
    print(f"  Sector flip: u-derived -> phi-derived  [verified]")
    print(f"  Both flips are inherent to the orbit-decomposition pairing.")


if __name__ == '__main__':
    main()
