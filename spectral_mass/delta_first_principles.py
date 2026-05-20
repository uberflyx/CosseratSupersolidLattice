#!/usr/bin/env python3
"""
delta_first_principles.py
=========================

First-principles spectral derivation of the Delta(1232) baryon mass on
the void-activated coordination shell, in the Cosserat-FCC framework.

Cluster
-------
The Delta sits on the proton's 13-node coordination shell (centre + 12
cuboctahedron vertices) augmented by four FCC tetrahedral voids of one
T_d orbit, giving 17 nodes total.  Each void is the centroid of a
{origin, three NN} tetrahedron, sitting at distance sqrt(6)/4 from each
of its four corners.

The void orbit choice breaks the cluster's symmetry from O_h to T_d:
inversion swaps the two T_d orbits of voids, so the cluster activated
for the Delta lacks inversion.  Its point group is T_d (order 24).

Cosserat dynamical matrix
-------------------------
Built with central NN springs of unit stiffness on the 42 d=1 bonds
(36 within-shell + 6 void-void) and on the 16 void-corner bonds at
d = sqrt(6)/4.  The Cosserat coupling alpha = 1 is used throughout, in
the same convention as proton_first_principles.py and
rho_first_principles.py.

T_d irrep classification
------------------------
All 102 eigenmodes are projected onto the five T_d irreps using the
canonical isotypic projectors P_R = (d_R/|G|) sum_g chi_R(g) D(g).
Multiplicities (per orthogonality and verified by trace):
    A_1: 4 copies,  A_2: 4 copies,  E: 8 copies,
    T_1: 13 copies, T_2: 13 copies   (total 4 + 4 + 16 + 39 + 39 = 102).

Delta identification
--------------------
The Delta is J^P = 3/2^+, isospin 3/2.  The mass mode lives in a
T_1 stiff phi-dominant eigenvector (T_1 inherits from T_{1g} of O_h,
the microrotation-block irrep on the parent shell).

T_1 has 13 distinct eigenvalues at alpha = 1.  Four are stiff
(lambda > 4) and phi-dominant (>= 90%), at lambda in {6.66, 9.05,
9.91, 18.26}.  Two structural criteria fix the rest-mass mode:
  (i)  shell-concentrated (void localisation < 5%) -- inheritance
       from the proton's bare-shell T_{1g} phi-block at lambda = 8.07.
  (ii) lowest such root -- mirrors the proton's selection of
       lambda^{A_2u}_+ = 8.303 over higher A_2u stiff partners.
Of the four candidates the void-localisations are {21.7, 0.1, 79.1,
2.5}%, so (i) keeps {9.05, 18.26} and (ii) picks

    lambda^{T_1}_+ = 9.0515   (phi 99.2%, void 0.1%).
Master formula:

    m_Delta = N m_0 - N (4 - lambda) m_e
            = 17 (m_e / alpha) - 17 (4 - 9.0515) m_e
            = 1234.31 MeV   against PDG 1232 +/- 2 MeV,
            residual +0.19%.
"""

import numpy as np
from scipy.linalg import eigh
import sys
from spectral_classifier import fcc_nn_vectors, ELL, A_LAT


# ============================================================================
# Cluster geometry
# ============================================================================

def void_positions_Td():
    """The 4 FCC tetrahedral voids of one T_d orbit, at the centroid of the
    tetrahedron (origin, three NN sites).  Signs have product +1."""
    signs = [(+1, +1, +1), (+1, -1, -1), (-1, +1, -1), (-1, -1, +1)]
    return np.array(signs, dtype=float) * (A_LAT / 4.0)


def cluster_delta():
    """17-node Delta cluster: centre + 12 cuboctahedron + 4 voids."""
    centre = np.array([[0., 0., 0.]])
    shell  = fcc_nn_vectors()
    voids  = void_positions_Td()
    return np.vstack([centre, shell, voids])


# ============================================================================
# Mixed-length Cosserat dynamical matrix
# ============================================================================
# The cluster has two distinct bond lengths:
#   d = 1            (intra-shell NN, including void-void at this distance)
#   d = sqrt(6)/4    (void-to-corner)
# Both bond types carry the same K_u = K_phi = 1 and alpha = 1.
# The Cosserat curl operator omega = (1/(2d)) sum_NN rhat x u picks up the
# 1/d factor per bond.

def build_cosserat_matrix_two_d(coords, K_u=1.0, K_phi=1.0, alpha=1.0,
                                 ell_nn=ELL, void_d=None, tol=1e-6):
    if void_d is None:
        void_d = np.sqrt(6.0) / 4.0
    n = len(coords)
    P = np.zeros((6 * n, 6 * n))
    I3 = np.eye(3)
    bonds = []
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(coords[i] - coords[j])
            if abs(d - ell_nn) < tol or abs(d - void_d) < tol:
                bonds.append((i, j, d))

    # Central-force spring on u, isotropic Laplacian on phi.
    Phi_uu = P[:3*n, :3*n]
    Phi_pp = P[3*n:, 3*n:]
    for (i, j, d) in bonds:
        rhat = (coords[j] - coords[i]) / d
        outer = np.outer(rhat, rhat)
        Phi_uu[3*i:3*i+3, 3*j:3*j+3] -= K_u * outer
        Phi_uu[3*j:3*j+3, 3*i:3*i+3] -= K_u * outer
        Phi_uu[3*i:3*i+3, 3*i:3*i+3] += K_u * outer
        Phi_uu[3*j:3*j+3, 3*j:3*j+3] += K_u * outer
        Phi_pp[3*i:3*i+3, 3*j:3*j+3] -= K_phi * I3
        Phi_pp[3*j:3*j+3, 3*i:3*i+3] -= K_phi * I3
        Phi_pp[3*i:3*i+3, 3*i:3*i+3] += K_phi * I3
        Phi_pp[3*j:3*j+3, 3*j:3*j+3] += K_phi * I3

    # Discrete-curl Cosserat coupling: alpha * |phi - omega|^2 with
    # omega_i = (1/(2 d)) sum_NN rhat_{ij} x u_j (per-bond contributions).
    C = np.zeros((3*n, 3*n))
    for (i, j, d) in bonds:
        rhat = (coords[j] - coords[i]) / d
        for (k, m, sign) in [(i, j, +1), (j, i, -1)]:
            rh = sign * rhat
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        eps = 0.0
                        if (a, b, c) in [(0,1,2), (1,2,0), (2,0,1)]: eps = +1.0
                        elif (a, b, c) in [(0,2,1), (2,1,0), (1,0,2)]: eps = -1.0
                        if eps != 0.0:
                            C[3*k+a, 3*m+c] += eps * rh[b] / (2.0 * d)

    P[:3*n, :3*n] += alpha * (C.T @ C)
    P[3*n:, 3*n:] += alpha * np.eye(3*n)
    P[:3*n, 3*n:]  = -alpha * C.T
    P[3*n:, :3*n]  = -alpha * C
    return P


# ============================================================================
# T_d group machinery
# ============================================================================

def generate_Td():
    """T_d as the subgroup of O_h preserving the positive-orbit void
    tetrahedron, found by direct search over the 48 signed permutations."""
    from spectral_classifier import generate_Oh
    pos_voids = void_positions_Td()
    Td = []
    for R in generate_Oh():
        new_voids = (R @ pos_voids.T).T
        if all(np.sign(np.prod(v)) > 0.5 for v in new_voids):
            Td.append(R)
    return Td


TD_CHARS = {  # T_d character table; key is class label, value is character
    'A_1': {'E': +1, '8C_3': +1, '3C_2': +1, '6sigma_d': +1, '6S_4': +1},
    'A_2': {'E': +1, '8C_3': +1, '3C_2': +1, '6sigma_d': -1, '6S_4': -1},
    'E':   {'E': +2, '8C_3': -1, '3C_2': +2, '6sigma_d':  0, '6S_4':  0},
    'T_1': {'E': +3, '8C_3':  0, '3C_2': -1, '6sigma_d': -1, '6S_4': +1},
    'T_2': {'E': +3, '8C_3':  0, '3C_2': -1, '6sigma_d': +1, '6S_4': -1},
}


def classify_Td(R, tol=1e-4):
    """Identity, 8C_3, 3C_2, 6sigma_d, 6S_4 by (det, trace)."""
    det = round(np.linalg.det(R))
    tr  = np.trace(R)
    if det == +1:
        if abs(tr - 3) < tol: return 'E'
        if abs(tr - 0) < tol: return '8C_3'
        if abs(tr + 1) < tol: return '3C_2'
    else:
        if abs(tr - 1) < tol: return '6sigma_d'
        if abs(tr + 1) < tol: return '6S_4'
    return 'unknown'


def vertex_perm(R, coords, tol=1e-6):
    n = len(coords)
    out = []
    for i in range(n):
        target = R @ coords[i]
        match = -1
        for j in range(n):
            if np.linalg.norm(coords[j] - target) < tol:
                match = j; break
        if match == -1:
            return None
        out.append(match)
    return out if len(set(out)) == n else None


def build_cosserat_rep(R, coords, tol=1e-6):
    """6n x 6n: u transforms as R, phi as det(R) R."""
    perm = vertex_perm(R, coords, tol)
    if perm is None:
        return None
    n = len(coords)
    det = round(np.linalg.det(R))
    D = np.zeros((6*n, 6*n))
    for i in range(n):
        j = perm[i]
        D[3*j:3*j+3, 3*i:3*i+3] = R
        D[3*n + 3*j:3*n + 3*j+3, 3*n + 3*i:3*n + 3*i+3] = det * R
    return D


def build_td_projector(coords, Td, chars):
    """Canonical isotypic projector P_R = (d_R/|G|) sum_g chi_R(g) D(g)."""
    n = len(coords)
    P = np.zeros((6*n, 6*n))
    d_R = chars['E']
    for R in Td:
        D = build_cosserat_rep(R, coords)
        if D is None:
            return None
        P += chars[classify_Td(R)] * D
    return (d_R / len(Td)) * P


def project_modes(coords, projector, alpha):
    """List eigenmodes that live in the irrep, with u/phi sector fractions."""
    M = build_cosserat_matrix_two_d(coords, K_u=1.0, K_phi=1.0, alpha=alpha)
    ea, va = eigh(M)
    n = len(coords)
    out = []
    i = 0
    while i < len(ea):
        j = i
        while j < len(ea) and abs(ea[j] - ea[i]) < 1e-6:
            j += 1
        V = va[:, i:j]
        S = V.T @ projector @ V
        u_eig, v_eig = eigh(S)
        for k in range(j - i):
            if u_eig[k] > 0.5:
                mode = V @ v_eig[:, k]
                mode = mode / np.linalg.norm(mode)
                u_c = float(np.sum(mode[:3*n]**2))
                p_c = float(np.sum(mode[3*n:]**2))
                out.append((float(ea[i]), u_c, p_c))
        i = j
    return out


# ============================================================================
# Main: derive m_Delta
# ============================================================================
def main():
    # CODATA 2022
    M_E   = 0.510_998_950_69      # MeV
    ALPHA = 1.0 / 137.035_999_177
    M_0   = M_E / ALPHA           # = 70.0254 MeV
    PDG_DELTA = 1232.0            # MeV, Breit-Wigner

    coords = cluster_delta()
    n = len(coords)
    Td = generate_Td()

    print("=" * 78)
    print("Delta(1232) spectral closure on the void-extended coordination shell")
    print("=" * 78)
    print(f"\nCluster:")
    print(f"  N = {n} nodes  (13 shell + 4 T_d-orbit voids)")
    print(f"  Point group: T_d  (|Td| = {len(Td)})")
    print(f"  Bond inventory:")
    bonds_short = sum(1 for i in range(n) for j in range(i+1, n)
                      if abs(np.linalg.norm(coords[i]-coords[j]) - 1.0) < 1e-6)
    bonds_void  = sum(1 for i in range(n) for j in range(i+1, n)
                      if abs(np.linalg.norm(coords[i]-coords[j]) - np.sqrt(6)/4) < 1e-6)
    print(f"    d = 1 (NN, intra-shell + void-void): {bonds_short} bonds")
    print(f"    d = sqrt(6)/4 (void-corner):         {bonds_void} bonds")
    print(f"  Cosserat coupling: K_u = K_phi = alpha = 1")

    # Verify T_d projector traces
    print(f"\nT_d irrep multiplicities (n_Γ = tr(P_Γ)/dim Γ):")
    total = 0
    for irrep, chars in TD_CHARS.items():
        P = build_td_projector(coords, Td, chars)
        tr = np.trace(P).real
        total += tr
        print(f"  {irrep:>4s}: tr(P) = {tr:5.1f}, n = {int(round(tr/chars['E']))} copies")
    print(f"  Total dim = {total:.0f} (= 6n = {6*n})")

    # ------------------------------------------------------------------
    # T_1 mode classification with void-localisation
    # ------------------------------------------------------------------
    # The Delta mass lives in a stiff phi-dominant T_1 mode.  Four such
    # modes exist (lambda > 4, phi >= 0.9).  Two structural criteria pin
    # the rest-mass mode uniquely:
    #
    #  (i)  shell-concentrated (void localisation < 5%): the mass partner
    #       must inherit from the proton's coordination shell, the same
    #       way A_2u on the shell inherits from a shell-mode of T_{1g}
    #       under the O_h -> T_d subduction.
    #  (ii) lowest stiff shell-concentrated root: mirrors the proton's
    #       selection of lambda^{A_2u}_+ = 8.303 over higher A_2u stiff
    #       roots on the closed shell.
    # ------------------------------------------------------------------
    P_T1 = build_td_projector(coords, Td, TD_CHARS['T_1'])
    M    = build_cosserat_matrix_two_d(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    ea, va = eigh(M)
    void_idx = list(range(13, 17))   # last 4 sites are the activated voids

    print(f"\nT_1 irrep, phi-dominant stiff modes (lambda > 4, phi >= 0.9):")
    print(f"  {'lambda':>10s}  {'phi %':>6s}  {'void-loc %':>10s}  {'shell-loc %':>11s}")
    print(f"  {'-'*10}  {'-'*6}  {'-'*10}  {'-'*11}")
    candidates = []
    seen = set()
    i = 0
    while i < len(ea):
        j = i
        while j < len(ea) and abs(ea[j] - ea[i]) < 1e-6:
            j += 1
        V = va[:, i:j]
        S = V.T @ P_T1 @ V
        eigs, vecs = eigh(S)
        for k in range(j - i):
            if eigs[k] > 0.5:
                mode = V @ vecs[:, k]
                mode /= np.linalg.norm(mode)
                u_part = mode[:3*n]
                p_part = mode[3*n:]
                p_c = float(np.sum(p_part**2))
                void_amp = sum(np.sum(u_part[3*v:3*v+3]**2) for v in void_idx) \
                         + sum(np.sum(p_part[3*v:3*v+3]**2) for v in void_idx)
                if ea[i] > 4 and p_c > 0.9:
                    key = round(ea[i], 4)
                    if key not in seen:
                        seen.add(key)
                        candidates.append((float(ea[i]), p_c, void_amp))
        i = j
    candidates.sort()
    for lam, p, void in candidates:
        shell = 1.0 - void
        m_pred = n * M_0 - n * (4 - lam) * M_E
        print(f"  {lam:8.4f}  {p*100:6.1f}  {void*100:10.3f}  {shell*100:11.3f}    m = {m_pred:7.2f} MeV")

    # Apply selection rule
    shell_cands = [(lam, p, void) for (lam, p, void) in candidates if void < 0.05]
    delta_lam, _, _ = sorted(shell_cands)[0]      # lowest shell-concentrated

    # ------------------------------------------------------------------
    # Structural inheritance check: shell-only T_1g stiff phi-block
    # ------------------------------------------------------------------
    # The Delta's T_1 mass mode descends from a T_{1g} phi-block on the
    # bare 13-node shell.  Verify the parent eigenvalue is at lambda ~ 8.07.
    from cosserat_classifier import build_cosserat_matrix
    from spectral_classifier import generate_Oh, OH_CHARACTERS, vertex_perm as vp_oh, classify_operation as cls_oh
    shell_coords = np.vstack([np.array([[0., 0., 0.]]), fcc_nn_vectors()])
    nS = len(shell_coords)
    Ms = build_cosserat_matrix(shell_coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    ea_s, va_s = eigh(Ms)
    # T_{1g} projector on O_h, with d_R = 3
    PT1g = np.zeros((6*nS, 6*nS))
    for R in generate_Oh():
        perm = vp_oh(R, shell_coords)
        if perm is None: continue
        det = round(np.linalg.det(R))
        D = np.zeros((6*nS, 6*nS))
        for ii in range(nS):
            jj = perm[ii]
            D[3*jj:3*jj+3, 3*ii:3*ii+3] = R
            D[3*nS+3*jj:3*nS+3*jj+3, 3*nS+3*ii:3*nS+3*ii+3] = det * R
        PT1g += OH_CHARACTERS['T_1g'][cls_oh(R)] * D
    PT1g *= 3.0 / 48.0
    # Pick stiff phi-dominant T_1g modes on the shell
    print(f"\nStructural inheritance: shell-only T_{{1g}} stiff phi-dominant modes")
    i = 0
    shell_t1g_phi = []
    while i < len(ea_s):
        j = i
        while j < len(ea_s) and abs(ea_s[j] - ea_s[i]) < 1e-6:
            j += 1
        V = va_s[:, i:j]
        S = V.T @ PT1g @ V
        eigs, vecs = eigh(S)
        for k in range(j - i):
            if eigs[k] > 0.5:
                mode = V @ vecs[:, k]
                mode /= np.linalg.norm(mode)
                p_c = float(np.sum(mode[3*nS:]**2))
                if p_c > 0.9 and ea_s[i] > 4 and round(ea_s[i], 4) not in {round(x, 4) for x, _, _ in shell_t1g_phi}:
                    shell_t1g_phi.append((float(ea_s[i]), p_c, None))
                    print(f"  shell lambda^{{T_1g}}_+ = {ea_s[i]:.4f}   phi = {p_c*100:.1f}%")
        i = j

    # ------------------------------------------------------------------
    # Master-formula closure
    # ------------------------------------------------------------------
    m_pred = n * M_0 - n * (4 - delta_lam) * M_E
    resid_mev = m_pred - PDG_DELTA
    resid_pct = resid_mev / PDG_DELTA * 100
    print(f"\nDelta assignment: lowest shell-concentrated stiff phi-dominant T_1 mode")
    print(f"  lambda^{{T_1}}_+ = {delta_lam:.4f}")
    print(f"  m_Delta = N m_0 - N (4 - lambda) m_e")
    print(f"          = {n} x {M_0:.4f} - {n} x (4 - {delta_lam:.4f}) x {M_E:.5f}")
    print(f"          = {m_pred:.2f} MeV")
    print(f"  PDG     = {PDG_DELTA} +/- 2 MeV (Breit-Wigner)")
    print(f"  Residual = {resid_mev:+.2f} MeV  ({resid_pct:+.3f}%)")

    # ------------------------------------------------------------------
    # Sanity: T_2 phi-dominant counterpart and singlet pair check
    # ------------------------------------------------------------------
    print(f"\nFirst stiff phi-dominant mode in each other irrep (for context):")
    for irrep in ['T_2', 'E', 'A_2', 'A_1']:
        P = build_td_projector(coords, Td, TD_CHARS[irrep])
        modes = sorted(project_modes(coords, P, alpha=1.0))
        for (lam, u, p) in modes:
            if lam > 4 and p > 0.5:
                m = n * M_0 - n * (4 - lam) * M_E
                resid = (m - PDG_DELTA) / PDG_DELTA * 100
                print(f"  {irrep:>4s}: lambda = {lam:.4f}  phi = {p*100:.1f}%  m = {m:.2f} MeV  ({resid:+.3f}%)")
                break


if __name__ == "__main__":
    main()
