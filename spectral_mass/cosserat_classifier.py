"""
cosserat_classifier.py
======================

Cosserat extension of the spectral classifier.

Builds the 6n x 6n dynamical matrix with both displacement (u) and
microrotation (phi) degrees of freedom at each site.  Decomposes the
eigenspaces under the cluster's symmetry group into integer-spin irreps
of O_h.

The Cosserat free energy on a discrete lattice, in the simplest isotropic
form with central NN couplings, is

    F = F_uu + F_phiphi + F_coup

where

    F_uu       = (K_u/2)  sum_{<ij>}  ( (u_i - u_j) . rhat_ij )^2
                                                        [longitudinal stretch]
    F_phiphi   = (K_phi/2) sum_{<ij>} | phi_i - phi_j |^2
                                                        [isotropic Laplacian]
    F_coup     = (alpha/2) sum_i | phi_i - omega_i |^2
                                                        [Cosserat coupling]

with omega_i = (1/2) (curl u)_i computed discretely.

This module decouples the analysis for clarity:
  - run with alpha = 0 to see the irrep content of each sector separately;
  - run with alpha > 0 to see how the coupling mixes irreps of matching type.

Authors: M. Cox, with Claude (Anthropic).  License: MIT.
"""

import numpy as np
from scipy.linalg import eigh
import sys
from spectral_classifier import (
    fcc_nn_vectors, ELL, A_LAT,
    cluster_hex_cap, cluster_cuboctahedron, cluster_coord_shell, cluster_born,
    generate_Oh, classify_operation, vertex_perm,
    OH_CHARACTERS, OH_CLASSES, OH_CLASS_SIZES,
    D3D_CHARACTERS, D3D_CLASSES, D3D_CLASS_SIZES,
    group_eigvalues,
)


# ============================================================================
# Build Cosserat dynamical matrix (block form)
# ============================================================================

def build_uu_block(coords, K=1.0, ell=ELL, tol=1e-6):
    """Displacement-displacement block: central NN springs.

    Phi_uu[i*3+a, j*3+b] = -K rhat^a rhat^b  for bonded i != j
    Phi_uu[i*3+a, i*3+b] = +K rhat^a rhat^b  summed over j in NN(i)
    """
    n = len(coords)
    P = np.zeros((3*n, 3*n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            r = coords[j] - coords[i]
            d = np.linalg.norm(r)
            if abs(d - ell) < tol:
                rhat = r / d
                M = np.outer(rhat, rhat)
                P[3*i:3*i+3, 3*j:3*j+3] -= K * M
                P[3*i:3*i+3, 3*i:3*i+3] += K * M
    return P


def build_phiphi_block(coords, K=1.0, ell=ELL, tol=1e-6):
    """Microrotation-microrotation block: isotropic NN springs.

    Phi_phiphi[i*3+a, j*3+b] = -K delta^{ab}  for bonded i != j
    Phi_phiphi[i*3+a, i*3+b] = +K delta^{ab}  summed over j in NN(i)

    (This is the graph Laplacian times identity-3.)
    """
    n = len(coords)
    P = np.zeros((3*n, 3*n))
    I3 = np.eye(3)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            r = coords[j] - coords[i]
            d = np.linalg.norm(r)
            if abs(d - ell) < tol:
                P[3*i:3*i+3, 3*j:3*j+3] -= K * I3
                P[3*i:3*i+3, 3*i:3*i+3] += K * I3
    return P


def build_coupling_block(coords, alpha=1.0, ell=ELL, tol=1e-6):
    """Cosserat displacement-microrotation coupling: alpha * (phi - omega)^2.

    omega_i = (1/2) sum_{j in NN(i)} (rhat_ij x (u_j - u_i)) / |r_ij|
            ~ discrete curl evaluated at site i.

    For simplicity we use omega_i = (1/(2*ell)) sum_{j in NN(i)} (rhat_ij x u_j).

    The block matrix entries are:
        Phi_uu coupling correction:  +alpha * (omega contributions)^2
        Phi_phiphi diagonal:         +alpha * I
        Phi_u_phi off-diagonal:      -alpha * (omega operator)
    """
    n = len(coords)
    # Build the curl operator C such that omega_i^a = C[i*3+a, j*3+b] u_j^b
    C = np.zeros((3*n, 3*n))
    for i in range(n):
        # Find NN of i
        nn_indices = []
        nn_rhats = []
        for j in range(n):
            if i == j:
                continue
            r = coords[j] - coords[i]
            d = np.linalg.norm(r)
            if abs(d - ell) < tol:
                nn_indices.append(j)
                nn_rhats.append(r / d)
        # omega_i^a = (1/(2*ell)) sum_{j in NN} epsilon^{abc} rhat^b u_j^c
        for j_idx, j in enumerate(nn_indices):
            rh = nn_rhats[j_idx]
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        # epsilon^{abc}
                        eps = 0.0
                        if (a, b, c) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
                            eps = 1.0
                        elif (a, b, c) in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
                            eps = -1.0
                        if eps != 0:
                            C[3*i + a, 3*j + c] += eps * rh[b] / (2.0 * ell)

    # The coupling adds:
    #   (alpha/2) (phi - C u)^T (phi - C u)
    #   = (alpha/2) phi^T phi - alpha phi^T C u + (alpha/2) u^T C^T C u
    # So contributions to full 6n x 6n matrix:
    #   uu block:    +alpha C^T C
    #   phiphi block:+alpha I
    #   u-phi block: -alpha C^T  (acting on phi to give u-side contribution)
    #   phi-u block: -alpha C
    uu_corr = alpha * (C.T @ C)
    phiphi_corr = alpha * np.eye(3*n)
    u_phi_corr = -alpha * C.T   # 3n x 3n, [u row, phi col]
    phi_u_corr = -alpha * C
    return uu_corr, phiphi_corr, u_phi_corr, phi_u_corr


def build_cosserat_matrix(coords, K_u=1.0, K_phi=1.0, alpha=0.0):
    """Build the full 6n x 6n Cosserat dynamical matrix."""
    n = len(coords)
    Phi = np.zeros((6*n, 6*n))
    # u-u block (top-left 3n x 3n)
    Phi[:3*n, :3*n] = build_uu_block(coords, K=K_u)
    # phi-phi block (bottom-right 3n x 3n)
    Phi[3*n:, 3*n:] = build_phiphi_block(coords, K=K_phi)
    if alpha > 0:
        uu_c, pp_c, up_c, pu_c = build_coupling_block(coords, alpha=alpha)
        Phi[:3*n, :3*n] += uu_c
        Phi[3*n:, 3*n:] += pp_c
        Phi[:3*n, 3*n:] += up_c
        Phi[3*n:, :3*n] += pu_c
    return Phi


# ============================================================================
# Cosserat displacement representation: action of O_h
# ============================================================================

def cosserat_rep_matrix(R, coords, tol=1e-6):
    """6n x 6n matrix representing R acting on (u, phi) space.

    u transforms as T_1u:  (g u)_i^a = R^{ab} u_{perm^{-1}(i)}^b  
    phi transforms as T_1g: (g phi)_i^a = (det R) R^{ab} phi_{perm^{-1}(i)}^b

    The factor (det R) distinguishes T_1g from T_1u: under inversion, u flips sign,
    phi does not.
    """
    perm = vertex_perm(R, coords, tol)
    if perm is None:
        return None
    n = len(coords)
    detR = np.linalg.det(R)
    M = np.zeros((6*n, 6*n))
    for i in range(n):
        j = perm[i]
        # u block
        M[3*j:3*j+3, 3*i:3*i+3] = R
        # phi block (T_1g: extra factor of det R)
        M[3*n + 3*j:3*n + 3*j+3, 3*n + 3*i:3*n + 3*i+3] = detR * R
    return M


def cluster_subgroup(coords, oh_elements):
    return [R for R in oh_elements if vertex_perm(R, coords) is not None]


def decompose_cosserat(eigvecs, indices, coords, subgroup,
                       characters, classes, class_sizes):
    """Apply character formula to the 6n-dim eigenspace using the Cosserat
    representation matrices."""
    V = eigvecs[:, indices]
    chi_by_class = {c: 0.0 for c in classes}
    for R in subgroup:
        M = cosserat_rep_matrix(R, coords)
        if M is None:
            continue
        chi_g = np.trace(V.T @ M @ V)
        cls = classify_operation(R)
        if cls in chi_by_class:
            chi_by_class[cls] += chi_g
    H = sum(class_sizes.values())
    out = {}
    for ir, ch in characters.items():
        n_ir = sum(ch[c] * chi_by_class[c] for c in classes) / H
        if abs(n_ir) > 0.01:
            out[ir] = round(n_ir)
    return out


# ============================================================================
# Driver
# ============================================================================

def report(builder, alpha, oh):
    coords, name = builder()
    sub = cluster_subgroup(coords, oh)
    H = len(sub)
    Phi = build_cosserat_matrix(coords, K_u=1.0, K_phi=1.0, alpha=alpha)
    ev, ec = eigh(Phi)
    groups = group_eigvalues(ev)

    if H == 48:
        chars, classes, sizes, sym = OH_CHARACTERS, OH_CLASSES, OH_CLASS_SIZES, 'O_h'
    elif H == 12:
        chars, classes, sizes, sym = D3D_CHARACTERS, D3D_CLASSES, D3D_CLASS_SIZES, 'D_3d'
    else:
        chars = None
        sym = f'|H|={H}'

    print(f'\n{name}: n={len(coords)}, |H|={H}, sym={sym}, alpha={alpha}')
    print(f'  Cosserat dim: 6n = {6 * len(coords)}, spectrum '
          f'[{ev.min():.4f}, {ev.max():.4f}]')
    print(f'  {"lambda":>10s}  {"mult":>4s}  irrep content')
    print(f'  {"-"*10}  {"-"*4}  {"-"*60}')
    irrep_totals = {}
    for val, idx in groups:
        if chars is not None:
            d = decompose_cosserat(ec, idx, coords, sub, chars, classes, sizes)
            dstr = ' + '.join([f'{n}*{ir}' if n > 1 else ir for ir, n in d.items()])
            for ir, n in d.items():
                irrep_totals[ir] = irrep_totals.get(ir, 0) + n
        else:
            dstr = f'(low-sym, |H|={H})'
        marker = ' [ZERO]' if abs(val) < 1e-5 else ''
        print(f'  {val:8.4f}{marker:7s}  {len(idx):4d}  {dstr}')

    if chars is not None:
        # Print total irrep content
        print(f'\n  Total irrep multiplicities (Cosserat = u + phi):')
        sorted_irreps = ['A_1g', 'A_2g', 'E_g', 'T_1g', 'T_2g',
                        'A_1u', 'A_2u', 'E_u', 'T_1u', 'T_2u']
        for ir in sorted_irreps:
            if ir in irrep_totals:
                print(f'     {ir}: {irrep_totals[ir]}')


def main():
    print('=' * 78)
    print('Cosserat extension: 6n x 6n dynamical matrix')
    print('Displacement + microrotation, with optional coupling alpha')
    print('=' * 78)

    oh = generate_Oh()
    print(f'\n|O_h| = {len(oh)}\n')

    # Decoupled (alpha = 0)
    print('\n*** DECOUPLED CASE: alpha = 0 ***')
    for builder in [cluster_cuboctahedron, cluster_coord_shell]:
        report(builder, alpha=0.0, oh=oh)

    print('\n\n*** COUPLED CASE: alpha = 1 (Cosserat) ***')
    report(cluster_coord_shell, alpha=1.0, oh=oh)


if __name__ == '__main__':
    main()
