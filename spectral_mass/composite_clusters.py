"""
composite_clusters.py
=====================

Spectral analysis of composite hadron clusters.  Builds the dynamical matrix
on each composite as a single graph, computes the spectrum, and decomposes
under whatever subgroup of O_h preserves the cluster.

Composites treated:
  - Rho meson: 11-vertex crossed fault (centre + 2 intersecting hex equators)
  - Phi meson: GMO via 14- or 15-vertex strange-strange construction
  - Lambda baryon: 16-vertex shell + hex cap
  - K* meson: 13-vertex strangeness-pinned coord shell (same as shell; no composition)

Authors: M. Cox, with Claude (Anthropic).  License: MIT.
"""

import numpy as np
from scipy.linalg import eigh
import sys
from spectral_classifier import (
    fcc_nn_vectors, ELL, A_LAT,
    generate_Oh, classify_operation, vertex_perm,
    OH_CHARACTERS, OH_CLASSES, OH_CLASS_SIZES,
    group_eigvalues, build_dynamical_matrix, disp_rep,
    cluster_subgroup,
)


# ============================================================================
# Cluster builders for composite hadrons
# ============================================================================

def cluster_crossed_fault():
    """The rho's 11-node crossed-fault cluster.

    Construction: centre + 2 intersecting {111} hex equators on the cuboctahedron.
    Two equators share 2 vertices, giving 6 + 6 - 2 = 10 ring vertices + 1
    centre = 11 vertices total.

    We choose the (111) and (-1,1,1) planes; the shared axis is along (0,-1,1).
    """
    centre = np.array([[0., 0., 0.]])
    co = fcc_nn_vectors()
    n1 = np.array([1., 1., 1.])
    n2 = np.array([-1., 1., 1.])
    tol = 1e-6
    on1 = [v for v in co if abs(n1 @ v) < tol]
    on2 = [v for v in co if abs(n2 @ v) < tol]
    # Deduplicated union
    union = list(on1)
    for v in on2:
        if not any(np.linalg.norm(v - u) < tol for u in union):
            union.append(v)
    coords = np.vstack([centre, *union])
    assert len(coords) == 11, f"Expected 11 vertices, got {len(coords)}"
    return coords, "Rho crossed fault (11 nodes)"


def cluster_lambda():
    """The Lambda baryon's 16-node cluster: coord shell + hex cap.

    The hex cap adds 3 new vertices to the coord shell (the 6 ring vertices
    of a (111) hex cap that lie on a parallel plane half a lattice constant
    away).  These 3 vertices are the second-shell {111}-stacking sites.

    For the spectral content we use the simplest construction: coord shell
    (13 vertices, full O_h) + 3 hex cap apex vertices at second-NN distance
    along (1,1,1)/sqrt(3) direction, displaced by ell.
    """
    centre = np.array([[0., 0., 0.]])
    co = fcc_nn_vectors()
    # The cuboctahedron has a (111) equator with 6 vertices.  The next layer up
    # (along +(1,1,1)/sqrt(3) by ell*sqrt(2/3)) has 3 vertices that complete
    # the hex cap.  These 3 vertices have positions a/2 * permutations of (1,1,0)
    # shifted by ell along (1,1,1).  Simpler construction: place 3 nodes at
    # the cuboctahedron's second-shell {111}-positions.
    a = A_LAT
    apex_offsets = (a/2) * np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]], dtype=float)
    # Shift these along +z to put them just above the cuboctahedron
    # Actually for proper hex-cap geometry on (111) plane: the 3 second-layer
    # vertices sit at (a/2)*(1,1,0) + ell*(1,1,1)/sqrt(3)... we take simpler form.
    # Use the (1,1,1)/sqrt(3) shift by ell:
    shift = ELL * np.array([1, 1, 1]) / np.sqrt(3)
    apex = apex_offsets + shift
    coords = np.vstack([centre, co, apex])
    return coords, "Lambda (shell + 3-apex hex cap, 16 nodes)"


def cluster_phi_gmo():
    """The phi's GMO 14-node construction.

    The framework's phi uses (29/2) m_0 = 14.5 m_0.  Spectrally, this is
    14 shell-shell nodes plus 1/2 from the J=1 spin self-energy on the shell.
    The composite is essentially the coord shell + 1 strange apex (analog of
    eta-prime coord shell + apex), with the apex doubled because phi has
    s-sbar content.

    For simplicity we use coord shell + 1 apex along +(0,0,1):
    """
    centre = np.array([[0., 0., 0.]])
    co = fcc_nn_vectors()
    a = A_LAT
    # Single apex along +z at second-NN distance
    apex = np.array([[0., 0., a]])
    coords = np.vstack([centre, co, apex])
    return coords, "Phi (coord shell + 1 apex, 14 nodes)"


# ============================================================================
# Spectral analysis with reduced symmetry decomposition
# ============================================================================

def report_composite(builder, oh_elements):
    """Compute the dynamical-matrix spectrum and decompose under the cluster's
    subgroup of O_h.  Highlights T_1u (vector), A_2u (pseudoscalar), and
    E_g/T_2g (tensor) content."""
    coords, name = builder()
    sub = cluster_subgroup(coords, oh_elements)
    H = len(sub)
    n = len(coords)
    Phi = build_dynamical_matrix(coords, K_nn=1.0)
    ev, ec = eigh(Phi)
    groups = group_eigvalues(ev)

    # Identify class composition of the subgroup
    cls_count = {}
    for R in sub:
        c = classify_operation(R)
        cls_count[c] = cls_count.get(c, 0) + 1

    print(f'\n{name}')
    print(f'  n = {n}, |H| = {H}')
    print(f'  Subgroup classes: ' + ', '.join(f'{c}:{k}' for c, k in cls_count.items()))
    print(f'  Dynamical-matrix dim: 3n = {3*n}')
    print(f'  Spectrum range: [{ev.min():.4f}, {ev.max():.4f}]')

    # Compute character of each eigenspace under the subgroup
    print(f'\n  {"lambda/K":>10s}  {"mult":>4s}  irrep content (where |H|=48 allows decomposition)')
    print(f'  {"-"*10}  {"-"*4}  {"-"*60}')
    for val, idx in groups:
        if H == 48:
            decomp = _decompose_oh(ec, idx, coords, sub)
            dstr = ' + '.join([f'{m}*{ir}' if m > 1 else ir for ir, m in decomp.items()])
        else:
            # Lower symmetry: report eigenspace multiplicity and parity content
            decomp_info = _eigenspace_parity_content(ec, idx, coords, sub)
            dstr = decomp_info
        marker = ' [ZERO]' if abs(val) < 1e-5 else ''
        print(f'  {val:8.4f}{marker:7s}  {len(idx):4d}  {dstr}')

    # Identify likely vector mode (J = 1^-)
    print(f'\n  Identifying T_1u (J^P = 1^-) vector modes:')
    _identify_vector_modes(ec, ev, groups, coords, sub)


def _decompose_oh(eigvecs, indices, coords, subgroup):
    """Full O_h decomposition."""
    V = eigvecs[:, indices]
    chi = {c: 0.0 for c in OH_CLASSES}
    for R in subgroup:
        M = disp_rep(R, coords)
        if M is None:
            continue
        chi_g = np.trace(V.T @ M @ V)
        cls = classify_operation(R)
        if cls in chi:
            chi[cls] += chi_g
    H = sum(OH_CLASS_SIZES.values())
    out = {}
    for ir, ch in OH_CHARACTERS.items():
        n_ir = sum(ch[c] * chi[c] for c in OH_CLASSES) / H
        if abs(n_ir) > 0.01:
            out[ir] = round(n_ir)
    return out


def _eigenspace_parity_content(eigvecs, indices, coords, subgroup):
    """For reduced-symmetry clusters, report g/u parity content and degeneracy."""
    V = eigvecs[:, indices]
    # Find inversion in subgroup, if present
    info_parts = []
    info_parts.append(f'dim={len(indices)}')
    for R in subgroup:
        cls = classify_operation(R)
        if cls == 'i':  # inversion
            M = disp_rep(R, coords)
            chi_i = np.trace(V.T @ M @ V)
            # chi_i = (n_g - n_u) for parity decomposition
            n_g = (len(indices) + chi_i) / 2
            n_u = (len(indices) - chi_i) / 2
            info_parts.append(f'parity g:{int(round(n_g))} u:{int(round(n_u))}')
            break
    return '  '.join(info_parts)


def _identify_vector_modes(eigvecs, eigvals, groups, coords, subgroup):
    """Find non-trivial vector modes by projecting test patterns onto the
    eigenspace AFTER orthogonalising against rigid translations.

    Two test patterns:
      (1) Uniform translation: every site displaced by +x_hat, +y_hat, +z_hat
          (the rigid-body translation, at lambda=0).
      (2) Centre-vs-ring antiphase: central node displaced by +x_hat times one
          factor, ring/shell vertices displaced by an opposite factor, weighted
          so the net momentum is zero.

    The second pattern is the natural "vector mode of the cluster":
    polar vector character, parity odd, but orthogonal to translation.
    """
    n = len(coords)
    if n < 2:
        return
    # Find central node (closest to origin)
    norms = np.array([np.linalg.norm(c) for c in coords])
    central = int(np.argmin(norms))
    # Build centre-vs-rest antiphase test patterns along each axis
    test_orth = np.zeros((3 * n, 3))
    n_rest = n - 1
    for a in range(3):
        # +1 on central, -1/n_rest on each other site (mean-zero, vector along axis a)
        for i in range(n):
            weight = 1.0 if i == central else -1.0 / n_rest
            test_orth[3 * i + a, a] = weight
        # Normalise
        norm = np.linalg.norm(test_orth[:, a])
        if norm > 1e-12:
            test_orth[:, a] /= norm
    # Compute overlap with each eigenmode
    overlap_orth = np.abs(eigvecs.T @ test_orth) ** 2
    overlap_per_mode = overlap_orth.sum(axis=1)
    val_overlap = []
    for val, idx in groups:
        if abs(val) < 1e-5:
            continue  # skip translations
        total = overlap_per_mode[idx].sum()
        val_overlap.append((val, len(idx), total))
    val_overlap.sort(key=lambda x: -x[2])
    print(f'    {"lambda":>10s}  {"mult":>4s}  {"centre-vs-rest overlap":>25s}')
    for val, mult, ov in val_overlap[:5]:
        if ov > 1e-3:
            print(f'    {val:8.4f}        {mult:4d}  {ov:>25.4f}')


# ============================================================================
# Driver
# ============================================================================

def main():
    print('=' * 78)
    print('Composite-cluster spectral analysis')
    print('=' * 78)

    oh = generate_Oh()
    print(f'\n|O_h| = {len(oh)}\n')

    for builder in [cluster_crossed_fault, cluster_phi_gmo, cluster_lambda]:
        report_composite(builder, oh)


if __name__ == '__main__':
    main()
