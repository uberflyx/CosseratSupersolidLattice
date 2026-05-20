#!/usr/bin/env python3
"""
hadron_spectral_mass.py
=======================

Algorithmic mass closure for hadrons and light nuclei via the spectral
expectation rule.

The rule
--------
For a hadron with cluster C and reference defect mode psi_ref defined
on a sub-cluster S subset of C:

    lambda_eff = <psi_ref^embed | H_C | psi_ref^embed>
               = sum_k lambda_k * |<psi_ref^embed | psi_k>|^2

    m = N * m_0 - N * (4 - lambda_eff) * m_e

where psi_ref^embed is psi_ref zero-extended to the full cluster's
degree-of-freedom space.  No picking, no thresholds, no judgement: one
Rayleigh quotient.

Reference modes
---------------
- Screw dislocation (baryons, nucleons): psi_p = A_(2u)^+ at lambda = 8.303
  on the bare 13-node cuboctahedral shell.

Constituent combinations for compounds
--------------------------------------
For B identical constituents at positions {r_c}, the natural reference
is a parity-matched combination:

    psi_ref = sum_c eps_c * embed(psi_p, r_c) / sqrt(sum eps_c^2)

where eps_c in {+1, -1} encodes the constituent-label permutation
representation.  Parity matching: the hadron's observed parity is the
product of the local A_(2u) parity (parity-odd) and the label permutation
parity (+/-).

For totally-symmetric label combinations (eps_c = +1 for all c), the
state has parity (-1)^B in the A_(2u) tower (since each constituent
contributes -1 from local parity, and the label part is even).

Practical: for each compound, try both totally-symmetric and 
parity-flipped label combinations; report both.
"""

import numpy as np
from scipy.linalg import eigh
import sys

from spectral_classifier import fcc_nn_vectors
from cosserat_classifier import build_cosserat_matrix
from proton_first_principles import (build_oh_elements, CHARS_A2U,
                                       build_irrep_projector)

# =========================================================================
# Constants
# =========================================================================

M_E = 0.51099895   # MeV, CODATA 2022
ALPHA = 1.0 / 137.0359992
M_0 = M_E / ALPHA  # 70.0253 MeV
ELL = 1.0          # FCC NN distance in working units; r_e in SI

PDG = {
    'proton':   938.272,
    'Delta':    1232.0,
    'deuteron': 1875.613,
    'triton':   2808.921,
    'helium4':  3727.379,
    'Sigma*':   1384.4,    # mean of Sigma*+ Sigma*0 Sigma*-
    'Xi*':      1533.4,
    'Omega-':   1672.45,
}

# =========================================================================
# Reference: the proton mass mode psi_p (A_(2u)^+ on bare 13-shell)
# =========================================================================

def get_psi_p():
    """Returns (lambda_p, psi_p) where psi_p is on the bare 13-shell.

    Layout: psi_p has 78 components = [u_0,...,u_12, phi_0,...,phi_12]
    with u_k and phi_k each being 3-vectors.  Atom 0 is the center; atoms
    1..12 are the FCC nearest neighbours in canonical order from
    fcc_nn_vectors().
    """
    centre = np.array([[0., 0., 0.]])
    shell  = fcc_nn_vectors()
    coords = np.vstack([centre, shell])
    M = build_cosserat_matrix(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    vals, vecs = eigh(M)
    oh = build_oh_elements()
    P_A2u = build_irrep_projector(coords, oh, CHARS_A2U)

    a2u = [(vals[k], vecs[:, k].copy())
           for k in range(len(vals))
           if vecs[:, k] @ (P_A2u @ vecs[:, k]) > 0.5]
    a2u.sort(key=lambda t: -t[0])
    return a2u[0]


# =========================================================================
# Cluster builders: FCC coordinates for each hadron
# =========================================================================

def cluster_proton_coords():
    """13-atom cuboctahedral shell."""
    centre = np.array([[0., 0., 0.]])
    shell  = fcc_nn_vectors()
    return np.vstack([centre, shell])


def cluster_delta_coords():
    """17-atom void-activated cluster.

    Inner shell + 4 T_d^+ void atoms at body-diagonal positions.
    The void atoms are at the centroids of the 4 tetrahedral voids on
    the (1,1,1)-type body diagonals.
    """
    shell = cluster_proton_coords()  # 13 atoms
    # T_d^+ void positions: face-centered at (1,1,1)/4 etc in conventional
    # cubic cell; in our FCC NN units the void corners are at
    # +/- (1,1,1)/(2 sqrt(2)) along the 4 body diagonals (the proper T_d^+
    # subset has all-positive or all-negative signs).
    v = 1.0 / (2.0 * np.sqrt(2.0))
    voids = np.array([
        [+v, +v, +v],
        [+v, -v, -v],
        [-v, +v, -v],
        [-v, -v, +v],
    ])
    return np.vstack([shell, voids])


def cluster_deuteron_coords():
    """20-atom Mode P deuteron compound.

    Two coordination shells centred at A and B, sharing 4 atoms on the
    rectangular face perpendicular to AB.  For FCC NN distance 1 between
    centres, the constituents at A=(0,0,0) and B=(1,1,0)/sqrt(2) share
    4 atoms.
    """
    # Constituent A: 13 atoms (center + 12 NN)
    A = np.array([0., 0., 0.])
    shell_A = fcc_nn_vectors() + A
    cluster_A = np.vstack([[A], shell_A])  # 13 atoms

    # Constituent B: 13 atoms shifted by one NN vector
    B = fcc_nn_vectors()[0]  # take first NN as B's centre
    shell_B = fcc_nn_vectors() + B
    cluster_B = np.vstack([[B], shell_B])  # 13 atoms

    # Union with deduplication
    return _dedup_atoms(np.vstack([cluster_A, cluster_B]))


def cluster_tritium_coords():
    """25-atom equilateral-triangle compound.

    Three shells at vertices of an FCC NN equilateral triangle.
    The three centres are at unit NN distance from each other.
    """
    nn = fcc_nn_vectors()
    centres = np.array([
        [0., 0., 0.],
        nn[0],   # (1,1,0)/sqrt(2)
        nn[1],   # (1,0,1)/sqrt(2)  -- forms equilateral triangle with origin
    ])
    # Verify it's equilateral
    d12 = np.linalg.norm(centres[1] - centres[0])
    d13 = np.linalg.norm(centres[2] - centres[0])
    d23 = np.linalg.norm(centres[2] - centres[1])
    assert abs(d12 - 1) < 1e-8 and abs(d13 - 1) < 1e-8 and abs(d23 - 1) < 1e-8, \
        f"Triangle not equilateral: {d12}, {d13}, {d23}"

    atoms = []
    for c in centres:
        atoms.append(c)
        for nn_v in nn:
            atoms.append(c + nn_v)
    return _dedup_atoms(np.vstack(atoms))


def cluster_helium4_coords():
    """28-atom T_d-void tetrahedron compound.

    Four shells at vertices of a tetrahedral void, with the void centroid
    NOT being an FCC site (so no quadruple-shared atom).
    """
    nn = fcc_nn_vectors()
    centres = np.array([
        [0., 0., 0.],
        nn[0],   # (1,1,0)/sqrt(2)
        nn[1],   # (1,0,1)/sqrt(2)
        nn[2],   # (0,1,1)/sqrt(2)
    ])
    # Verify it's a regular tetrahedron
    ds = []
    for i in range(4):
        for j in range(i+1, 4):
            ds.append(np.linalg.norm(centres[i] - centres[j]))
    assert all(abs(d - 1) < 1e-8 for d in ds), f"Not a regular tetrahedron: {ds}"

    atoms = []
    for c in centres:
        atoms.append(c)
        for nn_v in nn:
            atoms.append(c + nn_v)
    return _dedup_atoms(np.vstack(atoms))


def _dedup_atoms(coords, tol=1e-8):
    """Return atoms with duplicates removed (within tolerance)."""
    keep = []
    for a in coords:
        is_dup = False
        for b in keep:
            if np.linalg.norm(a - b) < tol:
                is_dup = True
                break
        if not is_dup:
            keep.append(a)
    return np.array(keep)


# =========================================================================
# Cosserat matrix builders (handle uniform NN and mixed-bond cases)
# =========================================================================

def build_cosserat_uniform(coords, ell=1.0, tol=1e-6):
    """NN bonds only at distance ell.  K_u = K_phi = 1, alpha = 1."""
    return build_cosserat_matrix(coords, K_u=1.0, K_phi=1.0, alpha=1.0)


def build_cosserat_with_voids(coords, n_shell, ell=1.0):
    """Mixed bonds: NN at d=ell, plus void-shell bonds at d=sqrt(6)/4 ell.

    n_shell: number of shell atoms (first n_shell atoms are the inner shell).
    The rest are void atoms.  Bonds:
      - NN bonds at d=ell among shell atoms
      - Void-shell bonds at d=sqrt(6)/4 ell between each void atom and its
        3 nearest shell atoms
    """
    # Use the existing delta_first_principles builder which handles this
    from delta_first_principles import build_cosserat_matrix_two_d
    return build_cosserat_matrix_two_d(coords, K_u=1.0, K_phi=1.0, alpha=1.0)


# =========================================================================
# Reference embedding
# =========================================================================

def find_constituent_shell_atoms(coords, centre, tol=1e-6):
    """Return atom indices in coords that belong to the shell at centre.

    Inner shell = centre + 12 NN, so 13 atoms.  Returns indices in the
    order [centre, NN_0, NN_1, ..., NN_11] where NN_k follows the
    canonical order of fcc_nn_vectors().
    """
    nn = fcc_nn_vectors()
    shell_positions = [centre] + [centre + v for v in nn]

    indices = []
    for sp in shell_positions:
        # Find the index in coords matching this position
        best_idx = None
        best_d = np.inf
        for i, c in enumerate(coords):
            d = np.linalg.norm(c - sp)
            if d < best_d:
                best_d = d
                best_idx = i
        if best_d > tol:
            raise ValueError(f"Shell atom at {sp} not found in coords (closest {best_d})")
        indices.append(best_idx)
    return indices


def embed_psi_p_at(coords, centre, psi_p, sign=+1):
    """Build a 6N-dim vector with psi_p placed at the shell around centre.

    psi_p layout: [u_0..u_12, phi_0..phi_12].
    Output layout: [u_0..u_{N-1}, phi_0..phi_{N-1}] with N=len(coords).
    Non-shell atoms have zero u and phi.
    """
    N = len(coords)
    shell_idx = find_constituent_shell_atoms(coords, centre)

    out = np.zeros(6 * N)
    # Place u components (psi_p[0:39] are u_0..u_12, 3 per atom)
    for k, idx in enumerate(shell_idx):
        out[3*idx:3*(idx+1)] = sign * psi_p[3*k:3*(k+1)]
    # Place phi components (psi_p[39:78] are phi_0..phi_12)
    for k, idx in enumerate(shell_idx):
        out[3*N + 3*idx : 3*N + 3*(idx+1)] = sign * psi_p[39 + 3*k : 39 + 3*(k+1)]

    return out


def build_compound_reference(coords, constituent_centres, signs, psi_p):
    """Sum the per-constituent embeddings, then normalise.

    constituent_centres: list of 3-vectors (each constituent's centre)
    signs: list of +/-1 for each constituent
    """
    total = np.zeros(6 * len(coords))
    for c, s in zip(constituent_centres, signs):
        total = total + embed_psi_p_at(coords, c, psi_p, sign=s)
    nrm = np.linalg.norm(total)
    if nrm < 1e-12:
        raise ValueError(f"Combination has zero norm; check signs {signs}")
    return total / nrm


# =========================================================================
# Spectral expectation kernel
# =========================================================================

def spectral_expectation(H, psi_ref):
    """Rayleigh quotient <psi|H|psi>."""
    nrm = psi_ref @ psi_ref
    if nrm < 1e-12:
        return float('nan')
    return (psi_ref @ H @ psi_ref) / nrm


def mass_from_lambda(N, lam):
    return N * M_0 - N * (4 - lam) * M_E


def report_spectral_decomposition(H, psi_ref, top=5):
    """Print the top-weight eigenstates in the decomposition of psi_ref."""
    eigvals, eigvecs = eigh(H)
    overlaps = eigvecs.T @ psi_ref
    weights = overlaps**2
    order = np.argsort(-weights)
    print(f"    Top {top} eigenstates contributing to reference:")
    print(f"    {'lambda':>9}  {'weight':>9}  {'lambda*weight':>15}")
    for i in range(min(top, len(order))):
        k = order[i]
        if weights[k] > 1e-6:
            print(f"    {eigvals[k]:9.4f}  {weights[k]:9.4f}  "
                  f"{eigvals[k]*weights[k]:15.4f}")


# =========================================================================
# Hadron specifications
# =========================================================================

def hadron_spec_proton():
    coords = cluster_proton_coords()
    H = build_cosserat_uniform(coords)
    return {'name': 'proton', 'N': 13, 'coords': coords, 'H': H,
            'constituents': [np.zeros(3)], 'signs_options': [[+1]]}


def hadron_spec_delta():
    coords = cluster_delta_coords()
    H = build_cosserat_with_voids(coords, n_shell=13)
    return {'name': 'Delta', 'N': 17, 'coords': coords, 'H': H,
            'constituents': [np.zeros(3)], 'signs_options': [[+1]]}


def hadron_spec_deuteron():
    coords = cluster_deuteron_coords()
    H = build_cosserat_uniform(coords)
    centres = [np.zeros(3), fcc_nn_vectors()[0]]
    return {'name': 'deuteron', 'N': len(coords), 'coords': coords, 'H': H,
            'constituents': centres,
            # parity-even (g) = antisymmetric label
            # parity-odd  (u) = symmetric label
            'signs_options': [[+1, -1], [+1, +1]],
            'signs_labels': ['parity-even (g)', 'parity-odd (u)']}


def hadron_spec_tritium():
    coords = cluster_tritium_coords()
    H = build_cosserat_uniform(coords)
    nn = fcc_nn_vectors()
    centres = [np.zeros(3), nn[0], nn[1]]
    return {'name': 'triton', 'N': len(coords), 'coords': coords, 'H': H,
            'constituents': centres,
            # totally symmetric A_1 of C_3v
            'signs_options': [[+1, +1, +1]],
            'signs_labels': ['totally symmetric A_1']}


def hadron_spec_helium4():
    coords = cluster_helium4_coords()
    H = build_cosserat_uniform(coords)
    nn = fcc_nn_vectors()
    centres = [np.zeros(3), nn[0], nn[1], nn[2]]
    return {'name': 'helium4', 'N': len(coords), 'coords': coords, 'H': H,
            'constituents': centres,
            'signs_options': [[+1, +1, +1, +1]],
            'signs_labels': ['totally symmetric A_1']}


# =========================================================================
# Driver
# =========================================================================

def run_hadron(spec, psi_p, verbose=True):
    """Run the spectral expectation rule on one hadron."""
    name = spec['name']
    N = spec['N']
    coords = spec['coords']
    H = spec['H']

    results = []
    for opt_idx, signs in enumerate(spec['signs_options']):
        label = spec.get('signs_labels', [f'option {opt_idx+1}'])[opt_idx] \
            if spec.get('signs_labels') else f'option {opt_idx+1}'

        psi_ref = build_compound_reference(coords, spec['constituents'], signs, psi_p)
        lam_eff = spectral_expectation(H, psi_ref)
        m_pred  = mass_from_lambda(N, lam_eff)
        pdg = PDG.get(name, np.nan)
        residual_mev = m_pred - pdg
        residual_pct = 100.0 * residual_mev / pdg if not np.isnan(pdg) else np.nan

        results.append({
            'label': label, 'lambda_eff': lam_eff, 'm_pred': m_pred,
            'residual_mev': residual_mev, 'residual_pct': residual_pct,
        })

    if verbose:
        print(f"\n=== {name}  (N = {N} atoms,  PDG = {PDG.get(name, '?')} MeV) ===")
        for r in results:
            print(f"  [{r['label']:>20}]  lambda_eff = {r['lambda_eff']:7.4f}  "
                  f"m = {r['m_pred']:8.2f} MeV  "
                  f"residual = {r['residual_mev']:+7.2f} MeV ({r['residual_pct']:+.3f}%)")

        # Spectral decomposition for the first option
        signs = spec['signs_options'][0]
        psi_ref = build_compound_reference(coords, spec['constituents'], signs, psi_p)
        report_spectral_decomposition(H, psi_ref, top=5)

    return results


def main():
    print("=" * 76)
    print("SPECTRAL EXPECTATION RULE: full-catalogue audit")
    print("=" * 76)
    lam_p, psi_p = get_psi_p()
    print(f"\nReference (bare 13-shell):  psi_p at lambda = {lam_p:.4f}")

    specs = [
        hadron_spec_proton(),
        hadron_spec_delta(),
        hadron_spec_deuteron(),
        hadron_spec_tritium(),
        hadron_spec_helium4(),
    ]

    all_results = {}
    for spec in specs:
        all_results[spec['name']] = run_hadron(spec, psi_p)

    # Summary table
    print()
    print("=" * 76)
    print("SUMMARY")
    print("=" * 76)
    print(f"{'hadron':>10}  {'N':>3}  {'lambda_eff':>10}  {'m_pred':>9}  "
          f"{'PDG':>9}  {'residual':>9}  {'% PDG':>8}")
    for name, results in all_results.items():
        # Pick the first option as the default
        r = results[0]
        N = next(s['N'] for s in specs if s['name'] == name)
        pdg = PDG.get(name, np.nan)
        print(f"{name:>10}  {N:>3}  {r['lambda_eff']:>10.4f}  "
              f"{r['m_pred']:>9.2f}  {pdg:>9.2f}  "
              f"{r['residual_mev']:>+9.2f}  {r['residual_pct']:>+7.3f}%")


if __name__ == '__main__':
    main()
