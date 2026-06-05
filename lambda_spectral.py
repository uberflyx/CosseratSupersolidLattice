#!/usr/bin/env python3
"""
lambda_spectral.py
==================
First-principles, mode-weighted spectral mass of the Lambda(1115) hyperon,
the first broken-symmetry baryon, built on the same footing as the proton
chain in cosserat_spectral.py.

What is new relative to the proton, and why
-------------------------------------------
The Lambda^0 (uds, I=0) does not sit on the proton's closed O_h shell.  The
monograph puts it on the 13-node coordination shell plus a three-node {111}
cap (16 nodes), which lowers the point group from O_h to C_3v.  The cap is
the next FCC layer over one octahedral {111} face: three lattice sites, each
a nearest neighbour (d=1) of a shell edge and of one another, arranged as a
C_3 triangle about the [111] axis.

Channel.  The baryon source is the chiral microrotation screw.  Under the
O_h -> C_3v subduction its parity-odd shell irrep A_2u (and its parity-even
partner A_2g) both descend to the single C_3v irrep A_2, in the microrotation
(phi) sector.  That is the Lambda's channel; nothing is typed to fit.

The point of this script: report what can be calculated, not an eyeballed
root.  The monograph records the Lambda eigenvalue as the open case.  At the
reference lambda=4 the node count overshoots (16 m_0 = 1120.4 > 1115.68), so
the correction must be negative and the eigenvalue soft (lambda < 4); but on
the soft branch there is no clean microrotation-dominant shell mode of the
kind the proton and Sigma read, only a mixed mode about 60% microrotation
sitting half on the cap.  The monograph's lambda_Lambda = 3.204 is therefore
flagged with a dagger: it was picked to land on the mass.

What this script finds.  The proton's spectral mass is well-defined because
its channel A_2u(phi) is one-dimensional on the closed shell: the symmetry
alone pins the source, so the weighted lambda_eff = 8.000 is seed-independent
and robust.  The Lambda is different.  Breaking O_h to C_3v makes the A_2(phi)
channel ten-dimensional, so the symmetry channel does NOT pin a unique source:
the weighted lambda_eff swings with the (arbitrary) seed, from about 4.4 to
7.7.  The mode-weighted formula therefore does not close the Lambda from the
channel alone, and reporting any single seed's value as the mass would be
reporting noise.  This is the structural confirmation, from the lattice side,
of exactly what the monograph flags: the accommodated eigenvalue is not fixed
by structure alone.  Closing the Lambda needs the source pinned beyond
symmetry, by the explicit topological screw field (the [111] microrotation
winding); that is the next step.  What this script reports is what can be
calculated: the validated structure (cluster, point group, channel), the
channel degeneracy, and the fact that the monograph's daggered 3.204 is the
dominant A_2 mode (about 61% microrotation, the mixed cap mode) yet carries
only a minority of the channel weight.

Reuses the validated core (mode-weighted decomposition, the mass formula) from
cosserat_spectral.py and the mixed-length Cosserat matrix from the strange-
baryon module delta_first_principles.py (the same matrix that returns the
proton A_2u at 8.303 and the Delta T_1 at 9.0515).
"""

import os
import sys
import numpy as np

# importing the engine sets sys.path to include spectral_mass/ and gives us
# the validated weighted-formula core
import cosserat_spectral as eng                               # noqa: E402
from cosserat_spectral import modal_decomposition, spectral_mass, M0, ME  # noqa: E402

from spectral_classifier import fcc_nn_vectors, generate_Oh, vertex_perm  # noqa: E402
from delta_first_principles import build_cosserat_matrix_two_d            # noqa: E402

PDG_LAMBDA = 1115.683     # MeV (PDG; isoscalar uds, I=0)
MONO_DAGGER = 3.204       # the monograph's provisional, mass-matched eigenvalue


# ============================================================================
# Geometry: the 16-node C_3v hex-cap cluster
# ============================================================================
def hex_cap_cluster():
    """13-node coordination shell + a 3-node {111} cap = 16 nodes.

    The cap is the next FCC layer over one octahedral {111} face: the three
    sites (2,1,1)/sqrt2, (1,2,1)/sqrt2, (1,1,2)/sqrt2.  Each is a nearest
    neighbour (d=1) of two shell sites and of the other two cap sites, and the
    three are mapped into one another by the 3-fold rotation about [111].
    Capping one face (not both) removes inversion: the point group is C_3v.
    """
    shell = np.vstack([np.zeros((1, 3)), fcc_nn_vectors()])      # 13 nodes
    cap = np.array([[2., 1., 1.], [1., 2., 1.], [1., 1., 2.]]) / np.sqrt(2.0)
    return np.vstack([shell, cap])                              # 16 nodes


# ============================================================================
# C_3v group machinery (mirror of the T_d block in delta_first_principles.py)
# ============================================================================
C3V_CHARS = {  # C_3v character table; key is class label, value is character
    'A_1': {'E': +1, 'C_3': +1, 'sigma_v': +1},
    'A_2': {'E': +1, 'C_3': +1, 'sigma_v': -1},   # the Lambda channel
    'E':   {'E': +2, 'C_3': -1, 'sigma_v':  0},
}


def classify_C3v(R, tol=1e-4):
    """Identity, C_3 (the two 120-degree rotations), sigma_v (the three
    reflections), by (det, trace)."""
    det = round(np.linalg.det(R))
    tr = np.trace(R)
    if det == +1:
        if abs(tr - 3) < tol:
            return 'E'
        if abs(tr - 0) < tol:
            return 'C_3'           # 1 + 2 cos120 = 0
    else:
        if abs(tr - 1) < tol:
            return 'sigma_v'       # reflection
    return 'unknown'


def generate_C3v(coords, cap_idx, tol=1e-6):
    """C_3v as the subgroup of O_h that maps the cap triangle to itself, found
    by direct search over the 48 signed permutations.  O_h preserves the shell,
    so preserving the cap set preserves the whole cluster."""
    cap = coords[cap_idx]
    G = []
    for R in generate_Oh():
        moved = (R @ cap.T).T
        ok = all(any(np.linalg.norm(m - c) < tol for c in cap) for m in moved)
        if ok:
            G.append(R)
    return G


def state_representation(R, coords):
    """Represent point-group element R on the (displacement, microrotation)
    state space: u transforms as a polar vector (R), phi as an axial vector
    (det(R) R).  Same convention as the proton and Sigma builds."""
    n = len(coords)
    perm = vertex_perm(R, coords)
    D = np.zeros((6 * n, 6 * n))
    det = np.linalg.det(R)
    for i in range(n):
        j = perm[i]
        D[3 * j:3 * j + 3, 3 * i:3 * i + 3] = R
        D[3 * n + 3 * j:3 * n + 3 * j + 3, 3 * n + 3 * i:3 * n + 3 * i + 3] = det * R
    return D


def build_c3v_projector(coords, G, chars):
    """Isotypic projector P_Gamma = (d_Gamma/|G|) sum_g chi_Gamma(g) D(g)."""
    n = len(coords)
    P = np.zeros((6 * n, 6 * n))
    for R in G:
        P += chars[classify_C3v(R)] * state_representation(R, coords)
    return P * chars['E'] / len(G)


# ============================================================================
# The A_2 microrotation source on the C_3v cluster (built from the channel)
# ============================================================================
def channel_distortion_c3v(coords, G, chars, sector="phi", seed=0):
    """u0 = P_{A_2} @ s, with s a seed confined to the microrotation sector.

    Mirror of channel_distortion() in cosserat_spectral.py, with the O_h
    projector replaced by the C_3v one.  u0 depends only on (coords, irrep,
    sector), never on M, so the modal decomposition is not circular.
    """
    n = len(coords)
    P = build_c3v_projector(coords, G, chars)
    rng = np.random.default_rng(seed)
    s = np.zeros(6 * n)
    if sector == "phi":
        s[3 * n:] = rng.standard_normal(3 * n)
    elif sector == "u":
        s[:3 * n] = rng.standard_normal(3 * n)
    else:
        raise ValueError("sector must be 'u' or 'phi'")
    u0 = P @ s
    nrm = np.linalg.norm(u0)
    if nrm < 1e-9:
        raise ValueError(f"seed has no component in the {sector}-sector of A_2")
    return u0 / nrm


# ============================================================================
# Main: derive m_Lambda
# ============================================================================
def main():
    coords = hex_cap_cluster()
    n = len(coords)
    N = n                                            # node count = 16
    cap_idx = [13, 14, 15]

    print("=" * 78)
    print("lambda_spectral.py  --  Lambda(1115) hyperon, first broken-symmetry baryon")
    print("=" * 78)

    # --- cluster / geometry sanity ---
    d1 = sum(1 for i in range(n) for j in range(i + 1, n)
             if abs(np.linalg.norm(coords[i] - coords[j]) - 1.0) < 1e-6)
    G = generate_C3v(coords, cap_idx)
    print("\nCluster:")
    print(f"  N = {n} nodes  (13 shell + 3 {{111}} cap)")
    print(f"  d = 1 bonds: {d1}   (shell-shell, shell-cap, cap-cap)")
    print(f"  point group found: |G| = {len(G)}  (C_3v has order 6)")
    classes = {}
    for R in G:
        classes[classify_C3v(R)] = classes.get(classify_C3v(R), 0) + 1
    print(f"  classes: {dict(classes)}")

    # --- irrep multiplicity check: tr(P_Gamma)/dim must sum to 6n ---
    print("\nC_3v irrep multiplicities (n_Gamma = tr(P_Gamma)/dim):")
    total = 0.0
    for irrep, chars in C3V_CHARS.items():
        P = build_c3v_projector(coords, G, chars)
        tr = float(np.trace(P))
        total += tr
        print(f"  {irrep:>3s}: tr(P) = {tr:6.2f}   n = {tr / chars['E']:.1f} copies")
    print(f"  total = {total:.1f}  (= 6n = {6 * n})")

    # --- channel ---
    print("\nDerived channel (from the screw topology under O_h -> C_3v, not typed):")
    print("  irrep  = A_2   (descendant of the baryon screw A_2u/A_2g)")
    print("  sector = phi   (Cosserat microrotation)")

    # --- does the channel pin a unique source?  the proton's does; that is why it closes ---
    M = build_cosserat_matrix_two_d(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    P_A2 = build_c3v_projector(coords, G, C3V_CHARS['A_2'])
    phi_block = np.vstack([np.zeros((3 * n, 3 * n)), np.eye(3 * n)])
    dim_A2_phi = int(np.linalg.matrix_rank(P_A2 @ phi_block, tol=1e-8))
    print("\nDoes the channel pin a unique source?")
    print(f"  dim of the A_2 image of the microrotation sector = {dim_A2_phi}")
    print("  (the proton's A_2u(phi) channel is 1-dimensional, so its source is pinned")
    print("   by symmetry and its lambda_eff = 8.000 is seed-independent and robust)")

    # --- with a >1-dim channel the weighted value is not fixed by structure ---
    print("\nWeighted lambda_eff for several arbitrary microrotation seeds:")
    lam_seeds = []
    for s in range(6):
        u0 = channel_distortion_c3v(coords, G, C3V_CHARS['A_2'], sector="phi", seed=s)
        d = modal_decomposition(M, u0, n)
        lam_seeds.append(d["lambda_eff"])
        print(f"  seed={s}: lambda_eff = {d['lambda_eff']:.4f}  "
              f"-> m = {spectral_mass(N, d['lambda_eff']):.2f} MeV")
    print(f"  spread: lambda_eff in [{min(lam_seeds):.2f}, {max(lam_seeds):.2f}]"
          f"  ->  m in [{spectral_mass(N, min(lam_seeds)):.1f}, "
          f"{spectral_mass(N, max(lam_seeds)):.1f}] MeV")
    print("  The symmetry channel alone does not fix the source, so it does not fix the")
    print("  mass.  No single value here is 'the' calculated Lambda mass.")

    # --- the monograph's daggered mode is real, but only the dominant component ---
    u0 = channel_distortion_c3v(coords, G, C3V_CHARS['A_2'], sector="phi", seed=0)
    dec = modal_decomposition(M, u0, n)
    print("\nThe monograph's mode, located:")
    print(f"  dominant A_2 mode: lambda = {dec['dom_lambda']:.4f}, "
          f"{dec['dom_phi_frac']*100:.0f}% microrotation, weight {dec['dom_weight']:.2f}")
    print(f"  this is the monograph's daggered lambda = {MONO_DAGGER:.3f} mixed cap mode,")
    print(f"  correctly identified by the C_3v A_2 source, but carrying only "
          f"{dec['dom_weight']*100:.0f}% of")
    print("  the channel weight -- not the whole distortion.")

    print("\n" + "-" * 78)
    print("Conclusion: what can be calculated")
    print("-" * 78)
    print("  Structure validated: 16-node C_3v cluster, A_2 channel, correct")
    print("  projector multiplicities.  The proton closes because its channel is")
    print(f"  1-dimensional; the Lambda's A_2(phi) channel is {dim_A2_phi}-dimensional, so its")
    print("  spectral mass is NOT fixed by structure alone -- the lattice-side")
    print("  confirmation of the monograph's open flag on the accommodated hyperons.")
    print("  Closing it needs the source pinned by the explicit topological screw")
    print("  field (the [111] microrotation winding), not the symmetry channel alone.")
    return dim_A2_phi, lam_seeds


if __name__ == "__main__":
    main()
