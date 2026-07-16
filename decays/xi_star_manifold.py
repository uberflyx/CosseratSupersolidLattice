#!/usr/bin/env python3
"""
xi_star_manifold.py
====================

Diagonalises the Xi*(1530) void-deactivation graph and keeps the
EIGENVECTORS of its two soft modes, not just the eigenvalues.

decays/cosserat_decay_engine.py already builds this graph and uses its
two lowest nonzero eigenvalues (lambda_2 = 1.076, lambda_3 = 1.839) in a
product formula for the decay width Xi* -> Xi + pi (Sec. decuplet_decays
of the dynamics chapter).  That formula misses the observed width by
about a quarter (-24% with the Clebsch-Gordan-weighted channel sum and
the corrected shadow factor sigma = 0.0087), and the chapter names the
reason: at |S|=2 there are two soft modes and a coupling built from
lambda_2 alone cannot see the second (Sec. fiedler_manifold).

This script answers a question the eigenvalues alone cannot: what do the
two soft eigenVECTORS actually look like, and does the leading formula
use a channel that a symmetry of the cluster forbids outright?

RESULT
------
The two active strange caps sit on {111} planes related by an exact
180-degree rotation of the cuboctahedral shell (the rotation about the
axis bisecting the two planes' normals swaps both the two caps and the
two remaining voids; verified below by direct application to the four
{111} plane normals).  Because the deactivation Laplacian commutes with
this exchange, every non-degenerate eigenvector must be purely symmetric
or purely antisymmetric under it -- the same statement that fixes the
normal modes of two identical pendulums joined by a weak spring to be
only the in-phase and out-of-phase swings.

Diagonalising the 21-node graph confirms it directly: lambda_2 (1.076)
is the ODD (antisymmetric) combination and lambda_3 (1.839) is the EVEN
(symmetric) one.  Neither sits on one cap alone, which is itself a
correction to the chapter's earlier reading.

Because pion emission is the same physical event at either cap (the
graph-democracy argument of Sec. weak_decays, the same one that fixes
g_piNN = 13), the source that represents it is symmetric under the cap
exchange by construction.  A symmetric source has exactly zero overlap
with an antisymmetric mode.  So the lambda_2 channel that both the
Sigma*-style single-ratio guess and the existing product formula include
is not merely small: it is forbidden by an exact selection rule.  Only
lambda_3 can carry the decay.

WHAT THIS DOES NOT DO
----------------------
It does not produce a corrected MeV width.  lambda_3 alone still needs
the correct partner on the Delta (|S|=0) side: at |S|=0 all four voids
are equivalent and the reference eigenvalue (2.4384) is threefold
degenerate, and that degenerate space needs the SAME symmetry
decomposition before lambda_3 can be read as a ratio against it.  That
decomposition is not attempted here.
"""

import sys
import os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)
from cosserat_graph_legacy import FCCLattice, Defect


def verify_cap_exchange_is_a_symmetry():
    """Confirm the 180-degree rotation that swaps {111} planes 0 and 1
    also swaps planes 2 and 3, so the deactivation graph as a whole (both
    active caps and both remaining voids) is exchanged onto itself, not
    just relabelled at the two caps."""
    normals = [np.array([1., 1, 1]), np.array([1., -1, -1]),
               np.array([-1., 1, -1]), np.array([-1., -1, 1])]
    axis = normals[0] + normals[1]
    axis = axis / np.linalg.norm(axis)

    def rotate_pi(v):
        # Rodrigues' formula specialised to a 180-degree rotation.
        return 2 * np.dot(v, axis) * axis - v

    mapping = {}
    for i, n in enumerate(normals):
        r = rotate_pi(n)
        j = next(k for k, m in enumerate(normals) if np.allclose(r, m, atol=1e-8))
        mapping[i] = j
    assert mapping == {0: 1, 1: 0, 2: 3, 3: 2}, mapping
    return mapping


def build_deactivation_graph(n_strange):
    """Rebuild the deactivation graph exactly as
    cosserat_decay_engine.py::_decuplet_deact_coupling does (correct T_d
    void positions, void-swap per strange step, void-void bonds excluded),
    but return the node list and role map too, so eigenvectors can be
    read back onto physical nodes."""
    lat = FCCLattice()
    d = Defect(lat)
    d.add_shell()
    for vi in range(4 - n_strange):
        nm = f'v{vi}'
        d.nodes.add(nm); d.roles[nm] = 'void'; d.pos[nm] = lat.void_positions[vi].copy()
    d.has_voids = (4 - n_strange) > 0
    for pi in range(n_strange):
        d.add_strange_ext(pi)

    nodes = sorted(d.nodes)
    idx = {n: i for i, n in enumerate(nodes)}
    _, _, L_full = d.graph_matrices()
    A = -np.array(L_full)
    np.fill_diagonal(A, 0)
    voids = [idx[n] for n in nodes if d.roles.get(n) == 'void']
    for i in voids:
        for j in voids:
            if i != j:
                A[i, j] = 0.0
    D = np.diag(A.sum(axis=1))
    L_deact = D - A
    return nodes, idx, d.roles, L_deact


def cap_projection(vec, idx, cap_nodes):
    """Overlap of an eigenvector with a unit-normalised indicator vector
    on one strange cap's extension nodes."""
    ii = [idx[n] for n in cap_nodes]
    unit = np.zeros(len(vec))
    unit[ii] = 1.0 / np.sqrt(len(ii))
    return float(np.dot(unit, vec))


def main():
    mapping = verify_cap_exchange_is_a_symmetry()
    print(f"cap-exchange rotation: plane 0<->{mapping[0]}, plane 2<->{mapping[2]}  (confirms an exact O_h symmetry)\n")

    nodes, idx, roles, L = build_deactivation_graph(n_strange=2)
    evals, evecs = np.linalg.eigh(L)
    print(f"Xi* (|S|=2) deactivation graph: {len(nodes)} nodes")
    print(f"lowest eigenvalues: {np.round(evals[:4], 4)}  (chapter: lambda_2=1.076, lambda_3=1.839)\n")

    lam2, v2 = evals[1], evecs[:, 1]
    lam3, v3 = evals[2], evecs[:, 2]
    cap0 = [n for n in nodes if roles.get(n) == 'extension' and n.startswith('e0_')]
    cap1 = [n for n in nodes if roles.get(n) == 'extension' and n.startswith('e1_')]

    for label, lam, v in [('lambda_2', lam2, v2), ('lambda_3', lam3, v3)]:
        c0 = cap_projection(v, idx, cap0)
        c1 = cap_projection(v, idx, cap1)
        parity = 'odd (antisymmetric)' if abs(c0 + c1) < 1e-6 else 'even (symmetric)'
        print(f"{label} = {lam:.4f}:  c0 = {c0:+.4f}, c1 = {c1:+.4f}  -> {parity}")

    print("\nDemocratic (symmetric) pion-emission source (c0, c1) = (+1, +1)/sqrt(2):")
    for label, v in [('lambda_2 (odd)', v2), ('lambda_3 (even)', v3)]:
        c0 = cap_projection(v, idx, cap0)
        c1 = cap_projection(v, idx, cap1)
        overlap = (c0 + c1) / np.sqrt(2)
        print(f"  overlap with {label} mode: {overlap:+.4f}  "
              f"{'(exactly forbidden)' if abs(overlap) < 1e-6 else ''}")


if __name__ == '__main__':
    main()


def delta_reference_decomposition():
    """The Delta-side decomposition the monograph's Fiedler-manifold section
    names as the missing step. The cap-exchange rotation splits the
    threefold-degenerate lambda_2 = 2.4384 space into 1 even + 2 odd; the
    even (symmetric) piece keeps the degenerate eigenvalue, so the reference
    does not move. The even-only reading of the Xi* product then overshoots
    (suppression 0.754 against the 0.441 the data requires, with the LO
    product at 0.333): the closure needs the vertex matrix element, not a
    further eigenvalue combination."""
    import numpy as np
    nodes, idx, roles, L = build_deactivation_graph(0)
    w, v = np.linalg.eigh(L)
    trip = [k for k in range(len(w)) if abs(w[k] - 2.4384) < 1e-3]
    lat = FCCLattice(); d = Defect(lat); d.add_shell()
    for vi in range(4):
        nm = f'v{vi}'; d.nodes.add(nm); d.roles[nm] = 'void'
        d.pos[nm] = lat.void_positions[vi].copy()
    n0 = lat.plane_normals[0]/np.linalg.norm(lat.plane_normals[0])
    n1 = lat.plane_normals[1]/np.linalg.norm(lat.plane_normals[1])
    ax = (n0 + n1); ax = ax/np.linalg.norm(ax)
    P = np.zeros((len(nodes), len(nodes)))
    for n in nodes:
        r = 2*np.dot(d.pos[n], ax)*ax - d.pos[n]
        m = next(k for k in nodes if np.allclose(d.pos[k], r, atol=1e-8))
        P[idx[m], idx[n]] = 1.0
    assert np.allclose(P @ L, L @ P)
    sub = v[:, trip]
    M = sub.T @ P @ sub
    ev = np.linalg.eigvalsh((M + M.T)/2)
    return ev  # parities of the triplet branches: [-1, -1, +1]
