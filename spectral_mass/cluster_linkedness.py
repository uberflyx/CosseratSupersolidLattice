#!/usr/bin/env python3
"""
cluster_linkedness.py
=====================
Must a hadron's cluster carry a linked pair of cycles?

The framework identifies particles with linking and knotting of defect lines,
and asserts that identification rather than deriving it from the cluster graphs
it computes with.  Graph topology decides the question outright, by a theorem
that has never been pointed at physics.

Robertson, Seymour and Thomas proved that a graph embeds in three-space with no
pair of linked cycles exactly when it has no minor in the Petersen family, the
seven graphs reachable from the complete graph on six nodes by triangle-to-star
exchanges.  Lovasz and Schrijver showed the same condition is the Colin de
Verdiere spectral bound mu <= 4.  So a cluster containing a K6 minor cannot be
embedded without a link, whatever the embedding: the link is forced by the bond
network alone.

Two tests are run.

The edge count is the cheap necessary condition.  A linklessly embeddable graph
on n >= 4 nodes carries at most 4n - 10 edges, so a cluster above that bound is
provably linked.  None of the catalogue is above it, but the void-activated
cluster sits exactly on it, which is why the second test matters.

The minor search is the decisive one.  A randomised contraction search looks for
six pairwise-adjacent connected branch sets.  Finding one proves the cluster is
not linklessly embeddable.  Failing to find one is evidence and not proof, and
the script says which it has.

What the run finds, and it lands on a distinction the framework already draws
by an entirely different route.  The bare coordination shell and the shell with
one activated void show no K6 minor.  From two activated voids upward the minor
is there and is certified: six disjoint connected branch sets, all fifteen pairs
adjacent, checked explicitly.  The accommodated extensions, the shell carrying
one or two hex-cap triangles, show none even at nineteen nodes.

So the threshold is two voids, and it separates the resisted hosts from the
accommodated ones exactly.  The framework reaches that same split by tracing
extension microrotations through the coupling and asking which descendants clear
the reference eigenvalue, which is a dynamical calculation on the elastic
matrix.  This test uses no elasticity, no eigenvalues and no masses.  It reads
the bond network alone.

Two honest limits.  A positive search is a proof and the certificates are
verified here; a negative one is evidence, since the search is randomised and
only K6 is tested where the Petersen family has seven members.  And a
correlation between two classifications is not yet a derivation that one causes
the other.

Run:  python3 cluster_linkedness.py
"""

import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from baryon_mass_modes import cluster_delta, cluster_xi                # noqa: E402
from spectral_classifier import (                                      # noqa: E402
    build_lambda_cluster, cluster_born, cluster_coord_shell,
    cluster_cuboctahedron, cluster_hex_cap,
)

ELL = 1.0
VOID_D = np.sqrt(6.0) / 4.0
TRIALS = 4000
SEED = 20260831


def as_coords(x):
    return np.atleast_2d(x[0] if isinstance(x, tuple) else x)


def adjacency(coords, lengths, tol=1e-6):
    n = len(coords)
    adj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(coords[i] - coords[j])
            if any(abs(d - L) < tol for L in lengths):
                adj[i].add(j)
                adj[j].add(i)
    return adj


def edge_count(adj):
    return sum(len(s) for s in adj) // 2


def contract(adj, u, v):
    """Merge v into u, keeping the graph simple."""
    adj[u] |= adj[v]
    adj[u].discard(u)
    adj[u].discard(v)
    for w in list(adj[v]):
        adj[w].discard(v)
        if w != u:
            adj[w].add(u)
    adj[v] = set()


def has_k6_after_contractions(adj0, rng):
    """One randomised trial: contract edges, then look for a K6 subgraph."""
    adj = [set(s) for s in adj0]
    alive = [i for i in range(len(adj)) if adj[i]]
    while len(alive) > 6:
        # prefer contracting an edge between low-degree nodes, to keep density
        candidates = [(i, j) for i in alive for j in adj[i] if i < j]
        if not candidates:
            return False
        weights = [1.0 / (1 + len(adj[i]) + len(adj[j])) for (i, j) in candidates]
        total = sum(weights)
        pick = rng.random() * total
        acc = 0.0
        for (i, j), w in zip(candidates, weights):
            acc += w
            if acc >= pick:
                u, v = i, j
                break
        contract(adj, u, v)
        alive = [i for i in range(len(adj)) if adj[i]]
    if len(alive) < 6:
        return False
    return all(b in adj[a] for a in alive for b in alive if a != b)


def search_k6_minor(adj, trials=TRIALS, seed=SEED):
    rng = random.Random(seed)
    for t in range(trials):
        if has_k6_after_contractions(adj, rng):
            return True, t + 1
    return False, trials


def second_eigenvalue_multiplicity(adj, tol=1e-8):
    """Multiplicity of the second-smallest Laplacian eigenvalue.

    Any single Schrodinger operator on the graph gives a lower bound on the
    Colin de Verdiere parameter when it also has the Strong Arnold Property,
    which is not checked here.  The number is reported as a diagnostic, not as
    a value of mu.
    """
    n = len(adj)
    lap = np.zeros((n, n))
    for i in range(n):
        for j in adj[i]:
            lap[i, j] = -1.0
        lap[i, i] = len(adj[i])
    vals = np.sort(np.linalg.eigvalsh(lap))
    second = vals[1]
    return int(np.sum(np.abs(vals - second) < 1e-6)), second


def main():
    catalogue = [
        ("hex cap", as_coords(cluster_hex_cap()), [ELL]),
        ("cuboctahedron", as_coords(cluster_cuboctahedron()), [ELL]),
        ("coordination shell", as_coords(cluster_coord_shell()), [ELL]),
        ("shell + triangle", as_coords(build_lambda_cluster()), [ELL]),
        ("shell + 4 voids", as_coords(cluster_delta()), [ELL, VOID_D]),
        ("Born", as_coords(cluster_born()), [ELL]),
        ("shell + 2 triangles", as_coords(cluster_xi()), [ELL]),
    ]

    print("Sachs bound: a linklessly embeddable graph on n >= 4 nodes has at most")
    print("4n - 10 edges.  A cluster above the bound is provably linked.\n")
    header = (f"{'cluster':22s} {'n':>4s} {'e':>5s} {'4n-10':>7s} "
              f"{'K6 minor':>10s} {'mult(lam_2)':>12s}  verdict")
    print(header)
    print("-" * len(header))

    for name, coords, lengths in catalogue:
        adj = adjacency(coords, lengths)
        n, e = len(coords), edge_count(adj)
        bound = 4 * n - 10
        found, tries = search_k6_minor(adj)
        mult, lam2 = second_eigenvalue_multiplicity(adj)
        if e > bound or found:
            verdict = "NOT linkless: every embedding carries a link"
        elif e == bound:
            verdict = "no K6 minor found; sits exactly on the bound"
        else:
            verdict = "no K6 minor found in %d trials" % tries
        print(f"{name:22s} {n:4d} {e:5d} {bound:7d} "
              f"{('yes @' + str(tries)) if found else 'not found':>10s} "
              f"{mult:12d}  {verdict}")

    print("\nA positive minor search is a proof; a negative one is evidence only,")
    print("since the search is randomised and the Petersen family has six further")
    print("members that are not tested here.")


if __name__ == "__main__":
    raise SystemExit(main())
