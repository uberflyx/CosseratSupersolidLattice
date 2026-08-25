"""Proper three-colourings of the FCC coordination shell, and why there are four.

The neutrino chapter needs the number of {111} slip-plane families to be a
derived quantity rather than an observation about crystals.  It is: the four
families are the four structurally distinct ways of properly three-colouring
the cuboctahedral coordination graph, the graph whose vertices are a node's
twelve nearest neighbours in FCC and whose edges join neighbours that are
themselves nearest neighbours.

A proper colouring assigns one of three colours to each vertex so that no edge
joins two vertices of the same colour.  Two colourings that differ only by
permuting the three colour labels are the same structure, so the count that
matters is P(G, 3) divided by 3! = 6.

Vertices: the twelve permutations of (+-1, +-1, 0), which are the midpoints of
a cube's edges.  Two are adjacent when they are a nearest-neighbour pair,
separated by sqrt(2) in these units.
"""

import itertools

import numpy as np


def coordination_graph():
    """Vertices and edges of the FCC coordination shell (cuboctahedron)."""
    verts = [v for v in itertools.product([-1, 0, 1], repeat=3)
             if sorted(map(abs, v)) == [0, 1, 1]]
    edges = [(i, j)
             for i in range(len(verts)) for j in range(i + 1, len(verts))
             if np.isclose(np.linalg.norm(np.subtract(verts[i], verts[j])),
                           np.sqrt(2))]
    return verts, edges


def proper_colourings(n_vertices, edges, n_colours=3):
    """Every proper colouring, by direct enumeration.

    3^12 = 531441 assignments, so brute force is honest and instant; no
    chromatic-polynomial machinery is needed and none is trusted.
    """
    return [c for c in itertools.product(range(n_colours), repeat=n_vertices)
            if all(c[i] != c[j] for i, j in edges)]


def canonical(colouring):
    """Colouring with the colour labels renamed in order of first appearance.

    Two colourings sharing a canonical form differ only by a relabelling and
    describe the same structure.
    """
    mapping, nxt = {}, 0
    out = []
    for c in colouring:
        if c not in mapping:
            mapping[c], nxt = nxt, nxt + 1
        out.append(mapping[c])
    return tuple(out)


def main():
    verts, edges = coordination_graph()
    degree = 2 * len(edges) / len(verts)
    print(f"coordination graph: {len(verts)} vertices, {len(edges)} edges, "
          f"degree {degree:.0f}")

    colourings = proper_colourings(len(verts), edges)
    distinct = {canonical(c) for c in colourings}
    print(f"  P(G, 3)            = {len(colourings)}")
    print(f"  P(G, 3) / 3!       = {len(colourings) // 6}")
    print(f"  distinct up to relabelling = {len(distinct)}")
    assert len(colourings) == 24 and len(distinct) == 4

    # Each colouring splits the twelve neighbours into three groups of four,
    # and the group a vertex lands in names the {111} family it belongs to.
    print("\n  the four colourings, as the partition each induces:")
    for k, c in enumerate(sorted(distinct), start=1):
        groups = [[verts[i] for i in range(len(verts)) if c[i] == col]
                  for col in range(3)]
        sizes = [len(g) for g in groups]
        print(f"    {k}: group sizes {sizes}")


if __name__ == "__main__":
    main()
