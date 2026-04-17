#!/usr/bin/env python3
"""
fcc_defect_catalogue.py -- Composition-based catalogue of FCC defect cores.

Instead of enumerating every FCC-realisable subgraph up to N vertices (which
explodes combinatorially), this script builds the catalogue the monograph
actually uses: explicit compositions of a handful of canonical building
blocks, deduplicated under D_3d / O_h symmetry, with multiplicities.

Building blocks of the monograph's hadron catalogue:
    cell_pair     (N=2)    unoriented NN bond
    hex_cap_i     (N=7)    origin + 6 shell sites on {111} plane i (i=0..3)
    coord_shell   (N=13)   origin + all 12 shell sites
    bilayer_i     (N=8)    hex_cap_i + 1 NN on the next {111} layer
    tet_voids     (N=4)    4 tetrahedral voids around the origin

Overlays (do not add vertices, but change Q):
    winding       baryon-number label on a coord_shell
    strange_arm_i add a hex_cap-style extension on plane i
    radial n_r    append n_r bilayer quanta

This script enumerates all realisable compositions up to a total-vertex
bound, canonicalises each under the chosen point group, and records the
multiplicity of equivalent (plane-permutation-related) compositions.

The output catalogue is at most a few hundred entries — tractable by
inspection, usable directly as the invariant selector in cosserat_graph.py.

Usage:
    python fcc_defect_catalogue.py --nmax 25 --symmetry oh --out catalogue.json

Requires: networkx, numpy, pynauty (same as fcc_enumerate.py).
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from hashlib import blake2b
from pathlib import Path

import networkx as nx
import numpy as np
import pynauty


# ============================================================================
# FCC GEOMETRY (shared with fcc_enumerate_v2)
# ============================================================================

NN_OFFSETS: tuple[tuple[int, int, int], ...] = tuple(
    (s1, s2, 0) for s1, s2 in itertools.product((-1, 1), repeat=2)
) + tuple(
    (s1, 0, s2) for s1, s2 in itertools.product((-1, 1), repeat=2)
) + tuple(
    (0, s1, s2) for s1, s2 in itertools.product((-1, 1), repeat=2)
)
assert len(NN_OFFSETS) == 12

N111: tuple[tuple[int, int, int], ...] = (
    (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1),
)


def edge_character(off: tuple[int, int, int]) -> tuple[str, str]:
    dx, dy, dz = off
    if dz == 0:
        return ('C', '+' if dx * dy > 0 else '-')
    if dy == 0:
        return ('B', '+' if dx * dz > 0 else '-')
    return ('A', '+' if dy * dz > 0 else '-')


def _oh_matrices():
    return [
        tuple(
            tuple(signs[i] if perm[i] == j else 0 for j in range(3))
            for i in range(3)
        )
        for perm in itertools.permutations(range(3))
        for signs in itertools.product((1, -1), repeat=3)
    ]


OH_GROUP = _oh_matrices()


def _d3d_matrices(axis=(1, 1, 1)):
    stab = []
    for R in OH_GROUP:
        img = (
            R[0][0] * axis[0] + R[0][1] * axis[1] + R[0][2] * axis[2],
            R[1][0] * axis[0] + R[1][1] * axis[1] + R[1][2] * axis[2],
            R[2][0] * axis[0] + R[2][1] * axis[1] + R[2][2] * axis[2],
        )
        if img == axis or img == (-axis[0], -axis[1], -axis[2]):
            stab.append(R)
    return stab


D3D_GROUP = _d3d_matrices()


def apply_rotation(R, s):
    return (
        R[0][0] * s[0] + R[0][1] * s[1] + R[0][2] * s[2],
        R[1][0] * s[0] + R[1][1] * s[1] + R[1][2] * s[2],
        R[2][0] * s[0] + R[2][1] * s[1] + R[2][2] * s[2],
    )


def canonical_under_group(sites, group):
    best = None
    for R in group:
        img = tuple(sorted(apply_rotation(R, s) for s in sites))
        if best is None or img < best:
            best = img
    return best


# ============================================================================
# BUILDING BLOCKS — the 5 canonical defect constituents
# ============================================================================

def _shell_sites_in_plane(plane_idx: int) -> list:
    n = N111[plane_idx]
    return [s for s in NN_OFFSETS
            if s[0] * n[0] + s[1] * n[1] + s[2] * n[2] == 0]


def _shell_sites_above_plane(plane_idx: int) -> list:
    n = N111[plane_idx]
    return [s for s in NN_OFFSETS
            if s[0] * n[0] + s[1] * n[1] + s[2] * n[2] == 2]


def _tet_void_sites() -> list:
    """Four tetrahedral void positions in an FCC primitive cell.

    In the monograph's integer representation these are fractional
    coordinates, but for the graph we just treat them as distinct lattice
    labels — the combinatorics only cares about connectivity to shell sites.
    """
    return [(1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1)]


BUILDING_BLOCKS: dict = {
    # name : sites (origin always implicit as (0,0,0))
    'cell_pair':  {'sites': [(0, 0, 0), NN_OFFSETS[0]],
                   'plane_dependent': False, 'N': 2},
    'coord_shell': {'sites': [(0, 0, 0)] + list(NN_OFFSETS),
                    'plane_dependent': False, 'N': 13},
    'tet_voids':  {'sites': _tet_void_sites(),
                   'plane_dependent': False, 'N': 4,
                   'requires': 'coord_shell'},
    # plane-dependent blocks, one variant per {111} plane:
    **{f'hex_cap_{i}': {'sites': [(0, 0, 0)] + _shell_sites_in_plane(i),
                        'plane_dependent': True, 'plane_idx': i, 'N': 7}
       for i in range(4)},
    **{f'bilayer_{i}': {'sites': [(0, 0, 0)] + _shell_sites_in_plane(i)
                                  + _shell_sites_above_plane(i)[:1],
                        'plane_dependent': True, 'plane_idx': i, 'N': 8}
       for i in range(4)},
}


# ============================================================================
# COMPOSITION
# ============================================================================

@dataclass(frozen=True)
class Composition:
    """A multiset of building blocks composed around the origin."""
    blocks: tuple[str, ...]          # sorted tuple of block names

    def sites(self) -> frozenset:
        s: set = set()
        for name in self.blocks:
            s.update(BUILDING_BLOCKS[name]['sites'])
        return frozenset(s)

    def N(self) -> int:
        return len(self.sites())

    @property
    def block_counter(self) -> Counter:
        return Counter(self.blocks)

    def canonical_form(self, group) -> tuple:
        return canonical_under_group(self.sites(), group)


def induced_graph(sites) -> nx.Graph:
    sites = sorted(sites)
    idx = {s: i for i, s in enumerate(sites)}
    G = nx.Graph()
    for s, i in idx.items():
        G.add_node(i, coord=s)
    for s, i in idx.items():
        for d in NN_OFFSETS:
            t = (s[0] + d[0], s[1] + d[1], s[2] + d[2])
            j = idx.get(t)
            if j is not None and j > i:
                G.add_edge(i, j, character=edge_character(d), offset=d)
    return G


# ============================================================================
# NAUTY CERTIFICATE (with edge colouring)
# ============================================================================

_EDGE_COLOUR_MAP = {
    ('A', '+'): 0, ('A', '-'): 1,
    ('B', '+'): 2, ('B', '-'): 3,
    ('C', '+'): 4, ('C', '-'): 5,
}


def nauty_cert(G: nx.Graph) -> bytes:
    n_orig = G.number_of_nodes()
    m = G.number_of_edges()
    pg = pynauty.Graph(n_orig + m)
    nodes = sorted(G.nodes())
    remap = {v: i for i, v in enumerate(nodes)}
    adj: dict = {i: [] for i in range(n_orig + m)}
    dummy_by_colour = [set() for _ in range(6)]
    next_dummy = n_orig
    for u, v, d in G.edges(data=True):
        a, b = remap[u], remap[v]
        dv = next_dummy
        next_dummy += 1
        adj[a].append(dv)
        adj[dv].append(a)
        adj[b].append(dv)
        adj[dv].append(b)
        dummy_by_colour[_EDGE_COLOUR_MAP[d['character']]].add(dv)
    for i, nbrs in adj.items():
        if nbrs:
            pg.connect_vertex(i, sorted(set(nbrs)))
    colouring = [set(range(n_orig))]
    for cs in dummy_by_colour:
        if cs:
            colouring.append(cs)
    pg.set_vertex_coloring(colouring)
    return pynauty.certificate(pg)


# ============================================================================
# ENUMERATOR — compositions up to max_N total vertices
# ============================================================================

def enumerate_compositions(
    max_n: int,
    *,
    max_blocks: int = 8,
    symmetry: str = 'oh',
    verbose: bool = True,
) -> list[dict]:
    """Enumerate all block compositions with |vertices| <= max_n.

    Deduplicates by graph iso-class (nauty cert with edge-character colouring,
    applied AFTER symmetry-canonicalisation of the site set) and records
    multiplicity = how many compositions share that canonical form.

    max_blocks: cap on total number of blocks in any composition.
    symmetry:   'oh' or 'd3d' — the point group under which compositions
                that differ by a rotation are considered equivalent.
                'oh' collapses all 4 {111}-plane orientations.
                'd3d' keeps them distinct (monograph convention).
    """
    group = OH_GROUP if symmetry == 'oh' else D3D_GROUP

    block_names = list(BUILDING_BLOCKS.keys())
    N_of: dict = {b: BUILDING_BLOCKS[b]['N'] for b in block_names}
    requires: dict = {
        b: BUILDING_BLOCKS[b].get('requires')
        for b in block_names
    }

    t0 = time.time()
    # Key by (canonical-form hash) -> list of (composition, full_sites)
    all_classes: dict = defaultdict(list)
    total_compositions = 0

    def composition_admissible(combo):
        names = set(combo)
        for b in combo:
            if requires[b] is not None and requires[b] not in names:
                return False
        if 'cell_pair' in names and any(
            n.startswith(('hex_cap', 'coord_shell', 'bilayer', 'tet'))
            for n in names
        ):
            return False
        return True

    for k in range(1, max_blocks + 1):
        # Compositions are SETS of blocks (no repeats — repeating a block
        # doesn't change the vertex set).
        for combo in itertools.combinations(block_names, k):
            if not composition_admissible(combo):
                continue
            comp = Composition(blocks=tuple(sorted(combo)))
            sites = comp.sites()
            if not (2 <= len(sites) <= max_n):
                continue
            # Canonicalise under group FIRST, then build graph
            canon_sites = canonical_under_group(sites, group)
            canon_G = induced_graph(canon_sites)
            if not nx.is_connected(canon_G):
                continue
            cert = nauty_cert(canon_G)
            all_classes[cert].append((comp, canon_sites))
            total_compositions += 1

    dt = time.time() - t0

    records = []
    for cert, items in all_classes.items():
        # Pick the MINIMAL composition (fewest blocks) as the representative.
        items.sort(key=lambda x: (len(x[0].blocks), x[0].blocks))
        rep_comp, rep_sites = items[0]
        G = induced_graph(rep_sites)
        deg_seq = sorted(d for _, d in G.degree())
        edge_hist: Counter = Counter()
        for _, _, d in G.edges(data=True):
            edge_hist[d['character']] += 1

        # All compositions (as block sets) that produce this iso class.
        distinct_block_sets = sorted({tuple(sorted(c.blocks)) for c, _ in items})
        # Minimal compositions only — those not strictly containing a simpler
        # composition on the same iso class.
        set_list = [set(bs) for bs in distinct_block_sets]
        minimal_idx = [
            i for i, si in enumerate(set_list)
            if not any(sj < si for j, sj in enumerate(set_list) if i != j)
        ]
        minimal_compositions = [list(distinct_block_sets[i]) for i in minimal_idx]

        # O_h orbit multiplicity: the {111}-orientation count this class covers.
        # We compute it as |G| / |stabiliser|.  Stabiliser = elements of group
        # that map the canonical-sites set to itself.
        stab_size = sum(
            1 for R in group
            if tuple(sorted(apply_rotation(R, s) for s in rep_sites)) == rep_sites
        )
        oh_orbit_size = len(group) // stab_size if stab_size > 0 else len(group)

        records.append({
            'iso_hash':             blake2b(cert, digest_size=8).hexdigest(),
            'cert_hex':             cert.hex(),
            'N':                    len(rep_sites),
            'E':                    G.number_of_edges(),
            'degree_sequence':      deg_seq,
            'edge_hist':            {f"{f}{s}": c for (f, s), c in edge_hist.items()},
            'minimal_compositions': minimal_compositions,
            'n_block_sets':         len(distinct_block_sets),
            f'{symmetry}_orbit_size': oh_orbit_size,
            'stabiliser_size':      stab_size,
            'is_bipartite':         nx.is_bipartite(G),
            'n_triangles':          sum(nx.triangles(G).values()) // 3,
            'canonical_sites':      [list(s) for s in rep_sites],
        })

    records.sort(key=lambda r: (r['N'], r['iso_hash']))

    if verbose:
        print(f"\nEnumerated {total_compositions} compositions "
              f"-> {len(records)} distinct iso classes "
              f"(N <= {max_n}, max_blocks={max_blocks}, sym={symmetry}) "
              f"in {dt:.2f}s", file=sys.stderr)
        by_N: dict = defaultdict(int)
        for r in records:
            by_N[r['N']] += 1
        print("Iso classes by N:", file=sys.stderr)
        for n in sorted(by_N):
            print(f"  N={n:3d}: {by_N[n]} classes", file=sys.stderr)

    return records


# ============================================================================
# CLI
# ============================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Composition-based catalogue of FCC defect cores.")
    ap.add_argument('--nmax', type=int, default=25,
                    help="Maximum total vertex count of any composition.")
    ap.add_argument('--max-blocks', type=int, default=8,
                    help="Cap on number of blocks per composition.")
    ap.add_argument('--symmetry', choices=['oh', 'd3d'], default='oh')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args(argv)

    records = enumerate_compositions(
        max_n=args.nmax,
        max_blocks=args.max_blocks,
        symmetry=args.symmetry,
        verbose=not args.quiet,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({
        'nmax': args.nmax,
        'max_blocks': args.max_blocks,
        'symmetry': args.symmetry,
        'n_building_blocks': len(BUILDING_BLOCKS),
        'building_blocks': {
            name: {'N': info['N'],
                   'plane_dependent': info['plane_dependent'],
                   'plane_idx': info.get('plane_idx')}
            for name, info in BUILDING_BLOCKS.items()
        },
        'n_iso_classes': len(records),
        'classes': records,
    }, indent=2))

    print(f"Catalogue ({len(records)} classes) -> {args.out}",
          file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
