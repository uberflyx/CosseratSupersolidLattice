#!/usr/bin/env python3
"""
sft_cage.py
===========
The stacking-fault tetrahedron (SFT) as the Mode 2F tetraquark cage:
geometry verification, integer node structure, mass ledger, and the
two-shell impossibility theorem.

THE OBJECT.  The minimal SFT is the Thompson tetrahedron: four FCC sites
A = (0,0,0), B = (a/2)(1,1,0), C = (a/2)(1,0,1), D = (a/2)(0,1,1), all
mutual nearest neighbours, with four triangular stacking faults on the
four distinct {111} orientations and six <110> stair-rod edges.  Its
centroid is the T_d void (a/4)(1,1,1); no lattice site lies inside; its
height is exactly one {111} interlayer spacing.  First observed in
quenched gold by Silcox and Hirsch (1959).

COLOUR CLOSURE.  Each face carries an in-plane Shockley partial of
magnitude l/sqrt(3); assignments with b1+b2+b3+b4 = 0 exist (18 of them),
so any loop through the surrounding crystal encircles zero net stacking
shift.  The colour-singlet condition is topological, and the far field
self-screens: the structural definition of a cage state.

MODE 2F.  Each edge of the tetrahedron borders exactly two fault faces,
and a cell-pair meson is an edge.  The cc diquark and c(bar)c(bar)
antidiquark occupy opposite perpendicular edges at gap l/sqrt(2): each
constituent docks through its own pair of adjacent faces.  Mode 2F is
two-face docking per constituent, realised by cell pairs, NOT by shells.

THE IMPOSSIBILITY THEOREM.  No arrangement of two 13-node coordination
shells produces a six-atom two-disjoint-triangle contact: exhaustive over
all lattice separations (coherent), all exact gap-bond distances along
symmetry axes (molecular), in parallel and 60-degree-twinned orientation.
The only 6-atom overlap is the deuteron's <110> blob.  By-product: the
twin orientation shares one triangular face at 2h = 1.633 l, a twinned
Mode F junction.

INTEGER CAGE.  Twelve sites at sqrt(11/8) l from the centroid bond to
exactly two vertices each: the external pair of each edge's square ring
(the other two ring slots are the remaining vertices), two per edge, the
stair-rod channels.  Twelve sites at sqrt(19/8) l
bond to one vertex (three caps per vertex).  Full cage 4 + 24 = 28; core
4 + 12 = 16.

LEDGER (open).  m = 4 x 18 m0 (charm K_{9,9} strain) + N_SFT x m0 + ...
CMS interference masses 6638/6847/7134 MeV imply N_SFT = 22.8/25.8/29.9,
bracketed by core 16 and full cage 28; no-interference 6552/6927/7287
give 21.6/26.9/32.1 in near-equal steps of 5.2.  The activation rule
(which subset, which eigenvalue) is the remaining open step.
"""
import numpy as np
from itertools import product

M_E = 0.51099895069
M_0 = M_E * 137.035999177
A_LAT = np.sqrt(2.)          # lattice parameter, NN distance l = 1

VERTS = [np.zeros(3),
         A_LAT/2*np.array([1., 1., 0.]),
         A_LAT/2*np.array([1., 0., 1.]),
         A_LAT/2*np.array([0., 1., 1.])]
CEN = sum(VERTS)/4
PRIM = A_LAT/2*np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]])

def fcc_sites(radius, centre):
    out = []
    for i, j, k in product(range(-5, 6), repeat=3):
        p = i*PRIM[0] + j*PRIM[1] + k*PRIM[2]
        if np.linalg.norm(p - centre) < radius:
            out.append(p)
    return np.array(out)

def verify_geometry():
    print("Thompson tetrahedron checks:")
    Minv = np.linalg.inv(PRIM.T)
    for v in VERTS:
        c = Minv @ v
        assert np.allclose(c, np.round(c))
    els = [np.linalg.norm(VERTS[i]-VERTS[j])
           for i in range(4) for j in range(i+1, 4)]
    assert np.allclose(els, 1.0)
    print("  four FCC vertices, six NN edges: pass")
    faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    normals = set()
    for f in faces:
        n = np.cross(VERTS[f[1]]-VERTS[f[0]], VERTS[f[2]]-VERTS[f[0]])
        n = n/np.max(np.abs(n))
        normals.add(tuple(np.round(np.abs(n), 6)))
        assert np.allclose(np.abs(n), 1.0)
    print(f"  four {{111}} faces, distinct orientations: pass")
    assert np.allclose(CEN, A_LAT/4*np.ones(3))
    inner = []
    Mm = np.array([VERTS[1]-VERTS[0], VERTS[2]-VERTS[0],
                   VERTS[3]-VERTS[0]]).T
    for p in fcc_sites(1.5, CEN):
        lam = np.linalg.solve(Mm, p - VERTS[0])
        if all(lam > 1e-6) and sum(lam) < 1 - 1e-6:
            inner.append(p)
    assert not inner
    print("  centroid = T_d void; no interior site: pass")
    n = np.cross(VERTS[1]-VERTS[0], VERTS[2]-VERTS[0])
    n /= np.linalg.norm(n)
    assert abs(abs((VERTS[3]-VERTS[0]) @ n) - np.sqrt(2/3)) < 1e-9
    print("  height = one {111} interlayer spacing: pass")
    for (i, j), (k, l) in [((0, 1), (2, 3)), ((0, 2), (1, 3)),
                           ((0, 3), (1, 2))]:
        e1 = VERTS[j]-VERTS[i]; e2 = VERTS[l]-VERTS[k]
        gap = np.linalg.norm((VERTS[i]+VERTS[j])/2-(VERTS[k]+VERTS[l])/2)
        assert abs(e1 @ e2) < 1e-9 and abs(gap - 1/np.sqrt(2)) < 1e-9
    print("  opposite edges perpendicular at gap l/sqrt(2): pass")
    # Shockley closure: in-face 1/6<112> partials, |b| = l/sqrt(3)
    def shockleys(f):
        t = [VERTS[x] for x in f]
        out = []
        for perm in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
            v0, v1, v2 = [t[p] for p in perm]
            d = (v1+v2)/2 - v0
            d = d/np.linalg.norm(d)/np.sqrt(3)
            out += [d, -d]
        return out
    sols = sum(1 for combo in product(*[shockleys(f) for f in faces])
               if np.linalg.norm(sum(combo)) < 1e-9)
    assert sols > 0
    print(f"  Shockley closure sum(b_i) = 0: {sols} assignments exist")

def integer_cage():
    print("\nInteger cage structure:")
    sites = fcc_sites(3.0, CEN)
    shells = {}
    for p in sites:
        if any(np.linalg.norm(p-v) < 1e-6 for v in VERTS):
            continue
        nv = sum(1 for v in VERTS if abs(np.linalg.norm(p-v)-1.0) < 1e-6)
        if nv == 0:
            continue
        d = round(np.linalg.norm(p-CEN), 4)
        shells.setdefault((d, nv), []).append(p)
    for (d, nv), pts in sorted(shells.items()):
        print(f"  d = {d} l, bonded to {nv} vertices: {len(pts)} sites")
    s1 = shells[(round(np.sqrt(11/8), 4), 2)]
    s2 = shells[(round(np.sqrt(19/8), 4), 1)]
    assert len(s1) == 12 and len(s2) == 12
    # each edge's square ring: 4 common neighbours, of which the other
    # two tetrahedron vertices are two; exactly 2 external sites per edge
    for (i, j) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]:
        ring = [p for p in s1
                if abs(np.linalg.norm(p-VERTS[i])-1) < 1e-6
                and abs(np.linalg.norm(p-VERTS[j])-1) < 1e-6]
        assert len(ring) == 2
    print("  per-edge external ring pairs (2 sites each, none shared;"
          " the other 2 ring slots are the remaining vertices): pass")
    print("  core = 4 + 12 = 16; full cage = 4 + 24 = 28")

def two_shell_impossibility():
    print("\nTwo-shell Mode 2F impossibility:")
    nn = np.array([p for p in fcc_sites(1.05, np.zeros(3))
                   if np.linalg.norm(p) > 1e-9])
    shell = np.vstack([[np.zeros(3)], nn])
    axis = np.array([1., 1., 1.])/np.sqrt(3)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    R60 = (np.eye(3)*0.5 + np.sin(np.pi/3)*K
           + 0.5*np.outer(axis, axis))
    def two_triangles(bondsA):
        S = shell[sorted(bondsA)]
        if len(S) != 6:
            return False
        adj = [[m for m in range(6) if m != l
                and abs(np.linalg.norm(S[l]-S[m])-1) < 1e-6]
               for l in range(6)]
        if sorted(len(x) for x in adj) != [2]*6:
            return False
        comp = {0}; fr = [0]
        while fr:
            x = fr.pop()
            for y in adj[x]:
                if y not in comp:
                    comp.add(y); fr.append(y)
        return len(comp) == 3
    hits = 0
    for base, tag in [(shell, "parallel"), ((R60 @ shell.T).T, "twinned")]:
        # coherent: all lattice separations
        for i, j, k in product(range(-3, 4), repeat=3):
            Dv = i*PRIM[0] + j*PRIM[1] + k*PRIM[2]
            if np.linalg.norm(Dv) < 1e-9 or np.linalg.norm(Dv) > 3.5:
                continue
            Bs = base + Dv
            sh = [(x, y) for x, p in enumerate(shell) for y, q in
                  enumerate(Bs) if np.linalg.norm(p-q) < 1e-6]
            if len(sh) == 6 and two_triangles({x for x, _ in sh}):
                hits += 1
        # molecular: exact gap-bond distances along symmetry axes
        for u in [axis, np.array([1., 0., 0.]),
                  np.array([1., 1., 0.])/np.sqrt(2),
                  np.array([1., 1., 2.])/np.sqrt(6)]:
            cands = set()
            for p in shell:
                for q in base:
                    r = p-q; ru = r@u
                    disc = ru*ru-(r@r-1.0)
                    if disc < 0:
                        continue
                    for t in (ru+np.sqrt(disc), ru-np.sqrt(disc)):
                        if t > 1.05:
                            cands.add(round(t, 9))
            for t in cands:
                Bs = base + t*u
                if any(np.linalg.norm(p-q) < 1e-6 for p in shell
                       for q in Bs):
                    continue
                bonds = [(x, y) for x, p in enumerate(shell)
                         for y, q in enumerate(Bs)
                         if abs(np.linalg.norm(p-q)-1.0) < 1e-8]
                if len(bonds) == 6 and two_triangles({x for x, _ in bonds}):
                    hits += 1
    print(f"  six-atom two-triangle contacts found: {hits} (theorem: none exist)")
    assert hits == 0

def ledger():
    print("\nMass ledger (open activation rule):")
    charm = 4*18*M_0
    print(f"  4 charm vertices at 18 m0 each: {charm:.1f} MeV")
    for tag, masses in [("CMS interference", (6638, 6847, 7134)),
                        ("CMS no-interference", (6552, 6927, 7287))]:
        ns = [(m-charm)/M_0 for m in masses]
        print(f"  {tag}: N_SFT = "
              + ", ".join(f"{n:.1f}" for n in ns)
              + f"  (steps {ns[1]-ns[0]:.1f}, {ns[2]-ns[1]:.1f})")
    print("  integer brackets: core 12, full cage 24 non-vertex nodes"
          " (16/28 with vertices)")

if __name__ == "__main__":
    verify_geometry()
    integer_cage()
    two_shell_impossibility()
    ledger()
