"""Enumerate the two-shell docking junctions of the FCC coordination shell.

A junction between two coordination shells is fixed by two independent
questions, and keeping them apart is what makes the count come out right:

  REGISTRY  coherent  = B is an FCC translate of A, so both defects sit in
                        one grain
            molecular = B is a reoriented grain meeting A at a boundary
  CONTACT   overlapping = the two shells share atoms
            gapped      = no shared atom; the junction is held by the
                          nearest-neighbour bonds crossing the interface

The criterion for "within reach" is where this enumeration can go wrong.
Sharing an atom is not what makes a junction; carrying a bond across the
interface is, and the two part company at |d|^2 = 8, beyond which the shells
still bond but no longer overlap.  Cutting at shared atoms stops the coherent
family four members short and pushes those four into the molecular family,
where they do not belong: <222> at |d|^2 = 12 carries six interface bonds
across twenty-six positioned nodes, which is the face-docking junction the
P_c pentaquarks occupy.

Reported for every translate out to loss of contact: shared atoms, interface
bonds, positioned nodes and the resulting contact type.

Companion to the compound-defect chapter of the monograph.
"""
import itertools

import numpy as np

NN2 = 2                                   # squared nearest-neighbour distance

NN = [np.array(v) for v in itertools.product((-1, 0, 1), repeat=3)
      if sum(abs(x) for x in v) == 2]
assert len(NN) == 12
shellA = [np.array((0, 0, 0))] + NN       # 13 atoms: centre plus 12 vertices


def key(v):
    return tuple(int(x) for x in v)


Aset = {key(a) for a in shellA}


def fcc_sites(rmax2):
    """Integer points with even coordinate sum, out to |d|^2 = rmax2."""
    R = int(np.ceil(np.sqrt(rmax2))) + 1
    return [(x, y, z)
            for x in range(-R, R + 1)
            for y in range(-R, R + 1)
            for z in range(-R, R + 1)
            if (x + y + z) % 2 == 0 and 0 < x * x + y * y + z * z <= rmax2]


def oh_orbit_key(v):
    """Canonical representative of a vector's O_h orbit, to group separations."""
    return min(tuple(s[i] * v[p[i]] for i in range(3))
               for p in itertools.permutations(range(3))
               for s in itertools.product((1, -1), repeat=3))


def shape(pts):
    """Classify the geometry of a shared-atom set."""
    n = len(pts)
    if n == 0:
        return "none (gapped)"
    if n == 1:
        return "1 node (bridge)"
    P = np.array(pts)
    d2 = sorted({int(round(np.dot(P[i] - P[j], P[i] - P[j])))
                 for i in range(n) for j in range(n) if i < j})
    if n == 2 and d2 == [2]:
        return "2 atoms (edge)"
    if n == 3 and d2 == [2]:
        return "3 atoms: triangular face"
    if n == 4 and d2 == [2, 4]:
        return "4 atoms: square face"
    if n == 6:
        return "6 atoms: rectangular face + 2 cores"
    return f"{n} atoms (dists^2 {d2})"


def interface_bonds(dv):
    """Nearest-neighbour bonds crossing the interface at separation dv."""
    return sum(1 for p in shellA for q in shellA
               if abs(int(np.dot(q + dv - p, q + dv - p)) - NN2) == 0)


MODE = {2: "Mode P (deuteron)",
        8: "Mode V (d*(2380), strange decuplet)",
        12: "Mode F (P_c, P_cs pentaquarks)"}

print("COHERENT junctions: B an FCC translate of A, enumerated by bonded")
print("contact rather than by shared atoms.\n")
print(f"{'d':14s}{'|d|^2':7s}{'shared':8s}{'bonds':7s}{'nodes':7s}"
      f"{'contact':13s}{'shape':38s}{'mode'}")

seen = set()
for d in sorted(fcc_sites(24), key=lambda v: (np.dot(v, v), v)):
    ok = oh_orbit_key(d)
    if ok in seen:
        continue
    dv = np.array(d)
    Bset = {key(a + dv) for a in shellA}
    shared = sorted(Aset & Bset)
    bonds = interface_bonds(dv)
    if not shared and bonds == 0:
        continue                          # contact lost: the family ends here
    seen.add(ok)
    nodes = len(Aset | Bset)
    contact = "overlapping" if shared else "gapped"
    d2 = int(np.dot(d, d))
    print(f"{str(d):14s}{d2:<7d}{len(shared):<8d}{bonds:<7d}{nodes:<7d}"
          f"{contact:13s}{shape(shared):38s}{MODE.get(d2, '')}")

print("\nThe two accountings are not yet a single number.  The severed-bond")
print("count n_bond that the bridge identity consumes is the shared-atom")
print("accounting for the overlapping members (Mode P reads 4 from its four")
print("shared face atoms) and the interface-bond count for the gapped ones")
print("(Mode F reads 6).  That they agree at Mode F, whose twin testbed reads")
print("6 as three shared atoms bonded to both centres, is what licenses the")
print("testbed; whether the agreement is structural is an open item.")

print("\nMOLECULAR junctions: B a reoriented grain.  Run double_cosets.py,")
print("which enumerates the orientations AND tests whether each mirrored")
print("grain actually leaves the lattice.  Only the {111} twin does, so the")
print("lattice admits exactly one genuine two-shell grain boundary.")
