"""Detours through the compact axis of D4: what a hop onto the other layer can
and cannot do to an excitation of the FCC slice.

The neutrino chapter reads the antineutrino as a hole in the filled band of the
medium's constituent fermions.  Half of a site's 24 nearest-neighbour bonds in
D4 lean 45 degrees into the compact axis, so a hole can step onto the other
layer and back.  This script checks the four statements the chapter makes
about such detours.

Units: the D4 integer unit s = l/sqrt(2), so the nearest-neighbour distance is
sqrt(2) s = l.  D4 is the set of integer points with even coordinate sum; the
compact axis x4 carries two layers per period L4 = 2 s, and the layer x4 = 2 is
the layer x4 = 0 again.

1. Where the other layer sits.  The layer x4 = 1 holds the points whose spatial
   coordinate sum is odd.  Seen from the slice these are the octahedral
   interstitial sites of FCC: each lies half-way between two neighbouring
   close-packed {111} planes, directly over the stacking position that neither
   plane occupies.

2. Ordinary hops already change plane.  Six of the slice's twelve nearest
   neighbours lie one {111} plane up or down, so the plane an excitation sits on
   is part of its position, not a label it carries.

3. Where a round trip lands.  Two crossing steps that return to the slice,
   including the return through the wrap-around x4 = 2 == 0, add to a lattice
   vector of the FCC slice (19 distinct landings, the origin included).  The
   Shockley vector a/6<112>, which changes a plane's stacking position, is not
   among them.

4. Interference.  With hop -t on the twelve slice bonds and
   h(n) = -t_e + i t_o n4 on the twelve crossing bonds (t_e, t_o the parts
   even and odd in the crossing direction n4 = +-1), the band is

       eps(k) = -4 t (c1 c2 + c2 c3 + c3 c1)
                - 4 (c1 + c2 + c3) [t_e cos(k4 s) - t_o sin(k4 s)],

   with c_i = cos(k_i s).  In the antiperiodic sector, k4 s = +-pi/2, the even
   part cancels (the two routes to the other layer are one circuit apart and
   carry opposite signs: an Aharonov-Bohm cage) and the odd part splits the
   two states by 8 t_o (c1 + c2 + c3).  The odd part changes sign under the
   compact mirror x4 -> -x4.

Every statement is asserted, and the script prints one line for each that holds.
"""

import itertools

import numpy as np

# D4 nearest neighbours: the 24 vectors +-e_i +- e_j.
NN = np.array([v for v in itertools.product((-1, 0, 1), repeat=4)
               if sorted(map(abs, v)) == [0, 0, 1, 1]], dtype=int)
IN_SLICE = NN[NN[:, 3] == 0]      # the 12 FCC nearest neighbours
CROSSING = NN[NN[:, 3] != 0]      # the 12 bonds leaning into x4


def in_fcc(v3):
    """True for a lattice point of the FCC slice (integer, even coordinate sum)."""
    return int(sum(v3)) % 2 == 0


def plane_index(v3):
    """Height along [111] in units of the plane spacing d111 = 2 s / sqrt(3)."""
    return sum(v3) / 2


def check_layer_geometry():
    """Statement 1: odd sites sit half-way between planes, over the empty position."""
    odd = [np.array(p) for p in itertools.product(range(-3, 4), repeat=3)
           if not in_fcc(p)]
    for p in odd:
        h = int(sum(p))
        assert abs(plane_index(p) % 1 - 0.5) < 1e-12          # half-way
        # The [111] column through p meets FCC sites p + t(1,1,1), t odd,
        # at heights (h + 3t)/2: one stacking label, A, B or C (mod 3).
        column = {((h + 3 * t) // 2) % 3 for t in (-3, -1, 1, 3)}
        assert len(column) == 1
        below, above = ((h - 1) // 2) % 3, ((h + 1) // 2) % 3
        assert column.pop() not in (below, above)
    print("1. layer x4 = 1: octahedral sites, half-way between {111} planes, "
          "over the stacking position neither neighbouring plane occupies")


def check_plane_changes():
    """Statement 2: half of the in-slice hops climb or descend one plane."""
    heights = [plane_index(v[:3]) for v in IN_SLICE]
    counts = (heights.count(0), heights.count(1), heights.count(-1))
    assert counts == (6, 3, 3)
    print("2. in-slice nearest neighbours: 6 in-plane, 3 one plane up, 3 one down")


def check_round_trips():
    """Statement 3: round trips land on FCC lattice vectors only."""
    landings = {tuple(int(x) for x in (a + b)[:3])
                for a in CROSSING for b in CROSSING
                if (a[3] + b[3]) % 2 == 0}
    assert len(landings) == 19 and all(in_fcc(v) for v in landings)
    shockley = np.array([2, -1, -1]) / 3.0     # a/6 [2 -1 -1], cubic edge a = 2 s
    assert np.isclose(np.linalg.norm(shockley), np.sqrt(2 / 3))   # = l/sqrt(3)
    assert not np.allclose(shockley, np.round(shockley))
    print("3. 19 round-trip landings, all FCC lattice vectors; "
          "the Shockley vector is not a lattice vector")


def band_sum(k, t=1.0, t_e=0.7, t_o=0.0):
    """eps(k) = sum_n h(n) exp(-i k.n), k in units of 1/s."""
    total = sum(-t * np.exp(-1j * k @ n) for n in IN_SLICE)
    total += sum((-t_e + 1j * t_o * n[3]) * np.exp(-1j * k @ n) for n in CROSSING)
    return total


def band_closed(k, t=1.0, t_e=0.7, t_o=0.0):
    """The closed form quoted in the module docstring."""
    c = np.cos(k[:3])
    return (-4 * t * (c[0] * c[1] + c[1] * c[2] + c[2] * c[0])
            - 4 * c.sum() * (t_e * np.cos(k[3]) - t_o * np.sin(k[3])))


def check_interference(trials=200, seed=20261010):
    """Statement 4: the closed form, the cancellation and the splitting."""
    rng = np.random.default_rng(seed)
    for _ in range(trials):
        k3 = rng.uniform(-np.pi, np.pi, 3)
        t_e, t_o = rng.normal(size=2)
        for k4 in (0.0, np.pi, np.pi / 2, -np.pi / 2):
            k = np.array([*k3, k4])
            z = band_sum(k, t_e=t_e, t_o=t_o)
            assert abs(z.imag) < 1e-12
            assert abs(z.real - band_closed(k, t_e=t_e, t_o=t_o)) < 1e-12
        kp, km = np.array([*k3, np.pi / 2]), np.array([*k3, -np.pi / 2])
        # even part drops out of the antiperiodic sector
        assert abs(band_closed(kp, t_e=t_e) - band_closed(kp, t_e=0.0)) < 1e-12
        # odd part splits the two antiperiodic states
        split = band_closed(kp, t_e=t_e, t_o=t_o) - band_closed(km, t_e=t_e, t_o=t_o)
        assert abs(split - 8 * t_o * np.cos(k3).sum()) < 1e-12
    # Mirror x4 -> -x4: hopping along the mirrored bonds with the same h(n)
    # gives the band of the original lattice with t_o reversed.
    for _ in range(20):
        k = np.array([*rng.uniform(-np.pi, np.pi, 3), np.pi / 2])
        t_e, t_o = rng.normal(size=2)
        mirrored = sum((-t_e + 1j * t_o * n[3]) * np.exp(-1j * k @ (n * [1, 1, 1, -1]))
                       for n in CROSSING)
        mirrored += sum(-np.exp(-1j * k @ n) for n in IN_SLICE)
        assert abs(mirrored - band_closed(k, t_e=t_e, t_o=-t_o)) < 1e-12
    print("4. band closed form verified; at k4 s = +-pi/2 the even part cancels "
          "and the odd part splits by 8 t_o sum_i cos(k_i s); the odd part is mirror-odd")


if __name__ == "__main__":
    check_layer_geometry()
    check_plane_changes()
    check_round_trips()
    check_interference()
