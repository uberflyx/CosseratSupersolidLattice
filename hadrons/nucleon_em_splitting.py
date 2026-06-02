#!/usr/bin/env python3
"""
nucleon_em_splitting.py
=======================
The electromagnetic part of the neutron-proton mass splitting, computed
directly on the FCC lattice as the Coulomb energy of the three quark
charges of a baryon.

Physical setting
----------------
A baryon is a Y-junction of three partial dislocations (the three colours),
collapsed onto the 13-node coordination shell of the FCC lattice. The three
quark charges sit on one <111> triangular face of the cuboctahedron: they are
the colour triple carried round by the threefold rotation about that axis, so
they lie at mutual nearest-neighbour spacing. A clean geometric fact, verified
below, is that each quark then sits exactly on a common-neighbour site of the
other two, which is where the "third-arm" enhancement of the pair Coulomb must
live.

The proton (uud) and neutron (udd) differ only in the charge bilinear
sum_{i<j} q_i q_j: 0 for the proton, -1/3 for the neutron, a difference odd in
the isospin projection I_3, so the splitting preserves the isoscalar average.

Three tiers for the nucleon Coulomb (all parameter-free except the bond length):
  1. bare point charge          -- the FCC lattice Coulomb between two
                                    nearest-neighbour fractional charges;
  2. cell-pair distributed       -- the charge spread over the five
                                    charge-active bonds of the cell pair,
                                    5 m_e per pair (the same count that fixes
                                    the up-down mass difference);
  3. observed-implied            -- the value the observed splitting wants,
                                    15/2 m_e per pair, reached only with the
                                    third-arm common-neighbour overlap.

The on-site lattice self-energy is the divergent piece that renormalises into
the mass; it does not enter the splitting and is reported only to show that
including it naively is nonsensical.

Accompanies: M. A. Cox, "The Cosserat Supersolid" (2026).
https://doi.org/10.5281/zenodo.18636501

Inputs: the electron mass m_e (CODATA 2022), the FCC geometry, and the
bond-length identity l = r_e, by which two unit charges one bond apart cost
e^2/(4 pi eps0 r_e) = m_e c^2 = 1 m_e. No fitted parameters.

Author: Mitchell A. Cox, University of the Witwatersrand
"""

import numpy as np
from itertools import combinations, product
from numba import njit, prange

# ── Physical constants ──
M_E_MEV = 0.51099895069       # electron rest energy [MeV], CODATA 2022

# Quark electric charges in units of e (screw content q = cos beta).
Q_UP, Q_DOWN = 2.0 / 3.0, -1.0 / 3.0

# Charge-active-bond count of a nearest-neighbour cell pair: the direct bond
# plus the four common nearest neighbours on the perpendicular bisector. This
# is the count that fixes m_d - m_u = 5 m_e (the strong-isospin part).
N_CB_CELL_PAIR = 5.0

# ── FCC geometry (lengths in units of the cubic lattice constant a) ──

def fcc_nearest_neighbours():
    """The 12 nearest neighbours of an FCC site (a cuboctahedron).

    Returned in units of a; the nearest-neighbour distance is a/sqrt(2).
    """
    return 0.5 * np.array(
        [v for v in product((-1, 0, 1), repeat=3) if sum(abs(c) for c in v) == 2],
        dtype=float,
    )


def baryon_face_sites():
    """The three quark sites: the C3 orbit about [111] nearest that axis.

    These are the vertices of one <111> triangular face of the cuboctahedron.
    """
    return {
        "A": np.array([0.5, 0.5, 0.0]),
        "B": np.array([0.0, 0.5, 0.5]),
        "C": np.array([0.5, 0.0, 0.5]),
    }


# FCC reciprocal primitive vectors (a = 1). Used for the Brillouin-zone sum.
B1 = 2.0 * np.pi * np.array([-1.0, 1.0, 1.0])
B2 = 2.0 * np.pi * np.array([1.0, -1.0, 1.0])
B3 = 2.0 * np.pi * np.array([1.0, 1.0, -1.0])

NN_DIST = 1.0 / np.sqrt(2.0)   # nearest-neighbour distance in units of a


def verify_geometry():
    """Check the load-bearing geometric claims; return (passes, fails)."""
    sites = baryon_face_sites()
    names = list(sites)
    passes, fails = 0, 0

    def check(name, cond):
        nonlocal passes, fails
        ok = bool(cond)
        passes, fails = passes + ok, fails + (not ok)
        print(f"  [{'PASS' if ok else 'FAIL'}]  {name}")

    def is_nn(p, q):
        return abs(np.linalg.norm(p - q) - NN_DIST) < 1e-9

    for n1, n2 in combinations(names, 2):
        check(f"{n1}-{n2} at nearest-neighbour spacing",
              is_nn(sites[n1], sites[n2]))

    # Each quark sits on a common-neighbour site of the other two.
    for i, j, k in [("A", "B", "C"), ("B", "C", "A"), ("C", "A", "B")]:
        check(f"{k} is a common neighbour of {i}-{j}",
              is_nn(sites[k], sites[i]) and is_nn(sites[k], sites[j]))

    return passes, fails


# ── Isovector charge sums on the triangular face ──

def isovector_sums(charges):
    """Return (sum q_i^2, sum_{i<j} q_i q_j, |dipole|^2) for three charges.

    The three sit on the equilateral triangle of vertices A, B, C above.
    """
    sites = list(baryon_face_sites().values())
    sum_q2 = sum(q * q for q in charges)
    sum_pairs = sum(charges[a] * charges[b] for a, b in combinations(range(3), 2))
    dipole = sum(q * r for q, r in zip(charges, sites))
    return sum_q2, sum_pairs, float(dipole @ dipole)


# ── FCC lattice Coulomb Green's function ──

@njit(parallel=True, cache=True)
def _green_bz_sum(targets, nn, b1, b2, b3, ngrid):
    """Brillouin-zone sum of the inverse lattice Laplacian.

    G(r) = (1/Nk) sum_k cos(k.r) / L(k),  L(k) = sum_delta (1 - cos k.delta),
    the sum running over the 12 nearest-neighbour bond vectors. A midpoint
    k-grid avoids the k = 0 singularity. Each outer index writes to its own
    row, so the parallel loop carries no race.
    """
    m = targets.shape[0]
    partial = np.zeros((ngrid, m))
    for i1 in prange(ngrid):
        u1 = (i1 + 0.5) / ngrid
        for i2 in range(ngrid):
            u2 = (i2 + 0.5) / ngrid
            for i3 in range(ngrid):
                u3 = (i3 + 0.5) / ngrid
                kx = u1 * b1[0] + u2 * b2[0] + u3 * b3[0]
                ky = u1 * b1[1] + u2 * b2[1] + u3 * b3[1]
                kz = u1 * b1[2] + u2 * b2[2] + u3 * b3[2]
                lap = 0.0
                for j in range(nn.shape[0]):
                    lap += 1.0 - np.cos(kx * nn[j, 0] + ky * nn[j, 1] + kz * nn[j, 2])
                if lap > 1e-9:
                    inv = 1.0 / lap
                    for t in range(m):
                        kr = kx * targets[t, 0] + ky * targets[t, 1] + kz * targets[t, 2]
                        partial[i1, t] += np.cos(kr) * inv
    out = np.zeros(m)
    for i1 in range(ngrid):
        for t in range(m):
            out[t] += partial[i1, t]
    return out / (ngrid ** 3)


def lattice_green(positions, ngrid=96):
    """Lattice Coulomb Green's function G(r) at each position in `positions`.

    `positions` is a list of length-3 arrays (units of a). `ngrid` is the
    linear resolution of the Brillouin-zone sum; the on-site value G(0)
    converges from above as ngrid grows, the ratios faster.
    """
    nn = fcc_nearest_neighbours()
    targets = np.ascontiguousarray(np.array(positions, dtype=float))
    vals = _green_bz_sum(targets, nn, B1, B2, B3, int(ngrid))
    return {tuple(np.round(p, 6)): v for p, v in zip(positions, vals)}


def continuum_factor(ngrid=96):
    """Lattice Coulomb at nearest-neighbour, in units of its continuum value.

    Fits G(r) ~ C/r along the <110> nearest-neighbour chain to extract the
    long-range 1/r coefficient, then divides the nearest-neighbour value by
    C/r_NN. Two unit charges one bond apart cost 1 m_e in the continuum, so
    this factor times 1 m_e is the bare point-charge cell-pair Coulomb.
    """
    chain = [np.array([0.5 * k, 0.5 * k, 0.0]) for k in range(1, 11)]
    g = lattice_green(chain, ngrid)
    r = np.array([0.5 * k * np.sqrt(2.0) for k in range(1, 11)])
    gv = np.array([g[tuple(np.round(p, 6))] for p in chain])
    c_inf = float(np.mean((gv * r)[-3:]))     # large-r limit of G(r)*r
    g_nn = gv[0]
    return (g_nn * r[0]) / c_inf


# ── Report ──

def main():
    print("=" * 70)
    print("Neutron-proton electromagnetic splitting on the FCC lattice")
    print("=" * 70)

    print("\nGeometry of the three quark charges (<111> triangular face):")
    passes, fails = verify_geometry()
    print(f"  -> {passes} passed, {fails} failed")

    print("\nIsovector charge content (charges on the triangle):")
    print(f"  {'state':10}{'sum q^2':>10}{'sum_<ij q_iq_j':>16}{'|dipole|^2':>12}")
    proton = isovector_sums([Q_UP, Q_UP, Q_DOWN])
    neutron = isovector_sums([Q_UP, Q_DOWN, Q_DOWN])
    for label, s in (("proton uud", proton), ("neutron udd", neutron)):
        print(f"  {label:10}{s[0]:10.4f}{s[1]:16.4f}{s[2]:12.4f}")
    d_sum_pairs = neutron[1] - proton[1]
    iso_avg = 0.5 * (proton[2] + neutron[2])
    print(f"  n-p:  d(sum_<) = {d_sum_pairs:+.4f}   "
          f"d|dipole|^2 = {neutron[2] - proton[2]:+.4f}")
    print(f"  |dipole|^2 isoscalar average = {iso_avg:.4f}; "
          f"proton and neutron sit symmetrically about it (odd in I_3).")

    print("\nFCC lattice Coulomb Green's function (Brillouin-zone sum):")
    onsite, nn_pos = np.array([0.0, 0.0, 0.0]), np.array([0.5, 0.5, 0.0])
    for ngrid in (64, 96, 128):
        g = lattice_green([onsite, nn_pos], ngrid)
        g0, gnn = g[(0.0, 0.0, 0.0)], g[(0.5, 0.5, 0.0)]
        print(f"  ngrid={ngrid:3d}:  G0={g0:.4f}  G_NN={gnn:.4f}  "
              f"G0/G_NN={g0 / gnn:.3f}")
    factor = continuum_factor(128)
    print(f"  lattice Coulomb at NN, in continuum units: {factor:.3f} "
          f"(only ~{100 * (factor - 1):.0f}% above the 1/r tail)")

    print("\nNucleon Coulomb, three tiers (n-p = N_CB m_e * d(sum_<)):")
    print(f"  {'tier':34}{'N_CB [m_e]':>12}{'n-p EM [MeV]':>14}{'total [MeV]':>13}")
    tiers = [
        ("bare point charge (lattice)", factor),
        ("cell-pair distributed (= 5)", N_CB_CELL_PAIR),
        ("observed-implied (= 15/2)", 7.5),
    ]
    strong = N_CB_CELL_PAIR * M_E_MEV          # m_d - m_u = 5 m_e, the strong part
    for label, ncb in tiers:
        em = ncb * d_sum_pairs * M_E_MEV
        print(f"  {label:34}{ncb:12.3f}{em:14.3f}{strong + em:13.3f}")
    print(f"  strong part m_d - m_u = 5 m_e = {strong:.3f} MeV; "
          f"observed n-p = 1.293 MeV; BMW EM part = -1.00 MeV.")

    print("\nWhy the on-site self-energy is excluded:")
    g = lattice_green([onsite, nn_pos], 128)
    g0_over_gnn = g[(0.0, 0.0, 0.0)] / g[(0.5, 0.5, 0.0)]
    self_term = 0.5 * (N_CB_CELL_PAIR * g0_over_gnn) * (neutron[0] - proton[0])
    inter_term = N_CB_CELL_PAIR * d_sum_pairs
    print(f"  adding the bare self-energy (G0/G_NN = {g0_over_gnn:.2f}) gives "
          f"n-p EM = {(self_term + inter_term) * M_E_MEV:+.3f} MeV,")
    print(f"  total {strong + (self_term + inter_term) * M_E_MEV:+.3f} MeV: "
          f"nonsensical, and it double-counts the screw self-energy already")
    print(f"  inside m_d - m_u. The renormalised splitting keeps the "
          f"inter-quark term only.")


if __name__ == "__main__":
    main()
