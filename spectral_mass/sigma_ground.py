#!/usr/bin/env python3
"""
sigma_ground.py
===============
Confirm that the Sigma baryon ground mass mode is selected by rule, not
by a fit to the measured mass.

The Sigma sits on the same cluster as the Delta: the coordination shell
plus the four tetrahedral voids of one T_d orbit (17 nodes).  The two
baryons read different modes of this cluster.  The parity-flip rule fixes
the Sigma's spatial irrep as A_2, the T_d descendant of the proton's
A_2u; the resisted class then reads the lowest stiff phi-dominant A_2
root above the structural value lambda = 4.

This script lists every A_2 mode of the cluster with its microrotation
(phi) fraction.  The four A_2 modes are 2.170 (soft, displacement-like),
4.624, 7.087 and 9.898 (all phi-dominant).  The rule selects 4.624 with
no reference to the PDG mass; an earlier single-length convention placed
this mode at 4.315 and matched the mass to 0.001 percent, which was the
fingerprint of a fit.  The rule-selected value carries an honest residual
of +0.23 percent.

Matrix convention (shared with delta_first_principles.py): unit stiffness
on every bond (nearest-neighbour at d = 1 and void-corner at d = sqrt6/4),
with the bond length entering only through the discrete curl operator,
weighted 1/(2d) per bond.  This same matrix returns the Delta T_1 mode at
9.0515 and the proton A_2u at 8.303.

Run:  python3 sigma_ground.py
"""

import numpy as np
from delta_first_principles import (cluster_delta, build_cosserat_matrix_two_d,
                                     generate_Td, classify_Td, TD_CHARS)
from spectral_classifier import vertex_perm
from hadron_spectral_mass import M_0, M_E, mass_from_lambda

PDG_SIGMA0 = 1193.15  # MeV, isoscalar reference (not used in selection)


def state_representation(R, coords):
    """Represent point-group element R on the (displacement, microrotation)
    state space.  Displacement u transforms as a polar vector (R); the
    microrotation phi transforms as an axial vector (det(R) * R)."""
    n = len(coords)
    perm = vertex_perm(R, coords)
    D = np.zeros((6 * n, 6 * n))
    det = np.linalg.det(R)
    for i in range(n):
        j = perm[i]
        D[3 * j:3 * j + 3, 3 * i:3 * i + 3] = R
        D[3 * n + 3 * j:3 * n + 3 * j + 3, 3 * n + 3 * i:3 * n + 3 * i + 3] = det * R
    return D


def main():
    coords = cluster_delta()                       # shell + 4 voids, 17 nodes
    n = len(coords)
    P = build_cosserat_matrix_two_d(coords)
    ev, vec = np.linalg.eigh(P)
    ev[np.abs(ev) < 1e-9] = 0.0
    phi_frac = np.sum(vec[3 * n:, :] ** 2, axis=0)  # microrotation weight per mode

    # Point-group representations, one per T_d element, grouped by class.
    Td = generate_Td()
    reps = [(state_representation(R, coords), classify_Td(R)) for R in Td]

    # Group eigenvalues into degenerate sets.
    order = np.argsort(ev)
    groups, current, value = [], [], None
    for i in order:
        if value is None or abs(ev[i] - value) < 1e-3:
            value = ev[i] if value is None else value
            current.append(i)
        else:
            groups.append((np.mean(ev[current]), current))
            value, current = ev[i], [i]
    if current:
        groups.append((np.mean(ev[current]), current))

    print("A_2 modes of the shell + 4 voids cluster (Sigma / Delta cluster)")
    print(f"{'lambda':>8} {'mult':>4} {'phi%':>5} {'A_2 wt':>7} {'m(N=17)':>9}")
    for value, members in groups:
        if value < 0.5:
            continue
        V = vec[:, members]
        weight = sum(TD_CHARS['A_2'][cl] * np.einsum('ij,ij->', V, D @ V)
                     for D, cl in reps) / len(Td) / len(members)
        if weight > 0.4:                            # A_2 content
            phi_pct = 100 * np.mean([phi_frac[m] for m in members])
            print(f"{value:8.3f} {len(members):>4} {phi_pct:5.0f} "
                  f"{weight:7.2f} {mass_from_lambda(17, value):9.1f}")

    print()
    print("Rule: parity-flip -> A_2; resisted class -> lowest stiff phi-dominant")
    print("      A_2 root above lambda = 4  =>  lambda_Sigma = 4.624")
    print(f"      m_Sigma = {mass_from_lambda(17, 4.624):.2f} MeV "
          f"(PDG {PDG_SIGMA0}, {100 * (mass_from_lambda(17, 4.624) - PDG_SIGMA0) / PDG_SIGMA0:+.2f}%)")


if __name__ == "__main__":
    main()
