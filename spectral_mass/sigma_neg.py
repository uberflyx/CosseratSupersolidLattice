#!/usr/bin/env python3
"""
sigma_neg.py
============
Show that the negative-parity Sigma resonances cannot be three-quark
clusters, by reading the gerade modes of the dual-orbit cluster.

A negative-parity baryon reads a gerade (inversion-even) spatial mode,
which only an inversion-symmetric cluster supplies.  The candidate is the
dual-orbit void cluster (coordination shell plus both T_d void orbits,
21 nodes, full O_h symmetry), the same cluster that carries the Delta(1600).

This script projects the spectrum onto the O_h irreps and reports the
gerade modes that the negative-parity Sigma would read: A_2g for the
J = 1/2 states and T_1g / T_2g for the J = 3/2 states.  The lowest stiff
phi-dominant gerade modes land near 1500 MeV:

    Sigma(1620) 1/2-  ->  lowest A_2g  ~ 1513 MeV   (PDG 1620, -6.6%)
    Sigma(1670) 3/2-  ->  lowest T_1g  ~ 1497 MeV   (PDG 1670, -10.3%)

Both fall well short.  The dual-orbit cluster tops out near 1540 MeV, while
the parity-flipping bilayer floor sits at 1699.5 MeV (see sigma_desert.py).
The observed states lie in the gap between the two, so no three-quark
assignment is available and they must be meson-baryon molecules.  This is
the negative result that motivates the molecular closure of the Sigma(1670)
as a Kbar-Delta bound state.

Matrix convention: unit stiffness on every nearest-neighbour and
void-corner bond, length entering only through the curl, weighted 1/(2d).

Run:  python3 sigma_neg.py
"""

import numpy as np
from delta1600_dual_orbit import (build_dual_orbit_cluster, build_oh_elements,
                                   build_cosserat_oh_action, project_onto_irrep,
                                   O_H_IRREPS)
from delta_first_principles import build_cosserat_matrix_two_d
from hadron_spectral_mass import mass_from_lambda

N_DUAL = 21
GERADE = ('A_2g', 'T_1g', 'T_2g')
PDG = {'A_2g': ('Sigma(1620) 1/2-', 1620),
       'T_1g': ('Sigma(1670) 3/2- (T_1g)', 1670),
       'T_2g': ('Sigma(1670) 3/2- (T_2g)', 1670)}


def main():
    coords = build_dual_orbit_cluster()              # shell + 8 voids, 21 nodes
    n = len(coords)
    P = build_cosserat_matrix_two_d(coords)
    ev, vec = np.linalg.eigh(P)
    ev[np.abs(ev) < 1e-9] = 0.0
    phi = np.sum(vec[3 * n:, :] ** 2, axis=0)

    elements, class_ids = build_oh_elements()
    actions = build_cosserat_oh_action(coords, elements)
    projectors = {name: project_onto_irrep(None, actions, class_ids, chars)
                  for name, chars in O_H_IRREPS.items()}

    # Group eigenvalues into degenerate sets.
    order = np.argsort(ev)
    groups, current, value = [], [], None
    for i in order:
        if value is None or abs(ev[i] - value) < 5e-3:
            value = ev[i] if value is None else value
            current.append(i)
        else:
            groups.append((np.mean(ev[current]), current))
            value, current = ev[i], [i]
    if current:
        groups.append((np.mean(ev[current]), current))

    def irrep_of(members):
        V = vec[:, members]
        best, best_w = None, -1.0
        for name, Pr in projectors.items():
            w = np.sum((Pr @ V) ** 2) / np.sum(V ** 2)
            if w > best_w:
                best, best_w = name, w
        return best

    print("Gerade modes of the dual-orbit cluster (N = 21, O_h)")
    print(f"{'lambda':>8} {'mult':>4} {'irrep':>6} {'phi%':>5} {'m(N=21)':>9}")
    found = {name: [] for name in GERADE}
    for value, members in groups:
        if value < 0.5:
            continue
        name = irrep_of(members)
        phi_pct = 100 * np.mean([phi[m] for m in members])
        if name in GERADE and phi_pct > 55 and value > 4:
            print(f"{value:8.3f} {len(members):>4} {name:>6} {phi_pct:5.0f} "
                  f"{mass_from_lambda(N_DUAL, value):9.1f}")
            found[name].append((value, mass_from_lambda(N_DUAL, value)))

    print("\nResisted rule = lowest stiff phi-dominant mode of each gerade irrep:")
    for name in GERADE:
        label, pdg = PDG[name]
        if found[name]:
            lam, m = min(found[name], key=lambda x: x[0])
            print(f"  {label:24} lambda = {lam:.3f}  m = {m:.1f}  "
                  f"PDG {pdg} ({100 * (m - pdg) / pdg:+.2f}%)")
        else:
            print(f"  {label:24} no stiff phi-dominant mode found")
    print("\nAll gerade modes fall short: the negative-parity Sigma are not")
    print("dual-orbit quark clusters.  See sigma_desert.py for the resolution.")


if __name__ == "__main__":
    main()
