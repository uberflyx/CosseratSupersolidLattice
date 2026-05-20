#!/usr/bin/env python3
"""
delta1600_irrep_check.py
========================

The Delta(1232) ground state's parent on the bare shell is T_1u
(stiff phi-dom at lambda=8.054), per the parity-flip rule confirmed by
adiabatic continuation.  Its radial excitation Delta(1600) on the
dual-orbit cluster (N=21, O_h restored) should therefore be the
T_1u descendant, not T_1g.

Verify: project the Delta(1600) mode (at lambda = 10.11, descended
from the parent T_1 mode at lambda = 9.05 on the N=17 cluster) onto
T_1g and T_1u and confirm which irrep it lives in.
"""

import numpy as np
import scipy.linalg
import sys

from delta1600_dual_orbit import (
    build_dual_orbit_cluster, build_oh_elements,
    build_cosserat_oh_action, project_onto_irrep, O_H_IRREPS
)
from delta_first_principles import build_cosserat_matrix_two_d
from hadron_spectral_mass import mass_from_lambda


def main():
    print("=" * 78)
    print("Delta(1600) irrep check: T_1g or T_1u?")
    print("=" * 78)

    coords = build_dual_orbit_cluster()
    N = len(coords)
    H = build_cosserat_matrix_two_d(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    vals, vecs = scipy.linalg.eigh(H)

    oh_elements, class_ids = build_oh_elements()
    oh_actions = build_cosserat_oh_action(coords, oh_elements)

    irrep_names = ['A_1g', 'A_2g', 'E_g', 'T_1g', 'T_2g',
                   'A_1u', 'A_2u', 'E_u', 'T_1u', 'T_2u']
    projectors = {}
    for name in irrep_names:
        projectors[name] = project_onto_irrep(
            vecs, oh_actions, class_ids, O_H_IRREPS[name])

    # Look at the spectrum near lambda = 10.11 (the Delta(1600) candidate
    # from my adiabatic trace).  Also list all T_1u and T_1g triplets.
    print(f"\nAll T_1g and T_1u triplets in the spectrum:")
    print(f"{'lambda':>8} {'m (MeV)':>9} {'irrep':>6} {'phi%':>5} "
          f"{'shell%':>7}")
    seen = set()
    for k in range(len(vals)):
        if vals[k] in seen:
            continue
        # Find degenerate group
        group = [k]
        for j in range(k + 1, len(vals)):
            if abs(vals[j] - vals[k]) < 1e-4:
                group.append(j)
        # Check irrep
        v = vecs[:, group[0]]
        max_irrep = None
        max_content = 0
        for name in irrep_names:
            Pv = projectors[name] @ v
            content = np.dot(v, Pv) / np.dot(v, v)
            if content > max_content:
                max_content = content
                max_irrep = name
        if max_irrep in ('T_1g', 'T_1u') and len(group) == 3:
            # Phi content
            u_tot = np.sum(v[:3*N]**2)
            phi_tot = np.sum(v[3*N:]**2)
            phi_frac = phi_tot / (u_tot + phi_tot)
            # Shell content
            amps = np.array([
                np.sum(v[3*i:3*(i+1)]**2) +
                np.sum(v[3*N+3*i:3*N+3*(i+1)]**2)
                for i in range(N)
            ])
            shell_frac = amps[:13].sum() / amps.sum()
            m = mass_from_lambda(N, vals[k])
            tag = ""
            if 10.0 < vals[k] < 10.2:
                tag = "  <-- Delta(1600) candidate range"
            print(f"{vals[k]:>8.4f} {m:>9.2f} {max_irrep:>6} "
                  f"{phi_frac*100:>4.1f}  {shell_frac*100:>5.1f}%{tag}")
        for j in group:
            seen.add(vals[j])


if __name__ == '__main__':
    main()
