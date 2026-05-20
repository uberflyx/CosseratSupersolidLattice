#!/usr/bin/env python3
"""
bare_shell_irreps.py
====================

Diagonalise the bare 13-shell Cosserat matrix and decompose every mode
by O_h irrep.  Identify the irrep content at lambda ~ 6.16, which is
the parent of the chapter's N(1535) (lambda = 8.046 on bilayer).

Also catalogue: where do the chapter's other "stiff phi-dom" modes
come from, irrep by irrep?
"""

import numpy as np
import scipy.linalg
import sys

from spectral_classifier import cluster_coord_shell
from cosserat_classifier import build_cosserat_matrix
from delta1600_dual_orbit import (
    build_oh_elements, build_cosserat_oh_action,
    project_onto_irrep, O_H_IRREPS, CLASS_ORDERS, G_ORDER
)


def main():
    print("=" * 78)
    print("Bare 13-shell Cosserat spectrum, decomposed by O_h irrep")
    print("=" * 78)

    shell, _ = cluster_coord_shell()
    N = len(shell)
    print(f"Bare shell: N = {N}")

    H = build_cosserat_matrix(shell, K_u=1.0, K_phi=1.0, alpha=1.0)
    vals, vecs = scipy.linalg.eigh(H)
    print(f"Cosserat spectrum: {len(vals)} eigenvalues, "
          f"range [{vals[0]:.4f}, {vals[-1]:.4f}]")

    # Build O_h elements and Cosserat action
    oh_elements, class_ids = build_oh_elements()
    oh_actions = build_cosserat_oh_action(shell, oh_elements)

    # Build projectors
    irrep_names = ['A_1g', 'A_2g', 'E_g', 'T_1g', 'T_2g',
                   'A_1u', 'A_2u', 'E_u', 'T_1u', 'T_2u']
    projectors = {}
    for name in irrep_names:
        projectors[name] = project_onto_irrep(
            vecs, oh_actions, class_ids, O_H_IRREPS[name])

    # For each eigenvalue, find dominant irrep
    print(f"\n{'k':>4} {'lambda':>8} {'phi%':>5}  {'dominant irrep':>15}  "
          f"{'(content)':>10}")
    for k in range(len(vals)):
        v = vecs[:, k]
        u_tot = np.sum(v[:3*N]**2)
        phi_tot = np.sum(v[3*N:]**2)
        phi_frac = phi_tot / (u_tot + phi_tot)

        # Compute irrep content
        max_irrep = None
        max_content = 0
        for name in irrep_names:
            Pv = projectors[name] @ v
            content = np.dot(v, Pv) / np.dot(v, v)
            if content > max_content:
                max_content = content
                max_irrep = name

        # Show all but skip zero modes
        if abs(vals[k]) > 1e-6:
            print(f"  {k:>3} {vals[k]:>8.4f} {phi_frac*100:>4.1f}  "
                  f"{max_irrep:>15}  {max_content*100:>8.1f}%")

    # Now report all stiff phi-dom modes (lambda > 4, phi > 90%)
    print(f"\n{'-' * 60}")
    print(f"Stiff phi-dominant modes (lambda > 4, phi > 90%):")
    print(f"{'k':>4} {'lambda':>8} {'phi%':>5}  {'irrep':>10}")
    for k in range(len(vals)):
        v = vecs[:, k]
        u_tot = np.sum(v[:3*N]**2)
        phi_tot = np.sum(v[3*N:]**2)
        phi_frac = phi_tot / (u_tot + phi_tot)
        if vals[k] > 4 and phi_frac > 0.90:
            max_irrep = None
            max_content = 0
            for name in irrep_names:
                Pv = projectors[name] @ v
                content = np.dot(v, Pv) / np.dot(v, v)
                if content > max_content:
                    max_content = content
                    max_irrep = name
            print(f"  {k:>3} {vals[k]:>8.4f} {phi_frac*100:>4.1f}  "
                  f"{max_irrep:>10}")


if __name__ == '__main__':
    main()
