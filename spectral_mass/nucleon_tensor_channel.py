#!/usr/bin/env python3
"""
nucleon_tensor_channel.py
=========================
Full O_h-irrep-resolved Cosserat spectrum of the nucleon's 13-node
coordination shell, at the same alpha = 1 that fixes the proton's
A_2u rest-mass mode at lambda = 8.3028.

Why this exists
---------------
The proton-neutron splitting is currently carried by a scalar character
cost (5 m_e, odd in I_3) minus a monopole Coulomb. That lands at 1.70 MeV,
about 0.80 m_e above the observed 1.2933 MeV. The Coulomb route to the
missing piece is exhausted (it is even in I_3, and the explicit bond-count
attempts are underdetermined fits).

This script tests a different channel. The charge break that distinguishes
the proton (uud) from the neutron (udd) transforms as the photon mode T_1u.
By the O_h selection rule

        A_2u  (x)  T_1u  =  T_2g,

that break mixes the rest-mass A_2u mode into the tensor modes T_2g. T_2g
(with E_g) is the couple-stress / microrotation tensor channel: the SAME
channel that carries the graviton and the f_2(1270) tensor meson in this
framework (sec:massive_tensor). So the question is whether the missing
piece of the nucleon splitting lives in the gravity/tensor sector rather
than in electromagnetism.

What is computed (no fitted parameters)
---------------------------------------
  1. The Cosserat eigenvalue of every O_h irrep on the coordination shell,
     at alpha = 1, with displacement (u) vs microrotation (phi) content.
  2. The selection-rule check A_2u (x) T_1u = T_2g by character arithmetic.
  3. The reference offsets (lambda - 4) and (lambda - 8.3028) for the tensor
     modes, which is what the modal mass sum weighs.

Place in the budget (see nucleon_magnetic_selfenergy.py)
-------------------------------------------------------
This channel is the spin-INDEPENDENT half of the couple-stress sector: the
static charge break feeding a microrotation pattern, worth about -0.075 MeV.
Its spin-DEPENDENT partner is the nucleon magnetic self-energy (the moment's
own field energy, computed in nucleon_magnetic_selfenergy.py), worth about
-0.11 MeV. The two together (~-0.18 MeV) are the microrotation-sector share;
the remaining residual (~-0.23 MeV) is the three-body charge overlap plus the
inelastic Compton part. So the tensor channel is one term of the budget, not
the whole of the missing piece.

Accompanies: M. A. Cox, "The Cosserat Supersolid" (2026).

Author: Mitchell A. Cox, with Claude (Anthropic).  License: MIT.
"""

import numpy as np
from cosserat_classifier import build_cosserat_matrix
from spectral_classifier import fcc_nn_vectors
from proton_first_principles import build_oh_elements, class_of, find_perm, build_rep

# ---------------------------------------------------------------------------
# Full O_h character table, keyed by the class names that class_of() returns.
# Classes: E, 8C_3, 6C_4, 6C_2, 3C_2, i, 8S_6, 6S_4, 6sigma_d, 3sigma_h
# ---------------------------------------------------------------------------
CLASS_KEYS = ['E', '8C_3', '6C_4', '6C_2', '3C_2', 'i', '8S_6', '6S_4', '6σ_d', '3σ_h']

CHARS = {
    'A_1g': [ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1],
    'A_2g': [ 1,  1, -1, -1,  1,  1,  1, -1,  1, -1],
    'E_g':  [ 2, -1,  0,  0,  2,  2, -1,  0,  0,  2],
    'T_1g': [ 3,  0,  1, -1, -1,  3,  0,  1, -1, -1],
    'T_2g': [ 3,  0, -1,  1, -1,  3,  0, -1,  1, -1],
    'A_1u': [ 1,  1,  1,  1,  1, -1, -1, -1, -1, -1],
    'A_2u': [ 1,  1, -1, -1,  1, -1, -1,  1,  1, -1],
    'E_u':  [ 2, -1,  0,  0,  2, -2,  1,  0,  0, -2],
    'T_1u': [ 3,  0,  1, -1, -1, -3,  0, -1,  1,  1],
    'T_2u': [ 3,  0, -1,  1, -1, -3,  0,  1, -1,  1],
}
DIM = {'A_1g': 1, 'A_2g': 1, 'E_g': 2, 'T_1g': 3, 'T_2g': 3,
       'A_1u': 1, 'A_2u': 1, 'E_u': 2, 'T_1u': 3, 'T_2u': 3}
CLASS_SIZE = {'E': 1, '8C_3': 8, '6C_4': 6, '6C_2': 6, '3C_2': 3,
              'i': 1, '8S_6': 8, '6S_4': 6, '6σ_d': 6, '3σ_h': 3}


def chars_dict(name):
    return {k: c for k, c in zip(CLASS_KEYS, CHARS[name])}


def build_projector(coords, oh_elements, name):
    """P_irrep = (dim/|G|) sum_g chi(g)* D(g)."""
    n = len(coords)
    chars = chars_dict(name)
    P = np.zeros((6 * n, 6 * n))
    for R in oh_elements:
        p = find_perm(R, coords)
        if p is None:
            return None
        P += chars[class_of(R)] * build_rep(R, p, n)
    return (DIM[name] / len(oh_elements)) * P


def irrep_modes(coords, oh_elements, name, alpha):
    """Eigenvalues of the Cosserat matrix lying in the given irrep, with
    u/phi content. Eigenvalues are reported once per distinct value (the
    irrep multiplet is degenerate)."""
    n = len(coords)
    P = build_projector(coords, oh_elements, name)
    if P is None:
        return None
    M = build_cosserat_matrix(coords, K_u=1.0, K_phi=1.0, alpha=alpha)
    ev, vec = np.linalg.eigh(M)
    out = []
    i = 0
    while i < len(ev):
        j = i
        while j < len(ev) and abs(ev[j] - ev[i]) < 1e-6:
            j += 1
        V = vec[:, i:j]
        S = V.T @ P @ V
        w = np.linalg.eigvalsh(S)
        weight = float(np.sum(w[w > 0.5]))  # how much of this eigenspace is in the irrep
        if weight > 0.5:
            # representative eigenvector for u/phi content
            s_eval, s_evec = np.linalg.eigh(S)
            mode = V @ s_evec[:, -1]
            mode /= np.linalg.norm(mode)
            u_c = float(np.sum(mode[:3 * n] ** 2))
            p_c = float(np.sum(mode[3 * n:] ** 2))
            mult = int(round(weight))
            out.append((float(ev[i]), u_c, p_c, mult))
        i = j
    return out


def selection_rule_check():
    """Decompose A_2u (x) T_1u into O_h irreps by character arithmetic."""
    prod = {k: chars_dict('A_2u')[k] * chars_dict('T_1u')[k] for k in CLASS_KEYS}
    print("  A_2u (x) T_1u, decomposed over O_h:")
    contains = []
    for name in CHARS:
        chars = chars_dict(name)
        n_copies = sum(CLASS_SIZE[k] * prod[k] * chars[k] for k in CLASS_KEYS) / 48
        n_copies = round(n_copies)
        if n_copies != 0:
            contains.append((name, n_copies))
            print(f"     {name}: x{n_copies}")
    return contains


def main():
    co = fcc_nn_vectors()
    centre = np.array([[0., 0., 0.]])
    coords = np.vstack([centre, co])  # 13-node coordination shell
    oh = build_oh_elements()

    M_E = 0.51099895069   # MeV
    LAM_REF = 4.0         # universal A_2u displacement reference
    LAM_A2U = 8.3028      # proton rest-mass mode (phi-derived A_2u)

    print("=" * 76)
    print("Selection rule for the proton-neutron charge break")
    print("=" * 76)
    selection_rule_check()
    print("\n  The charge break (uud -> udd) is a T_1u perturbation, so it feeds")
    print("  the T_2g tensor channel: the couple-stress / graviton / f_2 sector.")

    print("\n" + "=" * 76)
    print("Full O_h-irrep Cosserat spectrum, coordination shell, alpha = 1")
    print("=" * 76)
    order = ['A_1g', 'A_2g', 'E_g', 'T_1g', 'T_2g',
             'A_1u', 'A_2u', 'E_u', 'T_1u', 'T_2u']
    print(f"\n{'irrep':<6s} {'lambda':>9s}  {'u%':>5s} {'phi%':>5s}  {'mult':>4s}  {'lambda-4':>9s}  notes")
    print("-" * 76)
    tensor_eigs = {}
    for name in order:
        modes = irrep_modes(coords, oh, name, 1.0)
        if not modes:
            continue
        for lam, u_c, p_c, mult in sorted(modes):
            note = ""
            if name == 'A_2u' and abs(lam - LAM_A2U) < 1e-2:
                note = "proton rest mass"
            if name in ('T_2g', 'E_g'):
                note = (note + " " if note else "") + "TENSOR (couple-stress)"
            print(f"{name:<6s} {lam:>9.4f}  {u_c*100:>5.1f} {p_c*100:>5.1f}  "
                  f"{mult:>4d}  {lam-LAM_REF:>9.4f}  {note}")
        if name in ('T_2g', 'E_g'):
            tensor_eigs[name] = sorted(modes)

    print("\n" + "=" * 76)
    print("Tensor-channel offsets relevant to the modal mass sum")
    print("=" * 76)
    print("  The modal mass sum weighs (lambda - 4); the splitting is the change")
    print("  in tensor weight between the two nucleons, routed by A_2u (x) T_1u.")
    print("  Reference values: lambda_ref = 4, proton A_2u = 8.3028.\n")
    for name in ('T_2g', 'E_g'):
        for lam, u_c, p_c, mult in tensor_eigs.get(name, []):
            d4 = lam - LAM_REF
            dA = lam - LAM_A2U
            print(f"  {name} at lambda = {lam:.4f}: "
                  f"(lambda-4) = {d4:+.4f},  (lambda-8.3028) = {dA:+.4f},  "
                  f"u={u_c*100:.0f}% phi={p_c*100:.0f}%")

    print("\n  Target for the missing piece of the splitting:")
    need = 1.2933324 - (5 - 5.0/3.0) * M_E  # observed minus (strong - monopole Coulomb)
    print(f"     observed splitting 1.2933 MeV minus cell-pair (5 m_e - 5/3 m_e)")
    print(f"     = {need:+.4f} MeV = {need/M_E:+.4f} m_e   (additional reduction wanted)")


if __name__ == '__main__':
    main()
