#!/usr/bin/env python3
"""
docking_junction_binding.py
===========================
Binding of the coherent two-shell docking junctions by the weighted-mode
(geometric-mean local-consistency) closure, with the parity-sector step.

For two identical coordination shells whose centres differ by an FCC
translation, the junction is fixed by that separation.  Three overlapping
separations are tested here:

    <110> nearest neighbour   -> deuteron (Mode P), 6 shared atoms
    <100>                     -> square {100} face, 4 shared atoms
    <211>                     -> edge, 2 shared atoms

Method (see the mass-mode-identification section of the monograph):
  1. Build the compound Cosserat matrix on the deduplicated cluster.
  2. Classify each eigenmode by parity under inversion through the midpoint
     (displacement is polar, microrotation is axial).
  3. Within the parity-odd sector (the bound two-nucleon channel), rank modes
     by the geometric mean over constituents of the local squared cosine with
     the proton mass mode psi_p, and take the maximum.
  4. Mass from the spectral formula m = N_tot m_0 - N_tot (4 - lambda*) m_e
     with N_tot = N_A + N_B (not deduplicated).

Result: the deuteron binds (lambda* = 8.1995, 1876.5 MeV, ~1.4 MeV below
threshold); the square-face and edge junctions lie above threshold, so the
deuteron is the only bound coherent two-nucleon junction.

Authors: M. Cox, with Claude (Anthropic).  License: MIT.
"""
import numpy as np
from scipy.linalg import eigh
from spectral_classifier import fcc_nn_vectors
from cosserat_classifier import build_cosserat_matrix
from hadron_spectral_mass import get_psi_p, find_constituent_shell_atoms, M_0, M_E

LAM_P, PSI_P = get_psi_p()                       # proton mass mode on the bare 13-shell
NN = fcc_nn_vectors()                            # 12 unit nearest-neighbour vectors
M_P = 13 * M_0 - 13 * (4 - LAM_P) * M_E          # framework proton mass


def dedup(coords, tol=1e-8):
    keep = []
    for a in coords:
        if not any(np.linalg.norm(a - b) < tol for b in keep):
            keep.append(a)
    return np.array(keep)


def two_shell_cluster(separation):
    """Two coordination shells, centres at 0 and `separation`, deduplicated."""
    atoms = [np.zeros(3)] + list(NN) + [separation] + [separation + v for v in NN]
    return dedup(np.array(atoms))


def local_restriction(psi_k, coords, centre):
    """78-component restriction of compound mode psi_k to the shell at `centre`."""
    n = len(coords)
    idx = find_constituent_shell_atoms(coords, centre)
    out = np.zeros(78)
    for k, i in enumerate(idx):
        out[3 * k:3 * k + 3] = psi_k[3 * i:3 * i + 3]                  # displacement
        out[39 + 3 * k:39 + 3 * k + 3] = psi_k[3 * n + 3 * i:3 * n + 3 * i + 3]  # microrotation
    return out


def squared_cosine(a, b):
    na, nb = a @ a, b @ b
    return 0.0 if na < 1e-14 or nb < 1e-14 else (a @ b) ** 2 / (na * nb)


def inversion_operator(coords, centre):
    """Inversion through `centre`: displacement flips (polar), microrotation does not (axial)."""
    n = len(coords)
    perm = np.zeros(n, int)
    for i, r in enumerate(coords):
        ri = 2 * centre - r
        j = int(np.argmin([np.linalg.norm(ri - s) for s in coords]))
        if np.linalg.norm(ri - coords[j]) > 1e-6:
            return None                          # cluster not centrosymmetric
        perm[i] = j
    P = np.zeros((6 * n, 6 * n))
    for i in range(n):
        j = perm[i]
        for a in range(3):
            P[3 * i + a, 3 * j + a] = -1.0
            P[3 * n + 3 * i + a, 3 * n + 3 * j + a] = +1.0
    return P


def close_junction(separation, label):
    coords = two_shell_cluster(separation)
    centres = [np.zeros(3), separation]
    n_tot = 13 * len(centres)
    H = build_cosserat_matrix(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    vals, vecs = eigh(H)
    P = inversion_operator(coords, separation / 2.0)

    odd = []
    for k in range(len(vals)):
        if vals[k] < 1.0:                        # drop the six rigid modes
            continue
        gm = np.prod([squared_cosine(PSI_P, local_restriction(vecs[:, k], coords, c))
                      for c in centres]) ** (1.0 / len(centres))
        parity = vecs[:, k] @ (P @ vecs[:, k]) if P is not None else 0.0
        if parity < -0.5:                        # parity-odd: the bound two-nucleon channel
            odd.append((gm, vals[k]))
    odd.sort(reverse=True)
    lam = odd[0][1]
    m = n_tot * M_0 - n_tot * (4 - lam) * M_E
    binding = 2 * M_P - m
    print(f"{label:28s} shared {13*len(centres)-len(coords):2d}  "
          f"lambda* = {lam:7.4f}  m = {m:8.2f} MeV  "
          f"binding = {binding:+6.2f} MeV  [{'BOUND' if binding > 0 else 'unbound'}]")


if __name__ == "__main__":
    print(f"proton reference: lambda_p = {LAM_P:.4f}, m_p = {M_P:.2f} MeV "
          f"(threshold 2 m_p = {2*M_P:.2f} MeV)\n")
    close_junction(NN[0], "<110> NN (deuteron, Mode P)")     # validation
    close_junction(np.array([np.sqrt(2), 0., 0.]), "<100> square {100} face")
    close_junction(NN[0] + NN[1], "<211> edge")
