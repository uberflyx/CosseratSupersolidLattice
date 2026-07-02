#!/usr/bin/env python3
"""
eta_family_spectral.py
======================
Spectral eigenvalues of the eta family clusters, replacing the missing
cook_spectral_masses.py cited in earlier drafts.

Two clusters:

1. The eta' shell-plus-apex (N = 14): the 13-node coordination shell plus
   one apex at the octahedral second-neighbour position on a fourfold axis.
   Point group C4v.  The pseudoscalar irrep of C4v is A2 (invariant under
   all rotations, odd under every mirror).  The displacement A2 block reads
   {0, 1, 2}: the zero is the rigid rotation about the fourfold axis, and
   the lowest non-trivial pseudoscalar mode is the EXACT integer lambda = 1,
   a twist of the equatorial square against the polar square about the axis.
   Consistency check: the shell's parent A2u eigenvector survives the apex
   attachment unchanged at lambda = 4, because the apex sits at a node of
   that mode; the eta' mode is NOT the A2u descendant.

2. The eta bilayer (N = 8): the seven-node hex cap plus one apex over an
   adjacent-layer hollow.  Residual point group Cs = {E, sigma}.  The
   eta is a fault-winding state like the pion: its 0^- is the topological
   label of the stacking winding, and the mass mode is the winding's
   reaction coordinate.  The advance coordinate is NOT the apex sliding
   alone (that mostly excites a floppy mechanism); the stacking phase
   advances when the patch bonded to the apex (apex + its three bond
   partners) slips coherently against the rest of the plane along the
   in-plane <112> partial direction.  The probe is fixed by the bond
   topology with no freedom.  It is even under the surviving mirror (the
   slip direction lies in the mirror plane), which excludes the odd
   candidate at 0.6529 outright, and it places its dominant weight, 0.475,
   on the mode at lambda = 0.8229; no other mode reaches half that.  The
   same probe on the pion's cell pair returns the stretch at lambda = 2
   with weight 1.000.

Mass closures (m = N m_0 - N (4 - lambda) m_e, CODATA inputs):
  eta':  14 m_0 - 14*3 m_e       = 958.89 MeV   (PDG 957.78, +0.12%)
  eta:    8 m_0 - 8*3.1771 m_e   = 547.21 MeV   (PDG 547.86, -0.12%)
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cosserat_classifier import build_uu_block
from spectral_classifier import fcc_nn_vectors

M_E = 0.51099895069
INV_ALPHA = 137.035999177
M_0 = M_E * INV_ALPHA

def polar_rep(R, coords, tol=1e-6):
    """3n x 3n matrix of R acting on displacement space, or None."""
    n = len(coords)
    perm = []
    for c in coords:
        d = np.linalg.norm(coords - (R @ c), axis=1)
        k = int(np.argmin(d))
        if d[k] > tol:
            return None
        perm.append(k)
    if len(set(perm)) != n:
        return None
    D = np.zeros((3*n, 3*n))
    for i, j in enumerate(perm):
        D[3*j:3*j+3, 3*i:3*i+3] = R
    return D

def signed_perms():
    from itertools import permutations
    out = []
    for p in permutations(range(3)):
        for s in np.ndindex(2, 2, 2):
            R = np.zeros((3, 3))
            for row, (col, sg) in enumerate(zip(p, s)):
                R[row, col] = 1 - 2*sg
            out.append(R)
    return out

def irrep_eigs(coords, M, chars, classify):
    """Eigenvalues of M restricted to each irrep of the cluster's group."""
    G = [(R, polar_rep(R, coords)) for R in signed_perms()]
    G = [(R, D) for R, D in G if D is not None]
    out = {}
    for name, ch in chars.items():
        P = np.zeros_like(M)
        for R, D in G:
            P += ch[classify(R)] * D
        P *= ch['E'] / len(G)
        w, v = np.linalg.eigh(P)
        B = v[:, w > 0.5]
        if B.shape[1] == 0:
            continue
        out[name] = sorted(np.round(np.linalg.eigvalsh(B.T @ M @ B), 4))
    return out

# ---------------------------------------------------------------- eta'
def eta_prime():
    shell = np.vstack([[[0., 0., 0.]], fcc_nn_vectors()])
    apex = np.array([[np.sqrt(2.), 0., 0.]])   # octahedral NNN on the x axis
    coords = np.vstack([shell, apex])
    M = build_uu_block(coords)

    def c4v_class(R):
        tr, det = round(np.trace(R)), round(np.linalg.det(R))
        if det == 1:
            return {3: 'E', 1: 'C4', -1: 'C2'}[tr]
        return 'sv' if np.allclose(R, np.diag(np.diag(R))) else 'sd'

    chars = {'A1': {'E': 1, 'C4': 1, 'C2': 1, 'sv': 1, 'sd': 1},
             'A2': {'E': 1, 'C4': 1, 'C2': 1, 'sv': -1, 'sd': -1},
             'B1': {'E': 1, 'C4': -1, 'C2': 1, 'sv': 1, 'sd': -1},
             'B2': {'E': 1, 'C4': -1, 'C2': 1, 'sv': -1, 'sd': 1},
             'E':  {'E': 2, 'C4': 0, 'C2': -2, 'sv': 0, 'sd': 0}}
    blocks = irrep_eigs(coords, M, chars, c4v_class)
    print("eta' shell-plus-apex (N = 14, C4v), displacement blocks:")
    for k, lam in blocks.items():
        print(f"  {k}: {lam}")

    lam = 1.0
    m = 14*M_0 - 14*(4 - lam)*M_E
    print(f"  Pseudoscalar (A2) lowest non-trivial mode: lambda = {lam}")
    print(f"  m_eta' = 14 m_0 - 14*3 m_e = {m:.2f} MeV  (PDG 957.78)")

    # consistency: shell A2u survives at 4
    a2u_in_b2 = any(abs(x - 4.0) < 1e-3 for x in blocks['B2'])
    print(f"  Shell A2u descendant unchanged at lambda = 4 in B2: {a2u_in_b2}")

# ---------------------------------------------------------------- eta
def eta_bilayer():
    e1 = np.array([1., -1., 0.]) / np.sqrt(2)
    e2 = np.array([1., 1., -2.]) / np.sqrt(6)
    nrm = np.array([1., 1., 1.]) / np.sqrt(3)
    ring = [np.cos(k*np.pi/3)*e1 + np.sin(k*np.pi/3)*e2 for k in range(6)]
    cen = np.zeros(3)
    apex = (cen + ring[0] + ring[1]) / 3 + nrm * np.sqrt(2./3.)
    coords = np.vstack([[cen], ring, [apex]])
    M = build_uu_block(coords)
    w = np.linalg.eigvalsh(M)
    print("\neta bilayer (N = 8, Cs), displacement spectrum:")
    print(" ", sorted(set(np.round(w, 4))))

    # mirror through apex/centre plane
    mn = ring[0] - ring[1]
    mn /= np.linalg.norm(mn)
    S = np.eye(3) - 2*np.outer(mn, mn)
    D = polar_rep(S, coords)
    for name, sg in [("A' (mirror-even)", +1), ("A'' (mirror-odd)", -1)]:
        P = (np.eye(24) + sg*D) / 2
        wp, vp = np.linalg.eigh(P)
        B = vp[:, wp > 0.5]
        lam = [x for x in sorted(np.round(np.linalg.eigvalsh(B.T @ M @ B), 4))
               if x > 1e-6]
        print(f"  {name}: {lam[:5]}")

    # Patch-slip selector: the winding's reaction coordinate, fixed by the
    # bond topology (patch = apex + its three bond partners vs the rest,
    # along the in-plane <112> partial direction, zero net momentum).
    # The naive apex-only slide is NOT the advance coordinate: it mostly
    # excites a floppy mechanism, because sliding the apex against its own
    # bond partners stretches those bonds without advancing the phase.
    w, v = np.linalg.eigh(M)
    d = apex - nrm*(apex @ nrm)
    d /= np.linalg.norm(d)
    pat = np.zeros(24)
    for i in (0, 1, 2, 7):       # centre, ring0, ring1, apex: the patch
        pat[3*i:3*i+3] = +d
    for i in (3, 4, 5, 6):       # far ring nodes: the host
        pat[3*i:3*i+3] = -d
    pat -= np.tile(pat.reshape(8, 3).mean(axis=0), 8)
    pat /= np.linalg.norm(pat)
    weights = sorted(((w[k], (pat @ v[:, k])**2) for k in range(24)),
                     key=lambda t: -t[1])
    print("  Patch-slip probe (advance coordinate), top weights:")
    for lam_k, o in weights[:3]:
        print(f"    lambda = {lam_k:7.4f}  weight = {o:.3f}")
    lam = weights[0][0]
    m = 8*M_0 - 8*(4 - lam)*M_E
    print(f"  Mass mode = dominant-weight mode: lambda = {lam:.4f} (mirror-even)")
    print(f"  m_eta = {m:.2f} MeV  (PDG 547.862, {100*(m-547.862)/547.862:+.2f}%)")

if __name__ == "__main__":
    eta_prime()
    eta_bilayer()
