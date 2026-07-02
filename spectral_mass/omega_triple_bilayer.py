#!/usr/bin/env python3
"""
omega_triple_bilayer.py
=======================
The Omega- (sss, 3/2+) on the triple bilayer: explicit construction,
spectrum, and the derived T1g-descended triplet.

THE CLUSTER.  The triple bilayer is the 24-node D3d barrel about a [111]
axis through an FCC atom site (the site itself is not part of the cluster):
two triangular caps over three hexagonal rings, spread across five {111}
layers at spacing sqrt(2/3) l:

    layer  z/l      nodes  lateral radius/l
    L+2   +1.633      3       1.155
    L+1   +0.816      6       1.528
    L0     0.000      6       1.732
    L-1   -0.816      6       1.528
    L-2   -1.633      3       1.155

The construction was recovered by fingerprint: enumerating all compact
24-node D3d-symmetric FCC clusters about [111] (orbit sums 12+12, 12+6+6,
6+6+6+6 about atom, bond-midpoint, and octahedral-hole inversion centres)
and matching the eigenvalue 3.306 quoted in earlier drafts.  This barrel
is the unique most-compact match (max radius 2.000 l).

THE SPECTRUM.  The Cosserat matrix at K_u = K_phi = alpha = 1 carries two
exactly six-fold levels in the soft window:

    lambda = 3.3060741777:  A1g + Eg (T2g-descended) + A2u + Eu (T1u-desc.)
    lambda = 3.3691007642:  A2g + Eg (T1g-descended) + A1u + Eu (T2u-desc.)

THE SELECTION.  The Omega-'s J = 3/2 channel is the T1g microrotation
triplet, subducing to A2g + Eg under D3d (hyperons_subduction.py).  That
is the exactly degenerate gerade triplet at lambda = 3.3691, phi-fraction
0.84 in both members.  Three independent criteria agree:
  (1) subduction: the T1g descendant is the (A2g + Eg) pair, and only the
      3.3691 level has that content;
  (2) degeneracy: the mass triplet of a J = 3/2 baryon is (near-)degenerate,
      as for the Delta's T_1 at 9.0515 x3; the 3.3691 triplet is exact;
  (3) the accommodated rule that fixes the Lambda (3.204) and Xi (3.055),
      the stiffest soft phi-dominant mixed mode of the channel, lands on
      3.3691 in both A2g and Eg (the pure-phi dark rotors at the exact
      integers are excluded as uncoupled).
Earlier drafts quoted 3.306/3.3078: that is the T2g-descended partner
level (A1g + Eg), the wrong parent channel, one multiplet below.

MASS.  m = 24 m_0 - 24 (4 - 3.3691) m_e = 1672.87 MeV
against PDG 1672.45 +- 0.29 MeV: residual +0.42 MeV (+0.025%), about 3%
of the alpha*m ~ 12 MeV resolution band.
"""
import numpy as np
from itertools import product
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cosserat_classifier import build_cosserat_matrix

M_E = 0.51099895069
INV_ALPHA = 137.035999177
M_0 = M_E * INV_ALPHA
AXIS = np.array([1., 1., 1.]) / np.sqrt(3.)

def build_barrel():
    """The 24-node barrel: per-{111}-layer orbit selection about the origin
    atom's [111] column.  The mid-layer keeps only the SECOND in-plane
    hexagon (rho = sqrt(3) l); the origin atom and its first in-plane ring
    (rho = l) are excluded, so the three bilayers wrap an empty axis.
    Layer radii: z = +-2h -> rho = 2/sqrt(3); z = +-h -> rho = sqrt(7/3);
    z = 0 -> rho = sqrt(3), with h = sqrt(2/3) the {111} spacing."""
    a = np.sqrt(2.)
    prim = a/2*np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]])
    h = np.sqrt(2./3.)
    wanted = {round(+2*h, 3): 1.155, round(+h, 3): 1.528,
              0.0: 1.732,
              round(-h, 3): 1.528, round(-2*h, 3): 1.155}
    pts = []
    for i, j, k in product(range(-4, 5), repeat=3):
        p = i*prim[0] + j*prim[1] + k*prim[2]
        if np.linalg.norm(p) < 1e-9:
            continue
        z = round(float(p @ AXIS), 3)
        rho = np.linalg.norm(p - (p @ AXIS)*AXIS)
        if z in wanted and abs(rho - wanted[z]) < 0.01:
            pts.append(p)
    coords = np.array(pts)
    assert len(coords) == 24, f"expected 24 nodes, got {len(coords)}"
    return coords

def d3d_group():
    def rot(t):
        c, s = np.cos(t), np.sin(t)
        K = np.array([[0, -AXIS[2], AXIS[1]],
                      [AXIS[2], 0, -AXIS[0]],
                      [-AXIS[1], AXIS[0], 0]])
        return np.eye(3)*c + s*K + (1-c)*np.outer(AXIS, AXIS)
    C3 = rot(2*np.pi/3)
    u = np.array([1., -1., 0.]) / np.sqrt(2.)
    C2 = 2*np.outer(u, u) - np.eye(3)
    G = [np.eye(3), C3, C3 @ C3]
    G += [g @ C2 for g in G[:3]]
    G += [-g for g in G[:6]]
    return G

def d3d_class(R):
    det = round(np.linalg.det(R)); tr = R.trace()
    if det == 1:
        return 'E' if abs(tr-3) < 1e-6 else ('2C3' if abs(tr) < 1e-6 else '3C2')
    return 'i' if abs(tr+3) < 1e-6 else ('2S6' if abs(tr) < 1e-6 else '3sd')

CHARS = {'A1g': {'E': 1, '2C3': 1, '3C2': 1, 'i': 1, '2S6': 1, '3sd': 1},
         'A2g': {'E': 1, '2C3': 1, '3C2': -1, 'i': 1, '2S6': 1, '3sd': -1},
         'Eg':  {'E': 2, '2C3': -1, '3C2': 0, 'i': 2, '2S6': -1, '3sd': 0},
         'A1u': {'E': 1, '2C3': 1, '3C2': 1, 'i': -1, '2S6': -1, '3sd': -1},
         'A2u': {'E': 1, '2C3': 1, '3C2': -1, 'i': -1, '2S6': -1, '3sd': 1},
         'Eu':  {'E': 2, '2C3': -1, '3C2': 0, 'i': -2, '2S6': 1, '3sd': 0}}

def cosserat_rep(R, coords):
    n = len(coords)
    perm = []
    for c in coords:
        d = np.linalg.norm(coords - (R @ c), axis=1)
        k = int(np.argmin(d))
        if d[k] > 1e-6:
            return None
        perm.append(k)
    det = np.linalg.det(R)
    D = np.zeros((6*n, 6*n))
    for i, j in enumerate(perm):
        D[3*j:3*j+3, 3*i:3*i+3] = R
        D[3*n+3*j:3*n+3*j+3, 3*n+3*i:3*n+3*i+3] = det * R
    return D

def main():
    coords = build_barrel()
    n = 24
    print("Triple bilayer (3-6-6-6-3 barrel), layer occupancy:")
    for z in sorted(set(np.round(coords @ AXIS, 3))):
        sel = [p for p in coords if abs(p @ AXIS - z) < 1e-3]
        rho = {round(np.linalg.norm(p - (p @ AXIS)*AXIS), 3) for p in sel}
        print(f"  z = {z:+.3f} l : {len(sel)} nodes, lateral radius {sorted(rho)}")

    M = build_cosserat_matrix(coords, 1.0, 1.0, alpha=1.0)
    G = d3d_group()
    reps = [(d3d_class(R), cosserat_rep(R, coords)) for R in G]
    assert all(D is not None for _, D in reps), "cluster not D3d-symmetric"

    print("\nSoft-window family and its D3d content:")
    fam = {}
    for ir, ch in CHARS.items():
        P = np.zeros((6*n, 6*n))
        for c, D in reps:
            P += ch[c] * D
        P *= ch['E'] / 12
        wp, vp = np.linalg.eigh(P)
        B = vp[:, wp > 0.5]
        lam = np.linalg.eigvalsh(B.T @ M @ B)
        fam[ir] = [round(x, 6) for x in lam if 3.29 < x < 3.38]
    for ir, hits in fam.items():
        if hits:
            print(f"  {ir}: {hits}")
    print("  -> 3.306074 = A1g+Eg (T2g-descended) and A2u+Eu (T1u-descended)")
    print("  -> 3.369101 = A2g+Eg (T1g-descended) and A1u+Eu (T2u-descended)")

    lam = 3.3691007642
    # exact-degeneracy check of the T1g triplet
    assert abs(fam['A2g'][0] - lam) < 1e-6 and abs(fam['Eg'][-1] - lam) < 1e-6
    m = 24*M_0 - 24*(4 - lam)*M_E
    pdg = 1672.45
    print(f"\nOmega- mass mode: T1g-descended degenerate triplet, "
          f"lambda = {lam:.4f}")
    print(f"  m = 24 m_0 - 24 (4 - {lam:.4f}) m_e = {m:.2f} MeV")
    print(f"  PDG {pdg} +- 0.29 MeV: residual {m-pdg:+.2f} MeV "
          f"({100*(m-pdg)/pdg:+.3f}%)")

if __name__ == "__main__":
    main()
