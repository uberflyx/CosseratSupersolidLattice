#!/usr/bin/env python3
"""
gamma3_vertex.py
================

The cubic coupling Gamma^(3) of the FCC bond network, and the decay
amplitude as Gamma^(3) contracted with the spectral mass eigenvectors at
the graph cut.

Derivation (displacement sector, geometric nonlinearity).
A bond p-q stores V(dl) with the elongation, to cubic order in the
relative displacement du = u_q - u_p,

    dl = s + t^2/(2 l) - s t^2/(2 l^2),   s = rhat.du,  t^2 = |du|^2 - s^2.

Expanding V = (1/2) k dl^2 + (1/6) k3 dl^3 to cubic order in du gives the
per-bond cubic energy

    F3 = (k/2l) s t^2 + (1/6) k3 s^3.

Its third derivative w.r.t. du (the bare vertex in relative coords) is

    g3[a,b,c] = (k/2l)( rhat_a Pp_bc + rhat_b Pp_ac + rhat_c Pp_ab )
              + k3 rhat_a rhat_b rhat_c,        Pp = I - rhat rhat  (transverse).

The geometric piece (k/2l ...) needs no material input; the potential
piece (k3 ...) carries the spring anharmonicity.  We start with k3 = 0.

Mapping to nodes: du depends on u_q (+) and u_p (-), so each of the three
tensor indices is p (sign -) or q (sign +); 2^3 sign combinations per bond.

The decay amplitude for i -> f + pi is the vertex contracted with three
mode vectors,

    M = sum  Gamma^(3)_{(n1 a)(n2 b)(n3 c)} Psi_i Psi_f Psi_pi.

The point of this script: the bare overlap <Psi_f|Psi_i> = 0 by symmetry
(different O_h irreps); the vertex supplies the bridging irrep, so M != 0.
That is the synthesis -- the width lives entirely in the vertex, weighted
by the actual eigenvectors.

Author: M. Cox, with Claude (Anthropic).  License: MIT.
"""
# ============================================================================
# SUPERSEDED / EXPLORATORY -- not the authoritative decay path.
# The production engine is decays/cosserat_decay_engine.py, which derives
# g_piNN = N_H = 13 and the decuplet widths (Delta -3.7%, Sigma* -2.4%) from
# graph invariants with NO fitted coupling.  These scripts are session
# explorations kept for the record; their absolute scale is less accurate.
# See decays/README.md.
# ============================================================================

import sys as _sys, os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_os.path.dirname(_HERE))
_sys.path[:0] = [_ROOT, _os.path.join(_ROOT, 'spectral_mass')]
import numpy as np
from spectral_classifier import fcc_nn_vectors
from cosserat_classifier import build_cosserat_matrix
from delta_first_principles import build_cosserat_matrix_two_d
from composite_clusters import cluster_crossed_fault

ELL = np.linalg.norm(fcc_nn_vectors()[0])

# ----------------------------------------------------------------------
def build_gamma3_u(coords, k=1.0, k3=0.0, ell=ELL, tol=1e-2):
    """Rank-3 cubic tensor in the 3n displacement space."""
    n = len(coords)
    G = np.zeros((3 * n, 3 * n, 3 * n))
    I3 = np.eye(3)
    bonds = []
    for p in range(n):
        for q in range(p + 1, n):
            d = np.linalg.norm(coords[q] - coords[p])
            if abs(d - ell) < tol:
                bonds.append((p, q, (coords[q] - coords[p]) / d))
    for (p, q, rhat) in bonds:
        Pp = I3 - np.outer(rhat, rhat)
        # bare vertex in relative coords
        g = np.zeros((3, 3, 3))
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    g[a, b, c] = (k / (2 * ell)) * (
                        rhat[a] * Pp[b, c] + rhat[b] * Pp[a, c] + rhat[c] * Pp[a, b]
                    ) + k3 * rhat[a] * rhat[b] * rhat[c]
        # distribute over node indices with signs (p:-, q:+)
        for s1, n1 in ((-1, p), (1, q)):
            for s2, n2 in ((-1, p), (1, q)):
                for s3, n3 in ((-1, p), (1, q)):
                    sgn = s1 * s2 * s3
                    G[3*n1:3*n1+3, 3*n2:3*n2+3, 3*n3:3*n3+3] += sgn * g
    return G, len(bonds)

# ----------------------------------------------------------------------
def mass_mode(Phi, target, phi_min=0.0):
    w, V = np.linalg.eigh(Phi); n = Phi.shape[0] // 6
    best, bs = None, 1e9
    for kk in range(len(w)):
        if np.sum(V[3*n:, kk]**2) < phi_min:
            continue
        if abs(w[kk] - target) < bs:
            best, bs = kk, abs(w[kk] - target)
    return w[best], V[:, best]

def u_part(psi):
    """displacement (u) half of a 6n Cosserat eigenvector."""
    n = psi.size // 6
    return psi[:3*n]

def contract(G, a, b, c):
    """M = G_{ijk} a_i b_j c_k  for 3n-vectors a,b,c."""
    return np.einsum('ijk,i,j,k->', G, a, b, c)

# ======================================================================
if __name__ == '__main__':
    print("CUBIC VERTEX Gamma^(3): geometric (displacement) nonlinearity, k3 = 0\n")

    # ---- MESON: rho -> pi pi on the crossed fault ----------------------
    Cr, _ = cluster_crossed_fault(); nC = len(Cr)
    Phi_r = build_cosserat_matrix(Cr, 1, 1, 1.0)
    lam_r, psi_r = mass_mode(Phi_r, 4.891)
    G_r, nb_r = build_gamma3_u(Cr)
    ur = u_part(psi_r); ur = ur / np.linalg.norm(ur)

    def bond_stretch_vec(coords, p, q, n):
        v = np.zeros(3*n); rhat = (coords[q]-coords[p]); rhat /= np.linalg.norm(rhat)
        v[3*p:3*p+3] = -rhat/np.sqrt(2); v[3*q:3*q+3] = +rhat/np.sqrt(2)
        return v
    bonds = [(p,q) for p in range(nC) for q in range(p+1,nC)
             if abs(np.linalg.norm(Cr[q]-Cr[p])-ELL) < 1e-2]
    strain = sorted(((abs(bond_stretch_vec(Cr,p,q,nC)@ur), p, q) for p,q in bonds),
                    reverse=True)
    (_,p1,q1),(_,p2,q2) = strain[0], strain[1]
    pi1 = bond_stretch_vec(Cr,p1,q1,nC); pi2 = bond_stretch_vec(Cr,p2,q2,nC)
    s_wave = (pi1+pi2)/np.sqrt(2)        # symmetric  -> S-wave
    p_wave = (pi1-pi2)/np.sqrt(2)        # antisymmetric -> P-wave (rel. momentum)
    M_S = contract(G_r, ur, s_wave, s_wave)
    M_P = contract(G_r, ur, s_wave, p_wave)
    print(f"rho (1-): N={nC} lambda={lam_r:.3f}  u-content={1-np.sum(psi_r[3*nC:]**2):.2f}")
    print(f"   vertex reach ||G.u_rho|| = {np.linalg.norm(np.einsum('ijk,i->jk',G_r,ur)):.3f}")
    print(f"   S-wave amplitude = {M_S:+.5f}   (forbidden -> ~0)")
    print(f"   P-wave amplitude = {M_P:+.5f}   (allowed: the rho IS a P-wave resonance)\n")

    # ---- BARYON: Delta mass mode is microrotation-dominant -------------
    def shell(): return np.vstack([np.zeros((1,3)), fcc_nn_vectors()])
    def voids(nf):
        v = 1/(2*np.sqrt(2))
        return np.array([[+v,+v,+v],[+v,-v,-v],[-v,+v,-v],[-v,-v,+v]])[:nf]
    C = np.vstack([shell(), voids(4)]); n = len(C)
    Phi = build_cosserat_matrix_two_d(C,1,1,1.0)
    lamD, psiD = mass_mode(Phi, 9.05)
    G, _ = build_gamma3_u(C); uD = u_part(psiD); uD /= (np.linalg.norm(uD) or 1)
    print(f"Delta (3/2+): N={n} lambda={lamD:.3f}  u-content={1-np.sum(psiD[3*n:]**2):.3f}")
    print(f"   vertex reach ||G.u_Delta|| = {np.linalg.norm(np.einsum('ijk,i->jk',G,uD)):.4f}")
    print( "   The mass mode is ~99% microrotation, so the DISPLACEMENT vertex is")
    print( "   nearly blind to it.  Baryon widths need the Cosserat cross-coupling")
    print( "   vertex (phi-phi-u, the rotation nonlinearity) -- the next derivation.\n")
    print("Division of labour: mesons (u-rich modes) -> Gamma^(3)_uuu, validated by")
    print("the P-wave selection rule above; baryons (phi-rich) -> Cosserat vertex.")
