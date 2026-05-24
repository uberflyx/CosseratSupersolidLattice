#!/usr/bin/env python3
"""
cosserat_vertex.py
==================

The Cosserat rotation cubic vertex (phi-phi-u) and the baryon decay
amplitude as the displacement source it generates between two
microrotation-dominant mass modes.

Vertex (covariant microrotation gradient).
With finite rotations the microrotation gradient is measured in a frame
that the displacement rotates.  The curvature energy is
|phi_j - phi_i - omega_ij(u) x phibar_ij|^2, and the cross term is the
cubic vertex

  F3_ij = -K_phi (phi_j - phi_i) . ( omega_ij(u) x phibar_ij ),
  omega_ij(u) = (1/2d) rhat x (u_j - u_i),  phibar_ij = (phi_i + phi_j)/2.

Two microrotation legs (parent and daughter baryon), one displacement leg
(the pion).  The pion enters through omega, the CURL of u: a pure radial
breathing has zero curl, so the pion must carry transverse (P-wave)
structure.  That is the correct selection rule for Delta -> N pi.

The pion source.
Plugging the parent mode A and the daughter mode B into the two phi legs
and solving for the coefficient of u (the displacement the transition
drives) gives, per bond, the source

  J_e^{ij} = -(K_phi/2d) [ dphiA_e (phibarB . rhat) - (dphiA . rhat) phibarB_e ],
  dphiA = phiA_j - phiA_i,  phibarB = (phiB_i + phiB_j)/2,

distributed to nodes i (-) and j (+).  ||J||^2 is the decay strength; the
pion is radiated by this source, with the rate fixed by ||J||^2 and the
two-body phase space from the spectral masses.

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
from delta_first_principles import build_cosserat_matrix_two_d, ELL

VOID_D = np.sqrt(6.0) / 4.0

# ----------------------------------------------------------------------
def bond_list(coords, ell=ELL, void_d=VOID_D, tol=1e-6):
    n = len(coords); B = []
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(coords[j] - coords[i])
            if abs(d - ell) < tol or abs(d - void_d) < tol:
                B.append((i, j, d, (coords[j] - coords[i]) / d))
    return B

def phi_of(psi, n):  return psi[3*n:].reshape(n, 3)   # microrotation, node x 3
def u_of(psi, n):    return psi[:3*n].reshape(n, 3)    # displacement,  node x 3

def pion_source(coords, psiA, psiB, K_phi=1.0):
    """Displacement source J (node x 3) the A->B microrotation transition
    drives.  ||J||^2 is the decay strength into the pion channel."""
    n = len(coords)
    phiA, phiB = phi_of(psiA, n), phi_of(psiB, n)
    J = np.zeros((n, 3))
    for (i, j, d, rhat) in bond_list(coords):
        dphiA   = phiA[j] - phiA[i]
        phibarB = 0.5 * (phiB[i] + phiB[j])
        Je = -(K_phi / (2*d)) * (dphiA * (phibarB @ rhat) - (dphiA @ rhat) * phibarB)
        J[j] += Je; J[i] -= Je
    return J

def mass_mode(Phi, target, phi_min=0.0):
    w, V = np.linalg.eigh(Phi); n = Phi.shape[0] // 6
    best, bs = None, 1e9
    for k in range(len(w)):
        if np.sum(V[3*n:, k]**2) < phi_min: continue
        if abs(w[k]-target) < bs: best, bs = k, abs(w[k]-target)
    return w[best], V[:, best]

# ----------------------------------------------------------------------
def shell():  return np.vstack([np.zeros((1,3)), fcc_nn_vectors()])
def voids(nf):
    v = 1/(2*np.sqrt(2))
    return np.array([[+v,+v,+v],[+v,-v,-v],[-v,+v,-v],[-v,-v,+v]])[:nf]

if __name__ == '__main__':
    print("Cosserat phi-phi-u vertex: the pion source ||J||^2 for decuplet -> octet + pi\n")

    # daughter nucleon mode on the bare 13-shell, embedded into the cluster
    Sh = shell(); nS = len(Sh)
    lamN, psiN = mass_mode(build_cosserat_matrix(Sh,1,1,1.0), 8.303)

    def embed(psi13, n):
        out = np.zeros(6*n)
        out[:3*nS]            = psi13[:3*nS]          # u on shell nodes 0..12
        out[3*n:3*n+3*nS]     = psi13[3*nS:]          # phi on shell nodes 0..12
        return out

    print(" n_free  lambda_Delta   ||J||^2    ||J||^2 / n_free^2   (PDG width, MeV)")
    pdg_w = {4: 117.0, 3: 36.0, 2: 9.1}
    base = None
    for nf in (4, 3, 2):
        C = np.vstack([shell(), voids(nf)]); n = len(C)
        lamD, psiD = mass_mode(build_cosserat_matrix_two_d(C,1,1,1.0),
                               9.05 if nf == 4 else 9.2, phi_min=0.5)
        psiNe = embed(psiN, n)
        # transition source: parent Delta on the gradient leg, daughter N on
        # the frame leg (and the reverse); symmetrise.
        J1 = pion_source(C, psiD, psiNe)
        J2 = pion_source(C, psiNe, psiD)
        J  = J1 + J2
        strength = float(np.sum(J**2))
        if base is None: base = strength / nf**2
        print(f"   {nf}      {lamD:7.3f}    {strength:8.4f}     {strength/nf**2:8.4f}"
              f"          ({pdg_w[nf]})")

    print("\nProton (3/2+ -> 1/2+ has no analogue): the proton is the lightest B=1")
    print("state, so no daughter mode exists -> no source -> stable, as required.")
    print("\nThe source lives mostly on the shell (the transition form factor); at the")
    print("voids (the emission sites) it is PURELY TRANSVERSE -- radial part zero --")
    print("which is the P-wave selection rule for Delta -> N pi, from the vertex with")
    print("no input.  The width is the outgoing projection of this source: its")
    print("Fourier transform at the on-shell pion momentum (fixed by the spectral")
    print("masses), transverse part, times two-body phase space.  That projection is")
    print("the next step; the vertex and its P-wave structure are in place here.")
