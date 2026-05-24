"""
Synthesis test: the decay matrix element as a parent->daughter mode overlap
across the Fiedler cut.

  Gamma(i->f) = 2 pi | <psi_f (x) psi_pi | V_cut | psi_i> |^2 rho_f

We test the elastic core of this: the overlap of the parent's spectral mass
mode with the daughter's spectral mass mode, mediated by the coupling across
the cut (the void-face bonds the graph factorisation identifies).

Decuplet -> octet + pion:  Delta(17) -> N(13) + pi,  via void deactivation.
Parent = void-activated cluster mass mode.  Daughter = bare-shell nucleon mode.
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
from delta_first_principles import build_cosserat_matrix_two_d
from cosserat_classifier import build_cosserat_matrix

def shell():            return np.vstack([np.zeros((1,3)), fcc_nn_vectors()])
def voids(n):
    v = 1/(2*np.sqrt(2))
    a = np.array([[+v,+v,+v],[+v,-v,-v],[-v,+v,-v],[-v,-v,+v]])
    return a[:n]
def delta_cluster(n): return np.vstack([shell(), voids(n)])

def mass_mode(Phi, target, phi_min=0.5):
    """Eigenvector whose lambda is closest to target with phi-content >= phi_min."""
    w,V = np.linalg.eigh(Phi); n6 = Phi.shape[0]; n = n6//6
    best,bs = None,1e9
    for k in range(len(w)):
        if np.sum(V[3*n:,k]**2) < phi_min: continue
        if abs(w[k]-target) < bs: best,bs = k,abs(w[k]-target)
    return w[best], V[:,best]

# --- daughter nucleon mode on the bare 13-shell ---
Sh = shell(); nS = len(Sh)
Phi_N = build_cosserat_matrix(Sh, 1,1,1.0)
lamN, psiN = mass_mode(Phi_N, 8.303)
print(f"Daughter N: lambda={lamN:.4f}  (bare 13-shell A_2u mode)")

def embed_shell_into(coords):
    """Index map: the first 13 nodes of the delta cluster ARE the shell."""
    n = len(coords)
    P = np.zeros((6*n, 6*nS))     # maps a 6*13 shell vector into 6*n cluster space
    for i in range(nS):           # shell node i -> cluster node i (same ordering)
        P[3*i:3*i+3, 3*i:3*i+3] = np.eye(3)             # u
        P[3*n+3*i:3*n+3*i+3, 3*nS+3*i:3*nS+3*i+3] = np.eye(3)  # phi
    return P

print("\n n_free  lambda_parent  <psi_N|psi_Delta>  ||V_cut psi_Delta||  (cut=void-face block)")
for nf in (4,3,2):
    C = delta_cluster(nf); n = len(C)
    Phi = build_cosserat_matrix_two_d(C, 1,1,1.0)
    lamP, psiP = mass_mode(Phi, 9.05 if nf==4 else 9.2)
    # embed daughter
    P = embed_shell_into(C)
    psiN_emb = P @ psiN
    psiN_emb /= np.linalg.norm(psiN_emb)
    overlap = abs(psiN_emb @ psiP)
    # cut coupling operator: the off-diagonal block connecting voids to shell.
    # voids are cluster nodes 13..13+nf-1. Build a mask for void dofs.
    void_dof = np.zeros(6*n, bool)
    for vi in range(13, 13+nf):
        void_dof[3*vi:3*vi+3] = True
        void_dof[3*n+3*vi:3*n+3*vi+3] = True
    Vcut = Phi.copy()
    # keep only entries that cross the void/non-void partition (the cut)
    cross = np.outer(void_dof, ~void_dof) | np.outer(~void_dof, void_dof)
    Vcut = Phi * cross
    coupling = np.linalg.norm(Vcut @ psiP)
    print(f"   {nf}      {lamP:7.4f}        {overlap:7.4f}            {coupling:7.4f}")

# --- proton: ground-state baryon, no lighter B=1 daughter exists ---
print("\nProton: B=1 ground state. No lighter B=1 daughter cluster exists,")
print("so there is no psi_daughter to overlap with -> matrix element = 0 -> stable.")
