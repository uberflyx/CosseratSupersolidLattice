#!/usr/bin/env python3
"""
tier2_complex_mass.py
=====================

Mass AND width from a single complex diagonalisation of one cluster.

The closed-cluster picture (hadron_spectral_mass.py) diagonalises the
Hermitian Cosserat matrix Phi_Cos.  Its eigenvalues are real, so every
mass mode rings forever: the cluster has a sharp mass and no width.  A
real hadron sits in the infinite FCC medium, and its mass mode can leak
elastic energy across the cluster boundary into the surrounding lattice.
That leakage is the decay.

We open the cluster by adding a surface self-energy to the *energy*
operator (in MeV, not in dimensionless lambda units, because the width
is a leading strong effect at the m_0 scale while the mass mode-splitting
sits at the N*m_e scale):

    H      = E_0 * I + (N m_e) * Phi_Cos          (Hermitian, MeV)
    H_open = H - i W                              (non-Hermitian, MeV)

with E_0 = N (m_0 - 4 m_e) the mode-independent assembly offset, so the
mass-mode eigenvalue of H is exactly the spectral mass
m = N m_0 - N (4 - lambda) m_e.

The surface self-energy W is built from the framework's own radiation
quanta, not from a fit:

    q_phi = m_0 / pi   microrotation (node-rolling) radiation quantum,
                       the Cosserat rolling-coupled node energy that
                       already sets the node term in the factorisation
                       theorem Gamma = 2(|dE| m_e + |dV| m_0/pi).
    q_u   = m_e        displacement (bond-stretch) radiation quantum,
                       the electromagnetic Peach-Koehler bond energy
                       (one severed boundary bond = m_e, the molecular
                       Mode-F anchor: P_c -> 2*6*m_e = 6.13 MeV).

Each lattice node radiates in proportion to its coordination deficit
f_i = (12 - c_i)/12, the fraction of its bulk-FCC bonds that face the
medium rather than another cluster node.  In the 6n site-dof basis W is
diagonal:

    W[u_i]   = gate * f_i * q_u
    W[phi_i] = gate * f_i * q_phi

The complex eigenvalue tracking the mass mode is then m - i Gamma/2, with

    Gamma = 2 <psi | W | psi>
          = 2 sum_i gate * f_i * ( q_u |u_i(psi)|^2 + q_phi |phi_i(psi)|^2 ).

The 'gate' is the topological filter: 1 if the mass mode can radiate into
an open channel, 0 if a conserved topological charge (baryon winding for
the proton) seals every channel.  This is the one structural input that
the quadratic spectrum cannot see by itself; everything else is geometry
plus the two existing quanta.

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
from composite_clusters import cluster_crossed_fault

# ----------------------------------------------------------------------
# Constants (CODATA 2022)
# ----------------------------------------------------------------------
M_E   = 0.51099895          # MeV
ALPHA = 1.0 / 137.0359992
M_0   = M_E / ALPHA         # 70.0253 MeV  (node rest energy)
Q_PHI = M_0 / np.pi         # 22.29 MeV    (strong / rolling radiation quantum)
Q_U   = M_E                 # 0.511 MeV    (electromagnetic bond quantum)
ELL   = np.linalg.norm(fcc_nn_vectors()[0])   # FCC NN distance in working units

PDG = {
    'proton':   938.272,  'rho':    775.26,  'omega':  782.66,
    'Delta':    1232.0,   'Sigma*': 1384.4,  'Xi*':    1533.4,
}
PDG_WIDTH = {
    'proton':   0.0,      'rho':    147.4,   'omega':  8.68,
    'Delta':    117.0,    'Sigma*': 36.0,    'Xi*':    9.1,
}

# ----------------------------------------------------------------------
# Cluster builders (canonical geometry, reused from the mass scripts)
# ----------------------------------------------------------------------
def cluster_proton():
    return np.vstack([np.zeros((1, 3)), fcc_nn_vectors()])           # 13 nodes

def _voids_Td_plus():
    v = 1.0 / (2.0 * np.sqrt(2.0))
    return np.array([[+v, +v, +v], [+v, -v, -v],
                     [-v, +v, -v], [-v, -v, +v]])                    # T_d^+ orbit

def cluster_delta(n_voids=4):
    return np.vstack([cluster_proton(), _voids_Td_plus()[:n_voids]])

# ----------------------------------------------------------------------
# Surface exposure: coordination deficit relative to bulk FCC (12 NN)
# ----------------------------------------------------------------------
def coordination_deficit(coords, ell=ELL, tol=1e-2):
    """f_i = (12 - c_i)/12, the medium-facing bond fraction of each node."""
    n = len(coords)
    c = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if abs(np.linalg.norm(coords[j] - coords[i]) - ell) < tol:
                c[i] += 1
    return np.clip((12.0 - c) / 12.0, 0.0, 1.0)

# ----------------------------------------------------------------------
# The unified complex eigenvalue
# ----------------------------------------------------------------------
def complex_pole(coords, target_lambda, gate=1.0, want_phi_dominant=True,
                 voids=False):
    """Return (mass, width, lambda_re, phi_fraction) of the mass mode.

    Diagonalise the Hermitian Phi_Cos, lock onto the mass mode by matching
    target_lambda (and phi-dominance), build the energy operator and the
    surface self-energy, then diagonalise H - iW and track the same mode.

    voids=True selects the mixed-bond builder (NN bonds plus void-face bonds
    at sqrt(6)/4), required for the void-activated decuplet clusters.
    """
    n = len(coords)
    if voids:
        from delta_first_principles import build_cosserat_matrix_two_d
        Phi = build_cosserat_matrix_two_d(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    else:
        Phi = build_cosserat_matrix(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    w, V = np.linalg.eigh(Phi)

    # phi-fraction of each eigenvector (microrotation content)
    def phi_frac(vec):
        return float(np.sum(vec[3 * n:] ** 2))

    # lock onto the mass mode: closest lambda with the right sector character
    best, best_score = None, 1e9
    for k in range(len(w)):
        if want_phi_dominant and phi_frac(V[:, k]) < 0.5:
            continue
        score = abs(w[k] - target_lambda)
        if score < best_score:
            best, best_score = k, score
    psi = V[:, best]
    lam = w[best]
    pf  = phi_frac(psi)

    # energy operator (MeV): mass-mode eigenvalue = spectral mass
    E0 = n * (M_0 - 4.0 * M_E)
    H  = E0 * np.eye(6 * n) + (n * M_E) * Phi

    # surface self-energy (MeV), diagonal in the site-dof basis
    f = coordination_deficit(coords)
    Wdiag = np.zeros(6 * n)
    for i in range(n):
        Wdiag[3 * i:3 * i + 3]             = gate * f[i] * Q_U     # u dofs
        Wdiag[3 * n + 3 * i:3 * n + 3 * i + 3] = gate * f[i] * Q_PHI  # phi dofs
    Hopen = H - 1j * np.diag(Wdiag)

    # diagonalise the non-Hermitian operator; track the mass mode by overlap
    ev, evec = np.linalg.eig(Hopen)
    overlaps = np.abs(psi @ evec) ** 2
    m_idx = int(np.argmax(overlaps))
    pole = ev[m_idx]

    mass  = pole.real
    width = -2.0 * pole.imag
    return mass, width, lam, pf

# ----------------------------------------------------------------------
# Run the benchmark set
# ----------------------------------------------------------------------
def report(name, coords, target_lambda, gate=1.0, phi_dom=True, voids=False):
    m, g, lam, pf = complex_pole(coords, target_lambda, gate, phi_dom, voids)
    pdg_m, pdg_g = PDG[name], PDG_WIDTH[name]
    dm = (m - pdg_m) / pdg_m * 100
    print(f"{name:8s} N={len(coords):2d}  lam={lam:6.3f} phi={pf*100:4.0f}%  "
          f"m={m:8.2f} (PDG {pdg_m:7.1f}, {dm:+.2f}%)  "
          f"Gamma={g:7.2f} (PDG {pdg_g:6.1f})")
    return m, g

if __name__ == '__main__':
    print(f"Quanta:  m_0 = {M_0:.4f} MeV,  q_phi = m_0/pi = {Q_PHI:.4f} MeV,  "
          f"q_u = m_e = {Q_U:.4f} MeV")
    print(f"FCC NN distance in working units: ell = {ELL:.4f}\n")
    print("Mass = Re(pole), Width = -2 Im(pole), from ONE complex diagonalisation.")
    print("Width here uses a CHANNEL-BLIND surface self-energy (every medium-facing")
    print("dof radiates); see the notes on why this is only the first pass.\n")

    # proton: antenna is large but baryon winding seals every channel -> gate 0
    report('proton', cluster_proton(),      8.303, gate=0.0)
    m_open, g_open, _, _ = complex_pole(cluster_proton(), 8.303, gate=1.0)
    print(f"         (proton with gate forced open: Gamma_would-be = {g_open:.1f} MeV "
          f"-- sealed by B=1 winding)\n")

    report('rho',    cluster_crossed_fault()[0], 4.891, phi_dom=True)
    report('Delta',  cluster_delta(4),       9.052, voids=True)
    report('Sigma*', cluster_delta(3),       9.311, voids=True)
    report('Xi*',    cluster_delta(2),       9.331, voids=True)
