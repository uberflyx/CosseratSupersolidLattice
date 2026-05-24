#!/usr/bin/env python3
"""
decay_width.py
==============

The strong width as the outgoing projection of the pion source.

The source J(node) (cosserat_vertex.py) is the displacement the parent->
daughter microrotation transition drives.  The pion is radiated by this
source; by Fermi's golden rule the rate is the source's Fourier transform
at the on-shell pion momentum, projected on the pion polarisation,
integrated over emission directions, times the two-body phase space.

Lattice -> physical conversion (no new parameter).
The framework fixes ell = r_e through the Mott condition m_0 c^2 = hbar c/ell.
Hence a physical momentum p (MeV) is the dimensionless lattice phase

    p_hat = p / m_0,    m_0 = m_e/alpha = 70.025 MeV,

because  p . r_phys / hbar = (p/m_0) * r_lattice  with r_lattice in NN units.

On-shell pion momentum for i -> f + pi (rest frame of i):
    p* = sqrt( [m_i^2-(m_f+m_pi)^2][m_i^2-(m_f-m_pi)^2] ) / (2 m_i).

Relativistic two-body width:
    Gamma = g^2 * (p*/(8 pi m_i^2)) * <|A(p*)|^2>_dir,
where A is the pion-polarisation projection of the Fourier source and g is
the one overall coupling (calibrated once, then tested on Sigma*, Xi*).
The P-wave structure (Gamma ~ p*^3) must emerge from A(p*) ~ p* on its own.

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
from cosserat_vertex import pion_source, mass_mode, shell, voids
from cosserat_classifier import build_cosserat_matrix
from delta_first_principles import build_cosserat_matrix_two_d

M_E = 0.51099895
ALPHA = 1/137.0359992
M_0 = M_E/ALPHA                 # 70.025 MeV  -> sets p_hat = p/M_0

def p_star(mi, mf, mpi):
    return np.sqrt((mi**2-(mf+mpi)**2)*(mi**2-(mf-mpi)**2))/(2*mi)

def fibonacci_sphere(N):
    """N roughly-uniform directions on the unit sphere."""
    i = np.arange(N) + 0.5
    phi = np.arccos(1 - 2*i/N); th = np.pi*(1+5**0.5)*i
    return np.c_[np.sin(phi)*np.cos(th), np.sin(phi)*np.sin(th), np.cos(phi)]

def fourier_source(J, coords, p_vec):
    """tilde J(p) = sum_n J_n exp(-i p.r_n);  p_vec in lattice units."""
    ph = np.exp(-1j*(coords @ p_vec))
    return (J * ph[:, None]).sum(axis=0)        # complex 3-vector

def projected_strength(J, coords, p_hat, ndir=400):
    """Direction-averaged |A|^2 for longitudinal (scalar pion) and
    transverse (vector) polarisation projections of the Fourier source."""
    dirs = fibonacci_sphere(ndir)
    long2 = trans2 = 0.0
    for n in dirs:
        Jt = fourier_source(J, coords, p_hat*n)
        aL = n @ Jt                              # longitudinal (scalar)
        aT2 = np.sum(np.abs(Jt - n*aL)**2)       # two transverse pols
        long2 += abs(aL)**2; trans2 += aT2
    return long2/len(dirs), trans2/len(dirs)

# ----------------------------------------------------------------------
if __name__ == '__main__':
    print(f"m_0 = {M_0:.4f} MeV  ->  p_hat = p*/m_0\n")

    # --- Delta -> N pi -------------------------------------------------
    mD, mN, mpi = 1232.0, 938.92, 138.04
    ps = p_star(mD, mN, mpi); phat = ps/M_0
    print(f"Delta -> N pi:  p* = {ps:.1f} MeV  ->  p_hat = {phat:.3f} (BZ edge ~ pi)")

    nS = 13
    _, psiN = mass_mode(build_cosserat_matrix(shell(),1,1,1.0), 8.303)
    C = np.vstack([shell(), voids(4)]); n = len(C)
    _, psiD = mass_mode(build_cosserat_matrix_two_d(C,1,1,1.0), 9.05, phi_min=0.5)
    psiNe = np.zeros(6*n); psiNe[:3*nS]=psiN[:3*nS]; psiNe[3*n:3*n+3*nS]=psiN[3*nS:]
    J = pion_source(C, psiD, psiNe) + pion_source(C, psiNe, psiD)

    L2, T2 = projected_strength(J, C, phat)
    print(f"  direction-averaged |A|^2:  longitudinal {L2:.4e}   transverse {T2:.4e}")

    # P-wave check: how does the projected strength scale with p_hat?
    print("\n  P-wave scaling check (longitudinal A vs p_hat):")
    for f in (0.25, 0.5, 1.0):
        l2,_ = projected_strength(J, C, phat*f)
        print(f"    p_hat={phat*f:5.3f}:  |A|^2={l2:.4e}   |A|^2/p_hat^2={l2/(phat*f)**2:.4e}")

    # calibrate the one coupling g on the Delta width, report g
    L2use = max(L2, T2)
    pref = ps/(8*np.pi*mD**2)
    g2 = 117.0/(pref*L2use)
    print(f"\n  prefactor p*/(8 pi m^2) = {pref:.4e}")
    print(f"  to match Gamma_Delta = 117 MeV, need g^2 = {g2:.4e}, g = {np.sqrt(g2):.3f}")
    print(f"  compare framework quanta:  m_0 = {M_0:.1f},  m_0/pi = {M_0/np.pi:.1f},  "
          f"sqrt(g^2) in MeV-ish? {np.sqrt(g2):.2f}")
