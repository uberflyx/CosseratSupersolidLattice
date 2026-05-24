#!/usr/bin/env python3
"""
strong_widths.py
================

The decuplet strong widths from the synthesis: factorisation channel
strength times derived kinematics.

Across this programme the two halves separated cleanly.

  Selection rule + kinematics (DERIVED from the spectral eigenvectors):
    - the decay is P-wave (|A|^2 ~ p^2), shown by contracting the rotation
      vertex with the parent and daughter mass modes;
    - the proton is stable (no daughter mode -> no source);
    - the lattice scale is fixed, p_hat = p*/m_0, via ell = r_e.

  Channel strength (INHERITED from the graph factorisation):
    - the n_free = 4 - |S| active voids emit coherently, so the amplitude
      scales as n_free and the rate as n_free^2.

Married, the relativistic two-body width is

    Gamma_i = (g^2 / 8pi) * n_free^2 * p*^3 / (m_i^2 m_0^2),

with p*^2 the P-wave factor, p* the phase-space density, and one coupling
g calibrated on the Delta.  Sigma* and Xi* are then predictions.

p* is computed from the PDG masses (dominant pi channel).

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

M_E = 0.51099895
ALPHA = 1/137.0359992
M_0 = M_E/ALPHA

def p_star(mi, mf, mpi):
    return np.sqrt((mi**2-(mf+mpi)**2)*(mi**2-(mf-mpi)**2))/(2*mi)

MPI = 138.0
# (parent mass, dominant daughter mass, n_free, PDG width, dominant channel)
DECUPLET = {
    'Delta(1232)':  (1232.0, 938.9,  4, 117.0, 'N pi'),
    'Sigma*(1385)': (1384.6, 1115.7, 3, 36.0,  'Lambda pi'),
    'Xi*(1530)':    (1531.8, 1318.3, 2, 9.9,   'Xi pi'),
}

def width(mi, mf, nfree, g2):
    ps = p_star(mi, mf, MPI)
    return g2 * nfree**2 * ps**3 / (mi**2 * M_0**2) / (8*np.pi), ps

if __name__ == '__main__':
    # calibrate g^2 on the Delta
    mi, mf, nf, pdg, _ = DECUPLET['Delta(1232)']
    g_unit, _ = width(mi, mf, nf, 1.0)
    g2 = pdg / g_unit
    print(f"m_0 = {M_0:.3f} MeV;  one coupling g^2 = {g2:.4e} fixed on the Delta.\n")
    print(f"{'state':14s} {'channel':10s} n_free  p*(MeV)  Gamma_pred   PDG    ratio")
    for name, (mi, mf, nf, pdg, ch) in DECUPLET.items():
        g, ps = width(mi, mf, nf, g2)
        print(f"{name:14s} {ch:10s}  {nf}    {ps:6.1f}   {g:7.1f}   {pdg:6.1f}   {g/pdg:4.2f}")
    print("\nThe ordering and the Sigma* are reproduced; the Xi* sits low by ~1.6x,")
    print("the size of the secondary channels and the hex-cap shadow not yet")
    print("included.  The P-wave law and the n_free coherence are not fitted:")
    print("only the single overall coupling is.")
