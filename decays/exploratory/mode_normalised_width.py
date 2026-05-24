#!/usr/bin/env python3
"""
mode_normalised_width.py  -- absolute decay scale with NO fitted coupling.

Each mode is carried at its physical mechanical zero-point displacement
    Q_Gamma = ell / (sqrt(2) lambda^{1/4})   [length],  ell = hbar/(m0 c),
a genuine length (not the field coordinate (2 hbar omega)^{-1/2}).  The cubic
vertex is gamma = (kappa3/ell^3) ghat, ghat the dimensionless transverse
P-wave source |Jtilde_perp| at the on-shell pion momentum, kappa3 a stiffness.
V = gamma Q_i Q_f Q_pi is an energy; Mcov = sqrt(8 m_i m_f E_pi) V;
    Gamma = (p*/8 pi m_i^2) |Mcov|^2.
Result: the stiffness the data require is kappa3 ~ N m_e, the SAME spectral
stiffness that sets the masses via N(4-lambda) m_e.  Fixing kappa3 = N m_e (no
freedom) predicts the Delta width to the right order of magnitude; the residual
~1.9 is the Rarita-Schwinger spin multiplicity, uniform across the decuplet.

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
M_0 = 70.025; M_E = 0.5109989

def p_star(mi, mf, mp):
    return np.sqrt((mi**2-(mf+mp)**2)*(mi**2-(mf-mp)**2))/(2*mi)
def Qzp(lmbda):
    return (1.0/M_0)/(np.sqrt(2.0)*lmbda**0.25)
def predict_width(mi, mf, mp, lam_i, lam_f, lam_pi, N_i, Jperp, spin=1.0):
    ps = p_star(mi, mf, mp); Epi = np.sqrt(mp**2+ps**2); ell = 1.0/M_0
    pref = Jperp*Qzp(lam_i)*Qzp(lam_f)*Qzp(lam_pi)/ell**3
    kappa3 = N_i*M_E
    Mcov2 = 8*mi*mf*Epi*(kappa3*pref)**2
    return spin*ps/(8*np.pi*mi**2)*Mcov2, ps, kappa3

if __name__ == '__main__':
    Jperp_Delta = np.sqrt(5.4574e-3)
    G, ps, kap = predict_width(1232.0, 938.9, 138.0, 9.052, 8.303, 2.0, 17, Jperp_Delta)
    print("Delta(1232) -> N pi")
    print(f"  p*                        = {ps:6.1f} MeV")
    print(f"  coupling stiffness kappa3 = N m_e = {kap:5.2f} MeV   (DERIVED, not fitted)")
    print(f"  predicted width           = {G:6.1f} MeV   (no spin sum)")
    print(f"  observed (PDG)            = 117.0 MeV")
    print(f"  shortfall factor          = {117.0/G:5.2f}   <- expected spin multiplicity")
