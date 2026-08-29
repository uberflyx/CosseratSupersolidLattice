#!/usr/bin/env python3
"""The edge core width from its own Cosserat PN equilibrium, and what it settles.

THE CALCULATION (the neutrino sector's endgame, run blind of the demanded width)
Same machinery the alpha paper used for the screw, with the plane-strain
kernel. The screw's anti-plane kernel is certified exactly
(neutrinos/screw_vertex_calibration.py):
    K(k) = (mu_tot/2) |k| (k^2 + p^2)/(k^2 + q^2),
and the edge's in-plane kernel is stiffer by 1/(1-nu) with nu = 1/4 fixed by
the framework's mu_V/(3K) = 1/5. Since the misfit side is a property of the
interface and not of the dislocation character, the edge equilibrium is the
screw equilibrium with the right-hand side scaled by (1 - nu):
    Int_0^inf K(k)/k e^{-2kw} dk = (1-nu) mu_bar/(2 d_111).
Validation: the same code at (1-nu) -> 1 returns w/d = 0.78616, the alpha
paper's printed screw value, to all digits.

RESULT, AND IT IS AN INDEPENDENT CONFIRMATION
    w_edge = 0.59026 ell      (w/d = 1.02235)
    edge/screw width ratio = 1.3004  (classical elasticity: 1/(1-nu) = 1.3333;
                                      the Cosserat coupling trims it slightly)
The framework's brane geometry independently fixes the edge slip patch at
L_b = 0.589 ell. Two routes that share no machinery, continuum PN equilibrium
and D_4 brane compactification, agree to 0.2 per cent. The edge core width is
therefore settled, and it is NOT the 0.4584 ell that edge_width_prediction.py
inverted from the splittings.

WHAT THAT FALSIFIES: MY OWN COEFFICIENT, NOT THE MASS SCALE
At the settled width the coefficient formula of edge_width_prediction.py
returns 1.2869, which drives Dm2_21 to +51.5 sigma and Dm2_31 to +79.4 sigma.
The width cannot absorb that, so the formula is wrong. Auditing it against the
cancellation argument that preceded it exposes the defect: the two derivations
used DIFFERENT intermediate gaps for the same localised core libration.
    cancellation version : Delta = gamma kappa^2 ell^3       (gradient energy)
    coefficient version  : Delta = kappa_c ell^3             (contact restoring)
A localised mode confined at the Cosserat length ell_c = ell/2 carries both,
    H = gamma/(mu ell_c^2) + 2 kappa_c/mu = 4.000 + 3.504 = 7.504,
against the 1.752 the coefficient formula used. The honest statement is a
bracket over the gap treatment at the now-settled width:
    coefficient in [0.30, 1.29],  which straddles unity.
So the calculation is CONSISTENT with m_1 = theta_ch^2 m_0 and does not
confirm it. What it does settle is that the remaining uncertainty lives
entirely in the gap, a single localised-mode eigenvalue, and no longer in the
width, which two independent routes now fix.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

N2 = 1/np.pi
P2 = 4*N2*(1-N2)/(1-2*N2); Q2 = 4*N2/(1-2*N2)
MUT = 1/(1-2*N2); MUBAR = (1-N2)/(1-2*N2)
D111 = np.sqrt(2/3); DHOP = 1/np.sqrt(3); NU = 0.25
LC = 0.5

def kernel(k):
    return 0.5*MUT*(k*k + P2)/(k*k + Q2)

def equil(rhs_scale):
    f = lambda w: quad(lambda k: kernel(k)*np.exp(-2*k*w), 0, np.inf, limit=400)[0] \
                  - rhs_scale*MUBAR/(2*D111)
    return brentq(f, 0.15, 1.2, xtol=1e-13)

def main():
    w_s = equil(1.0); w_e = equil(1-NU)
    print(f"screw validation : w/d = {w_s/DHOP:.5f}  (alpha paper 0.78616)")
    print(f"edge             : w   = {w_e:.5f} ell,  w/d = {w_e/DHOP:.5f}")
    print(f"edge/screw ratio = {w_e/w_s:.4f}   (classical 1/(1-nu) = {1/(1-NU):.4f})")
    print(f"brane geometry   : L_b = 0.589 ell  ->  agree to {100*abs(w_e-0.589)/0.589:.1f} per cent")
    base = 6*w_e*(1-2*N2)/(np.pi*N2)
    kap = 2*N2/(1-2*N2)
    print("\ncoefficient bracket over the localised-gap treatment:")
    for lab, H in [("restoring only", kap), ("gradient only", 1/LC**2),
                   ("both (full mode)", 1/LC**2 + 2*kap)]:
        c = base*kap/H; s = c*c
        p21 = (74.711*s - 75.37)/(0.94 if 74.711*s > 75.37 else 1.00)
        p31 = (2523.4*s - 2511.0)/(21.0 if 2523.4*s > 2511.0 else 20.0)
        print(f"  H = {H:6.3f} [{lab:16s}]  coeff = {c:.4f}  "
              f"Dm21 {p21:+7.1f}s  Dm31 {p31:+7.1f}s")
    print("  bracket straddles unity: consistent with theta^2 m_0, not a confirmation.")

if __name__ == "__main__":
    main()
