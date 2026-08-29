#!/usr/bin/env python3
"""What alpha knows about the screw core width, and what the neutrino may not touch.

Committed as the record of the width-vs-alpha analysis (chat session
2026-08-28/29), before the first-principles misfit-potential calculation runs,
so the targets are on file and the later calculation can be checked for
blindness against them.

1. alpha^-1 is exponentially steep in w/d, so inverting CODATA pins the
   screw's Peierls-Nabarro width ratio to eleven digits:
       w/d in [0.785398163373, 0.785398163421]   (1 sigma)
   pi/4 = 0.785398163397 sits +0.02 sigma inside. pi/4 is underived but it
   is MEASURED, more precisely than anything else in the framework.
2. The framework's own variational PN equilibrium returns 0.78616, which
   gives alpha^-1 = 137.703 (+0.487%, +3.2e7 sigma). The 0.1% residual is
   real work still to do, attributed by the alpha paper to the single-
   harmonic (Frenkel) misfit fidelity and flagged there as an open
   first-principles calculation. That calculation is the next one run.
3. The neutrino-demanded edge width (edge_width_prediction.py), imposed as
   a ratio on the screw, gives alpha^-1 = 144.73 (+5.6%, +3.7e8 sigma):
   the two characters cannot share the ratio 0.79397. Dead branch.
4. Surviving reading: if w/d = pi/4 is universal across characters (the
   variational curve has a broad maximum in N^2, hence character-robust),
   the neutrino fixes the EDGE's misfit period instead:
       d_edge = 0.58365 ell   vs   d_screw = 1/sqrt(3) = 0.57735 ell,
   a 1.1% lattice-geometry question rather than a 33% elasticity violation
   (the classical edge width w_s/(1-nu) is rejected by the splittings at
   +58/+89 sigma).
5. Coincidence watch, flagged and NOT used: N^2 = 0.31699 would make the
   screw's own width give unit mass coefficient; 1/pi = 0.31831 (-0.42%).
"""
import numpy as np
from scipy.optimize import brentq

PI = np.pi; N2 = 1/PI
CODATA, SIG = 137.035999177, 2.1e-8

def alpha_inv(r):
    T = 1/137.0
    B = np.exp(2*PI*r)
    for _ in range(400):
        T = 1/(B - 2 - T - N2*T/(1-T) - 6*N2**2*T**3)
    return 1/T

def main():
    r0 = PI/4
    lo = brentq(lambda r: alpha_inv(r)-(CODATA-SIG), 0.78, 0.79)
    hi = brentq(lambda r: alpha_inv(r)-(CODATA+SIG), 0.78, 0.79)
    print(f"alpha^-1(pi/4) = {alpha_inv(r0):.9f}   CODATA {CODATA}")
    print(f"1-sigma band on w/d: [{lo:.12f}, {hi:.12f}]  (half-width {0.5*(hi-lo):.2e})")
    for lab, r in [("variational PN equilibrium", 0.78616),
                   ("neutrino ratio on the screw", 0.79397)]:
        a = alpha_inv(r)
        print(f"{lab:28s} r={r:.5f}: alpha^-1={a:.4f}  {(a-CODATA)/SIG:+.1e} sigma")
    print(f"shared-pi/4 reading: d_edge = {0.4584/r0:.5f} ell vs d_screw = {1/np.sqrt(3):.5f} ell")

if __name__ == "__main__":
    main()
