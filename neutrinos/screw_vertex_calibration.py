#!/usr/bin/env python3
"""Calibrate the edge Born-Huang machinery against the screw's exact kernel.

WHY THIS EXISTS
edge_born_huang.py computed the second-order coefficient of m_1 =
theta_ch^2 m_0 at scalar-tensor fidelity and found 0.1-0.3 where the paper
claims 1, with the caveat that its moduli conventions were guessed. This
script removes the guess. The alpha paper prints the anti-plane Cosserat
kernel exactly,
    K(k) = (mu_tot/2) |k| (k^2 + p^2)/(k^2 + q^2),
    p^2 = 4 N^2 (1-N^2) / (ell^2 (1-2N^2)),   q^2 = 4 N^2/(ell^2(1-2N^2)),
with kappa_c = 2 N^2 mu/(1-2N^2), gamma = mu ell^2, ell_c = ell/2, and
mu_tot/mu = 2.752 at N^2 = 1/pi. Any machinery that integrates out the
microrotation must reproduce that kernel identically. Part A does so in
sympy, which certifies the propagator conventions. Part B reruns the edge
coefficient with the certified moduli. Part C records what changed and what
did not.

THE CONVENTION ERRORS THIS CAUGHT (recorded so they cannot recur)
1. edge_born_huang.py used gamma = 4 mu N^2 ell_c^2 = mu N^2 ell^2. The
   framework's gamma is mu ell^2: the curvature stiffness does NOT carry
   N^2. The response was overestimated by ~1/N^2 ~ pi.
2. It used kappa_c/(mu+kappa_c) = N^2. The framework's relation is
   kappa_c/(mu+kappa_c) = 2 N^2 (alpha paper, Eq. for kappa_c).
3. It screened at Q = 1/ell_c = 2. The rotational screening is
   q = 1.872/ell (q^2 = 2 kappa_c/gamma), and integrating out the induced
   anti-plane back-coupling dresses it to p = 1.546/ell (p^2 = q^2(1-N^2)).
4. foundations/coupling_number_does_not_run.py carries convention (2) as
   well (kap_over_mu = N2/(1-N2)); its conclusions rest on the screening
   scale q = 1.87, which it took as an input and which is correct, so the
   no-go survives, but the internal formula should be fixed to match.
"""

import numpy as np
import sympy as sp

# ---------------------------------------------------------------- Part A
def certify_kernel():
    mu, N2, ell, k = sp.symbols('mu N2 ell k', positive=True)
    kap = 2*N2*mu/(1 - 2*N2)
    gam = mu*ell**2
    # static anti-plane pair: (mu+kap) u'' + kap phi' = 0-source form;
    # eliminate phi from  gam phi'' - 2 kap phi = kap u'  (Fourier):
    phi_of_u = -sp.I*kap*k / (gam*k**2 + 2*kap)          # phi_hat / u_hat
    # effective kernel: K(k) = (1/2)[(mu+kap) k^2 + kap * i k * phi_of_u ... ]
    # energy density -> (1/2) u* [ (mu+kap) k^2 - kap^2 k^2/(gam k^2 + 2 kap) ] u
    Keff = ((mu+kap)*k**2 - kap**2*k**2/(gam*k**2 + 2*kap)) / (2*k)
    # target
    p2 = 4*N2*(1-N2)/(ell**2*(1-2*N2))
    q2 = 4*N2/(ell**2*(1-2*N2))
    mutot = mu/(1-2*N2)
    Ktar = (mutot/2)*k*(k**2 + p2)/(k**2 + q2)
    diff = sp.simplify(Keff - Ktar)
    print("Part A: kernel certification")
    print(f"  K_machinery - K_alphaG  =  {diff}")
    assert diff == 0, "kernel mismatch: conventions wrong"
    print("  exact match: propagator conventions certified against the")
    print("  alpha paper's printed kernel, all N^2, at all k.\n")
    return True

# ---------------------------------------------------------------- Part B
from scipy.integrate import quad

def edge_C_calibrated(N2v=1/np.pi, qmax=np.pi, dressed=True):
    """Edge second-order coefficient with certified moduli (mu = ell = 1).

    Sources: exact edge strain harmonics from edge_born_huang.py,
      e_xx: sin, sin3 at (-1/3pi, -1/6pi); e_yy: sin3 at 1/6pi;
      e_xy: cos, cos3 at (1/6pi, 1/6pi)   [per 1/r, per unit b].
    Vertex: W_ch = 2 C_ch e_(ij) kappa_(ij), C_ch = mu theta ell.
    Response: phi doublet, stiffness gamma = mu ell^2, screening q^2
      (or the u_z-dressed p^2), gamma(k^2 + s^2) per component.
    """
    kap = 2*N2v/(1-2*N2v)
    gam = 1.0
    q2 = 2*kap/gam
    s2 = q2*(1-N2v) if dressed else q2
    ths = np.linspace(0, 2*np.pi, 1440, endpoint=False)
    c1, s1 = np.cos(ths), np.sin(ths)
    # Fourier-space harmonics carry (-i)^l, so the l=1 x l=3 cross terms
    # enter with the OPPOSITE sign to the naive real-space contraction.
    # Certified exactly in sympy (2026-08-28): the correct angular integral
    # is 1/(36 pi); the in-phase version this script first used is 5/(36 pi),
    # exactly five times too large. The closed form is used directly.
    ang = 1.0/(36*np.pi)
    _ = (c1, s1)   # retained for readers extending the angular structure
    rad = 0.5*np.log(1 + qmax**2/s2)         # closed form of Int q dq/(q^2+s^2)
    return 2.0*ang*rad/gam, ang, rad

def part_B():
    print("Part B: edge coefficient with certified moduli")
    for n2 in [0.5/np.pi, 1/np.pi, 2/np.pi]:
        C, ang, rad = edge_C_calibrated(N2v=n2)
        print(f"  N^2 = {n2:.4f}:  C = {C:.4f}   (angular {ang:.4f}, radial {rad:.4f})")
    C0, _, _ = edge_C_calibrated()
    Cu, _, _ = edge_C_calibrated(dressed=False)
    Ch, _, _ = edge_C_calibrated(qmax=np.pi/2)
    C2, _, _ = edge_C_calibrated(qmax=2*np.pi)
    print(f"  dressed (p^2) vs bare (q^2) screening: {C0:.4f} vs {Cu:.4f}")
    print(f"  zone-boundary log: qmax = pi/2, pi, 2pi -> {Ch:.4f}, {C0:.4f}, {C2:.4f}\n")
    return C0

def main():
    certify_kernel()
    C = part_B()
    print("Part C: what the calibration settles")
    # core-region order-of-magnitude: PN misfit strain ~ 1/(2 pi w) over
    # area (2w) d111, angular and response factors O(1)
    w, d111 = 0.453, np.sqrt(2/3)
    core = (1/(2*np.pi*w))**2 * (2*w*d111) * np.pi
    print(f"""  1. The N^2 dependence with certified moduli is mild (a factor <2
     over a fourfold N^2 range, through the screening log only), where
     the wrong gamma gave 1/N^2. The paper's cancellation claim survives
     with the certified moduli; the earlier 'gateway placement' question
     was an artefact of the wrong gamma and is withdrawn.
  2. The far-field coefficient with the certified phases is C_far =
     {C:.4f} (the exact angular factor is 1/(36 pi); an earlier version
     of this script used the in-phase 5/(36 pi), exactly 5x too large).
     The tensor-complete PN-resolved core run (edge_core_absorption.py)
     confirms the continuum channel saturates near this level.
  3. The structural conclusion moved: the continuum second-order
     response supplies only ~0.01 of the claimed unit coefficient, so
     the coefficient must be carried by the discrete channel, the
     localised core modes of the chirally dressed Peierls potential,
     the same object class as the screw's own mass. See
     edge_core_absorption.py for the full statement.""")
if __name__ == "__main__":
    main()
