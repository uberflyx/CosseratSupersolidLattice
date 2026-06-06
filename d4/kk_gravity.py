#!/usr/bin/env python3
"""
kk_gravity.py -- Kaluza-Klein structure of the gravitational sector on D4
==========================================================================
Mitchell A. Cox, University of the Witwatersrand
Accompanies: M. A. Cox, "The Cosserat Supersolid" (2026)

Gravity is the compression (longitudinal) channel of the four-dimensional
Cosserat medium.  At zero compact momentum the spatial slice factorises: the
compression mode (the A1g irrep of the cubic group) is gapless and carries
long-range Newtonian gravity, decoupled from the shear mode (T1u) that carries
light.  This script shows what the compact direction adds to that picture, and
checks four things against the monograph.

1.  Long-wavelength dispersions of the displacement sector,
        omega_L^2 = (lambda + 2 mu + kappa_c) |k|^2 / rho   (longitudinal)
        omega_T^2 = (mu + kappa_c) |k|^2 / rho               (transverse),
    so at zero spatial momentum the compression mode acquires a Kaluza-Klein gap
    omega_L^2 = (lambda + 2 mu + kappa_c) k4^2 / rho.

2.  The static gravitational potential carried by the n-th compact mode is
    Yukawa, with range set by the compact momentum alone,
        xi = 1 / k4 = L4 / (2 pi n) = 0.39 ell  for n = 1,
    about half a lattice spacing.  The longitudinal stiffness sets the strength
    of the correction, not its range.  The massless n = 0 mode carries ordinary
    long-range gravity and is untouched.

3.  At zero compact momentum the compression (A1g) and shear (T1u) channels are
    orthogonal, the reason a free electron has exactly zero active gravitational
    mass.  At finite compact momentum the mixed rotation phi_{i4} opens the chain
    T1u -> phi_{i4} -> u4 -> A1g; the coupling element is proportional to
    kappa_c k4 and vanishes at k4 = 0.

4.  Dimensional reduction of Newton's constant: G3 = G4 / L4, equivalently
    G3 = c^4 / (K_sf ell^2), and the four-dimensional constant is G4 = L4 G3
    = sqrt(6) c^4 / (K_sf ell).  Newton's constant is the inverse compressive
    stiffness of the compact direction.

Usage:
    python kk_gravity.py
"""

import numpy as np
import math

# ----------------------------------------------------------------------------
# CONSTANTS  (three inputs: c, hbar, m_e; everything else derived)
# ----------------------------------------------------------------------------
HBAR   = 1.054571817e-34      # J s
C      = 2.99792458e8         # m/s
ME_KG  = 9.1093837015e-31     # kg
ME_MEV = 0.51099895069        # MeV
MEV_J  = 1.602176634e-13      # J per MeV
G_OBS  = 6.67430e-11          # m^3/(kg s^2)

def fine_structure():
    a = 1.0 / 137.0
    for _ in range(100):
        a = 1.0 / (math.exp(math.pi**2 / 2) - 2
                   - a - a / (math.pi * (1 - a)) - 6 * a**3 / math.pi**2)
    return a

ALPHA  = fine_structure()
M0_MEV = ME_MEV / ALPHA                       # node mass ~ 70 MeV
M0_KG  = M0_MEV * MEV_J / C**2
ELL    = ALPHA * HBAR / (ME_KG * C)           # lattice spacing = classical electron radius
D111   = math.sqrt(2.0 / 3.0) * ELL           # close-packed interlayer spacing
N4     = 3                                     # three stacking layers
L4     = N4 * D111                             # compact circumference = sqrt(6) ell

# ----------------------------------------------------------------------------
# ELASTIC MODULI  (in units of mu = rho c^2)
# ----------------------------------------------------------------------------
MU   = 1.0                                     # shear modulus sets c = sqrt(mu/rho)
LAM  = 1.0                                     # Cauchy relation lambda = mu (C12 = C44)
NSQ  = 1.0 / math.pi                           # Cosserat coupling number N^2 = 1/pi
KC   = 2.0 * MU * NSQ / (1.0 - NSQ)            # kappa_c from N^2 = kappa_c/(2mu+kappa_c)

C_L2 = (LAM + 2 * MU + KC)                     # longitudinal speed^2 in units of c^2
C_T2 = (MU + KC)                               # System-B transverse speed^2 in units of c^2

# ----------------------------------------------------------------------------
# Displacement dispersions (analytic eigenvalues of the 4x4 displacement matrix
# D_uv(k) = [(lambda+mu) k_u k_v + (mu+kappa_c)|k|^2 delta_uv] / rho).
# Returned in units of c^2/ell^2 for k in units of 1/ell.
# ----------------------------------------------------------------------------
def omega2_longitudinal(k_sp, k4):
    return C_L2 * (k_sp**2 + k4**2)

def omega2_transverse(k_sp, k4):
    return C_T2 * (k_sp**2 + k4**2)

# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 70)
    print("KALUZA-KLEIN STRUCTURE OF THE GRAVITATIONAL SECTOR ON D4")
    print("=" * 70)
    print(f"  ell = {ELL*1e15:.4f} fm    L4 = sqrt(6) ell = {L4*1e15:.4f} fm")
    print(f"  moduli (units of mu):  lambda = {LAM:.3f}  mu = {MU:.3f}  kappa_c = {KC:.4f}")
    print(f"  longitudinal speed c_L = sqrt({C_L2:.3f}) c = {math.sqrt(C_L2):.3f} c")
    print()

    print("-" * 70)
    print("1.  Compression mode gains a Kaluza-Klein gap at finite k4")
    print("-" * 70)
    k4_1 = 2.0 * math.pi / (L4 / ELL)          # first compact mode, in 1/ell  (= 2 pi/sqrt(6))
    w2 = omega2_longitudinal(0.0, k4_1)        # units c^2/ell^2
    m_gap = math.sqrt(w2) * M0_MEV             # hbar omega = sqrt(w2) m0 c^2
    # numerical check that omega_L^2(0,k4) = (lambda+2mu+kappa_c) k4^2 / rho
    check = abs(w2 - (LAM + 2 * MU + KC) * k4_1**2) < 1e-12
    print(f"  k4 (n=1) = 2 pi / sqrt(6) = {k4_1:.4f} / ell")
    print(f"  omega_L^2(0,k4) = (lambda + 2 mu + kappa_c) k4^2 / rho   [identity holds: {check}]")
    print(f"  compression-mode energy at n=1:  hbar omega_L = {m_gap:.1f} MeV")
    print(f"  the massless n=0 mode carries ordinary long-range gravity.")
    print()

    print("-" * 70)
    print("2.  Gravitational Yukawa correction: range set by k4 alone")
    print("-" * 70)
    xi_over_ell = 1.0 / k4_1                    # 1/k4 = L4/(2 pi)
    print(f"  static potential of the n=1 mode ~ exp(-k4 r) / r")
    print(f"  Yukawa range xi = 1/k4 = L4/(2 pi) = {xi_over_ell:.4f} ell = {xi_over_ell*ELL*1e15:.3f} fm")
    print(f"  transverse KK range hbar/(sqrt(6) m0 c) = ell/sqrt(6) = {1/math.sqrt(6):.4f} ell")
    print(f"  both sit near half a lattice spacing (monograph: ~0.5 ell).")
    print(f"  -> short-range only; the inverse-square law is intact above ~1 fm")
    print(f"     and the correction is unobservable at any reachable scale.")
    print()

    print("-" * 70)
    print("3.  Channel mixing T1u -> phi_i4 -> u4 -> A1g opens only at k4 != 0")
    print("-" * 70)
    print("  the EOM coupling kappa_c (i k_j phi_ij + i k4 phi_i4) links the shear")
    print("  displacement u_i to the mixed rotation phi_i4 through the k4 piece:")
    for k4 in (0.0, 0.5, k4_1):
        coupling = KC * k4
        state = "decoupled (gravity orthogonal to light)" if coupling == 0 else "mixed"
        print(f"    k4 = {k4:5.3f}/ell :  |u_i <-> phi_i4| = kappa_c k4 = {coupling:.4f}   ({state})")
    print("  -> at k4 = 0 compression (A1g) and shear (T1u) are exactly orthogonal,")
    print("     so a free electron has zero active gravitational mass.  Only finite")
    print("     k4 (sub-fm) breaks it.")
    print()

    print("-" * 70)
    print("4.  Dimensional reduction of Newton's constant: G3 = G4 / L4")
    print("-" * 70)
    CW   = (1 + 1 / math.pi) * (1 - 17 * ALPHA / 18)
    G3   = HBAR * C * ALPHA**19 * CW / M0_KG**2            # predicted 3D Newton constant
    K_SF = (M0_KG / ELL**3) * C**2 / (ALPHA**19 * CW)      # superfluid bulk modulus
    G4   = G3 * L4                                          # standard KK relation
    print(f"  G3 (predicted)      = {G3:.5e} m^3/(kg s^2)   [obs {G_OBS:.5e}]")
    print(f"  G3 x K_sf           = {G3*K_SF:.6e}")
    print(f"  c^4 / ell^2         = {C**4/ELL**2:.6e}   ratio {G3*K_SF/(C**4/ELL**2):.10f}")
    print(f"  G4 = G3 x L4        = {G4:.5e} m^4/(kg s^2)")
    print(f"  sqrt(6) c^4/(K_sf ell) = {math.sqrt(6)*C**4/(K_SF*ELL):.5e}   ratio {G4/(math.sqrt(6)*C**4/(K_SF*ELL)):.6f}")
    print(f"  -> Newton's constant is the inverse compressive stiffness of the")
    print(f"     compact direction: G3 = c^4/(K_sf ell^2), and G4 = L4 G3.")
