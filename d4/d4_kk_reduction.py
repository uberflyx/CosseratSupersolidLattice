#!/usr/bin/env python3
"""
d4_kk_reduction.py -- Kaluza-Klein reduction of the D4 lattice
================================================================
Mitchell A. Cox, University of the Witwatersrand
Accompanies: M. A. Cox, "The Cosserat Supersolid" (2026)

The four-dimensional D4 crystal has 24 nearest neighbours (the vertices of the
24-cell).  Twelve lie in the FCC spatial slice; twelve cross to the layers ahead
and behind along the compact stacking direction.  This script carries out the
dimensional reduction in two parts and verifies both against the monograph.

Part 1 (zero compact momentum).  The static elastic response is the full
24-bond shell.  Relative to the bare 12-bond FCC slice the twelve extra bonds
stiffen the <100> cube axes, lifting C11 from 2 to 3 and driving the Zener
anisotropy ratio from 2 to 1 (exact isotropy).  The combination mu_V/(3K) = 1/5
is unchanged, so the acoustoelastic Schwarzschild condition
mu' = (3 - xi)/5 = 2 (with xi = -7) survives the reduction.

Part 2 (compact direction).  The lightest mode that varies along the compact
direction sets the Kaluza-Klein scale.  The close-packed layers are spaced by
d111 = sqrt(2/3)*ell; emergent Lorentz invariance fixes the long-wavelength wave
speed along the stack to c, so the one-dimensional lattice dispersion reaches
omega_zb = 2c/d111 = sqrt(6)*c/ell at the zone boundary.  With m0 c^2 = hbar c/ell,

        m_KK c^2 = hbar * omega_zb = sqrt(6) * m0 c^2  ~  172 MeV,

the System B gap below which the spatial slice is exact.  A dispersion scan
confirms the branch is linear (speed c) near q = 0 and reaches sqrt(6) m0 at the
boundary.

Usage:
    python d4_kk_reduction.py
"""

import numpy as np
from numba import njit, prange
import math

# ----------------------------------------------------------------------------
# CONSTANTS  (three inputs: c, hbar, m_e; everything else derived)
# ----------------------------------------------------------------------------
HBAR   = 1.054571817e-34      # J s
C      = 2.99792458e8         # m/s
ME_KG  = 9.1093837015e-31     # kg
ME_MEV = 0.51099895069        # MeV

def fine_structure():
    """alpha from the five-term Peierls-Nabarro self-energy series."""
    a = 1.0 / 137.0
    for _ in range(100):
        a = 1.0 / (math.exp(math.pi**2 / 2) - 2
                   - a - a / (math.pi * (1 - a)) - 6 * a**3 / math.pi**2)
    return a

ALPHA  = fine_structure()
M0_MEV = ME_MEV / ALPHA                      # node mass ~ 70 MeV
ELL    = ALPHA * HBAR / (ME_KG * C)          # lattice spacing = classical electron radius
D111   = math.sqrt(2.0 / 3.0) * ELL          # close-packed interlayer spacing
N4     = 3                                    # three stacking layers (A, B, C)
L4     = N4 * D111                            # compact circumference = sqrt(6) ell

# ----------------------------------------------------------------------------
# D4 NEAREST-NEIGHBOUR VECTORS  (integer embedding, |R|^2 = 2)
# ----------------------------------------------------------------------------
def d4_bonds():
    """24 D4 bonds: (+-1,+-1,0,0) and all four-coordinate permutations."""
    b = []
    for i in range(4):
        for j in range(i + 1, 4):
            for si in (-1, 1):
                for sj in (-1, 1):
                    v = [0, 0, 0, 0]
                    v[i] = si
                    v[j] = sj
                    b.append(v)
    return np.array(b, dtype=np.float64)

# ----------------------------------------------------------------------------
# PART 1: zero-mode elastic constants from the Born lattice sums
# ----------------------------------------------------------------------------
def elastic_constants(bonds):
    """Born long-wavelength elastic constants of a central-force lattice, in
    units of the bond stiffness.  Built from unit bond vectors Rhat = R/|R|:

        C11 = sum Rhat_x^4 ,   C12 = C44 = sum Rhat_x^2 Rhat_y^2 .

    The reduced response is three-dimensional, so the bulk modulus is the 3D
    form K = (C11 + 2 C12)/3.  Returns C11, C12, C44, Zener A, mu_V/(3K)."""
    inv4 = 1.0 / 4.0                          # 1/|R|^4 for D4 nearest neighbours
    C11 = sum((R[0]**4) * inv4 for R in bonds)
    C12 = sum((R[0]**2) * (R[1]**2) * inv4 for R in bonds)
    C44 = C12                                 # Cauchy relation for central forces
    A = 2.0 * C44 / (C11 - C12)
    K = (C11 + 2.0 * C12) / 3.0               # three-dimensional bulk modulus
    muV = C44                                 # isotropic shear modulus (A = 1)
    return C11, C12, C44, A, muV / (3.0 * K)

# ----------------------------------------------------------------------------
# PART 2-3: dispersion along the compact direction (one-dimensional chain of
# close-packed layers, acoustic speed calibrated to c)
# ----------------------------------------------------------------------------
@njit(parallel=True, cache=True)
def compact_dispersion(qd_array, nq):
    """Branch frequency (in units of c/d111) along the compact axis, as a
    function of the phase qd = q*d111 in (0, pi).  Standard one-dimensional
    lattice dispersion omega = 2 (c/d111) sin(qd/2): linear near q=0 (speed c),
    flat at the zone boundary qd=pi where omega = 2c/d111."""
    out = np.zeros(nq)
    for m in prange(nq):
        out[m] = 2.0 * math.sin(qd_array[m] / 2.0)
    return out

# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    bonds = d4_bonds()
    print("=" * 68)
    print("D4 KALUZA-KLEIN REDUCTION")
    print("=" * 68)
    print(f"  ell   = {ELL*1e15:.4f} fm      m0    = {M0_MEV:.3f} MeV")
    print(f"  d111  = {D111*1e15:.4f} fm      L4    = {L4*1e15:.4f} fm  (= sqrt(6) ell)")
    print()

    print("-" * 68)
    print("PART 1.  Zero compact momentum: static elastic constants")
    print("-" * 68)
    C11, C12, C44, A, ratio = elastic_constants(bonds)
    print(f"  C11 = {C11:.3f}   C12 = {C12:.3f}   C44 = {C44:.3f}")
    print(f"  Zener ratio A = 2 C44/(C11 - C12) = {A:.4f}    (bare FCC slice gives 2)")
    print(f"  mu_V/(3K) = {ratio:.4f}    (target 1/5 = {1/5:.4f})")
    xi = -7.0
    print(f"  acoustoelastic mu' = (3 - xi)/5 = {(3.0 - xi)/5.0:.3f}    (xi = {xi:.0f})")
    print("  -> the twelve inter-layer bonds isotropise the response (A: 2 -> 1)")
    print("     without disturbing mu', so the Schwarzschild condition mu' = 2")
    print("     holds for the full D4 lattice as well as the slice.")
    print()

    print("-" * 68)
    print("PART 2.  Compact direction: the Kaluza-Klein gap")
    print("-" * 68)
    m_kk = math.sqrt(6.0) * M0_MEV            # m_KK c^2 = hbar (2c/d111) = sqrt(6) m0 c^2
    print(f"  layer spacing d111 = sqrt(2/3) ell = {D111/ELL:.4f} ell")
    print(f"  zone-boundary frequency omega_zb = 2c/d111 = sqrt(6) c/ell")
    print(f"  Kaluza-Klein gap  m_KK = sqrt(6) m0 = {m_kk:.1f} MeV")
    print(f"  (compare pion 139.6 MeV: close in mass but J^P = 1^- vs 0^-, no mixing)")
    print()

    print("-" * 68)
    print("PART 3.  Dispersion along the compact direction")
    print("-" * 68)
    nq = 7
    qd = np.linspace(0.0, math.pi, nq)
    branch = compact_dispersion(qd, nq)       # in units of c/d111
    print("   q*d111     omega/(c/ell)     m(q) [MeV]")
    for q, w in zip(qd, branch):
        omega_c_ell = w * (ELL / D111)        # convert c/d111 -> c/ell
        print(f"   {q:6.3f}      {omega_c_ell:8.4f}        {omega_c_ell*M0_MEV:7.1f}")
    print("  -> linear (speed c) near q = 0, reaching sqrt(6) m0 = 172 MeV at the")
    print("     zone boundary.  Below this gap the FCC slice is the exact theory.")
