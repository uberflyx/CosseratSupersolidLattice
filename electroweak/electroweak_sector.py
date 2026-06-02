#!/usr/bin/env python3
"""
electroweak_sector.py -- the electroweak observables from FCC Cosserat mechanics
================================================================================
Mitchell A. Cox, University of the Witwatersrand

The whole electroweak sector reduces to two derived quantities plus the
order-parameter chain.  The only physical inputs are the electron mass m_e
and the fine-structure constant alpha; everything else is FCC geometry.

Two pillars
-----------
  (1) The W mass is the evanescent gap of the Cosserat optical branch.
      The weak channel is gapped, so a disturbance below the gap does not
      propagate: it tunnels in the frequency domain and decays as e^{-kappa r},
      which is why the weak force is short ranged.  The optical mode must
      tunnel through 2*pi/alpha Peierls-Nabarro barriers per coherent
      wavelength, each of energy C_F * m_0, with C_F = (Nc^2-1)/(2 Nc) = 4/3
      the colour Casimir of the cell pair.  Hence

          M_W = (2 pi C_F / alpha) * m_0 = (8 pi / 3) * M_EW.

  (2) The Weinberg angle is the channel-mixing fraction of the Cosserat
      constitutive tensor.  Summing bond stiffness over the 12 nearest
      neighbours splits into a force-stress (hypercharge) sector and a
      couple-stress (weak-isospin) sector with zero cross term.  The
      rank-2 cell pair contributes 2 hypercharge channels out of the
      Nc^2 = 9 colour-squared channels at the node, so

          sin^2(theta_W) = 2 / Nc^2 = 2/9.

The Z mass then follows from the Higgs mechanism with no further input:
M_Z = M_W / cos(theta_W).  This replaces the legacy "M_Z = 3 pi M_EW with a
13 alpha / 12 zone-boundary Coulomb correction": that correction is identically
8/(3 sqrt 7) - 1 to the precision of alpha, i.e. it was a roundabout way of
writing cos^2(theta_W) = 7/9.

The Higgs is the longitudinal optical mode (the amplitude mode of the
crystalline order parameter); its mass is the Landau curvature at the
minimum, set by the 13-node critical nucleation cluster, m_H = 13 M_EW (1+alpha/2).
The quartic coupling is the ratio of the cluster's repulsive core to the
coordination shell's relaxation phase space, lambda = 13/(32 pi); the VEV,
Fermi constant and top mass follow downstream.  The top is a broken
nearest-neighbour bond of the condensate, m_t = (v/sqrt2)(1-alpha), the
(1-alpha) being its fermionic self-energy.

No node-counting prefactor is used as a derivation: 8 pi / 3 is the colour
Casimir times the tunnelling phase 2 pi / alpha, and 2/9 is a channel count.

Authors: M. Cox, with Claude (Anthropic).  License: MIT.
"""

import numpy as np


# --- inputs: the only two physical constants -----------------------------
def solve_alpha():
    """Fine-structure constant from the PN fixed point (Chapter on alpha)."""
    a = 1.0 / 137.0
    for _ in range(100):
        a = 1.0 / (np.exp(np.pi**2 / 2) - 2 - a - a / (np.pi * (1 - a))
                   - 6 * a**3 / np.pi**2)
    return a


ALPHA = solve_alpha()          # 1/alpha = 137.035999...
M_E = 0.51099895069            # electron mass [MeV] (CODATA 2022)

# --- derived scales -------------------------------------------------------
M0 = M_E / ALPHA               # screened node mass m_0 = m_e/alpha  [MeV]
MEW = M0 / ALPHA               # unscreened node mass M_EW = m_e/alpha^2 [MeV]

# --- FCC / Cosserat structural integers -----------------------------------
NC = 3                         # colour = three FCC stacking phases A,B,C
CF = (NC**2 - 1) / (2 * NC)    # colour Casimir of the cell pair = 4/3
SIN2_THETA_W = 2 / NC**2       # Weinberg angle: 2 hypercharge / Nc^2 channels
COS_THETA_W = np.sqrt(1 - SIN2_THETA_W)
N_HIGGS = 13                   # central node + 12 nearest neighbours (cuboctahedron)
Z_COORD = 12                   # kissing number in 3D
V_PAIR_HAT = 8 * np.pi / 3     # breathing-mode volume of the cell pair (2 * 4pi/3)


def w_mass_GeV():
    """W mass as the evanescent optical gap: 2 pi C_F m_0 / alpha = (8pi/3) M_EW."""
    return (2 * np.pi * CF / ALPHA) * M0 / 1000.0


def z_mass_GeV(mw=None):
    """Z mass from the Higgs mechanism: M_W / cos(theta_W)."""
    mw = w_mass_GeV() if mw is None else mw
    return mw / COS_THETA_W


def higgs_mass_GeV():
    """Higgs (longitudinal optical / Landau curvature): 13 M_EW (1 + alpha/2)."""
    return N_HIGGS * (MEW / 1000.0) * (1 + ALPHA / 2)


def quartic_coupling():
    """Higgs quartic: repulsive core N_H over relaxation phase space Z * V_pair."""
    return N_HIGGS / (Z_COORD * V_PAIR_HAT)        # = 13/(32 pi)


def vev_GeV(mh=None, lam=None):
    """Higgs VEV downstream of m_H and lambda: v = m_H / sqrt(2 lambda)."""
    mh = higgs_mass_GeV() if mh is None else mh
    lam = quartic_coupling() if lam is None else lam
    return mh / np.sqrt(2 * lam)


def top_mass_GeV(v=None):
    """Top quark: broken NN bond of the condensate with fermionic self-energy."""
    v = vev_GeV() if v is None else v
    return (v / np.sqrt(2)) * (1 - ALPHA)


def top_yukawa():
    """Top Yukawa coupling: y_t = sqrt(2) m_t / v = 1 - alpha."""
    return 1 - ALPHA


def fermi_constant_GeV2(v=None):
    """Fermi constant downstream of the VEV: G_F = 1/(sqrt2 v^2)."""
    v = vev_GeV() if v is None else v
    return 1.0 / (np.sqrt(2) * v**2)


def weak_coupling():
    """SU(2) gauge coupling g_w = e/sin(theta_W) = 3 e / sqrt(2)."""
    e = np.sqrt(4 * np.pi * ALPHA)
    return e / np.sqrt(SIN2_THETA_W)


def summary():
    """Print every electroweak observable against PDG/CODATA 2022."""
    obs = {                         # central value, label, predicted
        "M_W   [GeV]": (80.369, w_mass_GeV()),
        "M_Z   [GeV]": (91.188, z_mass_GeV()),
        "m_H   [GeV]": (125.20, higgs_mass_GeV()),
        "lambda":      (0.1293, quartic_coupling()),
        "v     [GeV]": (246.22, vev_GeV()),
        "m_t   [GeV]": (172.57, top_mass_GeV()),
        "y_t":         (0.9910, top_yukawa()),
        "G_F e-5":     (1.16638, fermi_constant_GeV2() * 1e5),
        "sin^2_thW":   (0.22305, SIN2_THETA_W),
    }
    print(f"alpha^-1 = {1/ALPHA:.6f}   m_0 = {M0:.4f} MeV   M_EW = {MEW/1000:.4f} GeV")
    print(f"C_F = {CF:.5f}   sin^2(theta_W) = 2/9 = {SIN2_THETA_W:.5f}   "
          f"cos(theta_W) = sqrt(7/9) = {COS_THETA_W:.5f}")
    print("-" * 60)
    print(f"{'observable':<14}{'predicted':>12}{'measured':>12}{'resid %':>10}")
    print("-" * 60)
    for k, (meas, pred) in obs.items():
        resid = (pred / meas - 1) * 100
        print(f"{k:<14}{pred:>12.5f}{meas:>12.5f}{resid:>+10.3f}")
    print("-" * 60)
    print(f"g_w = e/sin(theta_W) = 3e/sqrt2 = {weak_coupling():.4f}")


if __name__ == "__main__":
    summary()
