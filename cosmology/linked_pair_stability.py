#!/usr/bin/env python3
r"""
linked_pair_stability.py
========================
The linked pair's stability gate, companion to the vortex-hadron-ladder
subsection. Three results:

  1. GEOMETRIC LEMMA (exact). In the symmetric Hopf configuration, two rings
     of radius R each threading the other through its centre, the squared
     strand separation is d^2 = R^2 [3 + 2 cos u - 2 cos t (1 + cos u)],
     minimised over t at cos t = 1 where it collapses to d^2 = R^2 for all u.
     The minimum strand separation IS the ring radius: the closest point of
     one ring to the other is the other's centre, which the first ring
     passes through. At the ground radius R0 = 1.68 ell against core reach
     xi ~ ell, the linked geometry itself holds the strands 0.68 ell apart.

  2. ENERGY LEDGER, at the adopted superfluid fraction f_s = 0.99932, where
     E1 = 2 pi^2 f_s m0 c^2 = 1.381 GeV. Reconnection conserves line length
     minus a one-core-length radiative toll (~131 MeV at the ground ring's
     line tension E1/(2 pi R0)). The merge channel (pair -> single loop of
     both circumferences, radius 2 R0) needs E(2R0) = 5.97 GeV against the
     2 E1 + 203 MeV = 2.97 GeV the pair owns: forbidden by 3.0 GeV. The
     exchange channel (pair -> two plain unlinked rings) releases the 203 MeV
     packing energy and pays the 131 MeV toll: net +72 MeV, EXOTHERMIC.
     Energy does not protect the pair. Every energy here scales as f_s, so
     every sign and ratio holds at any fraction in the window.

  3. DEFORMATION BARRIER. Closing the 0.68 ell gap over a chord of one to
     two core lengths adds delta^2/a of line, a tension cost of 30-60 MeV,
     against the 203 MeV internal (leapfrog) store. Barrier/store ~ 0.15-0.3:
     marginal, so survival is dynamical, not energetic. A Gross-Pitaevskii
     integration at ground-scale parameters decides it; the published
     untying simulations evolved links with decades of excess scale and do
     not transfer to a pair born at the bottom of its family.

Author: M. A. Cox.
"""

import numpy as np
from scipy.optimize import brentq
from scipy.special import ellipk, ellipe

E_CH = 1.602176634e-19            # elementary charge [C]
M0C2 = 0.51099895069e6 * 137.035999177   # node rest energy m_e c^2/alpha [eV]
F_S = 0.99932                     # the Leggett ceiling adopted as the working superfluid fraction
E1 = 2 * np.pi**2 * F_S * M0C2    # ground ring energy: line tension x one spacing [eV]
XI = 1.0                          # core reach [ell]
# Ground radius: E(R0) = E1 with core xi = ell reduces to x [ln(8x) - 2] = 1.
R0 = brentq(lambda x: x * (np.log(8 * x / XI) - 2) - 1, 1.0, 3.0)   # [ell]


def neumann_bracket(x1, x2, s):
    """sqrt(x1 x2) [(2/k - k) K - (2/k) E] for coaxial rings (scipy m = k^2)."""
    m = 4 * x1 * x2 / ((x1 + x2)**2 + s**2)
    k = np.sqrt(m)
    return np.sqrt(x1 * x2) * ((2 / k - k) * ellipk(m) - (2 / k) * ellipe(m))


# Pair packing energy at axial spacing 3 R0: rho_s kappa^2 ell = 2 E1 times
# the dimensionless Neumann bracket.
BOND = 2 * E1 * neumann_bracket(R0, R0, 3 * R0)   # [eV]
KT_CMB = 2.35e-4                  # CMB thermal quantum [eV]


def ring_energy_shape(R):
    """Thin-ring energy shape E(R) ~ R [ln(8R/xi) - 2], in units where E(R0)=E1."""
    return R * (np.log(8 * R / XI) - 2.0)


def geometric_lemma():
    """Numerically confirm d_min = R for the symmetric Hopf configuration."""
    print("=" * 72)
    print("1. GEOMETRIC LEMMA: minimum strand separation of the symmetric pair")
    print("=" * 72)
    R = 1.0
    t = np.linspace(0, 2 * np.pi, 3000)
    u = np.linspace(0, 2 * np.pi, 3000)
    T, U = np.meshgrid(t, u)
    d2 = ((R * np.cos(T) - R - R * np.cos(U)) ** 2
          + (R * np.sin(T)) ** 2 + (R * np.sin(U)) ** 2)
    print(f"  numerical min(d)/R = {np.sqrt(d2.min()):.8f}   (analytic: exactly 1)")
    print(f"  at R0 = {R0:.2f} ell vs core reach xi = {XI} ell: static margin "
          f"{R0 - XI:.2f} ell")
    print()


def energy_ledger():
    """The two unlinking channels: merge forbidden, exchange exothermic."""
    print("=" * 72)
    print("2. ENERGY LEDGER OF THE UNLINKING CHANNELS")
    print("=" * 72)
    t_line = E1 / (2 * np.pi * R0)          # line tension x ell [eV per core length]
    toll = t_line * XI                       # one core length radiated per event
    print(f"  f_s = {F_S}: E1 = {E1/1e9:.3f} GeV; line tension = {t_line/1e6:.0f} MeV/core length;"
          f" reconnection toll ~ {toll/1e6:.0f} MeV")
    e_merge = E1 * ring_energy_shape(2 * R0) / ring_energy_shape(R0)
    avail = 2 * E1 + BOND
    print(f"  MERGE  (pair -> one loop, radius 2R0): needs {e_merge/1e9:.2f} GeV, "
          f"owns {avail/1e9:.2f} GeV -> FORBIDDEN by {(e_merge-avail)/1e9:.2f} GeV")
    print(f"  EXCHANGE (pair -> two plain rings): releases {BOND/1e6:.0f} MeV, "
          f"pays {toll/1e6:.0f} MeV -> net +{(BOND-toll)/1e6:.0f} MeV EXOTHERMIC")
    print(f"  -> energy does not protect the pair.")
    print()
    return t_line


def deformation_barrier(t_line):
    """Line-tension cost of closing the geometric gap, vs the internal store."""
    print("=" * 72)
    print("3. DEFORMATION BARRIER vs INTERNAL STORE")
    print("=" * 72)
    delta = R0 - XI
    for a in (1.0, 2.0):
        dl = delta**2 / a
        print(f"  chord a = {a:.0f} ell: extra line ~ delta^2/a = {dl:.2f} ell "
              f"-> barrier ~ {t_line*dl/1e6:.0f} MeV")
    lo, hi = t_line * delta**2 / 2.0, t_line * delta**2 / 1.0
    print(f"  internal (leapfrog) store: {BOND/1e6:.0f} MeV; "
          f"barrier/store ~ {lo/BOND:.2f}-{hi/BOND:.2f} -> MARGINAL")
    print(f"  thermal activation: barrier/kT_CMB ~ {lo/KT_CMB:.0e} -> never;")
    print(f"  survival is decided by internal dynamics (GP at R0 = {R0:.2f} xi) or")
    print(f"  quantum tunnelling, not by temperature.")
    print()


if __name__ == "__main__":
    geometric_lemma()
    t_line = energy_ledger()
    deformation_barrier(t_line)
    print("Verdict: the pair is geometry-protected with an order-unity margin;")
    print("a ground-scale Gross-Pitaevskii run decides its lifetime.")
