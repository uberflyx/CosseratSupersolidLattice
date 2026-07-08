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

  2. ENERGY LEDGER. Reconnection conserves line length minus a one-core-
     length radiative toll (~105 MeV at the ground ring's line tension).
     The merge channel (pair -> single loop of both circumferences, radius
     2 R0) needs E(2R0) = 4.77 GeV against the 2 E1 + 162 MeV = 2.37 GeV
     the pair owns: forbidden by 2.4 GeV. The exchange channel (pair -> two
     plain unlinked rings) releases the 162 MeV packing energy and pays the
     105 MeV toll: net +57 MeV, EXOTHERMIC. Energy does not protect the pair.

  3. DEFORMATION BARRIER. Closing the 0.68 ell gap over a chord of one to
     two core lengths adds delta^2/a of line, a tension cost of 24-48 MeV,
     against the 162 MeV internal (leapfrog) store. Barrier/store ~ 0.15-0.3:
     marginal, so survival is dynamical, not energetic. A Gross-Pitaevskii
     integration at ground-scale parameters decides it; the published
     untying simulations evolved links with decades of excess scale and do
     not transfer to a pair born at the bottom of its family.

Author: M. A. Cox.
"""

import numpy as np

E_CH = 1.602176634e-19            # elementary charge [C]
M0C2 = 70.03e6                    # node rest energy [eV]
E1 = (8 * np.pi**2 / 5) * M0C2    # ground ring energy [eV]
R0 = 1.68                         # ground ring radius [ell]
XI = 1.0                          # core reach [ell]
BOND = 162e6                      # pair packing energy at s = 3 R0 [eV]
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
    print(f"  at R0 = {R0} ell vs core reach xi = {XI} ell: static margin "
          f"{R0 - XI:.2f} ell")
    print()


def energy_ledger():
    """The two unlinking channels: merge forbidden, exchange exothermic."""
    print("=" * 72)
    print("2. ENERGY LEDGER OF THE UNLINKING CHANNELS")
    print("=" * 72)
    t_line = E1 / (2 * np.pi * R0)          # line tension x ell [eV per core length]
    toll = t_line * XI                       # one core length radiated per event
    print(f"  E1 = {E1/1e9:.3f} GeV; line tension = {t_line/1e6:.0f} MeV/core length;"
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
    print(f"  internal (leapfrog) store: {BOND/1e6:.0f} MeV; "
          f"barrier/store ~ 0.15-0.3 -> MARGINAL")
    print(f"  thermal activation: barrier/kT_CMB ~ {30e6/KT_CMB:.0e} -> never;")
    print(f"  survival is decided by internal dynamics (GP at R0 = 1.68 xi) or")
    print(f"  quantum tunnelling, not by temperature.")
    print()


if __name__ == "__main__":
    geometric_lemma()
    t_line = energy_ledger()
    deformation_barrier(t_line)
    print("Verdict: the pair is geometry-protected with an order-unity margin;")
    print("a ground-scale Gross-Pitaevskii run decides its lifetime.")
