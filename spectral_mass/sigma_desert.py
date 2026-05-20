#!/usr/bin/env python3
"""
sigma_desert.py
===============
The FCC symmetry desert, and the Sigma(1670) as a forced Kbar-Delta molecule.

Clean parity (a gerade / ungerade label) requires an inversion centre,
hence an O_h-symmetric cluster.  The FCC graph supplies inversion-symmetric
clusters only at the node counts 13, 19, 21, 27 and 43: the centre plus
shell, the second-neighbour Born closure, the dual-orbit void set, and the
third-shell closures.  Between 21 and 27 there is none.

For the negative-parity Sigma this leaves a gap.  The gerade modes of the
dual-orbit cluster (N = 21) top out near 1540 MeV (see sigma_neg.py), while
the parity-flipping bilayer floor sits at 25 m_0 - 100 m_e = 1699.5 MeV and
its stiff mass mode near 1800 MeV.  A negative-parity Sigma in the window
1540 to 1700 MeV therefore has no three-quark assignment and must be a
meson-baryon molecule, exactly as the N = 20 desert forces the Lambda(1405)
to be a Kbar-N molecule.

This closes the Sigma(1670) (J^P = 3/2-) as a Kbar-Delta bound state:

    N = N_Kbar + N_Delta = 7 + 17 = 24
    m_LO = 24 m_0 = 1680.6 MeV   (PDG 1675, +0.33%)

bound about 45 MeV below the Kbar-Delta threshold.  Quantum numbers follow
from the constituents: Kbar(0-) Delta(3/2+) in S-wave gives 3/2-, isospin
1/2 x 3/2 = 1 selects the Sigma, and strangeness -1 + 0 = -1 is correct.

The Sigma(1620) (1/2-) sits in the same gap but has no clean molecular
cluster: Kbar-N (I = 1) at N = 20 gives 1400 MeV (too light) and eta-Sigma
at N = 25 gives 1751 MeV (too heavy).  It is a one-to-two-star state and
remains open.

The desert is a prediction: the coupled-channel hadrons are those whose node
count lands near N = 20 to 26, namely the Lambda(1405), the negative-parity
Sigma, and the scalar nonet, while every clean cluster closure lies outside it.

Run:  python3 sigma_desert.py
"""

from hadron_spectral_mass import M_0, M_E, mass_from_lambda

# Node counts (cluster sizes) for the relevant ground states and mesons.
N_NODES = {'N': 13, 'Lambda': 16, 'Sigma': 17, 'Delta': 17,
           'Kbar': 7, 'pi': 2, 'eta': 8}

# PDG masses of the constituents, MeV, for threshold checks.
M_PDG = {'Kbar': 493.7, 'N': 938.9, 'Lambda': 1115.7, 'Sigma': 1193.2,
         'Delta': 1232.0, 'pi': 138.0, 'eta': 547.9}


def mass_LO(*parts):
    """Leading-order molecular mass: node counts add, lambda = 4."""
    return sum(N_NODES[p] for p in parts) * M_0


def bilayer_floor(n_ground):
    """Lightest possible mass on the parity-flipping bilayer (lambda -> 0)."""
    n = n_ground + 8
    return n * M_0 - n * 4 * M_E


def main():
    print("Inversion-symmetric (parity-carrying) FCC node counts: 13, 19, 21, 27, 43")
    print(f"  dual-orbit gerade ceiling          ~ 1540 MeV  (sigma_neg.py)")
    print(f"  Sigma bilayer floor (lambda -> 0)  = {bilayer_floor(N_NODES['Sigma']):.1f} MeV")
    print(f"  => the window 1540-1700 MeV is a cluster gap\n")

    print("Sigma(1670) 3/2- = Kbar-Delta molecule:")
    m = mass_LO('Kbar', 'Delta')
    threshold = M_PDG['Kbar'] + M_PDG['Delta']
    print(f"  N = {N_NODES['Kbar']} + {N_NODES['Delta']} = "
          f"{N_NODES['Kbar'] + N_NODES['Delta']},  m_LO = {m:.1f} MeV  "
          f"(PDG 1675, {100 * (m - 1675) / 1675:+.2f}%)")
    print(f"  bound {threshold - m:.0f} MeV below the Kbar-Delta threshold "
          f"({threshold:.0f} MeV)")
    print(f"  forced: bilayer too heavy, dual-orbit too light\n")

    print("Sigma(1620) 1/2- (open): no clean molecular cluster at its node count")
    print(f"  Kbar-N (I=1) N=20  = {mass_LO('Kbar', 'N'):.0f} MeV  (too light)")
    print(f"  eta-Sigma    N=25  = {mass_LO('eta', 'Sigma'):.0f} MeV  (too heavy)")


if __name__ == "__main__":
    main()
