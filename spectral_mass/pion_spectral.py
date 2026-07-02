#!/usr/bin/env python3
"""
pion_spectral.py
================
The pion's spectral mass closure on the two-node cell pair, computed with
the same energy-derived Cosserat construction used for every other hadron
(cosserat_classifier.build_cosserat_matrix).

The coupling term is (alpha/2) |phi - (1/2) curl u|^2, a sum of squares, so
the dynamical matrix is positive semidefinite by construction: no eigenvalue
can be negative.  An earlier version of this script wired the u-phi coupling
in by hand as bare identity blocks; that construction is not derived from an
energy, is inconsistent with the matrix used for the proton, the rho, and
the Delta, and produced a spurious lambda = -2 triplet.  The consistent
construction removes it.

Results on the cell pair (two FCC nearest neighbours, one bond):
  - The full 12x12 spectrum is non-negative.
  - lambda = 2 is an exact eigenvector: the bond-stretch mode, the two
    nodes moving toward and away from each other along the bond.  This is
    the textbook two-body value 2K/mu at K = 1, mu = 1/2.
  - The antisymmetric relative microrotation along the bond (the naive
    pseudoscalar pattern) is an exact eigenvector at lambda = 3, one
    Cosserat unit above the stretch.
  - The pion's J^P = 0^- is carried by the Z_3 bond-field winding of the
    cell pair (a topological label), not by the parity of the mass mode.
    The mass mode is the stretch at lambda = 2.

Mass closure:
  m_pi = 2 m_0 - 2 (4 - 2) m_e = (2/alpha - 4) m_e = 138.007 MeV
against the PDG isospin average 138.039 MeV (residual -0.024%).
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cosserat_classifier import build_cosserat_matrix

# CODATA 2022
M_E = 0.51099895069          # MeV
INV_ALPHA = 137.035999177
M_0 = M_E * INV_ALPHA        # 70.0253 MeV

def main():
    # Two FCC nearest neighbours at unit separation (ELL = 1 convention:
    # displacements are measured in units of the bond length).
    coords = np.array([[0.0, 0.0, 0.0],
                       [1.0/np.sqrt(2.0), 1.0/np.sqrt(2.0), 0.0]])
    rhat = coords[1] - coords[0]
    rhat /= np.linalg.norm(rhat)

    M = build_cosserat_matrix(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    w, v = np.linalg.eigh(M)

    print("Cell-pair Cosserat spectrum (energy-derived construction):")
    for lam in sorted(set(np.round(w, 4))):
        mult = int(np.sum(np.abs(w - lam) < 5e-4))
        print(f"  lambda = {lam:7.4f}  (mult {mult})")
    print(f"  minimum eigenvalue = {w.min():.2e}  (positive semidefinite)")

    # Bond stretch: u antisymmetric along the bond, phi = 0.
    stretch = np.zeros(12)
    stretch[0:3] = rhat
    stretch[3:6] = -rhat
    stretch /= np.linalg.norm(stretch)

    # Antisymmetric relative microrotation along the bond, u = 0.
    twist = np.zeros(12)
    twist[6:9] = rhat
    twist[9:12] = -rhat
    twist /= np.linalg.norm(twist)

    lam_stretch = stretch @ M @ stretch
    lam_twist = twist @ M @ twist
    res_s = np.linalg.norm(M @ stretch - lam_stretch * stretch)
    res_t = np.linalg.norm(M @ twist - lam_twist * twist)
    print(f"\nBond stretch:      lambda = {lam_stretch:.4f}  "
          f"(eigenvector residual {res_s:.1e})")
    print(f"Antisym. rotation: lambda = {lam_twist:.4f}  "
          f"(eigenvector residual {res_t:.1e})")

    m_pi = 2*M_0 - 2*(4 - lam_stretch)*M_E
    pdg_iso = 138.0392
    print(f"\nMass closure on the stretch mode:")
    print(f"  m_pi = 2 m_0 - 2 (4 - {lam_stretch:.0f}) m_e = {m_pi:.4f} MeV")
    print(f"  PDG isospin average = {pdg_iso} MeV  "
          f"(residual {100*(m_pi-pdg_iso)/pdg_iso:+.3f}%)")

if __name__ == "__main__":
    main()
