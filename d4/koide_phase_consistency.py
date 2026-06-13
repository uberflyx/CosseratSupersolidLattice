#!/usr/bin/env python3
"""
koide_phase_consistency.py - the Koide Hamiltonian's flat diagonal is forced by
the same mirror that makes the two Shockley routes degenerate
=========================================================================

Backs the consistency note in the monograph's "Koide formula from a hopping
Hamiltonian" section.

The charged leptons are the three eigenstates of a Z3 tight-binding (Koide)
Hamiltonian on the three stacking sites A, B, C of a {111} glide plane:

        | eps      -t e^{+i phi}  -t e^{-i phi} |
    H = | -t e^{-i phi}   eps     -t e^{+i phi} |
        | -t e^{+i phi}  -t e^{-i phi}   eps    |

with equal on-site energies eps on the diagonal and a Berry phase phi on the
off-diagonal hopping.  This is a Hermitian circulant, with eigenvalues

    E_n = eps - 2 t cos(phi + 2 pi n / 3),   n = 0, 1, 2,

and the identification sqrt(m_n) = E_n (a common positive scale).  The physical
branch requires all three eigenvalues positive.

Three results:

  1. THE LEPTONS.  In this convention the phase that reproduces the measured
     e, mu, tau is phi = 47.3 deg.  Note 47.3 = 180 - 132.7: the value 132.7 deg
     quoted in the text is the supplement, the same physics in the +cos
     convention (hopping +t).  With the -t hopping written in the Hamiltonian,
     132.7 deg would send one eigenvalue negative, so 47.3 deg is the value
     consistent with the equation as written.

  2. Q IS FIXED BY t/eps ALONE.  On the physical branch (all E_n > 0),
     Q = sum(m)/(sum sqrt m)^2 = 1/3 + (2/3)(t/eps)^2, so Q = 2/3 exactly at
     t/eps = 1/sqrt(2), for EVERY phase.  The phase only selects which masses
     sit on the Q = 2/3 locus.

  3. A STATIC B/C SPLITTING WOULD BREAK KOIDE.  The mirror that exchanges the
     B and C routes (the swap a<->b, an exact D4 automorphism; see
     bc_route_mirror.py) forces the three on-site energies equal.  If instead a
     static splitting delta were allowed between the B and C sites, the diagonal
     would read (eps, eps+delta, eps-delta) and Q drifts off 2/3 with delta.  So
     the flat diagonal is not merely permitted by the mirror; it is required by
     the Koide relation.

The handedness therefore cannot be a static site-energy splitting (mirror-
forbidden, and Koide-breaking); it can only be the Berry phase phi, carried by
the imaginary, time-reversal-odd compact-direction coupling.

Inputs: the measured charged-lepton masses, used only to read off phi and to
check the observed Q.  Nothing else.

Usage:    python3 koide_phase_consistency.py
Requires: numpy.

Author: Mitchell A. Cox
Date:   June 2026
"""

import numpy as np

# ------------------------------------------------------------------
# measured charged-lepton masses (PDG 2024), MeV/c^2
# ------------------------------------------------------------------
M_E   = 0.51099895069
M_MU  = 105.6583755
M_TAU = 1776.86


def koide_Q(masses):
    """Koide ratio Q = sum(m) / (sum sqrt(m))^2."""
    s = np.sqrt(np.asarray(masses, float))
    return np.sum(masses) / np.sum(s)**2


def hamiltonian(eps, t, phi_rad, diag=None):
    """Z3 circulant Koide Hamiltonian with hopping -t e^{i phi} around
    A -> B -> C -> A.  'diag' overrides the three on-site energies."""
    d = np.full(3, eps, float) if diag is None else np.asarray(diag, float)
    H = np.diag(d).astype(complex)
    hop = -t*np.exp(1j*phi_rad)
    for k in range(3):
        H[k, (k+1) % 3] = hop
        H[k, (k-1) % 3] = np.conj(hop)
    return H


def spectrum(eps, t, phi_rad, diag=None):
    """Ascending real eigenvalues; sqrt(m_n) = E_n on the physical branch."""
    return np.sort(np.linalg.eigvalsh(hamiltonian(eps, t, phi_rad, diag)))


def Q_from_spectrum(E):
    """Q built from eigenvalues, with m_n = E_n^2 (sqrt m_n = |E_n|)."""
    return np.sum(E**2) / np.sum(np.abs(E))**2


def main():
    bar = "=" * 70
    print(bar)
    print("  Koide Hamiltonian: flat diagonal, Berry phase, and Q = 2/3")
    print(bar)

    sqrt_m = np.sqrt([M_E, M_MU, M_TAU])
    eps = sqrt_m.mean()                 # the three cosines sum to zero -> mean = eps
    t = eps/np.sqrt(2)                  # the Koide condition t/eps = 1/sqrt(2)
    Q_obs = koide_Q([M_E, M_MU, M_TAU])
    print(f"\n  measured masses (MeV): e={M_E:.6f}, mu={M_MU:.4f}, tau={M_TAU:.2f}")
    print(f"  observed Koide Q = {Q_obs:.6f}   (2/3 = {2/3:.6f})")
    print(f"  on-site energy eps = mean(sqrt m) = {eps:.4f} sqrt(MeV);  t = {t:.4f}")

    # ---- 1. the phase that reproduces the leptons ----
    # electron is the smallest sqrt(m); place it at n = 0 and solve for phi.
    phi = np.degrees(np.arccos((eps - np.sqrt(M_E))/(2*t)))
    print(f"\n1. Phase reproducing the leptons (hopping -t, E_n = eps - 2t cos):")
    print(f"   phi = {phi:.1f} deg   [supplement 180-phi = {180-phi:.1f} deg is the"
          f" +cos-convention value quoted in the text]")
    pr = np.radians(phi)
    E = np.array([eps - 2*t*np.cos(pr + 2*np.pi*n/3) for n in range(3)])
    print(f"   {'n':>3} {'E_n = sqrt(m_n)':>16} {'m_n (MeV)':>14}")
    for n in range(3):
        print(f"   {n:3d} {E[n]:16.4f} {E[n]**2:14.4f}")
    print(f"   reproduces e, tau, mu to {np.max(np.abs(np.sort(E)-np.sort(sqrt_m))):.1e}"
          f" sqrt(MeV); all eigenvalues positive: {np.all(E > 0)}")

    # ---- 2. Q = 2/3 for any phase on the physical branch ----
    print("\n2. On the physical branch, Q = 2/3 for every phase (t = eps/sqrt2):")
    print(f"   {'phi(deg)':>9} {'eigenvalues':>26} {'Q':>10} {'all>0':>6}")
    for ph in (47.3, 50, 55, 60, 65, 70):
        Es = spectrum(eps, t, np.radians(ph))
        print(f"   {ph:9.1f} {str(np.round(Es, 3)):>26} {Q_from_spectrum(Es):10.6f}"
              f" {str(bool(np.all(Es > 0))):>6}")
    print("   -> Q = 2/3 exactly, set by t/eps alone; the phase only moves the masses.")

    # ---- 3. a static B/C splitting breaks Koide ----
    print("\n3. A static B/C site-energy splitting drives Q off 2/3:")
    print("   diagonal = (eps, eps+delta, eps-delta), phi = 60 deg (a margin phase)")
    print(f"   {'delta/eps':>10} {'Q':>12} {'Q - 2/3':>12}")
    for fr in (0.0, 0.01, 0.02, 0.05, 0.10, 0.20):
        d = fr*eps
        Es = spectrum(eps, t, np.radians(60.0), diag=[eps, eps + d, eps - d])
        Q = Q_from_spectrum(Es)
        print(f"   {fr:10.2f} {Q:12.6f} {Q - 2/3:+12.6f}")
    print("   -> the deviation grows with the splitting; the flat diagonal the")
    print("      mirror enforces is exactly the condition Q = 2/3 needs.")

    print("\n" + bar)
    print("  Flat diagonal (mirror-forced) + Berry phase (dynamic handedness)")
    print("  => Koide Q = 2/3, with the phase setting the masses.  A static")
    print("     B/C splitting is both mirror-forbidden and Koide-breaking.")
    print(bar)


if __name__ == '__main__':
    main()
