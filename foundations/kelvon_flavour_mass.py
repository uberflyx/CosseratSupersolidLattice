"""The D_4 construction of the bound flavour mode: the muon as the n = +/-1
kelvon of the electron.

Outcome of the construction (supersedes the width-mode candidate):

STRUCTURE. The electron is a composite line: screw dislocation locked to a
condensate vortex, closed around the compact direction (circumference
L_4 = sqrt(6) ell, three sites). Its transverse translation branch at compact
wavenumber k_4 = 2 pi/L_4 is a helical Kelvin wave winding the ring once: it
carries compact quasi-momentum n = +/-1, the generation label. Three facts
make this branch the sharp candidate and retire the width/shape mode:
  (i)  the gyroscopic coefficient for line translation is rho_s * kappa
       EXACTLY (Thouless-Ao-Niu): no circulation-profile convention enters,
       unlike the shape-mode overlap that plagued the earlier candidate;
  (ii) the line tension is superfluid-dominated: rho_s kappa^2/4pi over
       mu_bar b^2/2pi = 16 pi^2 f_s / (something) = 31.6, so the dynamics is
       standard quantised-vortex hydrodynamics with kappa = 2 pi c ell, and
       crystal elasticity is a ~3% stiffness correction;
  (iii) at k_4 = 0 the same branch is the gapless Goldstone (no sharp e*),
       while at k_4 = +/-1 the crystal shear continuum is gapped at
       E_1 = 3 m_0 c^2/sqrt(2) = 148.546 MeV and condensate first sound is
       gapped far higher (v_p >> c), so a kelvon below E_1 is exactly stable
       and cannot decay by k_4 = 0 photon emission (mu -> e gamma forbidden).

DISPERSION. Hollow-core columnar vortex, m = 1 slow branch (Thomson 1880;
see e.g. Simula, "Vortex mass in a superfluid", Phys. Rev. A 97, 023609
(2018), arXiv:1704.08410, Eq. 14; retrograde character: Fetter,
cond-mat/0402208):
    omega_K = (kappa / 2 pi a^2) [ sqrt(1 + x K_0(x)/K_1(x)) - 1 ],  x = k a,
retrograde: each core element orbits against the circulation, so the bound
sense is locked to the charge sign. With kappa = 2 pi c ell and hbar = m_0 c
ell the mass is direct, no anchor freedom:
    hbar omega_K = m_0 c^2 (ell/a)^2 f(k_4 a),   k_4 ell = 2 pi/sqrt(6).

THE ONE REMAINING NUMBER. Everything is fixed except the effective core
parameter a. Derived candidates and results (MeV):
    a_L = (2 pi f_s)^{-1/2} ell = 0.4460  ->  124.0   (circulation column)
    w_par = 0.452 ell               ->  122.5   (PN width)
    one-site Wigner-Seitz disc, a = (sqrt(3)/2 pi)^{1/2} ell = 0.52504
                                    ->  106.48  (+0.77% on m_mu)
    a* = 0.52928 ell                ->  105.658 (exact m_mu)
All lie below the 148.5 MeV ceiling: bound and protected for every natural
radius. a*/a_L = 1.187, the size of the known Gross-Pitaevskii correction
by which a healing-core vortex's effective hydrodynamic core parameter
exceeds xi. Closure therefore needs ONE derivation: the effective core
parameter of the lattice-condensate vortex (Bogoliubov-de Gennes kelvon on
the actual condensate, or a discrete derivation of the one-site-cell
boundary condition). Percent-level corrections with known signs remain:
crystal elastic tension (+3% stiffness), PN valley pinning (+, small),
GP core softening (-).

DIVISION OF LABOUR WITH THE KOIDE MACHINERY. The kelvon sets the scale and
the protection (existence of one sharp ~100 MeV rung per compact winding);
the chiral Z_3 ring structure splits n = +1 from n = -1 (ABC vs ACB
handedness) and owns the precision triplet pattern. The reconciliation of
the two readings is the tau task.
"""

import numpy as np
from scipy.special import k0, k1
from scipy.optimize import brentq

def f_slow(x):
    """Thomson m=1 slow-branch form factor, x = k a."""
    return np.sqrt(1.0 + x*k0(x)/k1(x)) - 1.0

m0    = 0.51099895069*137.035999177     # MeV, CODATA 2022 (m_e c^2 / alpha)
m_mu  = 105.6583755                     # MeV, PDG (project file)
k4l   = 2*np.pi/np.sqrt(6.0)            # k_4 in units of 1/ell
E1    = 3*m0/np.sqrt(2.0)               # crystal shear rung at k_4 = 1 [MeV]

def kelvon_mass(a):
    """hbar omega_K in MeV for core parameter a in units of ell."""
    return m0*(1.0/a**2)*f_slow(k4l*a)

if __name__ == "__main__":
    print(f"k_4 = {k4l:.4f}/ell ; ceiling E_1 = {E1:.3f} MeV ; "
          f"target m_mu = {m_mu:.4f} MeV")
    a_star = brentq(lambda a: kelvon_mass(a)-m_mu, 0.3, 1.2)
    rows = [("a_L (circulation column)", (2*np.pi*0.8)**-0.5),
            ("w_par (PN width)",          0.452),
            ("one-site WS disc",          np.sqrt(np.sqrt(3)/(2*np.pi))),
            ("a* (solves m_mu)",          a_star),
            ("WS circumradius",           3**-0.5)]
    for name, a in rows:
        E = kelvon_mass(a)
        print(f"  {name:26s} a = {a:.5f} ell   "
              f"hbar omega_K = {E:8.3f} MeV   bound: {E < E1}")
    print(f"a*/a_L = {a_star*(2*np.pi*0.8)**0.5:.4f} "
          f"(GP core-parameter factor scale)")
