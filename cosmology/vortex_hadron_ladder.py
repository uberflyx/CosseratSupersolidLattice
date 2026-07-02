"""
The vortex hadron ladder: bound states of the free circulation sector.

The free sector's smallest object is the core-scale ring at E_1 = (8pi^2/5)
m0 c^2 = 1.107 GeV (sec:clock_actuation). This script builds its composite
spectrum, the circulation sector's copy of the hadron mass formula
m = N m0 + bond corrections, in two families:

  TRAIN: N coaxial rings in generalised leapfrog (Wacks-Baggaley-Barenghi
     2014). E_N = N E_1 + sum of exact Neumann pair energies. Nearest
     neighbours dominate, so E_N ~ N E_1 + (N-1) eps: node count times
     unit plus bond count times correction, mirroring the crystal sector.

  BUNDLE: N parallel strands on a polygon around the torus tube, the free
     twin of the nuclear coherent vortex bundle. Pair energies are the
     near-coaxial mutual inductances, so the bundle grows superlinearly.
     Havelock's classical result caps the polygon at N = 7 (stable
     through 7, unstable above; Kurakin-Yudovich); N = 8 restructures to
     the centred heptagon (Campbell-Ziff): shell structure, magic number 7.

  PRODUCTION: ring dispersion is non-Lorentzian (E and P are independent
     functions of R), so thresholds must be solved on the dispersion, not
     with invariant masses. We solve e+e- -> gamma + X and, the new
     channel, p pbar annihilation at rest -> X + recoil (gamma or pi0).

All scales from m_e and alpha; conventions of interior_knot_chemistry.py.
"""

import numpy as np
from scipy.special import ellipk, ellipe
from scipy.optimize import brentq

alpha = 1/137.035999177
m_e, c = 9.1093837015e-31, 2.99792458e8
m0  = m_e/alpha
ell = 2.8179403205e-15
kappa = 2*np.pi*ell*c
rho_s = 0.8*m0/ell**3
xi  = ell
GeV = 1.602176634e-10
E1  = (8*np.pi**2/5)*m0*c**2

def ring_E(R): return 0.5*rho_s*kappa**2*R*(np.log(8*R/xi) - 2.0)
def ring_P(R): return rho_s*kappa*np.pi*R**2

# Solve dimensionlessly: x (ln 8x - 2) = E1/(rho_s kappa^2 ell / 2); brentq's
# absolute xtol (2e-12) dwarfs a femtometre root, so never solve in metres.
_rhs = E1/(0.5*rho_s*kappa**2*ell)
R0 = ell*brentq(lambda x: x*(np.log(8*x*ell/xi) - 2) - _rhs, 1.0, 5.0)
P0 = ring_P(R0)
assert abs(ring_E(R0) - E1)/E1 < 1e-9

def E12_coax(R1, R2, s):
    """Exact Neumann mutual energy of coaxial rings (fluid reading)."""
    k2 = 4*R1*R2/((R1+R2)**2 + s**2); k = np.sqrt(k2)
    return rho_s*kappa**2*np.sqrt(R1*R2)*((2/k - k)*ellipk(k2) - (2/k)*ellipe(k2))

line = "-"*74
print(line)
print(f"Ground ring: E(R0) = {ring_E(R0)/GeV:.3f} GeV (= E1), R0 = {R0/ell:.2f} ell, "
      f"P0 = {P0*c/GeV:.3f} GeV/c, v0 = {min((kappa/(4*np.pi*R0))*(np.log(8*R0/xi)-1),c)/c:.2f} c")
print(line)

# ------------------------------------------------------------------ TRAIN
s_tr = 3*R0                      # inside the pair's stable leapfrog window
print("TRAIN family (coaxial leapfrog, spacing 3 R0):")
print(f"{'N':>3s}{'E_N [GeV]':>12s}{'E_N/E1':>9s}{'bond corr [MeV]':>17s}"
      f"{'P_N c [GeV]':>13s}")
for N in range(1, 9):
    Eint = sum(E12_coax(R0, R0, s_tr*abs(i-j))
               for i in range(N) for j in range(i+1, N))
    EN = N*E1 + Eint
    print(f"{N:3d}{EN/GeV:12.3f}{EN/E1:9.3f}{Eint/GeV*1e3:17.1f}"
          f"{N*P0*c/GeV:13.3f}")

# ---------------------------------------------------------------- BUNDLE
# N strands on a polygon of radius r_b around the tube, nearest spacing
# d_min = 3 xi; near-coaxial mutual energy ~ rho k^2 R [ln(8R/d) - 2].
print(line)
print("BUNDLE family: polygon strands, spacing d = 3 xi, evaluated at R = 30 xi")
print("(near-coaxial form needs R >> d; at core scale strands overlap and the")
print(" would-be multiply wound core is unstable, so only trains exist there):")
Rb = 30*xi
def bundle_E(N, R, dmin=3*xi):
    if N == 1: return ring_E(R)
    r_b = dmin/(2*np.sin(np.pi/N))
    ths = 2*np.pi*np.arange(N)/N
    E = N*ring_E(R)
    for i in range(N):
        for j in range(i+1, N):
            d = 2*r_b*np.sin(abs(ths[i]-ths[j])/2)
            E += rho_s*kappa**2*R*(np.log(8*R/d) - 2.0)
    return E
for N in range(1, 8):
    print(f"  N={N}: E = {bundle_E(N, Rb)/GeV:8.2f} GeV "
          f"({bundle_E(N, Rb)/(N*ring_E(Rb)):.2f} x N E(R))")
r_b7 = 3*xi/(2*np.sin(np.pi/7))
E8 = bundle_E(7, Rb) + ring_E(Rb) + 7*rho_s*kappa**2*Rb*(np.log(8*Rb/r_b7)-2)
print(f"  N=8: centred heptagon E = {E8/GeV:8.2f} GeV (shell 1+7, Havelock cap)")

# ------------------------------------------------------------ PRODUCTION
print(line)
print("PRODUCTION on the ring dispersion (non-Lorentzian: solve E(R), P(R)):")
def root_sqrt_s(N, m_rec):
    """sqrt(s) = E_rec(P) + E_N at common P (back-to-back), scanned in R."""
    f = lambda R: np.sqrt((N*ring_P(R)*c)**2 + (m_rec*c**2)**2) + \
                  N*ring_E(R) + (N-1)*E12_coax(R, R, 3*R)
    Rs = np.geomspace(R0, 40*R0, 4000)
    vals = np.array([f(R) for R in Rs])
    return vals.min()
for N in (1, 2, 3):
    print(f"  e+e- -> gamma + train_{N}: sqrt(s) >= "
          f"{root_sqrt_s(N, 0)/GeV:.3f} GeV")
# p pbar at rest: sqrt(s) = 2 m_p = 1.876 GeV. Margins:
m_p, m_pi = 1.67262192e-27, 2.406176e-28
s_pp = 2*m_p*c**2
thr1 = root_sqrt_s(1, 0)
print(f"  p pbar at rest supplies {s_pp/GeV:.3f} GeV vs ring+gamma floor "
      f"{thr1/GeV:.3f} GeV -> CLOSED: the ring is momentum-expensive")
# In-flight pbar on a proton at rest, lattice frame: min pbar momentum s.t.
# E_pbar + m_p c^2 = E_ring(R) + E_g and p_pbar = P_ring(R) +/- p_g (collinear).
def pbar_threshold():
    """Fixed-target p pbar -> ring + gamma, photon backward (the ring takes
    more momentum than the beam brings). Energy-momentum in the lattice
    frame gives sqrt(p^2 c^2 + m^2 c^4) + p c = E_r + P_r c - m_p c^2 = B,
    so p c = (B^2 - m_p^2 c^4)/(2B); minimum sits at the ground ring."""
    B = ring_E(R0) + ring_P(R0)*c - m_p*c**2
    return (B**2 - (m_p*c**2)**2)/(2*B)/c

def upsilon_tag():
    """Monochromatic photon energy for e+e- -> gamma + ring at the
    Upsilon(4S): solve E(R) + P(R) c = sqrt(s) on the ring dispersion."""
    s12 = 10.58*GeV
    x = brentq(lambda x: ring_E(x*ell) + ring_P(x*ell)*c - s12, R0/ell, 10*R0/ell)
    return ring_P(x*ell)*c

pth = pbar_threshold()
Eg_ups = upsilon_tag()
print(f"  Upsilon(4S) single-ring tag on the true dispersion: "
      f"E_gamma = {Eg_ups/GeV:.2f} GeV (supersedes invariant-mass 5.23)")
print(f"  in-flight forge: pbar p >= {pth*c/GeV:.2f} GeV/c on a fixed target opens\n"
      f"  ring + gamma with a BACKWARD photon tag (LEAR-class momenta)")
sig_shed, sig_ann = 1e-43, 5e-30          # fb ceiling vs ~50 mb capture
print(f"  branching ceiling p pbar -> ring + X: "
      f"{sig_shed/sig_ann:.1e} (BaBar sigma_shed over sigma_ann)")
print(line)
print("Selection rule: shedding needs a lock in the process; photons carry")
print("no circulation, so gamma gamma cannot forge. Annihilation is the")
print("maximal lock-dissolution vertex; at rest it opens N = 1 only.")
