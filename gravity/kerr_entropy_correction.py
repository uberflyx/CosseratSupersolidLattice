"""
Pricing the Kerr interior's entropy correction.

The annealed-interior argument gives S = 0 exactly for a Schwarzschild
interior. A Kerr interior holds the spin-protected Feynman array and its
chemistry (sec:interior_knot_chemistry), and this script prices what those
degrees of freedom add to the entropy-area law, from framework constants.

Three contributions:
  1. ZERO MODES. The array's ground state breaks translation and rotation
     symmetry; distinct placements at core resolution xi give
     S_zero ~ k_B [ ln(b/xi)^2 + ln(2 pi R / 3 xi) ]: a log-of-area-FORM
     term, present only when J != 0. Contrast: the universal c ln(A)
     corrections of quantum-gravity approaches (Kaul-Majumdar 2000) apply
     to all holes; the framework's log term is spin-tagged.
  2. THERMAL GASES at the horizon temperature: the 2D Tkachenko phonon gas
     of the array (speed c_T = sqrt(hbar Omega / 4 m0), Tkachenko 1969)
     and the 1D Kelvin gas on each line.
  3. CHEMISTRY. Bead states are gapped at E_1 = 1.107 GeV >> k_B T_H, so
     the thermal population is exp(-1e40)-suppressed: exactly empty unless
     accretion-driven.

Fiducial hole: 10 Msun, a* = 0.98.
"""

import numpy as np

alpha = 1/137.035999177
m_e, c, G = 9.1093837015e-31, 2.99792458e8, 6.6743e-11
hbar, k_B = 1.054571817e-34, 1.380649e-23
m0  = m_e/alpha
ell = 2.8179403205e-15
kappa = 2*np.pi*ell*c
xi = ell
lP2 = hbar*G/c**3
Msun = 1.989e30

M, astar = 10*Msun, 0.98
rg  = G*M/c**2
rp  = rg*(1 + np.sqrt(1 - astar**2))
rm  = rg*(1 - np.sqrt(1 - astar**2))
ag  = astar*rg
A   = 4*np.pi*(rp**2 + ag**2)                 # horizon area
Om  = astar*rg*c/(rp**2 + ag**2)              # Omega_H
T_H = hbar*c*(rp - rm)/(4*np.pi*k_B*(rp**2 + ag**2))

n_v = 2*Om/kappa
b   = 1/np.sqrt(n_v)
N_v = n_v*np.pi*rp**2
S_BH = A/(4*lP2)                              # in units of k_B

line = "-"*70
print(line)
print(f"Fiducial: M = 10 Msun, a* = {astar};  Omega_H = {Om:.2e} rad/s,")
print(f"T_H = {T_H:.2e} K,  N_v = {N_v:.2e} lines,  b = {b*1e6:.1f} um")
print(f"S_BH = A/(4 lP^2) = {S_BH:.2e} k_B")
print(line)

# 1. Zero modes: translational placements at core resolution + orientation.
S_trans = np.log((b/xi)**2)
S_rot   = np.log(2*np.pi*rp/(3*xi))
S_zero  = S_trans + S_rot
print(f"1. Zero modes: S_zero = {S_zero:.1f} k_B "
      f"(trans {S_trans:.1f} + orient {S_rot:.1f}); spin-tagged log term")

# 2a. Tkachenko gas (2D, linear dispersion, cutoff k < pi/b).
c_T = np.sqrt(hbar*Om/(4*m0))
kk  = np.linspace(1e-12, np.pi/b, 200000)
x   = hbar*c_T*kk/(k_B*T_H)
occ = np.where(x < 500, x/np.expm1(np.minimum(x, 500)) -
               np.log1p(-np.exp(-np.minimum(x, 500))), 0.0)
S_Tk = (np.pi*rp**2)/(2*np.pi)*np.trapezoid(occ*kk, kk)
# 2b. Kelvin gas: 1D quadratic dispersion per line, length ~ 2 rp.
lnf = np.log(b/xi)
kkK = np.geomspace(np.pi/(2*rp), np.pi/xi, 400000)  # IR cutoff: line length
xK  = hbar*(kappa*kkK**2/(4*np.pi))*lnf/(k_B*T_H)
xKc = np.clip(xK, 1e-300, 500)
occK = np.where(xK < 500,
                xKc/np.expm1(xKc) + np.where(xKc < 1e-6, 1 - np.log(xKc),
                                             -np.log1p(-np.exp(-xKc))), 0.0)
S_K = N_v*(2*rp/np.pi)*np.trapezoid(occK, kkK)
print(f"2. Thermal gases at T_H: Tkachenko {S_Tk:.2e} k_B, "
      f"Kelvin {S_K:.2e} k_B  (c_T = {c_T:.3f} m/s)")

# 3. Chemistry: gapped.
print(f"3. Chemistry: gap E1/k_B T_H = {1.774e-10/(k_B*T_H):.2e} "
      f"-> thermal population exactly negligible")

S_int = S_zero + S_Tk + S_K
print(line)
print(f"TOTAL interior term: {S_int:.2e} k_B "
      f"=  {S_int/S_BH:.1e} of S_BH")
print("Scaling: S_zero ~ k_B ln(A/ell^2) * Theta(J); thermal terms ~ T_H^2,")
print("vanishing at extremality. The area law gains its first interior")
print("term and is numerically untouched.")
