"""
Pre-crystallisation scoping: the temperature ladder of the freezing vacuum.

The medium before crystallisation is a pure superfluid (f_s = 1) of
fermion-antifermion pairs, each pair a node of mass m0 = m_e/alpha.
This script orders the temperature scales of that epoch from framework
constants alone, and records the census factor that the circulation
quantum carries.

  (1) Bose condensation of the pair gas at lattice density n = 1/l^3
      (ideal-gas estimate; the real system is dense and interacting, so
      this is an order-of-magnitude anchor, not a transition temperature):
          k_B T_BEC = (2 pi / zeta(3/2)^(2/3)) (hbar^2/m0) n^(2/3)
                    = 3.31 m0 c^2  ~ 232 MeV        [hbar = m0 c l]
  (2) Debye scale of the crystal that forms (transverse speed c):
          k_B Theta_D = hbar c (6 pi^2 n)^(1/3) = (6 pi^2)^(1/3) m0 c^2
                      ~ 273 MeV
  (3) Crystallisation / deconfinement: T_c = 156 MeV (monograph, from
      the compact-direction bond count).
  (4) Chirality window floor: T_geom = 28.6 MeV.

  Required ordering for the monograph's narrative (liquid is already
  superfluid when it freezes): T_BEC > T_c.  Verified: 232 > 156.

  Census factor: kappa = h/m0 counts the condensed PAIR (the node), the
  level-0 analogue of the flux quantum h/2e counting the Cooper pair.
  An unpaired-constituent condensate (mass ~ m0/2 each) would carry
  kappa' = 2 h/m0 instead.
"""

import numpy as np
from scipy.special import zeta

alpha = 7.2973525643e-3
m0_MeV = 0.51099895069 / alpha        # node mass [MeV]

# (1) ideal-Bose condensation at n = 1/l^3, hbar = m0 c l
coef_BEC = 2 * np.pi / zeta(1.5) ** (2.0 / 3.0)
T_BEC = coef_BEC * m0_MeV

# (2) Debye estimate with transverse speed c
coef_D = (6 * np.pi**2) ** (1.0 / 3.0)
Theta_D = coef_D * m0_MeV

T_c = 156.1          # crystallisation / deconfinement [MeV] (monograph)
T_geom = 28.6        # chirality window floor [MeV] (monograph)

print("=== the temperature ladder of the freezing vacuum [MeV] ===")
print(f"Theta_D (crystal Debye scale)    = {coef_D:.3f} m0 = {Theta_D:7.1f}")
print(f"T_BEC   (pair-gas condensation)  = {coef_BEC:.3f} m0 = {T_BEC:7.1f}")
print(f"T_c     (crystallisation)        =            {T_c:7.1f}")
print(f"T_geom  (chirality window floor) =            {T_geom:7.1f}")
print(f"\nordering check T_BEC > T_c : {T_BEC:.0f} > {T_c:.0f}  "
      f"(ratio {T_BEC/T_c:.2f})  -> liquid superfluid before it freezes: OK")

print("\n=== census factor ===")
h = 6.62607015e-34
m0_kg = 9.1093837139e-31 / alpha
print(f"kappa (pair, measured/bootstrap) = h/m0   = {h/m0_kg:.4e} m^2/s")
print(f"kappa' (unpaired constituents)   = 2 h/m0 = {2*h/m0_kg:.4e} m^2/s")
print("the bootstrap value kappa = 2 pi c l is the PAIR census, the")
print("level-0 analogue of Phi_0 = h/2e reading charge 2e for the Cooper pair")

print("\n=== what stays open ===")
print("constituent identity, the level-0 pairing glue, and the node's")
print("pair-breaking gap 2*Delta_node: the framework's effective-theory floor")
