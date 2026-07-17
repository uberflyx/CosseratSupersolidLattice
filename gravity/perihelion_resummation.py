"""Perihelion advance of Mercury: the trial of the second-order acoustoelastic law.

Three candidate completions of the acoustoelastic law share the first-order
index n = 1 + eps (eps = GM/rc^2) and part at second order:

  (ii) logarithmic (self-consistent): n = e^eps            -> beta = 1
       truncated index:               n = 1 + eps           -> beta = 3/2
  (i)  linear modulus law:            n = (1 - 2 eps)^-1/2  -> beta = 0

The orbit equation in the isotropic Gordon metric (A = 1/n^2, B = n^2) is
  u'' + u = (r_g c^2/L^2) + (C r_g^2 c^2/L^2) u,
with C read from the eps^2 coefficients of B/A and B:
  C = [eps^2 coeff of B/A] - [eps^2 coeff of B].
The perihelion advance per orbit is Delta phi = C pi r_g^2 c^2 / L^2
= (C/6) * 6 pi GM / (c^2 a (1-e^2)), and the PPN cross-check is
Delta ~ (2 + 2 gamma - beta)/3 with gamma = 1 for all three laws.

This script derives C symbolically for each law, evaluates Mercury, and
checks the PPN factors.
"""

import numpy as np
import sympy as sp

# --- symbolic: the orbit-equation coefficient for each law -------------------
eps = sp.symbols('epsilon', positive=True)
laws = {
    'exponential  (beta = 1)  ': sp.exp(eps),
    'truncated    (beta = 3/2)': 1 + eps,
    'linear law   (beta = 0)  ': (1 - 2*eps)**sp.Rational(-1, 2),
}
print("law                        :  [B/A]_2  [B]_2   C = difference   advance/GR")
results = {}
for name, n in laws.items():
    BA = sp.series(n**4, eps, 0, 3).removeO()       # B/A = n^4
    B  = sp.series(n**2, eps, 0, 3).removeO()       # B   = n^2
    c_BA = BA.coeff(eps, 2)
    c_B  = B.coeff(eps, 2)
    C = c_BA - c_B
    results[name] = C
    print(f"{name} :   {c_BA}       {c_B}      {C}            {sp.Rational(C,6)}")

# PPN cross-check: g00/c^2 = 1/n^2 = 1 - 2 eps + 2 beta eps^2 -> beta
print("\nPPN cross-check: beta from g00 = 1/n^2, advance factor (2 + 2*gamma - beta)/3, gamma = 1")
for name, n in laws.items():
    g00 = sp.series(n**-2, eps, 0, 3).removeO()
    beta = g00.coeff(eps, 2) / 2
    factor = sp.Rational(2 + 2 - 0, 1)  # placeholder to keep sympy exact below
    factor = (2 + 2*1 - beta) / 3
    print(f"{name} : beta = {beta},  (2+2g-b)/3 = {sp.nsimplify(factor)}  "
          f"(orbit route gave {sp.Rational(results[name],6)})")

# --- Mercury numbers ----------------------------------------------------------
c   = 2.99792458e8
GM  = 1.32712440018e20      # GM_sun, m^3/s^2
a   = 5.790905e10           # Mercury semi-major axis, m
e   = 0.205630
P   = 87.9691 * 86400.0     # orbital period, s
arcsec = np.degrees(1) * 3600

dphi = 6*np.pi*GM / (c**2 * a * (1 - e**2))          # rad per orbit, GR/exponential
orbits_per_century = 36525*86400.0 / P
per_century = dphi * orbits_per_century * arcsec
print(f"\nMercury: Delta phi = {dphi:.4e} rad/orbit; "
      f"{orbits_per_century:.1f} orbits/century")
print(f"exponential law : {per_century:.2f} arcsec/century  (observed anomalous advance)")
print(f"truncated index : {per_century*5/6:.2f} arcsec/century")
print(f"linear law      : {per_century*4/3:.2f} arcsec/century")

# exclusion level of the rejected laws from the LLR beta bound
beta_err = 5.6e-4    # Williams et al. 2009: beta - 1 = (-4.5 +/- 5.6) x 10^-4
print(f"\nLLR bound excludes beta = 3/2 at {0.5/beta_err:.0f} sigma "
      f"and beta = 0 at {1.0/beta_err:.0f} sigma")
