"""
Closing a_eff by spin quantisation of the winding's condensate flow.

The locked-faces identity of the superconductivity appendix requires the
electron winding's effective hydrodynamic ring radius to be

    a_req = (4 pi^3 f_s alpha)^(-1/4) l = 1.0841 l          (g = 2)
    a_1   = a_req / sqrt(2)            = 0.7665 l           (g = 1 shadow)

This script derives the radius from the OTHER end, with no reference to
the magnetic moment, in four steps:

  (1) LEMMA (exact): in a singly quantised vortex, every element of the
      condensate carries angular momentum hbar per node mass m0,
      independent of radius, because r * v_theta = kappa/(2 pi) = hbar/m0.

  (2) SPIN-HALF (exact): S = hbar/2 therefore fixes the effective
      circulating superfluid mass at M_circ = m0/2 exactly: half a node.
      The volume follows: V = M_circ/rho_s = l^3/(2 f_s) = 5 l^3/8.

  (3) SHAPE (geometric): a column one lattice spacing tall holding that
      volume has radius a_L = sqrt(1/(2 pi f_s)) l = 0.4460 l, landing on
      the Peierls-Nabarro core half-width w = 0.454 l to 1.8%. Two
      independent routes (elastic equilibrium; spin quantisation) to one
      core radius.

  (4) MOMENT (chain): one wrap of the compact circle L4 = sqrt(6) l
      crosses the three stacking layers of the ABC period, so the flow
      dipole collects three core discs, A -> 3A, radius x sqrt(3); the
      both-helicities coupling (the recorded g=1 -> g=2 doubling) gives
      moment x2, radius x sqrt(2). Total:

          a_eff(derived) = sqrt(6) a_L = sqrt(3/(pi f_s)) l = 1.0925 l

      against a_req = 1.0841 l: +0.78%, i.e. O(alpha) at the framework's
      stated resolution.

Also computed: the flow/microrotation share bracket, the equivalent
alpha-residual framing, and the height/width trade-off against w.
"""

import numpy as np
import sympy as sp

# ---------------------------------------------------------------- constants
alpha = 7.2973525643e-3
f_s = 4.0 / 5.0
w_PN = 0.454                      # PN core half-width [l]
d111 = np.sqrt(2.0 / 3.0)         # interplanar spacing [l]
L4 = np.sqrt(6.0)                 # compact circumference [l] (= 3 d111)

print("=== step 1-2: half a node (exact) ===")
print("angular momentum per unit condensate mass = kappa/2pi = hbar/m0")
print("S = hbar/2  =>  M_circ = m0/2 exactly")
Vol = 1.0 / (2.0 * f_s)           # [l^3]
print(f"V = M_circ/rho_s = l^3/(2 f_s) = {Vol:.4f} l^3  (= 5/8 at f_s = 4/5)")

print("\n=== step 3: the column and the PN core ===")
a_L = np.sqrt(1.0 / (2.0 * np.pi * f_s))
print(f"column height l  ->  a_L = sqrt(1/(2 pi f_s)) = {a_L:.5f} l")
print(f"PN core half-width w = {w_PN} l ; a_L/w = {a_L/w_PN:.4f}  ({(1-a_L/w_PN)*100:.1f}% below)")
h_at_w = Vol / (np.pi * w_PN**2)
print(f"inverted: column of radius w needs height h = {h_at_w:.4f} l")

print("\n=== step 4: the moment chain ===")
a1_derived = np.sqrt(3.0) * a_L                 # three stacked discs per wrap
a_eff_derived = np.sqrt(2.0) * a1_derived       # both-helicities doubling
a1_req = (16.0 * np.pi**3 * f_s * alpha) ** (-0.25)
a_req = (4.0 * np.pi**3 * f_s * alpha) ** (-0.25)
print(f"a_1 derived   = sqrt(3) a_L = {a1_derived:.5f} l   required {a1_req:.5f} l "
      f"({(a1_derived/a1_req-1)*100:+.2f}%)")
print(f"a_eff derived = sqrt(6) a_L = {a_eff_derived:.5f} l   required {a_req:.5f} l "
      f"({(a_eff_derived/a_req-1)*100:+.2f}%)")

# closed forms via sympy
fs_s, al_s = sp.symbols('f_s alpha', positive=True)
a_der_s = sp.sqrt(3 / (sp.pi * fs_s))
a_req_s = (4 * sp.pi**3 * fs_s * al_s) ** sp.Rational(-1, 4)
alpha_closure = sp.solve(sp.Eq(a_der_s, a_req_s), al_s)[0]
print(f"\nclosed form: a_eff(derived) = sqrt(3/(pi f_s)) l"
      f" = sqrt(15/(4 pi)) l at f_s = 4/5")
print(f"exact closure would need alpha = {sp.simplify(alpha_closure)}"
      f" = f_s/(36 pi) -> 1/(45 pi) at f_s = 4/5")
print(f"i.e. alpha^-1 = 45 pi = {45*np.pi:.3f} vs 137.036: "
      f"{(45*np.pi/137.036-1)*100:+.2f}% (the 0.78% radius residual, "
      f"quartically magnified)")

print("\n=== share bracket (flow vs microrotation) ===")
# If spin splits by energy share as in the proton partition, the flow's
# share is E_sf/(E_sf+E_cr) = 4 pi^2 f_s / (1 + 4 pi^2 f_s).
share = 4 * np.pi**2 * f_s / (1 + 4 * np.pi**2 * f_s)
a_eff_shared = a_eff_derived * np.sqrt(share)   # M_circ, hence a^2, scales
print(f"flow energy share = {share:.5f}")
print(f"a_eff (all spin in flow)   = {a_eff_derived:.5f} l  ({(a_eff_derived/a_req-1)*100:+.2f}%)")
print(f"a_eff (energy-share spin)  = {a_eff_shared:.5f} l  ({(a_eff_shared/a_req-1)*100:+.2f}%)")
print(f"required                   = {a_req:.5f} l  (inside the bracket, "
      f"{(a_req-a_eff_shared)/(a_eff_derived-a_eff_shared)*100:.0f}% of the way up)")

print("\n=== w/d cross-check against the alpha-fixing width ===")
d = 1.0 / np.sqrt(3.0)
print(f"spin route: w/d = a_L/d = {a_L/d:.5f} ; alpha route: 0.783 "
      f"({(a_L/d/0.783-1)*100:+.1f}%)")

print("\n=== L4 bookkeeping ===")
print(f"L4 = sqrt(6) l = {L4:.4f} l = 3 x d111 = {3*d111:.4f} l  (three layers per wrap)")
