"""
The electron winding's condensate ring, built from spin quantisation, against
the radius that reproduces the Bohr magneton.

The locked-faces identity of the superconductivity appendix asks for the
effective hydrodynamic ring radius at which the condensate flow dipole of the
electron's winding carries the same dipole energy as the magnetic moment,

    rho_s d^2 = mu0 mu_B^2,   d = kappa pi a^2
    ->  a_req = (4 pi^3 f_s alpha)^(-1/4) l = 1.0254 l     (f_s = 0.99932)

This script builds the radius from the OTHER end, with no reference to the
magnetic moment, in four steps, and reports how far it lands from a_req:

  (1) LEMMA (exact): in a singly quantised vortex every element of the
      condensate carries angular momentum hbar per node mass m0, independent
      of radius, because r * v_theta = kappa/(2 pi) = hbar/m0.

  (2) SPIN-HALF (exact): S = hbar/2 fixes the circulating superfluid mass at
      M_circ = m0/2, half a node, so V = M_circ/rho_s = l^3/(2 f_s).
      rho_s = f_s m0/l^3 takes m0/l^3 as the total density.

  (3) SHAPE (geometric): a column one lattice spacing tall holding V has
      radius a_L = l/sqrt(2 pi f_s) = 0.399 l, against the dressed
      Peierls-Nabarro half-width w = W d = 0.452 l (W = ln(1/alpha)/(2 pi),
      d = l/sqrt3): 11.7% below it.  The two routes put the core at one
      scale; a uniform column of assumed height does not resolve one radius.

  (4) MOMENT (chain): the three-layer wrap (area x3, radius x sqrt3) and the
      both-helicities doubling (moment x2, radius x sqrt2) give
          a_ring = sqrt(6) a_L = sqrt(3/(pi f_s)) l = 0.978 l,
      4.7% short of a_req in radius.  The moment grows as the ring's area,
      so the ring carries (a_ring/a_req)^2 mu_B = sqrt(36 pi alpha/f_s) mu_B
      = 0.909 mu_B.  Exact closure needs f_s = 36 pi alpha = 0.825.

Also computed: the flow/microrotation share bracket, the closure in coupling
language, and the census over the three charged leptons.  A lepton of mass m
carries e hbar/(2m); the ring's moment grows as a^2, so the radius it needs
falls as m^(-1/2), while the chain returns one radius for all three.
"""

import numpy as np
import sympy as sp

# ---------------------------------------------------------------- constants
alpha = 1 / 137.035999177          # CODATA 2022
mmu_me = 206.7682827               # CODATA 2022 muon-electron mass ratio
mtau_me = 1776.93 / 0.51099895069  # PDG tau mass over CODATA 2022 m_e c^2 [MeV]
f_s = 0.99932                      # adopted operating point: Leggett ceiling
d = 1.0 / np.sqrt(3.0)             # glide period [l]
W = np.log(1.0 / alpha) / (2.0 * np.pi)
w_PN = W * d                       # dressed PN half-width [l]
d111 = np.sqrt(2.0 / 3.0)          # interplanar spacing [l]
L4 = np.sqrt(6.0)                  # three-layer compact circumference [l] (= 3 d111)

print("=== step 1-2: half a node (exact) ===")
print("angular momentum per unit condensate mass = kappa/2pi = hbar/m0")
print("S = hbar/2  =>  M_circ = m0/2 exactly")
Vol = 1.0 / (2.0 * f_s)            # [l^3]
print(f"V = M_circ/rho_s = l^3/(2 f_s) = {Vol:.5f} l^3  (f_s = {f_s})")

print("\n=== step 3: the column and the PN core ===")
a_L = np.sqrt(1.0 / (2.0 * np.pi * f_s))
print(f"column height l  ->  a_L = sqrt(1/(2 pi f_s)) = {a_L:.5f} l")
print(f"dressed PN half-width w = W d = {w_PN:.5f} l ; a_L/w = {a_L/w_PN:.4f} "
      f"({(1 - a_L/w_PN)*100:.1f}% below)")
print(f"bare closed form (pi/4) d = {np.pi/4*d:.5f} l ; a_L/((pi/4)d) = {a_L/(np.pi/4*d):.4f}")
h_at_w = Vol / (np.pi * w_PN**2)
print(f"inverted: a column of radius w needs height h = {h_at_w:.4f} l")
print(f"a_L = (pi/4) d exactly would need f_s = 24/pi^3 = {24/np.pi**3:.6f}")

print("\n=== step 4: the moment chain ===")
a1_derived = np.sqrt(3.0) * a_L                 # three stacked discs per wrap
a_ring = np.sqrt(2.0) * a1_derived              # both-helicities doubling
a1_req = (16.0 * np.pi**3 * f_s * alpha) ** (-0.25)
a_req = (4.0 * np.pi**3 * f_s * alpha) ** (-0.25)
print(f"a_1    derived = sqrt(3) a_L = {a1_derived:.5f} l   required {a1_req:.5f} l "
      f"({(a1_derived/a1_req - 1)*100:+.2f}%)")
print(f"a_ring derived = sqrt(6) a_L = {a_ring:.5f} l   required {a_req:.5f} l "
      f"({(a_ring/a_req - 1)*100:+.2f}%)")
moment_e = (a_ring / a_req) ** 2
print(f"ring moment = (a_ring/a_req)^2 mu_B = {moment_e:.5f} mu_B "
      f"(closed form sqrt(36 pi alpha/f_s) = {np.sqrt(36*np.pi*alpha/f_s):.5f})")

# closed forms via sympy
fs_s, al_s = sp.symbols('f_s alpha', positive=True)
a_der_s = sp.sqrt(3 / (sp.pi * fs_s))
a_req_s = (4 * sp.pi**3 * fs_s * al_s) ** sp.Rational(-1, 4)
print(f"\nclosed form: (a_ring/a_req)^2 = {sp.simplify((a_der_s/a_req_s)**2)}")
alpha_closure = sp.solve(sp.Eq(a_der_s, a_req_s), al_s)[0]
print(f"exact closure needs alpha = {sp.simplify(alpha_closure)}, i.e. f_s = 36 pi alpha = {36*np.pi*alpha:.5f}")
print(f"at f_s = {f_s}: closure would want alpha^-1 = 36 pi/f_s = {36*np.pi/f_s:.2f} "
      f"against 137.036; the adopted f_s overshoots closure by {(f_s/(36*np.pi*alpha) - 1)*100:.1f}%")

print("\n=== share bracket (flow vs microrotation) ===")
# If spin splits by energy share, the flow's share is
# E_sf/(E_sf + E_cr) = 4 pi^2 f_s / (1 + 4 pi^2 f_s); M_circ, hence a^2, scales with it.
share = 4 * np.pi**2 * f_s / (1 + 4 * np.pi**2 * f_s)
a_shared = a_ring * np.sqrt(share)
print(f"flow energy share = {share:.5f}")
print(f"a_ring (all spin in flow)  = {a_ring:.5f} l  ({(a_ring/a_req - 1)*100:+.2f}%)")
print(f"a_ring (energy-share spin) = {a_shared:.5f} l  ({(a_shared/a_req - 1)*100:+.2f}%)")
print(f"required                   = {a_req:.5f} l  "
      f"({(a_req/a_ring - 1)*100:.1f}% above the pure-flow end: outside the bracket)")

print("\n=== census over the charged leptons ===")
# Required radius: rho_s (kappa pi a^2)^2 = mu0 (e hbar/2m)^2  ->  a^4 = alpha m0^2 l^4/(4 pi^3 f_s m^2)
m, m0 = sp.symbols('m m_0', positive=True)
a4 = al_s * m0**2 / (4 * sp.pi**3 * fs_s * m**2)
print(f"d ln a_req / d ln m = {sp.simplify(sp.diff(sp.log(a4)/4, m) * m)}")
print(f"{'lepton':8s} {'required a [l]':>15s} {'chain a [l]':>12s} {'pred/obs moment':>16s}")
for name, r in [("electron", 1.0), ("muon", mmu_me), ("tau", mtau_me)]:
    a_need = a_req / np.sqrt(r)
    print(f"{name:8s} {a_need:15.5f} {a_ring:12.5f} {(a_ring/a_need)**2:16.4g}")
print(f"muon ratio across the window 0.97 .. {f_s}: "
      f"{mmu_me*np.sqrt(36*np.pi*alpha/0.97):.1f} .. {mmu_me*np.sqrt(36*np.pi*alpha/f_s):.1f}")
print(f"f_s that would bring the muon to unity: 36 pi alpha (m_mu/m_e)^2 = {36*np.pi*alpha*mmu_me**2:.3g}")
print(f"column height ratio that would retune the muon (moment ~ 1/h at fixed V): {mmu_me:.1f}")

print("\n=== w/d cross-check against the alpha-fixing width ===")
print(f"spin route: a_L/d = {a_L/d:.5f} ; alpha route: W = {W:.5f} ({(a_L/d/W - 1)*100:+.1f}%)")

print("\n=== L4 bookkeeping ===")
print(f"L4 = sqrt(6) l = {L4:.4f} l = 3 x d111 = {3*d111:.4f} l  (three layers per wrap)")
