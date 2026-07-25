#!/usr/bin/env python3
r"""
line_anchor_density.py
========================
How much force can ordinary matter take from a vacuum vortex line, per
metre of line?

The naive estimate treats the line as straight. A straight line through a
solid passes within a core radius ell of only n*pi*ell^2 nuclei per metre,
which for iron is about two: matter is empty at the femtometre scale, so a
straight line is barely anchored at all.

That estimate is wrong, because the line is flexible. It carries a tension
T = E1/(2 pi R0) ~ 6 kN, and each nucleus is a pinning well of depth
U ~ T * d_N (the length of core the line does not have to pay for while it
lies inside already-spoiled order). A flexible string in a landscape of
wells does not run straight; it detours to visit them. The line therefore
threads far more nuclei than geometry alone would give.

The optimum follows from one trade-off. Let the line visit a nucleus every
s metres, reaching a transverse distance delta to do so.

  COST   a detour of transverse reach delta over a span s adds a path
         length ~2 delta^2/s, costing  2 T delta^2 / s.
  BUDGET the visit is worth making only if that is less than U, so
             delta^2 <= U s / (2 T).
  SUPPLY within a tube of radius delta and length s there are n pi delta^2 s
         nuclei, and at least one is needed, so  delta^2 >= 1/(pi n s).

Eliminating delta gives the closest spacing the line can sustain,

    s_min = sqrt( 2 T / (pi n U) ),

and the anchor the line can hand to the host is one chemical bond per
visit, so the force per unit length of line is F_bond / s_min.

Every number here is a scoping estimate at the order-of-magnitude level.
The well depth assumes a vortex core costs no tension inside a nucleus,
which is a framework premise rather than a derived result.
"""

import numpy as np

# ------------------------------------------------------------- primitives
alpha = 1/137.035999177
c     = 2.99792458e8
hbar  = 1.054571817e-34
h     = 2*np.pi*hbar
e     = 1.602176634e-19
ell   = 2.8179403205e-15               # lattice spacing = r_e [m]
m0c2  = 0.51099895069e6*e/alpha        # node rest energy [J]
m0    = m0c2/c**2
f_s   = 5/6
rho_s = f_s*m0/ell**3                  # superfluid density [kg/m^3]
kappa = h/m0                           # circulation quantum [m^2/s]

R0     = 1.68*ell                      # ground-state ring radius [m]
E1     = (8*np.pi**2/5)*m0c2           # ring energy [J]
T_line = E1/(2*np.pi*R0)               # line tension [N] = energy per length


def nuclear_diameter(A):
    """Nuclear diameter from the standard radius formula, r = 1.2 A^(1/3) fm."""
    return 2*1.2e-15*A**(1/3)


def anchor_density(n, A, F_bond, label):
    """Pin spacing, transverse reach and force per metre for one host material.

    n       nuclear number density [m^-3]
    A       mass number
    F_bond  force to tear one host nucleus out of its lattice [N]
    """
    U = T_line*nuclear_diameter(A)                 # pinning well depth [J]
    s = np.sqrt(2*T_line/(np.pi*n*U))              # pin spacing [m]
    delta = 1/np.sqrt(np.pi*n*s)                   # transverse reach [m]
    f_anchor = F_bond/s                            # holdable force [N/m]
    straight = n*np.pi*ell**2                      # naive straight-line count [1/m]
    print(f"  {label}")
    print(f"    well depth U          = {U/e/1e6:8.0f} MeV")
    print(f"    pin spacing s         = {s*1e9:8.1f} nm      ({1/s:.2e} sites/m)")
    print(f"    transverse reach delta= {delta*1e12:8.2f} pm")
    print(f"    straight-line estimate= {straight:8.1f} sites/m  "
          f"(low by {1/(s*straight):.1e}x)")
    print(f"    force per metre       = {f_anchor:8.3f} N/m")
    return f_anchor


print("="*72)
print("Line properties")
print("="*72)
print(f"  tension T = {T_line/1e3:.2f} kN, so a metre of vacuum vortex line")
print(f"  stores {T_line:.0f} J and must be written at nuclear energy density")
print(f"  ({m0c2/e/1e6:.0f} MeV per {ell*1e15:.2f} fm) along its whole length.")
print()

print("="*72)
print("Anchor density: a flexible line seeking its wells")
print("="*72)
f_fe = anchor_density(8.5e28, 56, 1.0e-8, "iron  (bond ~10 nN)")
print()
f_w = anchor_density(6.3e28, 184, 1.4e-8, "tungsten (bond ~14 nN)")
print()

print("="*72)
print("What that buys, and what it costs to make")
print("="*72)
print("  A loop must be only PARTLY gripped, or the closure identity")
print("  cancels the force, so take half the loop embedded and half free.")
for M, label in [(1.0, "1 kg"), (1e5, "100 t")]:
    F = M*9.8
    L_grip = F/f_fe
    L_loop = 2*L_grip
    print(f"    {label:>6}: {F:.2e} N needs {L_grip:.1f} m gripped,"
          f" a loop of {L_loop:.1f} m")
    print(f"            line energy {L_loop*T_line/1e6:.2f} MJ; drag speed cap"
          f" u = {f_fe/(rho_s*kappa)*1e12:.2f} pm/s")
print()
u_cap = f_fe/(rho_s*kappa)
print(f"  The drag speed is capped at u = f_anchor/(rho_s kappa) = {u_cap:.2e} m/s,")
print(f"  because the Magnus load per unit length is rho_s kappa u.")
print(f"  Sweeping such a line across a 10 m aperture then takes {10/u_cap/3.15e7:.1e} yr,")
print(f"  so the phase-pump reading stays dead even with an extended line.")
print(f"  What survives is the Magnus reading: thrust = {f_fe:.2f} N per metre")
print(f"  of gripped line, at a power F*u = {9.8*u_cap:.1e} W per 9.8 N.")
print("="*72)
