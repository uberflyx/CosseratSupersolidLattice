#!/usr/bin/env python3
r"""
condensate_drive_figure_of_merit.py
=====================================
The three ways to push on the vacuum, priced side by side.

The framework's clock-engineering chapter establishes that the condensate's
chemical potential mu = mu_int(rho) + (1/2) m0 v_s^2 + V_ext is the clock
rate, so a gradient of mu is a gravitational field. That gives three
propulsion routes, and this script prices each with the framework's own
numbers:

  ROUTE 1  (static compression / "warp"): hold a density slope so the craft
           free-falls down it. Cost: the condensate bulk modulus
           K_sf = c^4/(G ell^2) ~ 1.5e73 Pa makes the required pressure
           ~1e58 Pa. Shut. (This reproduces the monograph's own no-go.)

  ROUTE 2  (radiated second sound): emit momentum into the compression
           channel at v_2 = 3.65 c. Wave momentum flux F = P/v_2, i.e.
           3.65x WORSE thrust-per-watt than a photon rocket, before the
           transducer efficiency is even paid. A signalling channel, not
           a drive.

  ROUTE 3  (ambient-propellant jet / "swim"): pump condensate through the
           craft at slow speed w. The medium is the densest propellant in
           physics (rho_s = f_s m0/ell^3 ~ 4.5e15 kg/m^3), is on hand
           everywhere at zero carried mass, and below the Landau velocity
           flows without loss. Standard jet mechanics:
               thrust        F     = rho_s A w^2      (ambient intake)
               jet power     P_jet = (1/2) rho_s A w^3 = (1/2) F w
               figure of merit  F / P_jet = 2 / w
           Slow, massive jets are energy-cheap; this medium is the
           limiting case. All difficulty moves into ONE number: the
           actuation efficiency (the fluid-crystal cross-coupling), the
           chapter's named open quantity.

The point of the script is the comparison, not any one number: the flow
route is not an energy problem (as rocketry is) but a transducer problem.
"""

import numpy as np

# ----- framework constants ------------------------------------------------
alpha  = 1.0 / 137.035999
m_e    = 9.1093837015e-31          # kg
c      = 2.99792458e8              # m/s
G      = 6.67430e-11               # m^3 kg^-1 s^-2
ell    = 2.8179403262e-15          # m   (lattice spacing = r_e)
f_s    = 4.0 / 5.0                 # superfluid fraction
m0     = m_e / alpha               # node mass, kg
rho_s  = f_s * m0 / ell**3         # superfluid density, kg/m^3
K_sf   = c**4 / (G * ell**2)       # condensate bulk modulus, Pa
v2     = 3.65 * c                  # second sound speed

print("=" * 72)
print("PUSHING ON THE VACUUM: three routes, priced")
print("=" * 72)
print(f"  rho_s = f_s m0/ell^3 = {rho_s:.3e} kg/m^3   (lead: 1.1e4; ~4e11 x denser)")
print(f"  K_sf  = c^4/(G l^2)  = {K_sf:.3e} Pa")
print()

# ----- ROUTE 1: static compression (the no-go, monograph's own) -----------
a_want, L_craft = 9.8, 10.0
eps = a_want * L_craft / c**2               # fractional clock slope needed
dP  = K_sf * eps                            # pressure to hold it
u   = 0.5 * K_sf * eps**2                   # energy density stored
print("ROUTE 1 - static compression ('warp the medium'):")
print(f"  1 g over {L_craft:.0f} m needs delta-rho/rho = {eps:.1e}")
print(f"  pressure to hold it: K_sf * eps = {dP:.1e} Pa      -> SHUT")
print(f"  stored energy density: (1/2) K_sf eps^2 = {u:.1e} J/m^3")
print("  (Same stiffness that makes gravity weak makes it undrivable.)")
print()

# ----- ROUTE 2: radiated second sound --------------------------------------
print("ROUTE 2 - radiated second sound:")
print(f"  wave momentum: F/P = 1/v_2 = {1/v2:.2e} N/W")
print(f"  photon rocket:  F/P = 1/c  = {1/c:.2e} N/W  (2nd sound is 3.65x worse,")
print("  before transducer efficiency). A signalling channel, not a drive.")
print()

# ----- ROUTE 3: ambient-propellant jet -------------------------------------
print("ROUTE 3 - ambient-propellant jet ('swim in the medium'):")
print(f"  {'w [m/s]':>10s} {'F/P [N/W]':>12s} {'A for 1 MN [m^2]':>18s} {'P for 1 MN [W]':>16s}")
F_want = 1.0e6                              # hover ~100 t
for w in [1e-3, 1e-6, 1e-9]:
    FP = 2.0 / w
    A  = F_want / (rho_s * w**2)
    P  = 0.5 * F_want * w
    print(f"  {w:10.0e} {FP:12.1e} {A:18.3e} {P:16.2e}")
print(f"  At w = 1 um/s: hover a 100 t craft on {0.5*F_want*1e-6:.2f} W of jet power")
print(f"  through ~{F_want/(rho_s*1e-12):.0f} m^2 of effective intake.")
print(f"  Photon rocket for the same 1 MN: {F_want*c:.1e} W (three hundred TW).")
print(f"  Landau margin: w/v_2 = {1e-6/v2:.1e}  -> dissipationless by 15 orders.")
print()
print("  Mass processed at 1 um/s, 220 m^2: rho_s*A*w =",
      f"{rho_s*220*1e-6:.2e} kg/s of AMBIENT condensate (not carried).")
print()
print("=" * 72)
print("THE CATCH, STATED HONESTLY")
print("=" * 72)
print("  The half watt prices the jet once it exists. The price of gripping")
print("  the condensate to make the jet is the actuation efficiency, set by")
print("  the fluid-crystal cross-coupling: the chapter's single most")
print("  consequential open number. The comparison's value is what it says")
print("  about the KIND of problem: not an energy problem, a transducer one.")
print("=" * 72)
