#!/usr/bin/env python3
r"""
condensate_momentum_routes.py
===============================
The three channels for exchanging momentum with the vacuum, priced side by side.

The framework's clock-engineering appendix establishes that the condensate's
chemical potential mu = mu_int(rho) + (1/2) m0 v_s^2 + V_ext is the clock
rate, so a gradient of mu is a gravitational field. That gives three
candidate channels, and this script prices each with the framework's own
numbers.

  ROUTE 1  (static compression): hold a density slope so that a body
           free-falls down it. Gibbs-Duhem at T = 0 (dP = n dmu,
           K = n dP/dn) gives dP = n mu0 eps whatever the bulk modulus.
           With the chemical potential the framework adopts, reading (B),
           mu0 = m0 v_p^2, the condensate is a stiff fluid and
           n mu0 = rho v_p^2 = K_sf, so one standard gravity across ten
           metres asks ~1.7e58 Pa, some 1e25 times the crystal's shear
           modulus. Reading (A), mu0 = m0 c^2, would have given
           rho c^2 eps ~ 5.5e17 Pa. The wall fails as well: ordinary
           matter is transparent to the condensate below the Landau
           velocity, and a conservative well left to itself is erased at
           the speed of first sound (~1e20 c). See the monograph appendix,
           sec-clock-routes and sec-clock-energy-theorem.

  ROUTE 2  (radiated second sound): emit momentum into the compression
           channel at v_2 = 1.59 c. Wave momentum flux F = P/v_2, i.e.
           v_2/c = 1.59 times WORSE thrust-per-watt than a photon
           rocket, before the transducer efficiency is even paid. A
           signalling channel, not a thruster.

  ROUTE 3  (ambient-propellant jet / "swim"): pump condensate through the
           body at slow speed w. The medium is the densest reservoir in
           physics (rho_s = f_s m0/ell^3, near-unity f_s), is on hand
           everywhere at zero carried mass, and below the Landau velocity
           flows without loss. Standard jet mechanics:
               thrust        F     = rho_s A w^2      (ambient intake)
               jet power     P_jet = (1/2) rho_s A w^3 = (1/2) F w
               figure of merit  F / P_jet = 2 / w
           Slow, massive jets are energy-cheap; this medium is the
           limiting case. All difficulty moves into ONE number: the
           actuation efficiency (the fluid-crystal cross-coupling), the
           chapter's named open quantity, formalised as the closure
           identity in knot_closure_theorem.py.

The point of the script is the comparison, not any one number. Route 1 fails on price
and containment together, Route 2 is a poor figure of merit, and Route 3 is
a transducer problem.
"""

import numpy as np

# ----- framework constants ------------------------------------------------
alpha  = 1.0 / 137.035999
m_e    = 9.1093837015e-31          # kg
c      = 2.99792458e8              # m/s
G      = 6.67430e-11               # m^3 kg^-1 s^-2
ell    = 2.8179403262e-15          # m   (lattice spacing = r_e)
f_s    = 0.9999                    # superfluid fraction, Leggett-saturated
                                    # (sec-fs-bootstrap; f_n < 1e-4)
m0     = m_e / alpha               # node mass, kg
rho    = m0 / ell**3               # node density = total mass density, kg/m^3
rho_s  = f_s * rho                 # superfluid density, kg/m^3
mubar  = rho * c**2                # true shear modulus = crystal energy density, Pa
N2     = 1.0 / np.pi               # Cosserat coupling number
# The contacts violate the Cauchy relation by C_12 - C_44 = -kappa_c/2, so the
# D_4 Born sums give K_cr/mubar = (5-8N^2)/[3(1-N^2)] = 1.200, not 5/3, and
# C_11/C_44 = (3-4N^2)/(1-N^2) = 2.533, not 3.
K_cr   = (5.0 - 8.0*N2)/(3.0*(1.0 - N2)) * mubar   # crystal bulk modulus, Pa
K_sf   = c**4 / (G * ell**2)       # condensate bulk modulus, Pa
gamma_G = 7.0 / 6.0                # lattice Grueneisen parameter
v2     = np.sqrt((3.0 - 4.0*N2)/(1.0 - N2)) * c    # second sound (sec-second-sound-speed)

print("=" * 72)
print("PUSHING ON THE VACUUM: three routes, priced")
print("=" * 72)
print(f"  rho_s = f_s m0/ell^3 = {rho_s:.3e} kg/m^3   (lead: 1.13e4; {rho_s/1.134e4:.2e}x denser)")
print(f"  mubar = rho c^2      = {mubar:.3e} Pa   (Route 1's price under reading (A))")
print(f"  K_sf  = c^4/(G l^2)  = {K_sf:.3e} Pa   (= rho v_p^2, Route 1's price under reading (B))")
print()

# ----- ROUTE 1: static compression, correctly priced -----------------------
a_want, L_craft = 9.8, 10.0
eps = a_want * L_craft / c**2               # fractional clock slope needed
dP   = K_sf * eps                           # price, reading (B): n mu0 = rho v_p^2 = K_sf
dP_A = mubar * eps                          # price, reading (A): n mu0 = rho c^2
u_dep = dP / gamma_G                        # energy density if DEPOSITED (Mie-Grueneisen)
print("ROUTE 1 - static compression ('warp the medium'):")
print(f"  1 g over {L_craft:.0f} m needs a fractional clock shift eps = {eps:.2e}")
print(f"  price, Gibbs-Duhem with mu0 = m0 v_p^2 (reading B): {dP:.1e} Pa")
print(f"  (reading A, mu0 = m0 c^2, would give {dP_A:.1e} Pa)")
print(f"  pressure / crystal bulk modulus = {dP/K_cr:.1e}: no strain could hold it")
print(f"  energy density if deposited: P/gamma_G = {u_dep:.1e} J/m^3")
print("  First sound also erases an unwalled well at ~1e20 c, and the")
print("  condensate is transparent to ordinary matter below the Landau speed.")
print("  SHUT on price and on containment.")
print()

# ----- ROUTE 2: radiated second sound --------------------------------------
print("ROUTE 2 - radiated second sound:")
print(f"  wave momentum: F/P = 1/v_2 = {1/v2:.3e} N/W   (v_2 = {v2/c:.2f} c)")
print(f"  photon rocket:  F/P = 1/c  = {1/c:.3e} N/W  (2nd sound is {v2/c:.2f}x worse,")
print("  before transducer efficiency). A signalling channel, not a thruster.")
print()

# ----- ROUTE 3: ambient-propellant jet -------------------------------------
print("ROUTE 3 - ambient-propellant jet ('swim in the medium'):")
print(f"  {'w [m/s]':>10s} {'F/P [N/W]':>12s} {'A for 1 MN [m^2]':>18s} {'P for 1 MN [W]':>16s}")
F_want = 1.0e6                              # reference force, ~100 t weight
for w in [1e-3, 1e-6, 1e-9]:
    FP = 2.0 / w
    A  = F_want / (rho_s * w**2)
    P  = 0.5 * F_want * w
    print(f"  {w:10.0e} {FP:12.1e} {A:18.3e} {P:16.2e}")
w_ref, A_ref = 1e-6, F_want / (rho_s * (1e-6)**2)
print(f"  At w = 1 um/s: {F_want:.0e} N of momentum flux for {0.5*F_want*w_ref:.2f} W,")
print("  which is why energy was never the obstacle. The grip is (see")
print("  knot_closure_theorem.py: closure denies it entirely).")
print(f"  through {A_ref:.0f} m^2 of effective intake.")
print(f"  Photon rocket for the same 1 MN: {F_want*c:.1e} W (three hundred TW).")
print(f"  Landau margin: w/v_2 = {w_ref/v2:.1e}  -> dissipationless by fifteen orders.")
print()
print(f"  Mass processed at 1 um/s, {A_ref:.0f} m^2: rho_s*A*w =",
      f"{rho_s*A_ref*w_ref:.2e} kg/s of AMBIENT condensate (not carried).")
print()
print("=" * 72)
print("THE CATCH, STATED HONESTLY")
print("=" * 72)
print("  The half watt prices the jet once it exists. The price of gripping")
print("  the condensate to make the jet is the actuation efficiency, set by")
print("  the fluid-crystal cross-coupling: the chapter's single most")
print("  consequential open number. The comparison's value is what it says")
print("  about the KIND of problem each route is: Route 1 fails on price")
print("  and containment, Route 2 on figure of merit, Route 3 on grip.")
print("=" * 72)
