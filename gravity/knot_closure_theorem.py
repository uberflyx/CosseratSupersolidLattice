#!/usr/bin/env python3
r"""
knot_closure_theorem.py
=========================
PART 2 of the actuation question: can any topological object of the
condensate grip the medium, and what does the answer require of it?

Builds on vortex_transit_phase_slip.py, whose corollary was that pumping
the macroscopic vacuum phase needs a vacuum vortex line that is
macroscopically EXTENDED. This script asks whether the free sector's own
objects qualify, and finds that they do not. They are closed loops a few
lattice spacings across, so both of the following bite:

  A. THE CLOSURE IDENTITY (no net force). The Magnus force per unit
     length is f = rho_s kappa (v_L - v_s) x s_hat. Carry a CLOSED loop
     bodily through a flow that is uniform on the loop's own scale, so
     that (v_L - v_s) is one vector for the whole loop, and

         F = rho_s kappa (v_L - v_s) x (closed-loop integral of dl) = 0

     exactly, because the loop integral of dl vanishes. This is the same
     identity that shields ordinary matter from the 370 km/s condensate
     wind, where Lorentz invariance of uniform motion requires it. Read
     on a DRIVEN object it is just as exact: a closed loop dragged
     through a smooth condensate transmits no net momentum, at any speed
     and whatever it is made of. Only a SHEARED loop, one whose two
     sides move differently, delivers a force, and shearing needs two
     grips on an object larger than the gripping apparatus.

  B. THE APERTURE FACTOR (no macroscopic phase slip). The phase around a
     vortex loop is half the solid angle its spanning surface subtends,
     theta_s = Omega/2. A point passing THROUGH the loop's disc sees
     Omega swing from +2pi to -2pi and so gains exactly 2pi; a point
     passing outside sees Omega return to its start and gains nothing.
     Sweeping a loop of radius a across an aperture of area A therefore
     winds the aperture-averaged phase difference by

         delta-theta_eff = 2 pi (pi a^2 / A)

     per transit, not by 2pi. For an extended line spanning the aperture
     the ratio is one. For a core-scale loop across a square metre it is
     ~1e-29, the same wake-tube suppression that silences a transiting
     electron.

  C. WHAT THAT COSTS THE HOPEFUL DESIGN. Both factors are computed below
     for the framework's ground-state ring (R0 = 1.68 ell). Beating the
     aperture factor by number would take ~1e28 loops, tens of kilograms
     of forged GeV objects, and the closure identity would still deny
     each of them any grip. Growing one loop instead runs into the
     R-independent curvature pull 2 pi T ~ 37 kN: matter, gripping at a
     chemical bond per host lattice spacing, cannot supply it until the
     loop is already of order 100 m across.

  D. WHAT SURVIVES. The trapped-knot torque bound. A closed loop in a
     uniform wind feels no force but does feel a torque,
     tau = rho_s kappa A n_hat x V, exactly as a magnetic dipole does in
     a uniform field. For one minimal knot that is ~1e-12 N m, five
     orders above the floor of precision torsion instruments, swinging
     with the sidereal day. Nothing of the kind is seen. CAVEAT, and it
     matters: a dipole free to rotate ALIGNS and then feels no torque,
     so the bound is firm only for orientation-locked pinning and weak
     for a knot free to swivel on its pinning site.
"""

import numpy as np

# ---------------------------------------------------------------- constants
alpha = 1/137.035999177
c     = 2.99792458e8
hbar  = 1.054571817e-34
h     = 2*np.pi*hbar
e     = 1.602176634e-19
ell   = 2.8179403205e-15          # lattice spacing = r_e [m]
m0c2  = 0.51099895069e6*e/alpha   # node rest energy [J]  (70.03 MeV)
m0    = m0c2/c**2                 # node mass [kg]
f_s   = 5/6                       # bootstrap superfluid fraction
rho_s = f_s*m0/ell**3             # superfluid density [kg/m^3]
kappa = h/m0                      # circulation quantum [m^2/s]
V_cmb = 3.70e5                    # lab speed through the lattice frame [m/s]

# ground-state ring geometry from the vortex-hadron ladder
R0     = 1.68*ell                 # ring radius [m]
Lam    = 10*ell                   # minimal knot line length [m]
E1     = (8*np.pi**2/5)*m0c2      # ring energy [J]
T_line = E1/(2*np.pi*R0)          # line tension [N]

print("="*72)
print("Framework primitives")
print("="*72)
print(f"  m0 c^2      = {m0c2/e/1e6:.2f} MeV")
print(f"  rho_s       = {rho_s:.3e} kg/m^3      (f_s = 5/6)")
print(f"  kappa       = {kappa:.3e} m^2/s")
print(f"  R0          = {R0*1e15:.2f} fm ;  E1 = {E1/e/1e9:.3f} GeV")
print(f"  line tension T = {T_line/1e3:.2f} kN  ({T_line:.0f} J per metre of line)")
print()

print("="*72)
print("A. Closure identity: no net force on a closed loop")
print("="*72)
print("  F = rho_s kappa (v_L - v_s) x (loop integral of dl) = 0, exactly.")
print("  Same identity as a current loop in a uniform B: torque, no force.")
print("  Holds for the 370 km/s wind (passive) and for a carried loop")
print("  (driven). A sheared loop is the only exception:")
for dv in (1e-6, 1e-3, 1.0):
    print(f"    sheared segment, L = 1 m, dv = {dv:8.1e} m/s"
          f"  ->  F = {rho_s*kappa*1.0*dv:.3e} N")
print("  ...but shearing needs two grips on one loop, so the loop must be")
print("  larger than the apparatus. A 4.7 fm loop admits no such grip.")
print()

print("="*72)
print("B. Aperture factor: the phase slip a compact loop actually delivers")
print("="*72)
loop_area = np.pi*R0**2
for A, label in [(1e-4, "1 cm^2"), (1.0, "1 m^2"), (100.0, "100 m^2")]:
    print(f"  aperture {label:>7}: delta-theta_eff / 2pi = {loop_area/A:.2e}")
print(f"  (loop disc area pi R0^2 = {loop_area:.2e} m^2)")
print()

print("="*72)
print("C. What beating either factor would cost")
print("="*72)
A_target = 1.0
N_needed = A_target/loop_area
M_needed = N_needed*2.21e9*e/c**2
print(f"  pumping across {A_target:.0f} m^2 by number: N = {N_needed:.2e} loops")
print(f"     = {M_needed:.1f} kg of forged 2.21 GeV linked pairs,")
print(f"       and closure still denies every one of them a net grip.")
F_in = 2*np.pi*T_line
d_atom, F_bond = 2.5e-10, 1e-8            # host spacing [m], bond rupture [N]
R_cross = T_line*d_atom/F_bond
print(f"  growing one loop instead: curvature pull 2 pi T = {F_in/1e3:.1f} kN,")
print(f"     independent of R, against a matter ceiling (2 pi R/d) F_bond")
print(f"     that only overtakes it at R = T d / F_bond = {R_cross:.0f} m.")
print(f"     The path from the core scale is barred, not the destination.")
print()

print("="*72)
print("D. What survives: the trapped-knot sidereal torque, with its caveat")
print("="*72)
tau = rho_s*kappa*loop_area*V_cmb
print(f"  tau = rho_s kappa A |n_hat x V| = {tau:.2e} N m  (normal perpendicular to V)")
print(f"  torsion-instrument floor ~1e-17 N m  ->  margin {tau/1e-17:.0e}x")
print("  Not seen, so ordinary matter stores no orientation-locked knots.")
print("  CAVEAT: a loop free to swivel aligns n_hat with V and the torque")
print("  goes to zero, exactly as a magnetic dipole aligns with B. The null")
print("  therefore bounds locked storage firmly and free storage weakly.")
print("="*72)
