#!/usr/bin/env python3
r"""
knot_rotor_theorem.py
=======================
PART 2 of the actuation question: the dark-sector knot as the full-strength
transducer, what it can and cannot do, and the bound that says matter holds
none today.

Builds on vortex_transit_phase_slip.py, whose corollary was: pumping the
macroscopic vacuum phase needs a macroscopically EXTENDED vacuum vortex
line, and the only such objects are the dark-sector knots (closed vacuum
vortex loops, helicity-protected, from the Kibble-Zurek tangle).

  A. NO-SAIL THEOREM. A knot is a CLOSED loop, so in the uniform 370 km/s
     condensate wind the net Magnus force vanishes exactly:
         F = rho_s kappa (oint dl) x V = 0.
     Same mathematics as a current loop in a uniform magnetic field.
     The wind delivers torque and internal stress, never net thrust.
     -> there is no vacuum sail. Any drive must PUMP (make its own
        nonuniform flow); it cannot ride the ambient wind for free.

  B. WHAT A HELD KNOT FEELS ANYWAY. Per unit line length the local Magnus
     load in the wind is f = rho_s kappa V ~ 8.8e15 N/m: about 25 N per
     lattice site, against a pinning strength ~ m0 c^2/ell ~ 4 kN per
     site. A pinned knot survives the wind (stress ratio ~1e-2) and
     transmits pure torque.

  C. THE TRAPPED-KNOT BOUND. That torque swings with the sidereal day as
     the lab rotates in the wind. One minimal knot (line length ~10 ell)
     trapped in a test mass gives tau ~ rho_s kappa V Lambda R_knot
     ~ 1e-12 N m: five orders above the ~1e-17 N m torque floor of
     precision torsion instruments. Nothing of the kind is seen, so
     ordinary matter holds NO trapped vacuum vortices: passive capture
     does not stockpile knots. (Consistent with the topology: plain loops
     fuse with our dislocations and are destroyed; knots are helicity-
     protected against fusion and scatter instead of sticking.)

  D. THE ROTOR (the drive that survives). Carry trapped knots around a
     circle at frequency f_r: each circuit is one transit of every
     meridian, so the phase between the linked sides winds continuously,
         delta-mu = N_knots * h * f_r        (Josephson-Anderson),
     a sustained chemical-potential step: a clock-rate step: an
     artificial gravitational potential, made by machinery.
     One standard gravity across a 10 m craft needs
         N_knots * f_r ~ 1.9e7 /s,
     e.g. ~1200 knots on a rotor rim at 16 kHz (tip speed ~1 km/s,
     inside material limits), with a Magnus load of ~0.7 N per knot:
     ordinary engineering loads. The rotor also IS the pump of the
     ambient-propellant jet (F/P = 2/w economics, priced in
     condensate_drive_figure_of_merit.py).

  E. WHAT REMAINS OPEN, honestly: supply. Local knot abundance (the
     cascade fraction of the mixed dark sector) and a trap that holds a
     knot without the fusion channel are the two unsolved steps. The
     torsion-quiet of ordinary matter says nature does not do the
     trapping for us.
"""

import numpy as np

# framework constants
alpha = 1/137.035999
c     = 2.99792458e8
hbar  = 1.054571817e-34
h     = 2*np.pi*hbar
e     = 1.602176634e-19
ell   = 2.8179403262e-15
m0c2  = 0.51099895069e6*e/alpha
m0    = m0c2/c**2
f_s   = 4/5
rho_s = f_s*m0/ell**3
kappa = 2*np.pi*c*ell
V_cmb = 3.70e5                       # lab speed through the lattice frame

print("="*72)
print("A. No-sail theorem")
print("="*72)
print("  closed loop, uniform wind: F = rho_s kappa (oint dl) x V = 0.")
print("  (Identical to a current loop in uniform B: torque, no net force.)")
print("  -> no free ride on the 370 km/s wind. Drives must pump.")
print()

print("="*72)
print("B. Local loads on a held knot in the wind")
print("="*72)
f_line = rho_s*kappa*V_cmb
per_site = f_line*ell
pin_site = m0c2/ell
print(f"  Magnus per length  rho_s kappa V = {f_line:.2e} N/m")
print(f"  per lattice site   = {per_site:.1f} N   vs pinning ~ m0c^2/ell = {pin_site:.0f} N")
print(f"  stress ratio {per_site/pin_site:.1e}: a pinned knot survives, transmitting torque.")
print()

print("="*72)
print("C. Trapped-knot bound (why matter holds none)")
print("="*72)
Lam_min = 10*ell                      # minimal ring, circumference ~10 ell
R_knot  = 2*ell
tau = f_line*Lam_min*R_knot
print(f"  minimal knot: line {Lam_min:.1e} m, size {R_knot:.1e} m")
print(f"  sidereal torque tau ~ rho_s kappa V * Lambda * R = {tau:.1e} N m")
print(f"  torsion-instrument floor ~1e-17 N m  ->  margin {tau/1e-17:.0e}x")
print("  Not seen -> zero trapped knots per test mass. Passive capture fails;")
print("  loops fuse and die on our dislocations, knots bounce (helicity).")
print()

print("="*72)
print("D. The rotor: artificial potential by phase pumping")
print("="*72)
a_want, L_craft = 9.8, 10.0
dmu_need = (a_want*L_craft/c**2)*m0c2         # J across the craft
Nf_need  = dmu_need/h
print(f"  1 g across {L_craft:.0f} m: delta-mu = {dmu_need:.2e} J = {dmu_need/e:.2e} eV")
print(f"  requires N_knots * f_rotor = delta-mu/h = {Nf_need:.2e} /s")
r_rotor, u_tip = 1e-2, 1.0e3                  # 1 cm rim, 1 km/s tip
f_r = u_tip/(2*np.pi*r_rotor)
N_k = Nf_need/f_r
F_per_knot = rho_s*kappa*u_tip*Lam_min
print(f"  steel-class rotor: r = {r_rotor*100:.0f} cm, tip {u_tip:.0f} m/s -> f = {f_r:.2e} Hz")
print(f"  knots needed N = {N_k:.0f};  Magnus load per knot = {F_per_knot:.2f} N")
print("  Ordinary engineering loads. The same rotor is the pump of the")
print("  ambient-propellant jet (F/P = 2/w).")
print()
print("="*72)
print("E. Open: supply (local knot abundance; a fusion-free trap).")
print("="*72)
