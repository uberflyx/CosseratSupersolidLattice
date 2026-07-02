#!/usr/bin/env python3
r"""
vortex_transit_phase_slip.py
==============================
PART 1 of the actuation question: does coherent circulation stir the vacuum
condensate, at what strength, and why has no laboratory ever noticed?

Three results, each a short calculation:

  A. N or N^2. Circulation is topological and adds as N with NO coherence
     required (Stokes). What coherence buys is alignment and rigidity:
     coherently MOVING aligned vortices radiate as N^2 (phased array).
     But spin PAIRING anti-aligns vortex lines pairwise: two opposite-spin
     protons carry the same charge (same screw handedness) yet opposite
     line orientations, so their far-field flow patterns cancel. Coherent
     amplitude goes as the NET POLARISATION, not the particle count.
     -> An even-even magic nucleus (J=0, fully paired) is the WORST case,
        not the best. The amplifier needs polarised ensembles.

  B. The transit theorem (Josephson-Anderson applied to matter).
     The framework's clock law hbar d(theta)/dt = -mu has Anderson's
     corollary: every vortex transit across a path slips the phase
     difference along that path by 2pi [Anderson, RMP 38, 298 (1966)].
     Electrons ARE vacuum-condensate vortices, so a spin-polarised current
     I is a vortex transit rate P*I/e, and the NAIVE slip-rate reading is
         delta-mu = h * (P I / e)   per polarised ampere
                  = 25.8 keV/A  ->  fractional clock step 3.7e-4 per A.
     Catastrophically excluded: no clock near any current has ever moved.

  C. The compactness protection (why B fails, and why that is good news).
     A particle's winding is CLOSED and COMPACT, confined to its core
     scale ell ~ fm. Two consequences, one for force and one for phase:
       (i)  closed loop in uniform relative wind: net Magnus force
            F = rho_s kappa (closed-loop integral of dl) x V = 0 exactly,
            because the loop integral of dl vanishes. Torque only.
            -> no net force on matter from the 370 km/s condensate wind,
               which is required by (and consistent with) Lorentz
               invariance of uniform motion.
       (ii) a transiting closed compact winding slips the phase only for
            point pairs its ~ell wake tube separates. Macroscopic clock
            pairs at separation L pick up only the dipole tail,
            suppressed by ~(ell/L)^2 ~ 1e-25 at a centimetre.
     Net predicted clock shift near a polarised amp: 3.7e-4 * 8e-26
     ~ 3e-29. Nothing. The SAME compactness that renormalises the
     electron's field energy caps its phase leverage.

  COROLLARY: pumping the macroscopic vacuum phase requires a vacuum
  vortex line that is macroscopically EXTENDED. Particles are not.
  The dark-sector knots are. (Taken up in knot_rotor_theorem.py.)
"""

import numpy as np

# framework constants
alpha = 1/137.035999
c     = 2.99792458e8
hbar  = 1.054571817e-34
h     = 2*np.pi*hbar
e     = 1.602176634e-19
ell   = 2.8179403262e-15
m0c2  = 0.51099895069e6*e/alpha        # node rest energy, J (70.03 MeV)
m0    = m0c2/c**2
f_s   = 4/5
rho_s = f_s*m0/ell**3
kappa = 2*np.pi*c*ell                  # h/m0, with hbar = m0 c ell

print("="*72)
print("A. N vs N^2, and the pairing catch")
print("="*72)
print("  circulation of N same-charge vortices: N*kappa  (topology; no")
print("  coherence needed). Coherent radiated power of N aligned movers:")
print("  ~N^2 within a coherence patch (phased array). BUT spin pairing")
print("  anti-aligns lines pairwise -> far-field flow cancels in pairs.")
print("  coherent amplitude ~ (net polarisation) x N, coupling ~ (P N)^2.")
print("  Even-even magic nuclei: J = 0, P = 0 -> worst case, not best.")
print("  The transducer candidates are POLARISED: hyperpolarised nuclear")
print("  ensembles, half-metal spin currents, polarised itinerant bands.")
print()

print("="*72)
print("B. The transit slip, taken naively at full strength")
print("="*72)
P_pol, I = 1.0, 1.0
slip_rate = P_pol*I/e                  # transits per second
dmu = h*slip_rate                      # J
print(f"  polarised current {I:.0f} A: transit rate = {slip_rate:.2e} /s")
print(f"  naive delta-mu = h I/e = {dmu:.3e} J = {dmu/e/1e3:.1f} keV")
print(f"  fractional clock step  = {dmu/m0c2:.2e}   <- absurd; excluded")
print()

print("="*72)
print("C. The compactness protection")
print("="*72)
print("  (i) closed compact loop, uniform wind V: net Magnus force")
print("      F = rho_s kappa (oint dl) x V = 0     [oint dl = 0 exactly]")
print("      -> zero net wind force on matter; torque and core stress only.")
L = 1e-2                               # clock separation, 1 cm
supp = (ell/L)**2
print(f"  (ii) transit slip confined to ~ell wake tube; macroscopic pair")
print(f"       at L = 1 cm sees dipole tail ~ (ell/L)^2 = {supp:.1e}")
print(f"  net predicted clock shift near a polarised ampere:")
print(f"       {dmu/m0c2:.1e} x {supp:.0e} = {dmu/m0c2*supp:.1e}   -> unobservable, consistent")
print()
print("  The corollary is the whole point: macroscopic phase pumping needs")
print("  a macroscopically EXTENDED vacuum vortex line. No particle is one.")
print("  The dark-sector knots are the only such objects in the framework.")
print("="*72)
