#!/usr/bin/env python3
r"""
knot_trapping_candidates.py
=============================
PART 3 of the actuation question: can anything HOLD a dark-sector knot
without offering it the reconnection (fusion) channel that destroys plain
loops? The rotor drive of knot_rotor_theorem.py needs a trap; the torsion-
balance bound says nature does not provide one passively.

A quantised vortex is attracted to a region where the condensate is
DEPLETED, because its core is itself a density hole: sitting the core in a
pre-existing hole removes a length of core and lowers the energy. This is
the standard pinning of vortices to low-density paths and to impurities in
He-II and in BECs (electron bubbles, hydrogen tracers, droplet gaps). The
binding energy is ~ (core-energy density) x (overlap volume).

Two candidate "holes" exist in the vacuum's dictionary. We test both.

  CANDIDATE 1 - the VACANCY (dark-matter point defect). RULED OUT, twice:
    * WRONG SIGN. A vacancy's 12 neighbours relax INWARD (the medium is
      near-incompressible), so it is a compression MAXIMUM, not a density
      hole. A vortex is REPELLED, not trapped.
    * DECOUPLED. The framework already shows vacancies barely touch the
      condensate (vacancy-mediated superfluidity is ruled out; the
      condensate rearranges at FIXED number around a vacancy). No grip.

  CANDIDATE 2 - the STACKING FAULT (the ribbon a partial/quark drags).
    A fault is a genuine 2D sheet of DISRUPTED stacking, ~1 lattice
    spacing thick, where the local crystalline order (hence the sites the
    condensate coheres on) is spoiled. A vacuum vortex line lying IN the
    fault sheet overlaps a depleted slab along its whole length: a line
    trap, not a point trap. We estimate the binding per length and, more
    importantly, whether resting in the fault forces a reconnection.

THE CATCH (and why it is not fatal):
    A knot's helicity is protected against reconnection with ANOTHER
    vortex line. A stacking fault is NOT a vortex line; it is a crystal
    (Burgers-sector) defect and carries no circulation. So a knot laid on
    a fault sheet is NOT touching a second vortex, and there is no
    circulation-reconnection to unwind it. The fault grips the knot
    through core-overlap (a condensate-density effect), while the
    helicity, a distinct topological charge, is untouched. That is the
    exemption the plain-loop fusion argument does not have: plain loops
    die on our DISLOCATION LINES (which do carry circulation, being
    screws/partials), whereas a knot parked on the fault SHEET between
    partials meets no circulation to swap ends with.

RESULT: the fault sheet is a candidate line-trap that binds by core
overlap without triggering helicity-changing reconnection. The open
number is the binding-vs-thermal margin and whether a fabricated fault
(a controlled partial-dislocation loop) can be made to hold a knot on
demand.
"""

import numpy as np

# framework constants
alpha = 1/137.035999
c     = 2.99792458e8
e     = 1.602176634e-19
ell   = 2.8179403262e-15
m0c2J = 0.51099895069e6*e/alpha          # node rest energy, J
m0c2  = m0c2J/e                           # eV
f_s   = 4/5

# core energy density of a vacuum vortex ~ its line tension / core area.
# line tension T = (f_s pi) (m0 c^2/ell) * ln(...) ; take ln~1 (core scale).
# core area ~ pi ell^2. energy density u_core ~ T/(pi ell^2).
T_line_over_ln = f_s*np.pi*m0c2/ell        # eV/m per unit ln
u_core = T_line_over_ln/(np.pi*ell**2)     # eV/m^3 (ln~1)

print("="*72)
print("CANDIDATE 1 - VACANCY: ruled out")
print("="*72)
print("  * neighbours relax INWARD -> compression maximum -> vortex repelled")
print("  * condensate rearranges at fixed number around it -> no grip")
print("  (Same incompressibility that kills vacancy-mediated superfluidity.)")
print()

print("="*72)
print("CANDIDATE 2 - STACKING FAULT: a line-trap that spares helicity")
print("="*72)
# A knot segment of length L lying in the fault overlaps a slab of
# thickness ~ell (the fault's strain thickness) over which the condensate
# coherence is spoiled. Binding ~ fraction of core energy recovered.
# Depletion fraction in the fault, phi_dep, is O(gamma_SF-driven); take a
# conservative O(0.1) of the core sitting in already-disrupted medium.
phi_dep = 0.1
E_bind_per_length = phi_dep*T_line_over_ln          # eV/m
E_bind_per_site   = E_bind_per_length*ell           # eV per lattice length
print(f"  vortex core energy density  u_core ~ {u_core:.2e} eV/m^3 (ln~1)")
print(f"  fault depletion fraction (conservative)  phi_dep = {phi_dep}")
print(f"  binding per length  phi_dep * T ~ {E_bind_per_length:.2e} eV/m")
print(f"  binding per lattice length        ~ {E_bind_per_site:.3f} eV"
      f"  = {E_bind_per_site/m0c2*1e3:.1f} milli-(m0c^2)")
print()
print("  Helicity check: a stacking fault carries NO circulation (it is a")
print("  Burgers-sector crystal defect, the ribbon between partials). A knot")
print("  laid in the sheet touches no second vortex line, so there is no")
print("  circulation reconnection -> helicity is untouched. Plain loops die")
print("  on DISLOCATION LINES (which DO carry circulation); a knot on the")
print("  fault SHEET does not meet one. That is the exemption.")
print()

# Thermal margin: is the per-length binding above the ambient second-sound
# thermal scale? The relevant comparison is binding over a knot persistence
# length vs k_B T_universe; but the sharp lab statement is the binding over
# a minimal knot (line ~10 ell):
Lambda_min = 10*ell
E_bind_knot = E_bind_per_length*Lambda_min
print(f"  minimal knot (line {Lambda_min:.1e} m): total fault binding"
      f" ~ {E_bind_knot/1e6:.1f} MeV")
print(f"    = {E_bind_knot/m0c2:.2f} x (m0 c^2); a real well, several nodes deep")
print("    over the minimal knot, and it scales with the loop's length.")
print()
print("="*72)
print("VERDICT")
print("="*72)
print("  Vacancy: NO (wrong sign, no coupling).")
print("  Stacking fault: CANDIDATE. Binds a knot by core overlap along a")
print("  line, and being circulation-free it does not trigger the helicity-")
print("  changing reconnection that unwinds plain loops. A fabricated fault")
print("  (a pinned partial-dislocation loop) is therefore the natural knot")
print("  trap to try. Open: the binding-vs-thermal margin and controlled")
print("  loading. The transducer's 'trapping problem' has a concrete target.")
print("="*72)
