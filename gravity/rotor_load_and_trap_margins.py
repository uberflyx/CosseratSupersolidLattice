"""
The knot rotor's load path: three margins, one discovered bottleneck.

The actuation section prices the drive at ~1,200 trapped knots on a
centimetre rotor with 0.7 N of Magnus load each, and names the trap (a
stacking-fault sheet) while leaving "the binding-against-thermal margin
and the loading" open. This script computes the full mechanical chain,
knot -> pinning well -> fault -> ordinary matter -> rotor, and finds the
first two links absurdly strong and the last one the true design driver.

  MARGIN 1 (thermal): the pinning well is ~10^8 eV deep against a
  laboratory thermal quantum of 0.025 eV. Boltzmann escape is never.

  MARGIN 2 (depinning): the well delivers force up to ~E_bind/xi ~ 10^4 N
  per knot, four orders above the 0.7 N Magnus load. The vacuum-side grip
  is not the problem.

  MARGIN 3 (the anchor, DISCOVERED BOTTLENECK): the fault must hand its
  0.7 N to ordinary matter, and ordinary matter holds itself together
  with chemical bonds of ~1 nN each. One atom cannot anchor a knot; it
  would simply be dragged through the rotor. Every trap must therefore
  spread its load across >~ 10^9 lattice atoms: a mesoscopic engineered
  fault structure (tens of nanometres at least), not a single pinned
  hadron. Loading (steering a forged knot onto such a structure) inherits
  that scale. Total rotor load, ~840 N, is trivial once distributed.

  A DESIGN FORK, flagged: free rings self-propel at ~c/2, so a beam
  raceway could replace carriage entirely (one or two circulating rings
  through a metre circuit already give the required 2e7 transits/s). But
  self-speed makes the Magnus reaction on any guide ~25 kN, and the
  closure theorem that shields matter from the 370 km/s aether wind
  protects only NET force: a macroscopic free guide line would carry
  internal wind stress of ~1e16 N per metre and shred. The raceway is
  therefore not free; it trades the anchor problem for a guide problem.
"""
import numpy as np

alpha = 1.0/137.035999177
m0_J  = (0.51099895069/alpha) * 1.602176634e-13   # node mass energy [J]
l_m   = 2.8179403205e-15
c     = 2.99792458e8
kappa = 2*np.pi*l_m*c
rho_s = 0.8 * (m0_J/c**2) / l_m**3

line = "-"*78
print(__doc__.strip().splitlines()[0]); print(line)

# Margin 1: thermal
f_overlap = 0.5                     # MODEL INPUT: fault thins, not empties
L_knot = 10*l_m                     # minimal-knot line length ~ 10 l
E_bind = f_overlap * (rho_s*kappa**2/(8*np.pi)) * L_knot
kT_lab = 1.380649e-23 * 300
print(f"1. THERMAL: E_bind ~ {E_bind/1.602e-13:.0f} MeV vs kT(300K) = 0.025 eV")
print(f"   ratio ~ {E_bind/kT_lab:.0e}: Boltzmann escape never happens.")

# Margin 2: depinning force capacity
F_pin = E_bind / l_m
F_magnus = rho_s*kappa*np.sqrt(6)*l_m*1.0e3      # at 1 km/s tip speed
print(f"2. DEPINNING: well delivers up to E_bind/xi ~ {F_pin:.1e} N;")
print(f"   Magnus load at 1 km/s ~ {F_magnus:.2f} N: margin ~ {F_pin/F_magnus:.0e}.")

# Margin 3: the anchor
F_bond = 1e-9                        # MODEL INPUT: ~nN per chemical bond
N_anchor = F_magnus / F_bond
size_nm = (N_anchor**(1/3)) * 0.3    # ~0.3 nm per atom
print(f"3. ANCHOR (bottleneck): one bond holds ~1 nN, so each trap needs")
print(f"   >= {N_anchor:.0e} anchoring atoms, a structure of >= {size_nm:.0f} nm scale.")
print(f"   Total rotor load 1200 x {F_magnus:.1f} N ~ {1200*F_magnus:.0f} N: trivial")
print(f"   once distributed. The trap is mesoscopic engineering, not a hadron.")

# Design fork: the self-propelled raceway
v_ring = 0.5*c
F_guide = rho_s*kappa*np.sqrt(6)*l_m*v_ring
f_transit = v_ring/10.0              # one ring, 10 m circuit
wind = 3.7e5
stress_per_m = rho_s*kappa*wind
print(f"FORK: a self-propelled ring gives {f_transit:.1e} transits/s alone")
print(f"   (needs 2e7: one or two rings suffice), but its guide feels")
print(f"   {F_guide/1e3:.0f} kN of Magnus, and any macroscopic free guide line")
print(f"   carries wind stress ~{stress_per_m:.0e} N/m. Carriage keeps the")
print(f"   anchor problem; the raceway trades it for a guide problem.")
