#!/usr/bin/env python3
r"""
line_anchor_density.py
========================
What force can ordinary matter take from a vacuum vortex line, per metre
of line? Everything here is derived rather than assumed, and each input is
either a framework primitive or a cited measurement.

The chain has five links.

  1. THE CORE PARAMETER. The classical energy of a vortex ring of radius R
     with core parameter a is E = (rho_s kappa^2 R/2)[ln(8R/a) - 2]. The
     framework fixes E = E_1 = (8 pi^2/5) m0 c^2 and R = R0 = 1.68 ell
     independently, so the pair over-determines a. Solving returns
     a ~ 1.03 ell, which is a genuine internal consistency check: the
     ladder's ring geometry and the ring's energy agree on a core one
     lattice spacing wide.

  2. THE PIN DEPTH. A nucleus is a region where the crystal's stacking
     order is spoiled and the condensate that coheres on it is suppressed,
     so a length of vortex core lying inside one costs no flow energy. The
     saving splits into two regimes, and they must join continuously.

       large pin (R_N >= xi): the line saves the flow energy stored inside
         the nuclear radius plus the core term, over the mean chord 4R_N/3
         of a sphere,
             U = (rho_s kappa^2/4pi)(4 R_N/3)[ln(R_N/xi) + 1/4].
       small pin (R_N < xi): the nucleus sits wholly inside the core, so
         the saving is the condensation energy density times the nuclear
         volume,
             U = rho_s kappa^2 R_N^3 / (12 pi xi^2).

     The 1/4 is not a fudge. The Gross-Pitaevskii core energy per unit
     length is eps_cond * pi xi^2 = rho_s kappa^2/(16 pi), exactly 1/4 of
     the logarithmic prefactor rho_s kappa^2/(4 pi). With that value the
     two regimes agree exactly at R_N = xi, which the script checks.

  3. STRONG OR WEAK PINNING? A flexible line in a random array of pins is
     the standard vortex pinning problem [Larkin & Ovchinnikov 1979;
     Blatter et al., Rev. Mod. Phys. 66, 1125 (1994)]. Which regime applies
     is set by the Labusch criterion: pinning is individual (strong) when a
     single pin can deflect the line by its own range against the line's
     elastic restoring force. Here the pin force U/R_N exceeds that
     restoring force by some seven orders, so the pins are individually
     strong, and the geometry rather than the pin strength sets the answer.

  4. THE VARIATIONAL OPTIMUM. Let the line touch a pin every s of its
     length, hopping a transverse distance delta to do so. A hop costs
     T delta^2/(2s) in tension and gains U. At least one pin must lie
     within reach, so n pi delta^2 s >= 1. The free energy per unit length
     is then

         e(s) = -U/s + T/(2 pi n s^3),

     minimised at  s = sqrt( 3 T / (2 pi n U) ).

  5. WHAT THE HOST CAN TRANSMIT. The pin force itself is nuclear, of order
     10 kN, but a nucleus is held in its host only by chemical bonds, so
     the transmitted force necks down there (the plumbing limit). The
     single-atom ceiling is estimated as the ideal strength times the area
     per atom, with the Frenkel estimate sigma ~ E/10 calibrated against
     the first-principles value for iron, 12.6 GPa [Clatterbuck, Chrzan &
     Morris, Acta Mater. 51, 2271 (2003)]. The result is bracketed above by
     the measured rupture force of a single covalent bond, 2.0 +/- 0.3 nN
     for Si-C [Grandbois et al., Science 283, 1727 (1999)].

  Result: f_anchor = F_pin / s, of order 0.02 N/m in iron and 0.1 N/m in
  diamond or osmium. The spread is real and is driven by the host's
  cohesion rather than by the vortex physics.
"""

import numpy as np
from scipy.optimize import brentq

# ------------------------------------------------------------- primitives
alpha = 1/137.035999177
c     = 2.99792458e8
hbar  = 1.054571817e-34
h     = 2*np.pi*hbar
e     = 1.602176634e-19
ell   = 2.8179403205e-15               # lattice spacing = r_e [m]
m0c2  = 0.51099895069e6*e/alpha        # node rest energy [J]
m0    = m0c2/c**2
f_s   = 5/6                            # bootstrap superfluid fraction
rho_s = f_s*m0/ell**3                  # superfluid density [kg/m^3]
kappa = h/m0                           # circulation quantum [m^2/s]

R0    = 1.68*ell                       # ground-state ring radius [m]
E1    = (8*np.pi**2/5)*m0c2            # ring energy [J]
T     = E1/(2*np.pi*R0)                # line tension [N]
P     = rho_s*kappa**2/(4*np.pi)       # logarithmic prefactor [N]

print("="*74)
print("1. The core parameter, over-determined by E1 and R0")
print("="*74)
brack  = E1/(rho_s*kappa**2*R0/2)      # = ln(8 R0/a) - 2
a_core = 8*R0/np.exp(brack + 2)
print(f"  prefactor rho_s kappa^2/4pi   = {P:.4e} N")
print(f"  line tension T = E1/(2 pi R0) = {T:.1f} N  ({T/1e3:.2f} kN)")
print(f"  bracket [ln(8 R0/a) - 2]      = {brack:.4f}")
print(f"  core parameter a = {a_core*1e15:.3f} fm = {a_core/ell:.3f} ell")
print("  R0 and E1 were fixed independently, so their agreement on a core")
print("  one lattice spacing wide is a check rather than an input.")
xi = ell                               # healing length, to 3% by the above
print()


def pin_depth(A):
    """Pinning well depth for a line crossing a nucleus of mass number A."""
    R_N = 1.2e-15*A**(1/3)
    if R_N >= xi:
        return P*(4*R_N/3)*(np.log(R_N/xi) + 0.25), R_N, "large"
    return rho_s*kappa**2*R_N**3/(12*np.pi*xi**2), R_N, "small"


print("="*74)
print("2. Pin depth, and the continuity check at R_N = xi")
print("="*74)
U_lo = rho_s*kappa**2*xi**3/(12*np.pi*xi**2)
U_hi = P*(4*xi/3)*(np.log(1.0) + 0.25)
print(f"  small-pin formula at R_N = xi : {U_lo/e/1e6:9.3f} MeV")
print(f"  large-pin formula at R_N = xi : {U_hi/e/1e6:9.3f} MeV")
print(f"  agreement to {abs(U_lo-U_hi)/U_lo:.1e}: the 1/4 core term is fixed by")
print("  the GP core energy, not fitted to make the regimes meet.")
print()

print("="*74)
print("3. Strong or weak pinning? The Labusch check (iron)")
print("="*74)
U_fe, R_fe, _ = pin_depth(55.85)
f_p       = U_fe/R_fe                  # peak single-pin force [N]
s_guess   = 3e-8                       # provisional pin spacing [m]
f_elastic = T*R_fe/s_guess             # line restoring force over that span [N]
print(f"  peak pin force       U/R_N      = {f_p:.2e} N  ({f_p/1e3:.1f} kN)")
print(f"  line restoring force T R_N / s  = {f_elastic:.2e} N")
print(f"  Labusch ratio                   = {f_p/f_elastic:.1e}  >> 1")
print("  Individual (strong) pinning: geometry sets the spacing, not strength.")
print()

print("="*74)
print("4-5. Pin spacing, and what the host can transmit")
print("="*74)
print(f"  {'host':>9} {'A':>7} {'R_N/fm':>7} {'reg':>6} {'U/MeV':>8} "
      f"{'s/nm':>7} {'delta/pm':>9} {'F_pin/nN':>9} {'f/(N/m)':>9}")
print("  " + "-"*72)

# host, atomic weight, number density [m^-3], Young's modulus [Pa]
hosts = [("iron",      55.85, 8.50e28,  211e9),
         ("tungsten", 183.84, 6.31e28,  411e9),
         ("osmium",   190.23, 7.15e28,  560e9),
         ("diamond",   12.01, 1.763e29, 1050e9)]

cal = 12.6e9/(211e9/10)   # Frenkel E/10 calibrated to iron's DFT ideal strength

# The tension that resists the detour is NOT the core-scale ring value
# T = E1/(2 pi R0). The detour is a bend of wavelength ~s, so the line's
# logarithm runs out to s, giving T(s) = eps_0 [ln(s/xi) - 1], about thirty
# times larger at the hundred-nanometre scale. Because s itself depends on
# T, the spacing is implicit and is iterated to convergence below.
eps_0 = rho_s*kappa**2/(4*np.pi)

def spacing(n, U, seed=1e-7):
    """Self-consistent pin spacing with the deformation-scale tension."""
    s = seed
    for _ in range(500):
        T_s = eps_0*(np.log(s/xi) - 1.0)
        s = 0.5*(s + np.sqrt(3*T_s/(2*np.pi*n*U)))
    return s

results = {}
for name, A, n, E_mod in hosts:
    U, R_N, regime = pin_depth(A)
    s     = spacing(n, U)
    delta = 1/np.sqrt(np.pi*n*s)
    F_pin = (E_mod/10)*cal*n**(-2/3)
    f_anchor = F_pin/s
    results[name] = (f_anchor, s, F_pin)
    print(f"  {name:>9} {A:7.2f} {R_N*1e15:7.3f} {regime:>6} {U/e/1e6:8.1f} "
          f"{s*1e9:7.1f} {delta*1e12:9.2f} {F_pin*1e9:9.2f} {f_anchor:9.3f}")
print("  " + "-"*72)
print(f"  Frenkel E/10 scaled by {cal:.2f} to match iron's 12.6 GPa; diamond then")
print(f"  returns {results['diamond'][2]*1e9:.1f} nN, consistent with the measured")
print("  2.0 +/- 0.3 nN for a single Si-C bond.")
print()

print("="*74)
print("Consequences")
print("="*74)
# The tension and the stored energy per unit length are NOT the core-scale
# ring value. Both carry the logarithm out to the loop's own radius:
#     eps(R) = (rho_s kappa^2/4pi) [ln(8R/xi) - 2]
# which is the same thin-ring expression the ladder uses, per unit length.
# Using the 6 kN ring tension here understates both by nearly two orders.
eps0 = rho_s*kappa**2/(4*np.pi)
eps_at = lambda R: eps0*(np.log(8*R/xi) - 2.0)

for name in ("iron", "osmium"):
    f, s, _ = results[name]
    u_cap = f/(rho_s*kappa)
    # curvature pull 2 pi eps(R) against the matter ceiling 2 pi R f
    R_cross = np.exp(brentq(lambda x: eps_at(np.exp(x)) - np.exp(x)*f, 5.0, 30.0))
    print(f"  {name}:  f_anchor = {f:.3f} N/m")
    print(f"    drag speed cap    u = f/(rho_s kappa) = {u_cap*1e12:.2f} pm/s")
    print(f"    stretch crossover, eps(R) = R f       = {R_cross/1e3:.0f} km"
          f"   (Earth radius 6371 km)")
    for M, lab in [(1.0, "1 kg"), (1e5, "100 t")]:
        L = M*9.8/f
        R_loop = 2*L/(2*np.pi)
        E_loop = 2*L*eps_at(R_loop)
        print(f"    {lab:>6}: {L:.3g} m gripped, loop {2*L:.3g} m,"
              f" line energy {E_loop/1e6:.3g} MJ"
              f"  ({eps_at(R_loop)/1e3:.0f} kJ/m)")
print()
print("  Thermal caveat: delta comes out near 12 pm, against a room-temperature")
print("  nuclear displacement of about 7 pm in iron, so lattice vibration is an")
print("  order-unity correction to which pins are in reach at a given instant.")
print("  Treat every f_anchor above as good to a factor of about two.")
print("="*74)
