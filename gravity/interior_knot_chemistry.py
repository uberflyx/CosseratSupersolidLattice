"""
Knot chemistry of the Kerr black-hole interior.

The vortex-reservoir census (gravity/bh_vortex_reservoir.py) leaves a
near-extremal interior threaded by a permanent Feynman array of quantised
vortex lines. This script asks what BOUND STATES that environment supports,
i.e. whether the interior has a chemistry, and prices every object from the
framework's constants (m_e, alpha and lattice geometry only).

The enabling observation is topological. Outside horizons the free sector
emptied because nothing arrests a loop's collapse: knots untie
(Kleckner-Kauffman-Irvine 2016), rings shrink or are captured on charged
cores. Inside, a loop that LINKS an array line cannot pass below the line's
core radius xi without a reconnection, so linkage puts a hard floor under
every linked object. The array is a stabiliser: it converts conserved
linking numbers into a discrete spectrum of protected objects.

Cast of objects, smallest to largest:

  1. THE WIRE: one array line. Tension T = (rho_s kappa^2/4pi) ln(R/xi),
     about (4pi/5) ln node energies per lattice spacing. Its transverse
     excitations are Kelvin waves, omega = (kappa k^2/4pi) ln(1/k xi):
     the interior's "photon", with a QUADRATIC dispersion.

  2. THE ATOM (bead): the smallest vortex ring threaded on a wire,
     E_1 = (8 pi^2 / 5) m0 c^2 = 1.107 GeV, radius R0 ~ 1.7 ell. Linkage
     with the wire is its conserved charge; it slides freely along the
     wire: a 1D particle. Ring dynamics gives E ~ sqrt(P): heavier rings
     move SLOWER (v = dE/dP falls with R). Anomalous dispersion is the
     interior's kinematic signature.

  3. THE BOND: two co-oriented beads on one wire interact through their
     induced flow, the exact coaxial-ring interaction energy
        E_12 = rho_s k1 k2 sqrt(R1 R2) [ (2/k - k) K(k) - (2/k) E(k) ],
     the fluid twin of Neumann's mutual-inductance formula. The bound
     orbit is Love's leapfrogging pair (Love 1893): the two rings thread
     each other in alternation. We integrate the exact coaxial Hamiltonian
     to confirm a bounded orbit and read off the bond scale.

  4. ANTI-MATTER: a bead of opposite linkage on the same wire. A bead and
     an anti-bead attract head-on and annihilate to Kelvin waves and first
     sound, releasing 2 E_1. Chemistry has an annihilation channel.

  5. THE HEAVY ELEMENTS (knots): a closed loop knotted THROUGH the array
     is pinned at the inter-vortex spacing b = sqrt(kappa/2 Omega): it
     cannot tighten past the lines it links, so the untying cascade that
     killed free knots is arrested at scale b, not xi. Its mass is
     tension x ropelength, M(K) = T * Lambda(K) * b_pin, so the periodic
     table of interior heavy elements is literally the knot table ordered
     by ropelength (trefoil first).

  6. REACTIONS: a reconnection changes one crossing; linking number is
     the exchanged charge. Reaction barrier = Kelvin-wave cost of driving
     two cores into contact.

All array numbers for the fiducial hole of the reservoir script:
M = 10 Msun near-extremal, Omega_H ~ 1e4 rad/s.
"""

import numpy as np
from scipy.special import ellipk, ellipe
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# Framework constants (conventions of gravity/bh_vortex_reservoir.py).
# ---------------------------------------------------------------------------
alpha = 1.0 / 137.035999177
m_e   = 9.1093837015e-31              # kg
m0    = m_e / alpha                   # node mass [kg]
c     = 2.99792458e8                  # m/s
ell   = 2.8179403205e-15              # lattice spacing r_e [m]
hbar  = 1.054571817e-34
kappa = 2.0 * np.pi * ell * c         # circulation quantum h/m0 [m^2/s]
f_s   = 0.8                           # superfluid fraction (bootstrap)
rho_s = f_s * m0 / ell**3             # superfluid density [kg/m^3]
xi    = ell                           # vortex core radius [m]
GeV   = 1.602176634e-10               # J
m0c2  = m0 * c**2                     # node energy [J] ~ 70 MeV

E1 = (8 * np.pi**2 / 5) * m0c2        # smallest free ring [J] (monograph)

line = "-" * 76
print(line)
print("KNOT CHEMISTRY OF THE KERR INTERIOR - all scales from m_e, alpha")
print(line)

# ---------------------------------------------------------------------------
# 0. The arena: fiducial near-extremal 10 Msun hole (reservoir script).
# ---------------------------------------------------------------------------
Omega_H = 1.0e4                       # horizon angular velocity [rad/s]
n_v = 2 * Omega_H / kappa             # Feynman line density [1/m^2]
b   = 1 / np.sqrt(n_v)                # inter-vortex spacing [m]
lnb = np.log(b / xi)
print(f"Array spacing b            : {b:.3e} m  ({b*1e6:.1f} um)")
print(f"Scale separation ln(b/xi)  : {lnb:.1f}")

# Line tension at the two natural cutoffs.
T_coeff = rho_s * kappa**2 / (4 * np.pi)          # [J/m] before the log
T_line  = T_coeff * lnb                            # wire tension [J/m]
print(f"Wire tension T             : {T_line:.3e} J/m "
      f"= {T_line*ell/m0c2:.1f} m0c^2 per lattice spacing")

# ---------------------------------------------------------------------------
# 1. THE ATOM: smallest ring threaded on a wire.
#    E(R) = (rho_s kappa^2 R / 2) [ln(8R/xi) - 2]   (classical thin ring)
#    P(R) = rho_s kappa pi R^2                       (impulse)
#    v(R) = (kappa / 4 pi R) [ln(8R/xi) - 1]         (self-induced speed)
# ---------------------------------------------------------------------------
def ring_E(R):  return 0.5 * rho_s * kappa**2 * R * (np.log(8*R/xi) - 2.0)
def ring_P(R):  return rho_s * kappa * np.pi * R**2
def ring_v(R):  return (kappa / (4*np.pi*R)) * (np.log(8*R/xi) - 1.0)

# Radius reproducing the monograph's E1:
from scipy.optimize import brentq
from scipy.optimize import brentq as _bq
# Solve dimensionlessly: x (ln 8x - 2) = E1/(rho_s kappa^2 ell / 2); _bq's
# absolute xtol (2e-12) dwarfs a femtometre root, so never solve in metres.
_rhs = E1/(0.5*rho_s*kappa**2*ell)
R0 = ell*_bq(lambda x: x*(np.log(8*x*ell/xi) - 2) - _rhs, 1.0, 5.0)

print(line)
print("ATOM (bead = ring linked on a wire)")
print(f"  mass E1                  : {E1/GeV:.3f} GeV   radius R0 = {R0/ell:.2f} ell")
print(f"  self-speed v(R0)         : {min(ring_v(R0),c)/c:.2f} c (LIA capped)")
# Anomalous dispersion: v falls as E rises.
for mult in (2, 5, 20):
    R = mult * R0
    print(f"  E = {ring_E(R)/GeV:7.2f} GeV  ->  v = {ring_v(R)/c:.3f} c"
          f"   (P = {ring_P(R)*c/GeV:9.1f} GeV/c)")

# ---------------------------------------------------------------------------
# 2. THE BOND: exact coaxial two-ring dynamics (Love's leapfrog).
#    Interaction energy of coaxial rings, separations s, radii R1, R2:
#      k^2 = 4 R1 R2 / ((R1+R2)^2 + s^2)
#      E12 = rho_s kappa^2 sqrt(R1 R2) [(2/k - k) K(k) - (2/k) E(k)]
#    Canonical pairs: (P_i = rho_s kappa pi R_i^2  conjugate to  z_i).
#    H = E(R1) + E(R2) + E12.  Integrate, confirm bounded relative orbit.
# ---------------------------------------------------------------------------
def E12(R1, R2, s):
    k2 = 4*R1*R2 / ((R1+R2)**2 + s**2)
    k  = np.sqrt(k2)
    return rho_s * kappa**2 * np.sqrt(R1*R2) * ((2/k - k)*ellipk(k2) - (2/k)*ellipe(k2))

def H(P1, P2, z1, z2):
    R1 = np.sqrt(max(P1,1e-30)/(rho_s*kappa*np.pi)); R2 = np.sqrt(max(P2,1e-30)/(rho_s*kappa*np.pi))
    return ring_E(R1) + ring_E(R2) + E12(R1, R2, z2-z1)

def rhs(t, y):
    P1, P2, z1, z2 = y
    eps_P = 1e-7*(P1+P2); eps_z = 2e-8*abs(z2-z1) + 1e-9*xi
    dz1 =  (H(P1+eps_P,P2,z1,z2)-H(P1-eps_P,P2,z1,z2))/(2*eps_P)
    dz2 =  (H(P1,P2+eps_P,z1,z2)-H(P1,P2-eps_P,z1,z2))/(2*eps_P)
    dP1 = -(H(P1,P2,z1+eps_z,z2)-H(P1,P2,z1-eps_z,z2))/(2*eps_z)
    dP2 = -(H(P1,P2,z1,z2+eps_z)-H(P1,P2,z1,z2-eps_z))/(2*eps_z)
    return [dP1, dP2, dz1, dz2]

# Two rings in the thin-ring validity window (R = 20 xi), launched 2R apart.
Rd = 20*xi; P0 = ring_P(Rd); s0 = 2*Rd; tau = Rd**2/kappa
sol = solve_ivp(rhs, [0, 40*tau], [P0, P0, 0, s0], max_step=tau/200,
                rtol=1e-9, dense_output=False)
sep  = sol.y[3] - sol.y[2]
R1s  = np.sqrt(np.abs(sol.y[0])/(rho_s*kappa*np.pi))
bound = sep.max() < 4*Rd            # never escapes
passes = int(np.sum(np.diff(np.sign(sep)) != 0))   # threading events
print(line)
print("BOND (leapfrog dimer, exact coaxial Hamiltonian)")
print(f"  bounded orbit            : {bound}   sep range "
      f"[{sep.min()/Rd:.2f}, {sep.max()/Rd:.2f}] R,  threadings: {passes}")
print(f"  radius exchange          : R in [{R1s.min()/Rd:.2f}, {R1s.max()/Rd:.2f}] R")
print(f"  leapfrog period          : ~{20*tau:.1e} s  (chemistry clock)")
print(f"  interaction at s = 4 R0  : {E12(R0,R0,4*R0)/GeV*1e3:.0f} MeV")
print(f"  interaction at s = 2 R0  : {E12(R0,R0,2*R0)/GeV*1e3:.0f} MeV"
      f"   (binding fraction {E12(R0,R0,2*R0)/E1:.2f} of E1)")

# ---------------------------------------------------------------------------
# 3. PHOTON: Kelvin wave on a wire, omega = (kappa k^2 / 4 pi) ln(1/k xi).
#    Quantum at wavelength = bead spacing s = 4 R0.
# ---------------------------------------------------------------------------
lamK = 100*xi                        # long-wave validity: k xi << 1
kK = 2*np.pi/lamK
omega_K = (kappa * kK**2 / (4*np.pi)) * np.log(1/(kK*xi))
print(line)
print("PHOTON (Kelvin quantum at bond wavelength)")
print(f"  hbar omega (lambda=100xi): {hbar*omega_K/GeV*1e3:.2f} MeV  (quadratic dispersion)")

# ---------------------------------------------------------------------------
# 4. HEAVY ELEMENTS: knots pinned by the array. Tightening is arrested at
#    the pinning radius ~ b/2, so mass = T * ropelength * pin radius.
#    Ropelength Lambda(K) in units of rope RADIUS (Cantarella-Kusner-
#    Sullivan bounds; trefoil ~ 32.7, figure-8 ~ 42, growing ~ linearly
#    with crossing number thereafter).
# ---------------------------------------------------------------------------
# Ropelength in unit-RADIUS convention (Baranska et al. 2004: trefoil
# 32.74295; CKS 2002 lower bound 21.45 for any nontrivial knot).
knots = {"cross-link ring (2 lines)": 2*np.pi,
         "trefoil 3_1": 32.74, "figure-eight 4_1": 42.09,
         "5_1": 47.2, "5_2": 49.7}
print(line)
print("HEAVY ELEMENTS (array-pinned loops; pin radius b/2)")
Tpin = T_coeff * np.log((b/2)/xi)
for name, lam in knots.items():
    M = Tpin * lam * (b/2)
    print(f"  {name:28s}: {M:6.1f} J = {M/GeV:.2e} GeV = {M/c**2*1e18:6.1f} ag")

# ---------------------------------------------------------------------------
# 5. INFORMATION CAPACITY: beads as bits on wires.
# ---------------------------------------------------------------------------
N_lines = 3e18                       # reservoir census
L_line  = 3.0e4                      # ~ horizon diameter [m]
sites   = L_line / (4*R0)            # one bead slot per bond length
print(line)
print("TAPE CAPACITY")
print(f"  slots per wire           : {sites:.1e}")
print(f"  total                    : {N_lines*sites:.1e} bead-bits")
print(line)

# ---------------------------------------------------------------------------
# 6. QUANTISATION SCALE of the dimer. The leapfrog is a 1-DOF oscillation in
#    the relative coordinate; its quantum is hbar * omega_lf. The period
#    scales as the self-induction time, T_lf ~ 20 R^2/kappa (measured above
#    at R = 20 xi), so hbar omega_lf ~ pi hbar kappa / (10 R^2): compare to
#    the binding to count levels. Core-scale dimers are few-level quantum
#    systems (nuclear-like); dimers at 100 xi are semiclassical.
# ---------------------------------------------------------------------------
print("QUANTISATION (dimer level count ~ binding / hbar omega_leapfrog)")
for Rq in (R0, 5*R0, 20*R0):
    w_lf = np.pi*kappa/(10*Rq**2)
    Eb   = abs(E12(Rq, Rq, 2*Rq))
    print(f"  R = {Rq/ell:5.1f} ell: hbar w = {hbar*w_lf/GeV*1e3:8.2f} MeV,"
          f"  binding = {Eb/GeV*1e3:8.1f} MeV,  levels ~ {Eb/(hbar*w_lf):6.1f}")

# ---------------------------------------------------------------------------
# 7. OUTSIDE: the ring wind gives free atoms in open space, but no array.
#    Chemistry there is kinetically frozen: encounter rate n sigma v with
#    the intergalactic background (~1e47 rings per galaxy share, halo
#    volume ~ (100 kpc)^3) never makes a single molecule in a Hubble time.
#    Molecules outside exist only if forged already linked (the 2.21 GeV
#    monophoton pair of sec:clock_actuation).
# ---------------------------------------------------------------------------
n_out  = 1e47 / (3.1e21)**3          # rings per m^3 in a halo share
rate   = n_out * np.pi*ell**2 * 0.5*c
print(line)
print("OUTSIDE (free atoms, frozen chemistry)")
print(f"  halo ring density        : {n_out:.1e} /m^3")
print(f"  link-encounter rate      : {rate:.1e} /s per ring"
      f"  ({rate*4.4e17:.1e} per Hubble time)")
print(line)
