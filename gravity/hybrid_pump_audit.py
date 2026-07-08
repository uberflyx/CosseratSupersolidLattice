#!/usr/bin/env python3
r"""
hybrid_pump_audit.py
====================
Audit of the "pipe-array condensate pump" propulsion candidate for the
Cosserat-supersolid vacuum, and the transducer selection theorem it proves.

The candidate: a superconducting asymmetric pipe drives the chirality coupling
at its open end, tilting the vacuum's compact fourth axis by an angle phi at a
focal point. Via the boost lemma (cos phi = 1/gamma(c sin phi), exact), that
tilt is a chemical-potential well delta_mu = mu0 (1 - cos phi). The hope was to
pump the superfluid condensate with a travelling train of such wells, using
present-day magnets and no vortex-knot forge.

This script shows the hope fails on three independent grounds, each a clean
piece of hydrodynamics:

  Blade 1 (pressure screening): a conservative well in a near-incompressible
          medium drives no flow; pressure rises to cancel it.
  Blade 2 (energy clamp): the solenoidal (circulation) part of the tilt cannot
          be screened, so if the phase followed it the stored energy would be
          rho_s c^2 phi^2 / 2, ~1e15x the available EM energy density. Back-
          reaction clamps the achievable tilt to phi_eq << phi_pipe.
  Blade 3 (Landau slip / vacuum wind): a well fixed to the craft must trap the
          condensate against the 370 km/s aether wind; below that threshold a
          superfluid slips through frictionlessly and transmits no momentum.

The rotor of the actuation chapter is exempt from all three, because its grip is
a topological winding (Anderson phase slip, 2 pi per transit) rather than a
potential. That exemption is the theorem: only a topological grip can pump the
vacuum condensate.

Reproduces the appendix numbers. Framework constants from the monograph; pipe
tilt formula from the engineering note (kappa_EM = sqrt(rho/eps0)).

Author: M. A. Cox.  Companion to Appendix "Driving the vacuum".
"""

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants (CODATA 2018) and framework scales
# ---------------------------------------------------------------------------
ALPHA   = 7.2973525693e-3        # fine-structure constant
HBAR    = 1.054571817e-34        # reduced Planck constant [J s]
H       = 2 * np.pi * HBAR       # Planck constant [J s]
C       = 2.99792458e8           # speed of light [m/s]
E_CH    = 1.602176634e-19        # elementary charge [C]
MU0_MAG = 4e-7 * np.pi           # magnetic constant [H/m]
KB      = 1.380649e-23           # Boltzmann constant [J/K]
T_CMB   = 2.725                  # CMB temperature [K]

# Framework: node rest energy sets the master scale mu0 = m0 c^2.
M0C2_MEV = 70.3                  # node rest energy [MeV] (E1 = 8pi^2/5 m0c2 = 1.11 GeV)
M0C2     = M0C2_MEV * 1e6 * E_CH # [J]
M0       = M0C2 / C**2           # node mass [kg]
ELL      = HBAR * C / M0C2       # lattice scale = hbar/(m0 c) [m]
F_S      = 4 / 5                 # superfluid fraction (bootstrap)
RHO_S    = F_S * M0 / ELL**3     # superfluid density [kg/m^3]
MU0_NODE = M0C2                  # condensate chemical potential per node [J]
THETA_CH = ALPHA**2 / (2 * np.pi)      # chirality parameter
KAPPA_EM = 2.510e13              # electromechanical conversion sqrt(rho/eps0) [T]
KAPPA_CIRC = H / M0              # circulation quantum [m^2/s]
XI       = ELL                  # vortex core / healing scale ~ ell [m]
K_SF     = 1.5e73               # condensate bulk modulus [Pa]
V_WIND   = 370e3                # lab motion through condensate frame [m/s]

# Two pipe configurations bracketing today's magnet technology:
# (name, B [T], bore R [m], read distance r [m])
CONFIGS = {
    "A (0.1 m bore, 30 T, r=0.5 m)": (30.0, 0.05, 0.5),
    "B (1 m bore,   30 T, r=1 m)":   (30.0, 0.50, 1.0),
}


def phi_pipe(B, R, r):
    """Compact-axis tilt at distance r from an asymmetric pipe mouth [rad].

    Engineering-note result: phi = 2 theta_ch B R^2 / (kappa_EM ell r).
    """
    return 2 * THETA_CH * B * R**2 / (KAPPA_EM * ELL * r)


def well_depth(phi):
    """Chemical-potential well depth from a tilt phi, delta_mu = mu0 (1-cos phi) [J].

    Uses the versine 2 sin^2(phi/2) form, which is numerically stable for the
    very small phi (~1e-12) produced by the Blade 2 energy clamp, where
    1 - cos(phi) would underflow to zero in double precision.
    """
    return MU0_NODE * 2.0 * np.sin(phi / 2.0) ** 2


def ring_energy(R):
    """Vortex-ring self-energy [J] (thin-core approximation)."""
    return 0.5 * RHO_S * KAPPA_CIRC**2 * R * (np.log(8 * R / XI) - 2.0)


def ring_impulse(R):
    """Vortex-ring impulse [kg m/s]."""
    return RHO_S * KAPPA_CIRC * np.pi * R**2


def boost_lemma_check():
    """Verify cos(phi) = 1/gamma(c sin phi) to machine precision."""
    print("=" * 70)
    print("BOOST LEMMA: a compact-axis tilt is a frozen boost")
    print("=" * 70)
    for name, (B, R, r) in CONFIGS.items():
        phi = phi_pipe(B, R, r)
        v = C * np.sin(phi)
        gamma_inv = np.sqrt(1 - (v / C) ** 2)
        err = abs(np.cos(phi) - gamma_inv)
        print(f"  {name}")
        print(f"    phi = {phi:.3e} rad,  v = c sin(phi) = {v:.3e} m/s")
        print(f"    |cos(phi) - 1/gamma| = {err:.1e}  (exact identity)")
    print()


def tilt_and_well():
    """Report tilt, well depth, and frozen flow speed per configuration."""
    print("=" * 70)
    print("CANDIDATE DRIVE: tilt, well depth, frozen flow speed")
    print("=" * 70)
    for name, (B, R, r) in CONFIGS.items():
        phi = phi_pipe(B, R, r)
        dmu = well_depth(phi)
        print(f"  {name}")
        print(f"    phi = {phi:.3e} rad")
        print(f"    well depth = {dmu / E_CH * 1e3:.2f} meV/node "
              f"({dmu / (KB * T_CMB):.0f} x CMB thermal)")
        print(f"    frozen flow v = c phi = {C * phi:.3e} m/s")
    print()


def blade1_pressure_screen():
    """Blade 1: a conservative well is screened by pressure; no flow results."""
    print("=" * 70)
    print("BLADE 1: pressure screening of the conservative well")
    print("=" * 70)
    B, R, r = CONFIGS["A (0.1 m bore, 30 T, r=0.5 m)"]
    phi = phi_pipe(B, R, r)
    t_screen = r / C                       # first sound crosses the well
    drho_over_rho = (1 - np.cos(phi))       # fractional density shift ~ phi^2/2
    print(f"  Config A: screening time ~ r/c1 = {t_screen:.1e} s")
    print(f"  fractional density shift to cancel well: drho/rho ~ {drho_over_rho:.1e}")
    print(f"  -> mu_loc uniform in steady state; v_s = 0; no pumped flow.")
    print(f"  (holding the residual density shift needs ~K_sf drho/rho "
          f"= {K_SF * drho_over_rho:.0e} Pa: the shut compression route.)")
    print()


def blade2_energy_clamp():
    """Blade 2: the solenoidal part cannot be screened; energy clamps the tilt."""
    print("=" * 70)
    print("BLADE 2: energy clamp on the solenoidal tilt")
    print("=" * 70)
    B, R, r = CONFIGS["A (0.1 m bore, 30 T, r=0.5 m)"]
    phi = phi_pipe(B, R, r)
    u_flow = 0.5 * RHO_S * (C * phi) ** 2   # flow/strain energy if phase follows tilt
    u_em = B**2 / (2 * MU0_MAG)             # available EM energy density
    phi_eq = np.sqrt(2 * u_em / (RHO_S * C**2))
    print(f"  if phase follows tilt: u_flow = rho_s c^2 phi^2 / 2 = {u_flow:.3e} J/m^3")
    print(f"  available EM energy density (30 T): u_EM = {u_em:.3e} J/m^3")
    print(f"  shortfall factor u_flow/u_EM = {u_flow / u_em:.1e}")
    print(f"  -> back-reaction clamps tilt to phi_eq = {phi_eq:.3e} rad "
          f"({phi / phi_eq:.1e}x below pipe prediction)")
    print(f"  clamped well depth = {well_depth(phi_eq) / E_CH:.2e} eV/node")
    print(f"  (phi_eq is B-independent: both u_flow and u_EM scale as B^2)")
    print()
    return phi_eq


def blade3_wind_gate(phi_eq):
    """Blade 3: the well must trap against the vacuum wind; else frictionless slip."""
    print("=" * 70)
    print("BLADE 3: vacuum wind gate and the Landau slip")
    print("=" * 70)
    dmu_gate = 0.5 * M0 * V_WIND**2
    phi_gate = V_WIND / C
    print(f"  trap threshold: delta_mu > (1/2) m0 v_wind^2 = {dmu_gate / E_CH:.1f} eV")
    print(f"  equivalently phi > v_wind/c = {phi_gate:.3e}")
    for name, (B, R, r) in CONFIGS.items():
        phi = phi_pipe(B, R, r)
        dmu = well_depth(phi)
        verdict = "clears (ideal)" if dmu > dmu_gate else "FAILS"
        print(f"    {name}: ideal well {dmu / E_CH:.3e} eV -> {verdict}")
    print(f"  clamped well (Blade 2) = {well_depth(phi_eq) / E_CH:.2e} eV "
          f"-> short of gate by {dmu_gate / well_depth(phi_eq):.1e}x")
    print(f"  swept-pattern escape: sweep wells at v_wind; residual gate phi > w/c,")
    print(f"    at 370 km/s across metre wells this is a ~370 kHz modulation.")
    print()


def nucleation_barrier():
    """Metastability: vortex-ring nucleation barrier in the texture's counterflow."""
    print("=" * 70)
    print("METASTABILITY: vortex-ring nucleation barrier")
    print("=" * 70)
    B, R, r = CONFIGS["A (0.1 m bore, 30 T, r=0.5 m)"]
    v = C * phi_pipe(B, R, r)
    Rgrid = np.logspace(np.log10(2 * XI), -6, 4000)
    barrier = ring_energy(Rgrid) - v * ring_impulse(Rgrid)
    i = int(np.argmax(barrier))
    print(f"  counterflow v = {v:.2e} m/s")
    print(f"  barrier max = {barrier[i]:.2e} J at R* = {Rgrid[i]:.2e} m")
    print(f"  barrier / kT_CMB = {barrier[i] / (KB * T_CMB):.1e} "
          f"-> homogeneous nucleation rate negligible")
    print()


def rotor_exemption():
    """The topological rotor is immune to all three blades."""
    print("=" * 70)
    print("ROTOR EXEMPTION: the topological grip survives")
    print("=" * 70)
    L, g = 10.0, 9.81
    Nfr = M0 * g * L / H                 # dmu = N h f_r needed for 1 g across L
    print(f"  Anderson slip: 2 pi per vortex transit, independent of pressure,")
    print(f"    energy budget, and wind speed.")
    print(f"  1 g across {L:.0f} m needs N f_r = m0 g L / h = {Nfr:.2e} /s")
    print(f"  met by ~1200 knots on a cm rotor at 1 km/s tip speed.")
    print(f"  -> forge-gated, not transducer-gated: the theorem forces the forge.")
    print()


if __name__ == "__main__":
    print(f"\nlattice scale ell = {ELL:.3e} m")
    print(f"superfluid density rho_s = {RHO_S:.3e} kg/m^3")
    print(f"chirality parameter theta_ch = {THETA_CH:.3e}\n")
    boost_lemma_check()
    tilt_and_well()
    blade1_pressure_screen()
    phi_eq = blade2_energy_clamp()
    blade3_wind_gate(phi_eq)
    nucleation_barrier()
    rotor_exemption()
    print("Verdict: the smooth-well pump fails on three independent grounds;")
    print("only the topological rotor pumps the condensate. QED.")
