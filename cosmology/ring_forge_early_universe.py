#!/usr/bin/env python3
r"""
ring_forge_early_universe.py
============================
Why the early universe leaves no free vortex rings, and why the rings an
astrophysical or laboratory forge makes later survive.

THE QUESTION
------------
The great annihilation ran ~1e9 matter-antimatter annihilations per surviving
particle. A collider can in principle forge a free vortex ring from a
collision (the monophoton channel), so did the annihilation era forge rings
in bulk? And if not, why do the rings a modern forge makes survive when the
primordial Kibble-Zurek population did not?

THE ANSWER, IN THREE PARTS
--------------------------
1. THE FLOOR CLOSES EVERY ANNIHILATION AT REST.
   A ring has no rest mass: its energy E(R) and impulse P(R) are independent
   functions of its radius R (cosmology/vortex_hadron_ladder.py carries the
   dispersion). At the adopted superfluid fraction f_s = 0.99932 the ground
   ring carries ~3.9 GeV/c of impulse on E_1 = 2 pi^2 f_s m0 c^2 = 1.38 GeV,
   and the balancing momentum must go to light quanta, which need
   E_gamma >= |p_gamma| c. So any final state containing a ring needs
       E + |P| c >= min_R [E(R) + P(R) c] = E_1 [1 + (R0/ell)^2] = 5.27 GeV
   in the condensate frame, where E and P are the initial state's total
   energy and momentum. Every annihilation at rest, the nucleon pair
   included (2 m_N c^2 = 1.88 GeV), sits below that floor.

2. IN FLIGHT, KINEMATICS ALONE DOES NOT CLOSE THE CHANNEL.
   A nucleon and an antinucleon moving together at speed beta c have
       E + |P| c = 2 gamma m_N c^2 (1 + beta) = 2 m_N c^2 sqrt((1+beta)/(1-beta)),
   so they clear the floor at sqrt((1+beta)/(1-beta)) = floor / (2 m_N c^2),
   i.e. beta = 0.774, gamma = 1.58, with only 2 m_N c^2 (gamma - 1) = 1.09 GeV
   of combined kinetic energy. The co-moving pair is the cheapest way to
   reach the floor: writing E_i + p_i c = m_N c^2 e^{y_i} (y_i the rapidity),
   the floor fixes sum e^{y_i} while the kinetic energy sum (cosh y_i - 1)
   is convex, so the minimum sits at y_1 = y_2. (An antinucleon striking a
   nucleon at rest needs 2.06 GeV/c, 1.33 GeV of kinetic energy.) Below
   k_B T_geom = m0 c^2 / sqrt(2) = 49.5 MeV, the warmest epoch at which a
   ring exists (the compact direction must be frozen for a closed loop of
   circulation to be defined), nucleons and antinucleons are still in
   chemical equilibrium, annihilating and re-forming far faster than the
   universe expands, so in-flight forging is not shut off by a Boltzmann
   tail. The equilibrium antinucleon abundance is printed for scale, with
   degeneracy g = 4 (antiprotons and antineutrons, two spin states each).

3. THE SINK CLOSES IT: WHATEVER FORMS DIES.
   A plain ring fuses on the circulation-carrying core of any charged
   particle (the condensate-vortex half of every electron and proton). The
   destruction rate Gamma = n_pm sigma P_fuse c is linear in the ring number
   and fed by the whole e+e- bath, so it cannot self-quench the way a
   WIMP's quadratic annihilation does. At even the minimal fusion
   probability of 1e-15 per contact, the ring lifetime is under a
   microsecond against a millisecond Hubble time. While the e+e- bath is
   relativistic and g_* fixed, Gamma/H is proportional to T, so the optical
   depth a ring born at T accumulates down to temperature T_f is
   (Gamma/H)(T) - (Gamma/H)(T_f) ~ (Gamma/H)(T): thousands of e-folds for
   any ring born while antinucleons are abundant.

BORN HOT VERSUS BORN COLD
-------------------------
A forged ring today is the same unprotected object. What changed is the
density of executioners. The cores never vanished; their number density fell
from ~4e31 per cm^3 (the e+e- pair bath at 1 MeV) to ~1 per cm^3
(interstellar medium) to ~1e-7 (intergalactic medium). The mean free path
grows from a nanometre to a megaparsec to beyond the observable universe.
Same weapon, same fragile ring; the room emptied.

Light arithmetic only.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

# ---- constants -------------------------------------------------------------
hbarc  = 197.3269804       # MeV fm
hbar_s = 6.582119569e-22   # MeV s
c_cm   = 2.99792458e10     # cm/s
me     = 0.51099895        # MeV
mN     = 938.918           # MeV (nucleon average)
m0     = me * 137.035999   # node mass ~ 70 MeV
Tgeom  = m0 / np.sqrt(2.0) # k_B T_geom = hbar c/L4, L4 = sqrt(2) ell: ~49.5 MeV
zeta3  = 1.2020569
M_PL   = 1.22091e22        # Planck mass [MeV]
ell_cm = 2.8179403262e-13  # lattice spacing r_e [cm]

# ---- ladder inputs (closed forms of cosmology/vortex_hadron_ladder.py) ------
F_S = 0.99932              # the Leggett ceiling adopted as the working superfluid fraction
E1_GEV = 2 * np.pi**2 * F_S * m0 / 1e3            # ground-ring energy [GeV]
X0 = brentq(lambda x: x * (np.log(8 * x) - 2) - 1, 1.0, 3.0)   # R0/ell, f_s-free
P0_GEV = E1_GEV * X0**2                            # ground-ring impulse x c [GeV]
FLOOR_GEV = E1_GEV + P0_GEV                        # min_R [E + P c], at R0
MN_GEV = mN / 1e3


def n_over_ngamma_nr(m, T, g):
    """Equilibrium (particle)/(photon) number ratio, non-relativistic species
    of internal degeneracy g at zero chemical potential:
    n = g (m T / 2 pi)^{3/2} e^{-m/T},  n_gamma = (2 zeta3 / pi^2) T^3."""
    return g * (m / (2 * np.pi * T)) ** 1.5 * np.exp(-m / T) * np.pi**2 / (2 * zeta3)


def n_epm_cm3(T):
    """e+ plus e- number density [cm^-3] at temperature T [MeV], exact
    Fermi-Dirac with the electron mass, zero chemical potential, g = 2 each."""
    a = me / T
    I = quad(lambda x: x**2 / (np.exp(np.sqrt(x * x + a * a)) + 1.0), 0, 60)[0]
    n_fm3 = 2 * (2 / (2 * np.pi**2)) * (T / hbarc) ** 3 * I
    return n_fm3 * 1e39


def hubble_rate(T, gstar):
    """Radiation-era expansion rate H = 1.66 sqrt(g_*) T^2 / M_Pl [1/s]."""
    return 1.66 * np.sqrt(gstar) * T**2 / M_PL / hbar_s


def report():
    W = 72
    print("=" * W)
    print("The early-universe ring forge: the floor, the flight, the sink")
    print("=" * W)

    print("\n[1] The floor: E + |P|c >= min_R [E(R) + P(R)c] in the condensate frame")
    print(f"    E1 = 2 pi^2 f_s m0 c^2 = {E1_GEV:.3f} GeV,  R0 = {X0:.4f} ell,"
          f"  P0 c = E1 (R0/ell)^2 = {P0_GEV:.2f} GeV")
    print(f"    floor = E1 [1 + (R0/ell)^2] = {FLOOR_GEV:.3f} GeV")
    for name, s in [("e+e-", 2 * me), ("mu pair", 211.3), ("pi pair", 279.1),
                    ("K pair", 987.4), ("N-Nbar", 2 * mN)]:
        print(f"    {name:8s} at rest: E + |P|c = {s/1e3:6.3f} GeV  "
              f"vs floor {FLOOR_GEV:.2f} GeV  -> below")

    print("\n[2] In flight: a co-moving nucleon-antinucleon pair")
    r = FLOOR_GEV / (2 * MN_GEV)          # required sqrt((1+b)/(1-b))
    beta = (r**2 - 1) / (r**2 + 1)
    gamma = (r**2 + 1) / (2 * r)
    KE_pair = 2 * MN_GEV * (gamma - 1)
    print(f"    sqrt((1+beta)/(1-beta)) = {FLOOR_GEV:.3f} / {2*MN_GEV:.3f} = {r:.3f}")
    print(f"    -> beta = {beta:.4f}, gamma = {gamma:.4f}, "
          f"combined kinetic energy 2 m_N c^2 (gamma-1) = {KE_pair:.3f} GeV")
    # fixed target for comparison: sqrt(p^2+m^2) + m + p = floor
    B = FLOOR_GEV - MN_GEV
    p_ft = (B**2 - MN_GEV**2) / (2 * B)
    KE_ft = np.sqrt(p_ft**2 + MN_GEV**2) - MN_GEV
    print(f"    (antinucleon on a nucleon at rest: p = {p_ft:.2f} GeV/c, "
          f"kinetic energy {KE_ft:.2f} GeV; co-moving is the minimum)")
    nbar = n_over_ngamma_nr(mN, Tgeom, g=4)
    print(f"    equilibrium antinucleons per photon at T_geom = {Tgeom:.1f} MeV"
          f" (g = 4): {nbar:.1e}")
    print("    Nucleons and antinucleons are still in chemical equilibrium here,")
    print("    so kinematics alone does not close the channel. The sink does.")

    print("\n[3] The sink (linear, bath-fed): ring lifetime against the Hubble time")
    sig_geo = np.pi * ell_cm**2            # cm^2, core overlap
    P_fuse = 1e-15                         # minimal fusion probability per contact
    gstar = 10.75                          # photons, e+e-, three neutrinos
    print(f"    sigma = pi ell^2 = {sig_geo:.2e} cm^2, P_fuse = {P_fuse:.0e},"
          f" g_* = {gstar}")
    for T in (Tgeom, 25.0, 10.0):
        n_pm = n_epm_cm3(T)
        G = n_pm * sig_geo * P_fuse * c_cm
        H = hubble_rate(T, gstar)
        print(f"    T = {T:5.1f} MeV: n_pm = {n_pm:.2e} cm^-3, lifetime {1/G:.1e} s"
              f" vs Hubble time {1/H:.1e} s  (Gamma/H = {G/H:.1e})")
    H_cons = hubble_rate(Tgeom, 17.25)     # muons and pions fully relativistic
    G_geom = n_epm_cm3(Tgeom) * sig_geo * P_fuse * c_cm
    print(f"    conservative g_* = 17.25 at T_geom: Gamma/H = {G_geom/H_cons:.1e}")
    print("    A WIMP outlives freeze-out because it must meet another WIMP;")
    print("    a ring dies on anything charged, and charge never runs out.")

    print("\n[4] Born hot versus born cold: the executioner density")
    rows = [("e+e- bath at 1 MeV", n_epm_cm3(1.0)),
            ("interstellar medium", 1.0),
            ("intergalactic medium", 1e-7)]
    for label, n in rows:
        mfp = 1 / (n * sig_geo)
        print(f"    {label:22s} n = {n:9.2e} cm^-3   "
              f"mean free path = {mfp:.1e} cm ({mfp/3.086e24:.1e} Mpc)")
    print("    The cores are the condensate-vortex halves of electrons and")
    print("    protons. They never disappeared; their density fell by ~32")
    print("    orders to the interstellar medium and ~39 to the intergalactic,")
    print("    and the ring's fate is set by the room it is born into.")
    print("=" * W)


if __name__ == "__main__":
    report()
