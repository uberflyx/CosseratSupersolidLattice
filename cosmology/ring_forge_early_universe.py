#!/usr/bin/env python3
r"""
ring_forge_early_universe.py
============================
Why the early universe forges no free vortex rings, and why the rings an
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
1. THE MOMENTUM BILL CLOSES THE ANNIHILATION FORGE OUTRIGHT.
   A ring has no rest mass: its energy E(R) and impulse P(R) are independent
   functions of its radius R (see cosmology/vortex_hadron_ladder.py, which
   carries the dispersion). At the adopted superfluid fraction f_s = 0.99932
   the ground ring carries ~3.9 GeV/c of impulse on E_1 = 2 pi^2 f_s m0 c^2
   = 1.38 GeV of energy, so any final state containing a ring costs at least
   sqrt(s) >= min_R [E(R) + P(R) c] = E_1 [1 + (R0/ell)^2] = 5.27 GeV,
   with the balancing momentum carried by light quanta. Nucleon-antinucleon
   annihilation at rest supplies 1.88 GeV: closed, for every species, at any
   temperature. This is kinematics, not suppression.

2. THE IN-FLIGHT ROUTE NEEDS BEAM-GRADE MOMENTA THE BATH DOES NOT HAVE.
   In flight the channel opens at antiproton lab momentum 2.06 GeV/c
   (ladder script). Below k_B T_geom = hbar c / L4 = m0 c^2 / sqrt(2)
   = 49.5 MeV, the warmest epoch at which a ring exists (the compact
   direction must be frozen for a closed circulation loop to be defined),
   thermal nucleons carry ~sqrt(2 m T) ~ 0.3 GeV/c. The threshold momentum
   exceeds m_N c, so its kinetic energy is taken relativistically,
   E_k = sqrt(p^2 c^2 + m_N^2 c^4) - m_N c^2 = 1.33 GeV. The Boltzmann tail
   above it is exp(-E_k/T) ~ 2e-12 per antinucleon, on an antinucleon
   abundance already ~1e-7 per photon, times a branching ceiling of 2e-14:
   of order 1e-32 rings per photon, ever.

3. THE SINK NEVER FREEZES, SO EVEN THOSE DIE.
   A plain ring fuses on the circulation-carrying core of any charged
   particle (the condensate-vortex half of every electron and proton). The
   destruction rate is linear in the ring number and fed by the whole e+-
   bath, so it cannot self-quench the way a WIMP's quadratic annihilation
   does. Lifetime at 25 MeV: under a microsecond against a millisecond
   Hubble time, at even the minimal fusion probability of 1e-15 per contact.

BORN HOT VERSUS BORN COLD
-------------------------
A forged ring today is the same unprotected object. What changed is the
density of executioners. The cores never vanished; their number density fell
from ~3e31 per cm^3 (the e+- pair bath at 1 MeV) to ~1 per cm^3 (interstellar
medium) to ~1e-7 (intergalactic medium). The mean free path grows from a
nanometre to a megaparsec to beyond the observable universe. Same weapon,
same fragile ring; the room emptied.

Light arithmetic only; numba would add nothing.
"""

import numpy as np
from scipy.optimize import brentq

# ---- constants -------------------------------------------------------------
hbarc = 197.3269804        # MeV fm
me    = 0.51099895         # MeV
mN    = 938.918            # MeV (nucleon average)
m0    = me * 137.035999    # node mass ~ 70 MeV
Tgeom = m0 / np.sqrt(2.0)  # k_B T_geom = hbar c/L4, L4 = sqrt(2) ell: ~49.5 MeV
zeta3 = 1.2020569

# ---- ladder inputs (closed forms of cosmology/vortex_hadron_ladder.py) ------
F_S = 0.99932              # the Leggett ceiling adopted as the working superfluid fraction
E1_GEV = 2 * np.pi**2 * F_S * m0 / 1e3            # ground-ring energy [GeV]
X0 = brentq(lambda x: x * (np.log(8 * x) - 2) - 1, 1.0, 3.0)   # R0/ell, f_s-free
FLOOR_GEV = E1_GEV * (1 + X0**2)                   # sqrt(s) floor, E1 + P0 c
# In flight, fixed target: sqrt(p^2 + m^2) + p = FLOOR - m_N (photon thrown
# backward), so p = (B^2 - m_N^2)/(2B) with B = FLOOR - m_N.
_B = FLOOR_GEV - mN / 1e3
PBAR_THR = (_B**2 - (mN / 1e3)**2) / (2 * _B)      # GeV/c, in-flight opening
BRANCH_CEIL = 2.0e-14      # forge branching ceiling (BaBar sigma_shed / sigma_ann)

def n_over_ngamma_nr(m, T):
    """Equilibrium (particle)/(photon) number ratio, non-relativistic species."""
    return 4.106 * (m / (2 * np.pi * T)) ** 1.5 * np.exp(-m / T)

def report():
    W = 72
    print("=" * W)
    print("The early-universe ring forge: closed three ways")
    print("=" * W)

    print("\n[1] The momentum bill (kinematic closure, all species)")
    for name, s in [("e+e-", 2 * me), ("mu pair", 211.3), ("pi pair", 279.1),
                    ("K pair", 987.4), ("N-Nbar", 2 * mN)]:
        print(f"    {name:8s} at rest: sqrt(s) = {s/1e3:6.3f} GeV  "
              f"vs floor {FLOOR_GEV:.2f} GeV  -> closed")
    print("    The floor is final-state independent: the ground ring's own")
    print("    impulse must be balanced, and light balancers cost |P|c.")

    print("\n[2] The in-flight tail below T_geom = %.1f MeV" % Tgeom)
    # relativistic: the threshold momentum exceeds m_N c
    Ek_thr = np.sqrt((PBAR_THR * 1e3)**2 + mN**2) - mN   # MeV
    tail   = np.exp(-Ek_thr / Tgeom)
    nbar   = n_over_ngamma_nr(mN, Tgeom)
    print(f"    beam threshold {PBAR_THR:.2f} GeV/c -> kinetic energy {Ek_thr:.0f} MeV")
    print(f"    Boltzmann tail above it at T_geom: exp(-{Ek_thr/Tgeom:.1f}) = {tail:.1e}")
    print(f"    antinucleons per photon at T_geom: {nbar:.1e}")
    print(f"    rings ever forged per photon <= tail x nbar x branching")
    print(f"      = {tail*nbar*BRANCH_CEIL:.1e}")

    print("\n[3] The sink (linear, bath-fed, never freezes)")
    sig_geo = np.pi * (2.818e-13) ** 2             # cm^2, core overlap
    for T, P in [(25.0, 1e-15)]:
        n_pos = 0.5 * (3 * zeta3 / np.pi**2) * (T / hbarc) ** 3 * 2 * 1e39
        G = 2 * n_pos * sig_geo * P * 3e10
        H = 1.66 * np.sqrt(10.75) * T**2 / 1.22091e22 / 6.582e-22
        print(f"    T={T:.0f} MeV, P_fuse={P:.0e}: ring lifetime {1/G:.0e} s "
              f"vs Hubble {1/H:.0e} s")
    print("    A WIMP outlives freeze-out because it must meet another WIMP;")
    print("    a ring dies on anything charged, and charge never runs out.")

    print("\n[4] Born hot versus born cold: the executioner density")
    def ngamma_cm3(T_MeV):
        return 2 * zeta3 / np.pi**2 * (T_MeV / (hbarc * 1e-13)) ** 3
    rows = [("pair bath at 1 MeV", ngamma_cm3(1.0)),
            ("interstellar medium", 1.0),
            ("intergalactic medium", 1e-7)]
    for label, n in rows:
        mfp = 1 / (n * sig_geo)
        print(f"    {label:22s} n = {n:9.2e} cm^-3   "
              f"mean free path = {mfp:.1e} cm ({mfp/3.086e24:.1e} Mpc)")
    print("    The cores are the condensate-vortex halves of electrons and")
    print("    protons. They never disappeared; their density fell by ~31")
    print("    orders, and the ring's fate is set by the room it is born into.")
    print("=" * W)

if __name__ == "__main__":
    report()
