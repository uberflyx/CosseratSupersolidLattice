#!/usr/bin/env python3
r"""
neff_longitudinal_modes.py
==========================
Dark radiation from the vacuum supersolid's second sound: thermal contact,
entropy sharing during electron-positron annihilation, the relic N_eff, and
(optionally) the light-element abundances.

THE QUESTION
------------
Early-universe annihilation is a monopole source of second sound, and its
time-reverse (pair creation) absorbs it, so the mode is driven towards the
plasma temperature. How much radiation does it leave behind, and what does
it do to N_eff and to nucleosynthesis?

THE PHYSICS, IN FOUR STEPS
--------------------------
1. CONTACT. The exchange rate with the plasma is Gamma = f n_e+ sigma*v,
   with f the per-event branching into second sound. Contact (Gamma/H > 1)
   is lost only when the positron density crashes, at 14 keV for f = 1e-2,
   22 keV for the positronium-limited f ~ 1e-8 and 82 keV for f = 1e-15.

2. CAPACITY. A gapless mode with omega = v k holds the entropy and energy
   of g = (c/v)^3 photon polarisations, because a faster mode packs fewer
   states into a band of wavenumbers. Second sound runs at
   v2 = sqrt(C11/C44) c = 1.5916 c, so g = 0.2480.

3. ENTROPY SHARING. After neutrino decoupling (~2 MeV) the electromagnetic
   plasma conserves its entropy on its own. While second sound is in
   contact it belongs to that plasma, so it takes part of the heat the
   pairs release (between about 500 and 50 keV):
       (T_gamma/T_nu)^3 = (2 + 7/2 + g)/(2 + g) = 2.557   (standard 2.75).
   The photons are heated less, so the neutrinos end up warmer relative to
   them. (Same separate-entropy bookkeeping as Boehm, Dolan & McCabe,
   JCAP 2013(08) 041, eqs. 1 and 6.)

4. RELIC. With N_eff defined by rho_rad = rho_gamma [1 + (7/8)(4/11)^{4/3} N_eff],
       N_eff = N_nu (11/4)^{4/3} (T_nu/T_gamma)^4 + (4/7)(11/4)^{4/3} g
             = 3.354 + 0.546 = 3.900,    Delta N_eff = 0.856.
   The second term alone (0.546) is the mode's own energy; the first term's
   excess (0.310) is the neutrinos' gain. If the mode decoupled before
   annihilation it would sit at T_nu and give only (4/7) g = 0.142.

LITERATURE INPUTS (checked against the journal or arXiv text)
-------------------------------------------------------------
  sigma*v (e+e- -> 2 gamma, slow pairs) = pi r_e^2 c         Dirac (1930)
  N_eff^SM = 3.044                                           Bennett+ 2021
  Planck 2018:   N_eff = 2.99 +/- 0.17 (68%)                 arXiv:1807.06209
  ACT DR6 P-ACT-LB: N_eff = 2.86 +/- 0.13 (68%), eq. 30;
                 Delta N_eff < 0.17 (95%), eq. 32;
                 N_idr < 0.114 (95%, P-ACT), eq. 55          JCAP 2025(11) 063,
                                                             arXiv:2503.14454
  Y_p = 0.2449 +/- 0.0040                                    Aver+ 2015
  D/H = (2.527 +/- 0.030)e-5                                 Cooke+ 2018

NUCLEOSYNTHESIS (optional, needs PRyMordial)
--------------------------------------------
PRyMordial (Burns, Tait & Valli, EPJC 84, 86 (2024); github.com/vallima/PRyMordial)
follows neutrino decoupling and the full nuclear network. Second sound enters
as a massless species in equilibrium with the photons and pairs (NP_e_flag).
One setting matters: the code fixes the late-time scale factor by assuming
today's entropy is in photons alone (s0bar), which is right for a species
that annihilates but wrong for one that stays massless. Setting s0bar to the
photons-plus-mode value restores a T_gamma = T0 after nucleosynthesis, which is
what the paper's own eq. (2.8) gives. Results at Omega_b h^2 = 0.02230:
    PRIMAT rates:   SM  N_eff 3.044, Y_p 0.2469, D/H 2.45e-5
                    2nd N_eff 3.898, Y_p 0.2516, D/H 2.71e-5
    NACRE II rates: SM  Y_p 0.2468, D/H 2.51e-5
                    2nd Y_p 0.2515, D/H 2.78e-5
    The D/H shift, +0.26e-5 (about +10.5%), is the same for both compilations.
Run with:  python neff_longitudinal_modes.py --prym /path/to/PRyMordial

Plain arithmetic and 1-D quadrature; numba would add nothing here.
"""

import argparse
import os
import sys

import numpy as np
from scipy.integrate import quad

# ---------------- constants (natural units: MeV, fm, s) --------------------
hbar_MeV_s = 6.582119569e-22          # hbar [MeV s]
hbarc      = 197.3269804              # hbar c [MeV fm]
c_fm_s     = 2.99792458e23            # speed of light [fm/s]
me         = 0.51099895               # electron mass [MeV]
re         = 2.8179403262             # classical electron radius [fm]
Mpl        = 1.22091e22               # Planck mass [MeV]
N_NU_SM    = 3.044                    # Standard Model N_eff

N2_COUPLING = 1.0 / np.pi             # Cosserat coupling number
V2_OVER_C   = ((3 - 4 * N2_COUPLING) / (1 - N2_COUPLING)) ** 0.5
                                      # second sound / c = sqrt(C11/C44) = 1.5916;
                                      # C11/C44 = 2.533 rather than 3 because the
                                      # contact stiffness breaks the Cauchy relation
G_MODE      = V2_OVER_C ** -3         # entropy/energy dof of the mode, 0.2480

SIGMA_V = np.pi * re**2 * c_fm_s      # slow-pair annihilation coefficient [fm^3/s]


# ---------------- electron-positron thermodynamics -------------------------
def _fermi_moments(T):
    """Number density of positrons [fm^-3] and the e+e- entropy in units of
    one bosonic degree of freedom, s_e / (2 pi^2 T^3 / 45), at temperature T [MeV]."""
    x = me / T

    def n_int(y):          # y = E/T
        return y * np.sqrt(y * y - x * x) / (np.exp(y) + 1.0)

    def s_int(y):          # (rho + P)/T^4 integrand: sqrt(y^2-x^2) (y^2 + (y^2-x^2)/3)
        q2 = y * y - x * x
        return np.sqrt(q2) * (y * y + q2 / 3.0) / (np.exp(y) + 1.0)

    n_pos = 2.0 / (2 * np.pi**2) * T**3 * quad(n_int, x, x + 80, limit=200)[0] / hbarc**3
    s_e = 4.0 / (2 * np.pi**2) * quad(s_int, x, x + 80, limit=200)[0]
    return n_pos, s_e / (2 * np.pi**2 / 45)


def g_s_pairs(T):
    """e+e- entropy in boson-dof units: 3.5 for T >> me, 0 for T << me."""
    return _fermi_moments(T)[1]


def hubble(T):
    """Expansion rate H [1/s] in radiation domination (Standard Model content;
    the mode's own small contribution does not move the decoupling temperature)."""
    ge = g_s_pairs(T)                                   # entropy ~ energy dof here
    tnu = ((2.0 + ge) / 5.5) ** (1.0 / 3.0)             # T_nu / T in the SM
    gstar = 2.0 + ge + (7.0 / 8.0) * 6.0 * tnu**4
    return 1.66 * np.sqrt(gstar) * T**2 / Mpl / hbar_MeV_s


# ---------------- 1. contact and decoupling ---------------------------------
def contact_ratio(T, f_br):
    """Gamma/H for exchange of second-sound quanta through the positron bath."""
    n_pos, _ = _fermi_moments(T)
    return f_br * n_pos * SIGMA_V / hubble(T)


def T_decouple(f_br):
    """Temperature [MeV] at which Gamma/H falls through 1 (log bisection)."""
    lo, hi = np.log(1e-3), np.log(5.0)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if contact_ratio(np.exp(mid), f_br) > 1.0:
            hi = mid
        else:
            lo = mid
    return np.exp(0.5 * (lo + hi))


# ---------------- 2-4. entropy sharing and the relic ------------------------
def neff_with_mode(T_dec, g=G_MODE, n_nu=N_NU_SM):
    """N_eff when the mode shares the plasma's entropy until T_dec [MeV].

    Comoving entropies, with a T = 1 before annihilation (T_nu stays at a T_nu = 1):
      contact phase:  (2 + g_e(T) + g) (aT)^3 = 5.5 + g
      after T_dec:    the mode keeps (aT_dec)^3; photons take the rest of the
                      pairs' entropy, 2 (aT_gamma)^3 = (2 + g_e(T_dec)) (aT_dec)^3.
    """
    ge = g_s_pairs(T_dec)
    aTd3 = (5.5 + g) / (2.0 + ge + g)
    aTg3 = (2.0 + ge) * aTd3 / 2.0
    tx_over_tg4 = (aTd3 / aTg3) ** (4.0 / 3.0)
    nu_part = n_nu * (11.0 / 4.0) ** (4.0 / 3.0) * aTg3 ** (-4.0 / 3.0)
    mode_part = (4.0 / 7.0) * (11.0 / 4.0) ** (4.0 / 3.0) * g * tx_over_tg4
    return nu_part, mode_part


# ---------------- optional: nucleosynthesis with PRyMordial -----------------
def run_prymordial(path):
    """SM and second-sound runs of PRyMordial; returns two result arrays
    [N_eff, Omega_nu h2 1e6, ..., Y_p(BBN), D/H 1e5, He3/H 1e5, Li7/H 1e10]."""
    cwd = os.getcwd()
    os.chdir(path)                       # PRyMordial reads its tables relative to cwd
    sys.path.insert(0, path)
    try:
        import PRyM.PRyM_init as PRyMini
        PRyMini.numba_flag = True
        PRyMini.aTid_flag = True
        PRyMini.compute_bckg_flag = True
        PRyMini.compute_nTOp_flag = True
        import PRyM.PRyM_main as PRyMmain
        g = G_MODE
        sm = PRyMmain.PRyMclass().PRyMresults()
        PRyMini.NP_e_flag = True
        # today's entropy holds the photons and the still-massless mode
        PRyMini.s0bar = 4 * np.pi**2 / 45 * (2.0 + g) / 2.0
        rho = lambda T: g * np.pi**2 / 30 * T**4
        pres = lambda T: rho(T) / 3.0
        drho = lambda T: 4.0 * rho(T) / T
        ss = PRyMmain.PRyMclass(rho, pres, drho).PRyMresults()
    finally:
        os.chdir(cwd)
    return sm, ss


# ---------------- report ----------------------------------------------------
def report(prym_path=None):
    W = 72
    print("=" * W)
    print("Second sound as thermal dark radiation")
    print("=" * W)
    print(f"v2/c = {V2_OVER_C:.4f},  g = (c/v2)^3 = {G_MODE:.4f}")

    print("\n[1] Contact and decoupling  (per-event branching f)")
    print(f"    {'f':>8} {'T_dec [keV]':>12} {'N_eff':>8} {'Delta':>8}")
    for f in [1e-2, 1e-8, 1e-15]:
        Td = T_decouple(f)
        nu, mo = neff_with_mode(Td)
        print(f"    {f:8.0e} {1e3 * Td:12.1f} {nu + mo:8.4f} {nu + mo - N_NU_SM:8.4f}")

    print("\n[2] The relic, decoupling after annihilation")
    ratio = (5.5 + G_MODE) / (2.0 + G_MODE)
    nu, mo = neff_with_mode(1e-3)
    print(f"    (T_gamma/T_nu)^3 = {ratio:.4f}   (standard 2.75)")
    print(f"    neutrinos {nu:.4f} + mode {mo:.4f} = N_eff {nu + mo:.4f}")
    print(f"    mode decoupled before annihilation: Delta N_eff = {4 / 7 * G_MODE:.4f}")

    print("\n[3] Against the data")
    neff = nu + mo
    for name, mu, sig in [("Planck 2018", 2.99, 0.17),
                          ("ACT DR6 P-ACT-LB", 2.86, 0.13)]:
        print(f"    {name:<18} {mu:.2f} +/- {sig:.2f}  ->  {(neff - mu) / sig:4.1f} sigma high")
    print(f"    ACT DR6 bound on an added species: Delta N_eff < 0.17 (95%);")
    print(f"    the prediction {neff - N_NU_SM:.2f} exceeds it {(neff - N_NU_SM) / 0.17:.1f}-fold.")

    print("\n[4] Nucleosynthesis")
    print(f"    expansion-rate share at n/p freeze-out (T_x = T_nu): {4 / 7 * G_MODE:.3f}")
    if prym_path:
        sm, ss = run_prymordial(prym_path)
        for lab, r in [("SM", sm), ("second sound", ss)]:
            print(f"    {lab:<13} N_eff {r[0]:.4f}  Y_p {r[4]:.5f}  D/H {r[5]:.4f}e-5")
        print(f"    shifts: Y_p {ss[4] - sm[4]:+.5f} (obs error 0.0040),"
              f" D/H {ss[5] - sm[5]:+.4f}e-5 (obs error 0.030e-5)")
    else:
        print("    (pass --prym /path/to/PRyMordial for the network calculation)")
    print("=" * W)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    ap.add_argument("--prym", default=None, help="path to a PRyMordial checkout")
    report(ap.parse_args().prym)
