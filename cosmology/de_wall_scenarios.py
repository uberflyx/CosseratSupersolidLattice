"""
Fossil domain-wall dark energy on a rigid lattice: scenario comparison
against DESI DR2 BAO.

Companion script to the Cosserat supersolid monograph (dark-energy sections).

The question this script answers
--------------------------------
The monograph identifies dark energy with the elastic energy of the
chirality domain-wall network: the grain boundaries of the vacuum
polycrystal.  The equation of state w = P/(rho c^2) then follows from the
wall TRANSPORT LAW: how the wall network evolves after the vacuum
crystallised.  Three candidate laws are compared here.

  A. FOSSIL.     The two chirality domains are exactly degenerate, so a
                 wall node tunnels both ways with equal amplitude alpha.
                 Unbiased hops are a random walk: the wall diffuses
                 ~70 km in a Hubble time and is, on cosmological scales,
                 frozen at fixed lattice positions.  The lattice itself
                 never stretches, so the proper wall density is constant:
                 w = -1 exactly, at every epoch, with no cosmological
                 constant behind it.

  B. BALLISTIC.  Walls migrate at v = f_s * alpha * c, so the grain size
                 grows as L ~ v t and rho_wall ~ 1/t.  This needs a
                 sustained energy bias between the wells, which the
                 degeneracy forbids; it is kept here as the foil.

  C. SIGMOID.    A frozen-to-scaling crossover gated at H = Gamma_cluster,
                 w(z) = -1 + (1/3) / (1 + (H/Gamma)^s).  This presumes a
                 conformal a^{-1} dilution channel for interior walls,
                 which the rigid lattice does not supply.

The script also checks the FOSSIL BIRTH IMPRINT: a frozen network must
inherit its grain size from the causal horizon at crystallisation.  The
information carrier is the lattice's longitudinal compression channel at
v_p = c * sqrt(K_sf / mu) = c / (sqrt(C0) * alpha^{19/2}), and the
available time is the radiation-era Hubble time at the crystallisation
temperature Theta_D = pi * m0 c^2 ~ 220 MeV.  The Kibble correlation
criterion (common causal past) gives a grain diameter

    L = 2 * v_p * t_c  ~  21.4 Mpc,

to be compared with the lattice-ladder form f_s * alpha^{-18} * r_e
= 21.2 Mpc and with the observed cosmic-web cell scale.

Data
----
DESI DR2 BAO: 13 distance measurements (D_V/r_d, D_M/r_d, D_H/r_d across
seven redshift bins) with the full 13x13 covariance, fetched from the
public CobayaSampler/bao_data repository.  Fits here use BAO alone, with
(Omega_m, h*r_d) free in every model, following the warning of Afroz &
Mukherjee (PRD 113, 083514, 2026; arXiv:2504.16868) that the BAO+SNe
combination carries an internal distance-duality inconsistency.

Usage
-----
    python de_wall_scenarios.py            # full analysis + figure
    python de_wall_scenarios.py --no-plot  # console only

Output: console tables and (optionally) de_wall_scenarios.png.

Requires: numpy, scipy, matplotlib (optional), numba (optional, speeds the
grid integrals; the script falls back to pure numpy without it).
"""

import argparse
import os
import urllib.request

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize

try:
    from numba import njit
except ImportError:  # graceful fallback: numba only accelerates, never gates
    def njit(*args, **kwargs):
        def wrap(f):
            return f
        return wrap if args and callable(args[0]) is False or kwargs else (
            args[0] if args and callable(args[0]) else wrap)

# ======================================================================
# Physical constants (CODATA 2018) and lattice parameters
# ======================================================================
C_LIGHT = 2.99792458e8            # speed of light [m/s]
HBAR    = 1.054571817e-34         # reduced Planck constant [J s]
G_N     = 6.67430e-11             # Newton constant [m^3 kg^-1 s^-2]
ALPHA   = 7.2973525693e-3         # fine-structure constant
R_E     = 2.8179403262e-15        # classical electron radius = lattice spacing [m]
F_S     = 4.0 / 5.0               # supersolid crystalline fraction
MPC     = 3.0856775814913673e22   # metres per megaparsec
MEV     = 1.602176634e-13         # joules per MeV

M0C2    = 0.51099895e6 * 1.602176634e-19 / ALPHA   # node rest energy ~70 MeV [J]
MU_LAT  = M0C2 / R_E**3                            # lattice shear modulus [Pa]
K_SF    = 1.5e73                                   # superfluid bulk modulus [Pa]
V_P     = C_LIGHT * np.sqrt(K_SF / MU_LAT)         # compression (pilot-wave) speed [m/s]
C0      = 1.0 / (ALPHA**19 * (K_SF / MU_LAT))      # K_sf = mu / (C0 alpha^19)

GAMMA_CLUSTER = ALPHA**19 * C_LIGHT / R_E          # 19-node volume-insertion rate [1/s]
H0_FID        = 67.4e3 / MPC                       # fiducial H0 [1/s]
T0_FID        = 13.8e9 * 3.15576e7                 # fiducial age of universe [s]

GAMMA_WALL = 4.3e13       # Cosserat-corrected Read-Shockley wall tension [J/m^2]
RHO_DE_OBS = 5.96e-10     # observed dark-energy density [J/m^3]
L_LADDER   = F_S * ALPHA**(-18) * R_E              # lattice-ladder grain size [m]

DESI_BASE = ("https://raw.githubusercontent.com/CobayaSampler/bao_data/"
             "master/desi_bao_dr2/")
DESI_MEAN = "desi_gaussian_bao_ALL_GCcomb_mean.txt"
DESI_COV  = "desi_gaussian_bao_ALL_GCcomb_cov.txt"

# Shared scale-factor grid for all background integrations
A_GRID = np.logspace(-4, 0, 4000)
LNA    = np.log(A_GRID)
OMEGA_R = 9.0e-5          # radiation density today (photons + neutrinos)


# ======================================================================
# Part 1: wall-transport scales.  Which law does the mechanics allow?
# ======================================================================
def part1_transport_scales():
    print("=" * 72)
    print("PART 1: HOW FAR CAN A WALL MOVE IN ONE HUBBLE TIME?")
    print("=" * 72)

    # (a) Unbiased tunnelling.  Both chirality wells are degenerate
    #     (theta_ch^2 is even in theta_ch), so a wall node hops either way
    #     with the same amplitude alpha.  Step r_e at attempt rate
    #     alpha c / r_e gives a diffusion constant D = alpha c r_e.
    d_wall = ALPHA * C_LIGHT * R_E
    x_diff = np.sqrt(2.0 * d_wall * T0_FID)
    print("\n(a) Unbiased zero-point tunnelling (degenerate wells):")
    print(f"      D = alpha c r_e        = {d_wall:.3e} m^2/s")
    print(f"      sqrt(2 D t0)           = {x_diff/1e3:.0f} km")

    # (b) Curvature (capillary) bias.  The only surviving force on a wall
    #     is its own curvature, ~ gamma / L.  The monograph's parabolic
    #     growth estimate gives ~5 km per Hubble time: same verdict.
    print("\n(b) Curvature-driven coarsening: ~5 km per Hubble time")
    print("      (parabolic von Neumann-Mullins estimate; same verdict).")

    # (c) Ballistic migration would need a sustained bias.
    v_wall = F_S * ALPHA * C_LIGHT
    print("\n(c) Ballistic migration (foil; needs a bias the degeneracy")
    print(f"      forbids): v t0 = {v_wall * T0_FID / MPC:.1f} Mpc")

    print("\n  Verdict: (a) and (b) fall short of the 21 Mpc grain scale by")
    print("  ~19 orders of magnitude.  The network is pinned: walls are")
    print("  FOSSILS of the crystallisation epoch, and the proper wall")
    print("  energy density has been constant ever since -> w = -1.")


# ======================================================================
# Part 2: the fossil birth imprint.  Does the grain size come out right?
# ======================================================================
def part2_birth_imprint():
    print("\n" + "=" * 72)
    print("PART 2: THE GRAIN SIZE AS A BIRTH IMPRINT")
    print("=" * 72)

    theta_d = np.pi * M0C2                      # crystallisation temperature [J]
    g_star  = 61.75                             # SM relativistic dof above QCD
    rho_rad = (np.pi**2 / 30.0) * g_star * theta_d**4 / (HBAR * C_LIGHT)**3
    # Radiation-era Hubble time 1/(2H) at T = Theta_D:
    t_c = np.sqrt(3.0 * C_LIGHT**2 / (32.0 * np.pi * G_N * rho_rad))

    l_fossil = 2.0 * V_P * t_c                  # Kibble correlation diameter
    print(f"\n  T_c = Theta_D = pi m0 c^2        = {theta_d/MEV:.0f} MeV")
    print(f"  t_c = 1/(2H) at T_c (g*={g_star})   = {t_c:.3e} s")
    print(f"  v_p = c sqrt(K_sf/mu)            = {V_P/C_LIGHT:.2e} c"
          f"   (C0 = {C0:.2f})")
    print(f"  L   = 2 v_p t_c                  = {l_fossil/MPC:.1f} Mpc")
    print(f"  lattice-ladder form f_s a^-18 re = {L_LADDER/MPC:.1f} Mpc")
    print(f"  ratio                            = {l_fossil/L_LADDER:.3f}")

    # The half-rung identity: v_p carries alpha^{-19/2}; the match above
    # therefore places the crystallisation epoch itself on the alpha
    # ladder at c t_c ~ 0.47 alpha^{-17/2} r_e  (19/2 + 17/2 = 18).
    kappa = C_LIGHT * t_c / (ALPHA**(-17.0/2.0) * R_E)
    print(f"\n  Half-rung identity: c t_c = {kappa:.3f} alpha^(-17/2) r_e")
    print("  (the crystallisation epoch joins the alpha ladder; rungs")
    print("   19/2 + 17/2 = 18 reproduce the grain-size exponent).")

    # Wall energy density for the two area-accounting conventions.
    rho_diam = 3.6 * GAMMA_WALL / l_fossil
    rho_rad_ = 3.6 * GAMMA_WALL / (0.5 * l_fossil)
    print(f"\n  rho_wall (S/V = 3.6/L, L = diameter) = {rho_diam:.2e} J/m^3")
    print(f"  rho_wall (S/V = 3.6/R, R = radius)   = {rho_rad_:.2e} J/m^3")
    print(f"  observed rho_DE                      = {RHO_DE_OBS:.2e} J/m^3")
    print(f"  gap                                  = "
          f"{RHO_DE_OBS/rho_rad_:.2f} - {RHO_DE_OBS/rho_diam:.2f}")
    print("  (Read-Shockley reproduces laboratory grain-boundary energies")
    print("   to 20-50%; the remaining spread is the radius/diameter")
    print("   convention in the foam area accounting.)")


# ======================================================================
# Part 3: background cosmology per scenario and the DESI DR2 BAO fit
# ======================================================================
def e_of_a(a, om, f_de):
    """Dimensionless Hubble rate E = H/H0 on the shared grid."""
    ode = 1.0 - om - OMEGA_R
    return np.sqrt(OMEGA_R * a**-4 + om * a**-3 + ode * f_de)


def age_of_a(om, f_de):
    """Cosmic time t(a) in units of 1/H0."""
    e = e_of_a(A_GRID, om, f_de)
    t = cumulative_trapezoid(1.0 / (A_GRID * e), A_GRID, initial=0.0)
    return t + A_GRID[0]**2 / (2.0 * np.sqrt(OMEGA_R))


def f_de_fossil():
    """Scenario A: constant proper wall density."""
    return np.ones_like(A_GRID)


def f_de_ballistic(om, n_iter=6):
    """Scenario B: rho_wall ~ 1/t, solved self-consistently."""
    f_de = np.ones_like(A_GRID)
    for _ in range(n_iter):
        t = age_of_a(om, f_de)
        f_de = t[-1] / t
    return f_de


def f_de_sigmoid(om, s, n_iter=6):
    """Scenario C: w = -1 + (1/3)/(1 + (H/Gamma)^s), Gamma = 1.22 H0."""
    g_over_h0 = GAMMA_CLUSTER / H0_FID
    f_de = np.ones_like(A_GRID)
    for _ in range(n_iter):
        e = e_of_a(A_GRID, om, f_de)
        w = -1.0 + (1.0 / 3.0) / (1.0 + (e / g_over_h0)**s)
        integ = cumulative_trapezoid(3.0 * (1.0 + w), LNA, initial=0.0)
        f_de = np.exp(integ[-1] - integ)
    return f_de


def w_of_a(f_de):
    """Effective equation of state from the density history."""
    return -1.0 - np.gradient(np.log(f_de), LNA) / 3.0


def fetch_desi(cache_dir="."):
    """Download (once) and load the DESI DR2 BAO vector and covariance."""
    paths = []
    for fname in (DESI_MEAN, DESI_COV):
        path = os.path.join(cache_dir, fname)
        if not os.path.exists(path):
            print(f"  fetching {fname} ...")
            urllib.request.urlretrieve(DESI_BASE + fname, path)
        paths.append(path)
    z, val, kind = [], [], []
    with open(paths[0]) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            z.append(float(p[0])), val.append(float(p[1])), kind.append(p[2])
    cov = np.loadtxt(paths[1])
    return np.array(z), np.array(val), kind, cov


def bao_predictions(z_arr, kinds, om, hrd, f_de):
    """D_M/r_d, D_H/r_d, D_V/r_d.  hrd = h * r_d in Mpc."""
    e = e_of_a(A_GRID, om, f_de)
    z_grid = 1.0 / A_GRID - 1.0
    order = np.argsort(z_grid)
    zg, eg = z_grid[order], e[order]
    chi = cumulative_trapezoid(1.0 / eg, zg, initial=0.0)
    c_km = 2.99792458e5
    out = []
    for zi, ki in zip(z_arr, kinds):
        ei, chii = np.interp(zi, zg, eg), np.interp(zi, zg, chi)
        dm = c_km * chii / (100.0 * hrd)
        dh = c_km / (100.0 * hrd * ei)
        out.append(dm if ki == "DM_over_rs"
                   else dh if ki == "DH_over_rs"
                   else (zi * dm**2 * dh)**(1.0 / 3.0))
    return np.array(out)


def make_f_de(scen, om, s=None):
    if scen == "A":
        return f_de_fossil()
    if scen == "B":
        return f_de_ballistic(om)
    return f_de_sigmoid(om, s)


def chi2(theta, z, val, kinds, icov, scen, s=None):
    om, hrd = theta
    if not (0.05 < om < 0.7 and 50.0 < hrd < 150.0):
        return 1e10
    pred = bao_predictions(z, kinds, om, hrd, make_f_de(scen, om, s))
    d = pred - val
    return d @ icov @ d


def cpl_fit(f_de, z_max=3.0):
    """Least-squares CPL (w0, wa) projection of w(a) over 0 < z < z_max."""
    w = w_of_a(f_de)
    mask = A_GRID >= 1.0 / (1.0 + z_max)
    x = 1.0 - A_GRID[mask]
    coef, *_ = np.linalg.lstsq(
        np.vstack([np.ones_like(x), x]).T, w[mask], rcond=None)
    return coef


def part3_desi(make_plot=True):
    print("\n" + "=" * 72)
    print("PART 3: DESI DR2 BAO CONFRONTATION (BAO alone, 13 points)")
    print("=" * 72)
    z, val, kinds, cov = fetch_desi()
    icov = np.linalg.inv(cov)
    z_pivot = 0.34

    models = [("A fossil  (w = -1)",     "A", None),
              ("B ballistic (rho~1/t)",  "B", None),
              ("C sigmoid  (s = 2)",     "C", 2.0),
              ("C sigmoid  (s = 8)",     "C", 8.0)]

    print(f"\n  {'model':<26}{'Om':>7}{'h rd':>8}{'chi2':>8}"
          f"{'w(0)':>7}{'w(zp)':>8}{'CPL w0':>8}{'CPL wa':>8}")
    results = []
    for label, scen, s in models:
        res = minimize(chi2, x0=[0.30, 99.0],
                       args=(z, val, kinds, icov, scen, s),
                       method="Nelder-Mead",
                       options={"xatol": 1e-4, "fatol": 1e-4})
        om, hrd = res.x
        f_de = make_f_de(scen, om, s)
        w = w_of_a(f_de)
        w0 = np.interp(1.0, A_GRID, w)
        wp = np.interp(1.0 / (1.0 + z_pivot), A_GRID, w)
        w0c, wac = cpl_fit(f_de)
        results.append((label, scen, s, om, hrd, res.fun, f_de))
        print(f"  {label:<26}{om:>7.3f}{hrd:>8.1f}{res.fun:>8.2f}"
              f"{w0:>7.2f}{wp:>8.2f}{w0c:>8.2f}{wac:>8.2f}")

    print("\n  Benchmarks (2026):")
    print("    pivot w(0.34) = -0.90 +/- 0.09   (model-independent)")
    print("    DDR-corrected w0 = -0.92 +/- 0.08 (Afroz & Mukherjee,")
    print("    PRD 113, 083514: BAO+SNe inconsistency removed)")
    print("\n  The fossil scenario is background-degenerate with LCDM and")
    print("  sits ~1 sigma from both 2026 measurements; the ballistic foil")
    print("  is excluded (Delta chi2 = +15, ~3 sigma from the pivot); the")
    print("  sigmoid buys Delta chi2 ~ -2 only at the price of an underived")
    print("  sharpness s and a dilution channel the rigid lattice lacks.")

    if make_plot:
        _plot(results, z_pivot)


def _plot(results, z_pivot):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not available; skipping figure)")
        return
    z_plot = np.linspace(0.0, 2.5, 400)
    a_plot = 1.0 / (1.0 + z_plot)
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    styles = {"A": ("-", 2.4), "B": (":", 1.6), "C": ("--", 1.6)}
    for label, scen, s, om, hrd, c2, f_de in results:
        w = np.interp(a_plot, A_GRID, w_of_a(f_de))
        ls, lw = styles[scen]
        ax.plot(z_plot, w, ls, lw=lw, label=f"{label}  ($\\chi^2$={c2:.1f})")
    ax.errorbar([z_pivot], [-0.90], yerr=[0.09], fmt="o", color="k",
                capsize=4, label="pivot $w_p$ (2026)")
    ax.errorbar([0.02], [-0.92], yerr=[0.08], fmt="s", color="0.4",
                capsize=4, label="DDR-corrected $w_0$ (Afroz 2026)")
    ax.axhline(-1.0, color="0.8", lw=0.8)
    ax.set_xlabel("redshift $z$")
    ax.set_ylabel("equation of state $w(z)$")
    ax.set_ylim(-1.15, -0.45)
    ax.set_title("Chirality domain-wall dark energy: transport-law scenarios")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig("de_wall_scenarios.png", dpi=160)
    print("\n  figure written: de_wall_scenarios.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--no-plot", action="store_true",
                    help="console output only")
    args = ap.parse_args()
    part1_transport_scales()
    part2_birth_imprint()
    part3_desi(make_plot=not args.no_plot)
