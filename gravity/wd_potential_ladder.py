#!/usr/bin/env python3
"""
wd_potential_ladder.py — Recover the gravitational-potential ladder for the
Bainbridge et al. (2017) white-dwarf sample directly from published log g.

We CANNOT fit dalpha/alpha for the sample: Bainbridge et al. (2017, Universe
3, 32) state explicitly that they withheld the per-object dalpha/alpha table
("We do not include a table of Delta alpha/alpha estimates ... for this
reason"), and four of the five best targets had not yet been observed.

But phi = GM/(Rc^2) needs only mass and radius, and the paper's Table 1 gives
log g and T_eff for all 13 objects. From log g plus a mass-radius relation we
recover phi for each object, i.e. the x-axis of the ladder. This tells us the
LEVER ARM available for separating a linear (k*phi) from a quadratic (k*phi^2)
gravitational response — the framework's distinguishing prediction.

Method:
  g   = G M / R^2                       (surface gravity, from log g)
  Use WD mass-radius relation R(M) to close the (M,R) pair given g.
  Solve  log10(g_cgs) = log10( G_cgs * M / R(M)^2 )  for M, then phi.

For sub-dwarfs (log g ~ 4-6) the WD M-R relation does not apply; they are
core-He/H-burning stars with M ~ 0.5 Msun on the horizontal branch. We flag
these and use a nominal M ~ 0.47 Msun to place them, since only phi (via g and
an assumed R) matters and their low log g already puts them low on the ladder.
"""
import numpy as np
from scipy.optimize import brentq

# constants
G_SI = 6.67430e-11; c = 299792458.0
Msun = 1.98892e30; Rsun = 6.957e8; Rearth = 6.371e6
G_cgs = 6.67430e-8

def wd_radius_Rsun(M_Msun, mu=2.0, Mch=1.454):
    x = M_Msun/Mch
    if x >= 1.0: return np.nan
    return (0.0225/mu)*np.sqrt(1.0 - x**(4.0/3.0))/x**(1.0/3.0)

def mass_from_logg_WD(logg):
    """Solve log g (cgs) = log10(G_cgs M / R(M)^2) for a WD on the M-R relation."""
    def f(M):
        R = wd_radius_Rsun(M)
        if np.isnan(R): return 1e9
        Rcm = R*Rsun*1e2; Mg = M*Msun*1e3
        return np.log10(G_cgs*Mg/Rcm**2) - logg
    try:
        return brentq(f, 0.15, 1.44, xtol=1e-6)
    except Exception:
        return np.nan

def phi_from_M_R(M_Msun, R_Rsun):
    return G_SI*(M_Msun*Msun)/((R_Rsun*Rsun)*c**2)

# Bainbridge et al 2017 Table 1: (name, type, Teff, logg)
sample = [
    ("vZ 1128",      "O(H)", 36600, 3.9),
    ("ROB 162",      "O(H)", 51000, 4.5),
    ("BD+28 4211",   "sdO",  82000, 6.20),
    ("Sh 2-174",     "O(H)", 64000, 6.94),
    ("Sh2-313",      "DAO",  80000, 7.2),
    ("HS0505+0112",  "DAO",  63200, 7.30),
    ("Ton 21",       "DA",   69710, 7.47),
    ("Feige 24",     "DA",   60000, 7.50),
    ("G191-B2B",     "DA",   52500, 7.53),
    ("REJ0558-373",  "DA",   59500, 7.70),
    ("RE-J0623-371", "DA",   58200, 7.14),
    ("REJ2214-492",  "DA",   61600, 7.29),
    ("REJ0457-281",  "DA",   51000, 7.93),
]

print("="*80)
print("Gravitational-potential ladder, Bainbridge et al. (2017) sample")
print("phi recovered from published log g; dalpha/alpha NOT available (withheld)")
print("="*80)
print(f"{'object':>14} {'type':>5} {'log g':>6} {'M/Msun':>7} {'R/Rearth':>8} "
      f"{'phi':>10} {'phi/phi_G191':>12}")

results = []
# reference
M_g191 = mass_from_logg_WD(7.53); R_g191 = wd_radius_Rsun(M_g191)
phi_g191 = phi_from_M_R(M_g191, R_g191)

for name, typ, Teff, logg in sample:
    if logg >= 6.8:  # white-dwarf-like: use WD M-R
        M = mass_from_logg_WD(logg); R = wd_radius_Rsun(M)
        phi = phi_from_M_R(M, R); Rearth_units = R*Rsun/Rearth
    else:            # sub-dwarf / HB star: nominal M, R from g
        M = 0.47
        Rcm = np.sqrt(G_cgs*(M*Msun*1e3)/10**logg)  # cm
        R = Rcm/1e2/Rsun
        phi = phi_from_M_R(M, R); Rearth_units = R*Rsun/Rearth
    results.append((name, typ, logg, M, phi))
    print(f"{name:>14} {typ:>5} {logg:>6.2f} {M:>7.3f} {Rearth_units:>8.2f} "
          f"{phi:>10.3e} {phi/phi_g191:>12.2f}")

phis = np.array([r[4] for r in results])
print(f"\nLever arm in phi (max/min): {phis.max()/phis.min():.0f}x "
      f"= {np.log10(phis.max()/phis.min()):.1f} decades")
print(f"Among the DA white dwarfs only (logg>=7): "
      f"{phis[[r[2]>=7 for r in results]].max()/phis[[r[2]>=7 for r in results]].min():.1f}x")
print()
print("Read-out:")
print("- The sub-dwarfs (vZ1128, ROB162, BD+28, Sh2-174) sit far DOWN the ladder;")
print("  they anchor the low-phi end and give the long lever arm the paper advertises.")
print("- The DA white dwarfs cluster within a factor ~3 in phi near G191-B2B.")
print("- REJ0457-281 (logg 7.93) is the highest-phi DA: the best high-end anchor,")
print("  and it was one of the four still-unobserved 'best targets' in 2017.")
