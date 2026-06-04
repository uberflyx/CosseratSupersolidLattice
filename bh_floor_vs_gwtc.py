#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bh_floor_vs_gwtc.py
===================

Confront the lattice primordial-pocket mass floor with the observed
gravitational-wave black-hole mass distribution, and show how the framework's
two formation channels (ordinary stellar collapse plus primordial pockets)
together account for the shape of that distribution.

The framework predicts a *survival floor* for primordial pockets: the smallest
region of uncrystallised vacuum whose own gravity holds its phase boundary
against the advancing crystallisation front. Below it the front engulfs the
pocket; at and above it the pocket locks in as a permanent black hole. The
floor is fixed by lattice constants alone, with no astrophysical input.

Part 1 derives the floor from first principles (CODATA constants and the
lattice relations m_0 = m_e/alpha, ell = r_e). It reproduces
M_min ~ 28 M_sun (fluid density 3 rho_crystal) up to ~48 M_sun (rho_crystal),
with Schwarzschild radius R_S(M_min) ~ 83 km.

Part 2 builds a two-panel figure:
  - Top: the most informative individual GWTC events, with 90% mass intervals.
  - Bottom: the framework's prediction as the sum of two components,
        * stellar collapse: the dominant ~10 M_sun population, suppressed
          above the pair-instability gap edge near ~45-50 M_sun;
        * primordial pockets: a pile-up at the survival floor (28-48 M_sun)
          plus a declining "void-tracking" tail that fills the gap and the
          high-mass range.
    Their sum reproduces the documented features of the GWTC-4.0 distribution:
    over-densities near 10, 20 and 35 M_sun and a continuum steepening above
    35 M_sun (Abac et al. 2025, arXiv:2508.18083).

Both component curves are schematic: they reproduce the *locations* of the
observed features and the framework's predicted shape (pile-up plus declining
tail). The exact pocket birth-mass spectrum is set by the nucleation statistics
of the crystallisation transition and is an open calculation; the curves are
labelled accordingly. The scientific comparison rests on the floor location,
which is parameter-free, and on the robust peak locations, not on curve detail.

Observational inputs (source-frame masses, 90% credible intervals):
  - GW150914: m1 = 35.6 (+4.7/-3.1) M_sun (Abbott et al. 2016).
  - GW170608: m1 = 11.0 (+5.5/-1.7) M_sun (Abbott et al. 2017).
  - GW190521: m1 = 85 (+21/-14),  m2 = 66 (+17/-18) M_sun (arXiv:2009.01075).
  - GW231123: m1 = 137 (+22/-17), m2 = 103 (+20/-52) M_sun (arXiv:2507.08219).
  - GW250114: m1 = 33.6 (+0.8/-0.8) M_sun (arXiv:2509.08099); the loudest event.
  - Mass-function features and the steepening above 35 M_sun: GWTC-4.0
    population analysis, Abac et al. (2025), arXiv:2508.18083.
  - Pair-instability gap ~50-130 M_sun (lower edge uncertain, ~45-65):
    Farmer et al. (2019), arXiv:1910.12874; Woosley & Heger (2021).
"""

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Physical and astronomical constants (CODATA 2022)
# ----------------------------------------------------------------------
c      = 299_792_458.0          # speed of light            [m/s]
G      = 6.674_30e-11           # Newton's constant         [m^3 kg^-1 s^-2]
m_e    = 9.109_383_7139e-31     # electron mass             [kg]
alpha  = 7.297_352_5643e-3      # fine-structure constant   [-]
r_e    = 2.817_940_3205e-15     # classical electron radius [m]
M_sun  = 1.988_409_8e30         # solar mass                [kg]
MeV    = 1.602_176_634e-13      # 1 MeV in joules           [J]

# ----------------------------------------------------------------------
# 1. The lattice pocket-mass floor, from first principles
# ----------------------------------------------------------------------
m0   = m_e / alpha                          # node mass m_0 = m_e/alpha  [kg]
ell  = r_e                                  # lattice spacing ell = r_e  [m]
a_fcc       = np.sqrt(2.0) * ell            # FCC conventional cube edge [m]
rho_crystal = 4.0 * m0 / a_fcc**3           # post-crystallisation density [kg/m^3]


def pocket_floor(rho_fluid):
    """Minimum surviving pocket mass and its Schwarzschild radius."""
    M_min = np.sqrt(3.0 * c**6 / (32.0 * np.pi * G**3 * rho_fluid))
    R_S   = 2.0 * G * M_min / c**2
    return M_min, R_S


M_lo, RS_lo = pocket_floor(3.0 * rho_crystal)   # dense fluid  -> lighter floor (~28)
M_hi, RS_hi = pocket_floor(1.0 * rho_crystal)   # relaxed      -> heavier floor (~48)
floor_lo = M_lo / M_sun
floor_hi = M_hi / M_sun

print("=" * 64)
print("Lattice primordial-pocket mass floor (parameter-free)")
print("=" * 64)
print(f"  node mass  m_0 = m_e/alpha = {m0*c**2/MeV:.3f} MeV")
print(f"  spacing    ell = r_e       = {ell:.4e} m")
print(f"  rho_crystal                = {rho_crystal:.3e} kg/m^3")
print(f"  floor (3 rho_crystal): {floor_lo:5.1f} M_sun,  R_S = {RS_lo/1e3:5.1f} km")
print(f"  floor (1 rho_crystal): {floor_hi:5.1f} M_sun,  R_S = {RS_hi/1e3:5.1f} km")
print(f"  => floor band: {floor_lo:.0f}-{floor_hi:.0f} M_sun")
print("=" * 64)

# ----------------------------------------------------------------------
# 2. Two-component prediction (schematic) on a primary-mass grid
# ----------------------------------------------------------------------
m = np.logspace(np.log10(5.0), np.log10(200.0), 2000)   # primary mass [M_sun]


def lognormal(x, mu, sig):
    return np.exp(-0.5 * ((np.log(x) - np.log(mu)) / sig) ** 2)


def gaussian(x, mu, sig):
    return np.exp(-0.5 * ((x - mu) / sig) ** 2)


def sigmoid(x, x0, w):
    return 0.5 * (1.0 + np.tanh((x - x0) / w))


# Channel A: stellar collapse. Dominant ~9.5 M_sun peak, a ~20 M_sun shoulder,
# a falling continuum, and a hard suppression above the PISN gap edge (~47).
stellar = (1.00 * lognormal(m, 9.5, 0.30)
           + 0.16 * gaussian(m, 20.0, 3.0)
           + 0.45 * (m / 9.5) ** (-3.2))
stellar *= sigmoid(m, 6.8, 0.9)             # detector turn-on at low mass
stellar *= (1.0 - sigmoid(m, 47.0, 3.5))    # PISN suppression above the gap edge

# Channel B: primordial pockets. A pile-up at the survival floor (centred ~35,
# inside 28-48) plus a declining void-tracking tail that fills the gap and the
# high-mass range. Both turn on at the floor.
pile_up = 0.085 * gaussian(m, 35.0, 6.0)
tail    = 0.045 * (m / 35.0) ** (-1.5) * sigmoid(m, 30.0, 3.0)
pocket  = (pile_up + tail) * sigmoid(m, 27.0, 2.0)

total = stellar + pocket
norm  = total.max()
stellar /= norm
pocket  /= norm
total   /= norm

# ----------------------------------------------------------------------
# 3. Notable individual events (source-frame primary mass, 90% interval)
# ----------------------------------------------------------------------
# name, m1, +err, -err, channel tag for colour
events = [
    ("GW170608",  11.0,  5.5,  1.7, "stellar"),
    ("GW250114",  33.6,  0.8,  0.8, "floor"),
    ("GW150914",  35.6,  4.7,  3.1, "floor"),
    ("GW190521",  85.0, 21.0, 14.0, "gap"),
    ("GW231123", 137.0, 22.0, 17.0, "gap"),
]
col_map = {"stellar": "#2a7f7f", "floor": "#c0392b", "gap": "#34495e"}

# ----------------------------------------------------------------------
# 4. Figure: two stacked panels sharing the mass axis
# ----------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 11,
    "font.family": "serif",
    "axes.linewidth": 0.8,
    "mathtext.fontset": "cm",
})

fig, (axE, axM) = plt.subplots(
    2, 1, figsize=(9.4, 6.8), sharex=True,
    gridspec_kw={"height_ratios": [1.0, 2.6], "hspace": 0.06},
)

# Shared shaded regions (drawn in both panels).
for ax in (axE, axM):
    ax.axvspan(50, 130, color="0.86", zorder=0)                           # PISN gap
    ax.axvspan(floor_lo, floor_hi, color="#c0392b", alpha=0.12, zorder=0)  # floor band

# ---- Top panel: observed events --------------------------------------
for i, (name, m1, ep, em, tag) in enumerate(events):
    y = i + 1
    axE.errorbar(m1, y, xerr=[[em], [ep]], fmt="o", color=col_map[tag],
                 ms=7, capsize=3, lw=1.4, zorder=5)
    axE.annotate(name, xy=(m1, y), xytext=(m1, y + 0.30), ha="center",
                 va="bottom", fontsize=8.4, color=col_map[tag])
axE.set_ylim(0.4, len(events) + 1.3)
axE.set_yticks([])
axE.set_ylabel("observed\nevents", fontsize=9)
axE.set_title("Lattice primordial-pocket floor and the gravitational-wave "
              "black-hole mass function", fontsize=12, pad=8)
axE.text(np.sqrt(floor_lo * floor_hi), len(events) + 0.95, "floor",
         ha="center", va="center", fontsize=8.6, color="#c0392b", weight="bold")
axE.text(np.sqrt(50 * 130), len(events) + 0.95, "PISN gap",
         ha="center", va="center", fontsize=8.6, color="0.4", weight="bold")

# ---- Bottom panel: framework prediction ------------------------------
axM.plot(m, stellar, color="#2a7f7f", lw=1.8, ls="--", zorder=4,
         label="stellar collapse (channel A)")
axM.plot(m, pocket, color="#c0392b", lw=2.0, ls="-", zorder=4,
         label="primordial pockets (channel B)")
axM.plot(m, total, color="#1b2a4a", lw=2.6, ls="-", zorder=5,
         label="total prediction")
axM.fill_between(m, 1e-3, pocket, color="#c0392b", alpha=0.10, zorder=2)

# Floor-edge guides.
axM.axvline(floor_lo, color="#c0392b", lw=1.2, ls="--", zorder=3)
axM.axvline(floor_hi, color="#c0392b", lw=1.0, ls=":", zorder=3)

# Annotate the observed 35 M_sun over-density on the total curve.
i35 = np.argmin(np.abs(m - 35.0))
axM.plot(35.0, total[i35], "v", color="#1b2a4a", ms=9, zorder=6)
axM.annotate(r"observed $35\,M_\odot$ peak", xy=(35.0, total[i35]),
             xytext=(60, 0.34), fontsize=9, color="#1b2a4a", ha="left",
             arrowprops=dict(arrowstyle="->", color="#1b2a4a", lw=1.0))

# Story annotations.
axM.annotate("stellar peak\n(channel A)", xy=(9.5, 0.95), xytext=(5.6, 0.30),
             fontsize=8.4, color="#2a7f7f", ha="left",
             arrowprops=dict(arrowstyle="->", color="#2a7f7f", lw=0.9))
axM.annotate(r"survival floor $M_{\min}=28$–$48\,M_\odot$" "\n(parameter-free)",
             xy=(floor_lo, 0.52), xytext=(15.0, 0.52), fontsize=8.6,
             color="#c0392b", va="center",
             arrowprops=dict(arrowstyle="->", color="#c0392b", lw=0.9))
axM.annotate("void-tracking tail\nfills the gap",
             xy=(100, total[np.argmin(np.abs(m - 100))]),
             xytext=(58, 0.028), fontsize=8.4, color="#c0392b", ha="left",
             arrowprops=dict(arrowstyle="->", color="#c0392b", lw=0.9))

axM.set_xscale("log")
axM.set_yscale("log")
axM.set_xlim(5, 200)
axM.set_ylim(4e-3, 1.7)
axM.set_xlabel(r"primary black-hole mass  $m_1\ [M_\odot]$")
axM.set_ylabel(r"relative merger-rate density  $dR/dm_1$  (arb.)")
axM.set_xticks([5, 10, 20, 35, 50, 100, 200])
axM.set_xticklabels(["5", "10", "20", "35", "50", "100", "200"])
axM.legend(loc="upper right", fontsize=8.4, framealpha=0.96)

fig.savefig("bh_floor_vs_gwtc.pdf", bbox_inches="tight")
fig.savefig("bh_floor_vs_gwtc.png", dpi=175, bbox_inches="tight")
print("\nFigure written: bh_floor_vs_gwtc.pdf / .png")
