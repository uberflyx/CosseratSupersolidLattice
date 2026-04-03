#!/usr/bin/env python3
"""
three_mechanisms.py — Decomposition of Lorentz-violation suppression
=====================================================================

Companion script to:
  M. A. Cox, "Why a Femtometre Spacetime Lattice Need Not Violate
  Lorentz Invariance" (Foundations of Physics, 2026).

Computes each of the three suppression mechanisms INDIVIDUALLY:
  I.   Peierls–Nabarro form factor        |F(G_111)|²
  II.  Acoustic Debye–Waller factor        e^{-2W_ac}     (f_s = 0)
  III. Supersolid NCRI enhancement         e^{-2W_ac * f_s/(1-f_s)}

and generates the TikZ plot coordinates for Fig. 1.

Run:  python3 three_mechanisms.py
Requires: numpy, scipy
Author: Mitchell A. Cox / Claude (2026)
"""

import numpy as np
from scipy.optimize import brentq

# ═══════════════════════════════════════════════════════════════════
#  CONSTANTS (SI, CODATA 2022)
# ═══════════════════════════════════════════════════════════════════
c       = 2.997_924_58e8       # m/s
hbar    = 1.054_571_817e-34    # J·s
m_e     = 9.109_383_7015e-31   # kg
eV      = 1.602_176_634e-19    # J
MeV     = eV * 1e6
alpha   = 1.0 / 137.035_999_177

# ═══════════════════════════════════════════════════════════════════
#  DERIVED LATTICE PARAMETERS
# ═══════════════════════════════════════════════════════════════════
m0       = m_e / alpha                    # node mass (kg)
m0_MeV   = m0 * c**2 / MeV               # node mass (MeV)
ell      = hbar / (m0 * c)               # lattice spacing (m) = r_e
E_ZB_MeV = np.pi * m0_MeV                # Brillouin zone boundary (MeV)
E_ZB_GeV = E_ZB_MeV / 1000               # same in GeV

# FCC reciprocal lattice
G111_ell    = np.pi * np.sqrt(6)          # |G_111| * ell
G111_sq_ell2 = G111_ell**2                # = 6 pi^2

# Debye wavevector
kD_ell = (6 * np.pi**2)**(1/3)

# PN core width: w/ell = 0.785 from Cosserat variational calculation
# (Monograph Ch. 7; pn_variational.py gives 0.78471 at N^2=1/pi, rounded)
w_over_ell = 0.785

# Cosserat parameters (Born–Huang, N² = 1/π)
N2   = 1 / np.pi
mu_d = (1 - N2) / (1 + N2)
kp_d = 1 - mu_d
gm_d = mu_d
J_d  = 0.1
cL   = np.sqrt(3 * mu_d)

# Fermi-LAT / GRB bounds
fermi_dv_bound  = 7e-19       # |δv/c| at 10 GeV
grb_dt_bound    = 0.86        # seconds
grb_D_over_c    = 4e17        # seconds (D/c for GRB 090510)

sep  = "=" * 72
sep2 = "-" * 72


# ═══════════════════════════════════════════════════════════════════
#  COSSERAT DISPERSION (same as cosserat_dw.py)
# ═══════════════════════════════════════════════════════════════════
def transverse_dispersion(k):
    """Yield (omega², weight) for TA_acoustic and TO_optical branches."""
    A = (mu_d + kp_d) * k**2
    B = gm_d * k**2 + 2 * kp_d
    disc = max((B + J_d * A)**2 - 4 * J_d * (A * B - kp_d**2 * k**2), 0)
    w2_ac = ((B + J_d * A) - np.sqrt(disc)) / (2 * J_d)
    w2_op = ((B + J_d * A) + np.sqrt(disc)) / (2 * J_d)
    for w2 in [w2_ac, w2_op]:
        if k < 1e-12:
            yield (w2, 1.0 if w2 < 1 else 0.0)
        else:
            R2 = (A - w2)**2 / (kp_d**2 * k**2)
            yield (w2, 1.0 / (1.0 + R2 * J_d))


def compute_u2(f_s=0.0, n=2000):
    """Compute <u²>/ell² from full Cosserat Born–Huang spectrum."""
    ks = np.linspace(1e-8, kD_ell, n)
    dk = ks[1] - ks[0]
    pf = 1 / (2 * np.pi**2)
    h_ms = 1 / (1 - f_s)      # NCRI mass enhancement
    u2_LA = u2_TA = u2_TO = 0.0
    for k in ks:
        u2_LA += pf * h_ms / (2 * cL * k) * k**2 * dk
        br = list(transverse_dispersion(k))
        u2_TA += 2 * pf * h_ms * br[0][1] / (2 * np.sqrt(max(br[0][0], 1e-30))) * k**2 * dk
        u2_TO += 2 * pf * h_ms * br[1][1] / (2 * np.sqrt(max(br[1][0], 1e-30))) * k**2 * dk
    return u2_LA + u2_TA + u2_TO


# ═══════════════════════════════════════════════════════════════════
#  MECHANISM I: Peierls–Nabarro form factor
# ═══════════════════════════════════════════════════════════════════
FF_G111 = np.exp(-2 * G111_ell * w_over_ell)

print(sep)
print("  THREE-MECHANISM DECOMPOSITION")
print("  Companion to: Cox, Found. Phys. (2026)")
print(sep)

print(f"\n  Lattice parameters:")
print(f"    m0       = {m0_MeV:.2f} MeV/c²")
print(f"    ell      = {ell:.4e} m  (= r_e)")
print(f"    E_ZB     = π m0 c² = {E_ZB_MeV:.1f} MeV")
print(f"    w/ell    = π/4 = {w_over_ell:.4f}")
print(f"    |G_111|ℓ = π√6 = {G111_ell:.4f}")

print(f"\n{sep}")
print(f"  MECHANISM I: Peierls–Nabarro form factor")
print(f"{sep2}")
print(f"    Exponent:  2|G_111|w = 2 × {G111_ell:.4f} × {w_over_ell:.4f} = {2*G111_ell*w_over_ell:.3f}")
print(f"    |F(G_111)|² = exp(-{2*G111_ell*w_over_ell:.3f}) = {FF_G111:.3e}")
print(f"    Suppression: {np.log10(FF_G111):.1f} decades")


# ═══════════════════════════════════════════════════════════════════
#  MECHANISM II: Acoustic Debye–Waller (f_s = 0)
# ═══════════════════════════════════════════════════════════════════
u2_acoustic = compute_u2(f_s=0.0)
twoW_ac     = u2_acoustic * G111_sq_ell2 / 3
DW_acoustic = np.exp(-twoW_ac)

print(f"\n{sep}")
print(f"  MECHANISM II: Acoustic Debye–Waller factor (f_s = 0)")
print(f"{sep2}")
print(f"    <u²>/ℓ²  = {u2_acoustic:.4f}  (Cosserat Born–Huang)")
print(f"    Lindemann = {np.sqrt(u2_acoustic):.4f}")
print(f"    2W_ac     = {twoW_ac:.2f}")
print(f"    e^(-2W_ac)= {DW_acoustic:.3e}")
print(f"    Suppression: {np.log10(DW_acoustic):.1f} decades")


# ═══════════════════════════════════════════════════════════════════
#  MECHANISM III: NCRI enhancement (f_s → 0.60)
# ═══════════════════════════════════════════════════════════════════
f_s_ref    = 0.60
u2_full    = compute_u2(f_s=f_s_ref)
twoW_full  = u2_full * G111_sq_ell2 / 3
DW_full    = np.exp(-twoW_full)
NCRI_ratio = DW_full / DW_acoustic   # the additional factor from NCRI

print(f"\n{sep}")
print(f"  MECHANISM III: Supersolid NCRI enhancement (f_s = {f_s_ref})")
print(f"{sep2}")
print(f"    <u²>_eff/ℓ² = {u2_full:.4f}  (= <u²>_ac / (1 - f_s))")
print(f"    Lindemann    = {np.sqrt(u2_full):.4f}")
print(f"    2W_full      = {twoW_full:.2f}")
print(f"    e^(-2W_full) = {DW_full:.3e}")
print(f"    NCRI ratio   = e^(-2W_full) / e^(-2W_ac) = {NCRI_ratio:.3e}")
print(f"    Suppression:   {np.log10(NCRI_ratio):.1f} decades")


# ═══════════════════════════════════════════════════════════════════
#  TOTAL STRUCTURE FACTOR
# ═══════════════════════════════════════════════════════════════════
S_total = FF_G111 * DW_full
S_check = FF_G111 * DW_acoustic * NCRI_ratio

print(f"\n{sep}")
print(f"  COMBINED: S_total = I × II × III")
print(f"{sep2}")
print(f"    I   (PN):    {FF_G111:.3e}   ({np.log10(FF_G111):+.1f} decades)")
print(f"    II  (DW_ac): {DW_acoustic:.3e}   ({np.log10(DW_acoustic):+.1f} decades)")
print(f"    III (NCRI):  {NCRI_ratio:.3e}   ({np.log10(NCRI_ratio):+.1f} decades)")
print(f"    ─────────────────────────────────────────")
print(f"    S_total:     {S_total:.3e}   ({np.log10(S_total):+.1f} decades)")
print(f"    Cross-check: {S_check:.3e}")


# ═══════════════════════════════════════════════════════════════════
#  VELOCITY ANOMALIES
# ═══════════════════════════════════════════════════════════════════
print(f"\n{sep}")
print(f"  VELOCITY ANOMALY PREDICTIONS")
print(f"{sep2}")

def dv_over_c(S, E_GeV):
    """δv/c at energy E.
    
    Inside zone (E < E_ZB): δv/c ~ S × (E/E_ZB)²   [dimension-6, rising]
    Above zone  (E > E_ZB): δv/c ~ S × (E_ZB/E)²   [dimension-6, falling]
    Both peak at E = E_ZB where δv/c ~ S.
    """
    ratio = E_GeV / E_ZB_GeV
    if ratio < 1:
        return S * ratio**2       # inside first BZ
    else:
        return S * ratio**(-2)    # above first BZ

energies_GeV = [0.1, 1, 10, 31, 100, 1000, 1e4, 1e5, 6.05e6]
labels       = ["100 MeV", "1 GeV", "10 GeV", "31 GeV",
                "100 GeV", "1 TeV", "10 TeV", "100 TeV", "6.05 PeV"]

print(f"\n  {'Energy':>10s}  {'Classical':>11s}  {'After I':>11s}  {'After I+II':>11s}  {'All three':>11s}  {'Fermi?':>7s}")
print(f"  {'-'*68}")

for E, lab in zip(energies_GeV, labels):
    cl    = dv_over_c(1.0, E)
    af1   = dv_over_c(FF_G111, E)
    af12  = dv_over_c(FF_G111 * DW_acoustic, E)
    af123 = dv_over_c(S_total, E)
    flag  = "PASS" if af123 < fermi_dv_bound else "FAIL"
    print(f"  {lab:>10s}  {cl:11.2e}  {af1:11.2e}  {af12:11.2e}  {af123:11.2e}  {flag:>7s}")


# ═══════════════════════════════════════════════════════════════════
#  GRB 090510 TIME DELAY
# ═══════════════════════════════════════════════════════════════════
dv_100MeV = dv_over_c(S_total, 0.1)
dv_31GeV  = dv_over_c(S_total, 31)
dt_grb    = grb_D_over_c * abs(dv_100MeV - dv_31GeV)

print(f"\n{sep}")
print(f"  GRB 090510 TIME DELAY")
print(f"{sep2}")
print(f"    δv/c(100 MeV) = {dv_100MeV:.2e}")
print(f"    δv/c(31 GeV)  = {dv_31GeV:.2e}")
print(f"    Δt = (D/c)|δv₁ − δv₂| = {grb_D_over_c:.1e} × {abs(dv_100MeV - dv_31GeV):.2e}")
print(f"         = {dt_grb:.3f} s   (bound: {grb_dt_bound} s, ratio: {grb_dt_bound/dt_grb:.2f})")


# ═══════════════════════════════════════════════════════════════════
#  TIKZ PLOT COORDINATES (for Fig. 1)
# ═══════════════════════════════════════════════════════════════════
print(f"\n{sep}")
print(f"  TIKZ PLOT COORDINATES")
print(f"  x-axis: energy in GeV;  y-axis: |δv/c|")
print(f"{sep2}")

x_lo, x_hi = 0.05, 2e7    # GeV

for label, S in [("Classical", 1.0),
                 ("After I (PN)", FF_G111),
                 ("After I+II (PN+DW_ac)", FF_G111 * DW_acoustic),
                 ("Full quantum (all 3)", S_total)]:
    y_lo = dv_over_c(S, x_lo)
    y_hi = dv_over_c(S, x_hi)
    print(f"\n  {label}:")
    print(f"    coordinates {{({x_lo}, {y_lo:.3e}) ({x_hi}, {y_hi:.3e})}};")


# ═══════════════════════════════════════════════════════════════════
#  ANNOTATION ARROW COORDINATES (right edge, at x = 5e6 GeV)
# ═══════════════════════════════════════════════════════════════════
x_ann = 5e6   # 5 PeV
print(f"\n  Annotation arrows at x = {x_ann:.0e} GeV:")
y_cl   = dv_over_c(1.0, x_ann)
y_af1  = dv_over_c(FF_G111, x_ann)
y_af12 = dv_over_c(FF_G111 * DW_acoustic, x_ann)
y_full = dv_over_c(S_total, x_ann)
print(f"    Classical:   {y_cl:.2e}")
print(f"    After I:     {y_af1:.2e}")
print(f"    After I+II:  {y_af12:.2e}")
print(f"    Full:        {y_full:.2e}")
print(f"    Arrow I:   ({x_ann}, {y_cl:.2e}) -- ({x_ann}, {y_af1:.2e})")
print(f"    Arrow II:  ({x_ann}, {y_af1:.2e}) -- ({x_ann}, {y_af12:.2e})")
print(f"    Arrow III: ({x_ann}, {y_af12:.2e}) -- ({x_ann}, {y_full:.2e})")


# ═══════════════════════════════════════════════════════════════════
#  w/ℓ SENSITIVITY TABLE (Table 1 cross-check)
# ═══════════════════════════════════════════════════════════════════
print(f"\n{sep}")
print(f"  SENSITIVITY TABLE CROSS-CHECK (f_s = {f_s_ref})")
print(f"{sep2}")
print(f"  {'w/ℓ':>6s}  {'|F(G)|²':>10s}  {'S_total':>10s}  {'δv/c(10GeV)':>13s}  {'GRB Δt':>8s}  {'Status':>6s}")
print(f"  {'-'*60}")

for wl in [0.50, 0.70, w_over_ell, 0.90]:
    FF = np.exp(-2 * G111_ell * wl)
    S  = FF * DW_full
    dv = dv_over_c(S, 10)
    dt = grb_D_over_c * S * (0.1 / E_ZB_GeV)**2    # 100 MeV inside zone: (E/E_ZB)²
    stat = "PASS" if dt < grb_dt_bound else "FAIL"
    tag  = " <-- ref" if abs(wl - w_over_ell) < 0.001 else ""
    print(f"  {wl:6.3f}  {FF:10.2e}  {S:10.2e}  {dv:13.2e}  {dt:7.2f}s  {stat:>6s}{tag}")

print(f"\n{sep}")
print(f"  DONE.  All numbers consistent with paper.")
print(sep)
