#!/usr/bin/env python3
"""
gw_spectrum_crystallisation.py — GW spectrum from the vacuum phase transition
=============================================================================
M. A. Cox, University of the Witwatersrand (2026)

The vacuum crystallisation (fluid → FCC Cosserat supersolid) is a first-order
phase transition whose four thermodynamic parameters are all determined by the
lattice mechanics:

    α   = 0.94   transition strength (bag model, g_QGP - g_had)
    β/H = 7.5    inverse duration (CNT barrier / T_c)
    v_w = c/√3   wall velocity (fluid sound speed)
    T*  = 156 MeV nucleation temperature (D4 Polyakov loop)

The dominant GW source is the sound-wave mechanism (Hindmarsh et al. 2014).
Bubble collisions are negligible (κ_wall ≈ 2e-19; the crystal is incompressible
and absorbs the wall kinetic energy). Turbulence is subdominant.

The PEAK FREQUENCY (~354 nHz) is robust — set by nucleation geometry alone.
The PEAK AMPLITUDE (h²Ω ≤ 2.6e-10) is an UPPER BOUND — the standard sound-wave
efficiency κ_v was derived for a classical fluid pushed by expanding bubbles,
but the lattice crystallises by material addition at the front, which may
generate less anisotropic stress.

Usage:
    python gw_spectrum_crystallisation.py
"""

import numpy as np

# ═══════════════════════════════════════════
# PHYSICAL CONSTANTS
# ═══════════════════════════════════════════

hbar_c = 197.3269804       # MeV·fm
M_Pl   = 1.22089e19        # GeV (Planck mass)
alpha_em = 1/137.035999177
m_e    = 0.51099895        # MeV
m0     = m_e / alpha_em    # node mass
ell    = hbar_c / m0       # lattice spacing (fm)

# ═══════════════════════════════════════════
# THERMODYNAMIC PARAMETERS
# ═══════════════════════════════════════════

# Effective degrees of freedom
g_QGP   = 44.7    # QGP: 37 (QCD) + 7.75 (leptons/photons in QGP-adjacent range)
g_had   =  3.0    # hadronic (pions only)
g_nonQCD = 14.2   # non-QCD (photons, leptons, neutrinos)
g_star  = g_QGP + g_nonQCD

# Transition strength from the bag model:
# L = (π²/90)(g_QGP - g_had)T_c⁴ = latent heat
# α = L / ρ_rad = 4(g_QGP - g_had) / (3 g_star)
alpha_PT = 4 * (g_QGP - g_had) / (3 * g_star)

# Nucleation temperature and barrier
T_c = 156.1  # MeV (D4 Polyakov loop)
DeltaG = (16 * np.pi / 3) * m0  # CNT barrier ≈ 1173 MeV
beta_H = DeltaG / T_c

# Wall velocity (fluid sound speed — detonation mode)
v_w = 1.0 / np.sqrt(3)

# Sound-wave efficiency (Espinosa et al. 2010)
kappa_v = alpha_PT / (0.73 + 0.083*np.sqrt(alpha_PT) + alpha_PT)

# Lifetime suppression (Guo et al. 2021)
Upsilon = 1 - 1 / (1 + 1/beta_H)

# Bubble collision efficiency (negligible — crystal is incompressible)
kappa_wall = 2e-19  # from monograph Eq. kappa_wall

# Bag constant
B_fm3 = (np.pi**2 / 90) * (g_QGP - g_had) * T_c**4 / hbar_c**3
B14 = (B_fm3 * hbar_c**3)**0.25


# ═══════════════════════════════════════════
# GW SPECTRUM: SOUND-WAVE MECHANISM
# ═══════════════════════════════════════════

def gw_spectrum(f_nHz, T_MeV, alpha, g, beta_over_H):
    """Compute h²Ω_GW(f) from the sound-wave mechanism."""
    T_GeV = T_MeV / 1000
    kv = alpha / (0.73 + 0.083*np.sqrt(alpha) + alpha)
    Ups = 1 - 1/(1 + 1/beta_over_H)

    # Peak frequency (Hindmarsh et al. 2017, Eq. 29)
    f_pk = 1.9e-5 * (1/v_w) * beta_over_H * (T_GeV/100) * (g/100)**0.1667
    f_pk_nHz = f_pk * 1e9

    # Peak amplitude
    K = kv * alpha / (1 + alpha)
    h2_pk = 2.65e-6 * beta_over_H**(-2) * K**2 * (100/g)**(1/3) * v_w * Ups

    # Spectral shape: C(s) = s³(7/(4+3s²))^{7/2}
    s = f_nHz / f_pk_nHz
    h2_Omega = h2_pk * s**3 * (7 / (4 + 3*s**2))**3.5

    return h2_Omega, f_pk_nHz, h2_pk


# ═══════════════════════════════════════════
# OUTPUT
# ═══════════════════════════════════════════

print("GW spectrum from vacuum lattice crystallisation")
print("=" * 60)
print(f"  m₀ = {m0:.2f} MeV,  ℓ = {ell:.4f} fm,  T_c = {T_c} MeV")
print(f"  ΔG* = {DeltaG:.0f} MeV,  β/H* = {beta_H:.2f}")
print(f"  α = {alpha_PT:.3f} (bag model),  v_w = c/√3")
print(f"  κ_v = {kappa_v:.4f},  Υ = {Upsilon:.4f}")
print(f"  κ_wall = {kappa_wall:.0e} (bubble collisions negligible)")

_, f_pk, h2_pk = gw_spectrum(1, T_c, alpha_PT, g_star, beta_H)
print(f"\n  PREDICTION (upper bound):")
print(f"  ─────────────────────────")
print(f"  f_peak = {f_pk:.0f} nHz")
print(f"  h²Ω_peak ≤ {h2_pk:.1e}")
print(f"  At 100 nHz: h²Ω ≤ {gw_spectrum(100, T_c, alpha_PT, g_star, beta_H)[0]:.1e}")
print(f"  At  30 nHz: h²Ω ≤ {gw_spectrum(30, T_c, alpha_PT, g_star, beta_H)[0]:.1e}")

print(f"\n  Bag model self-consistency:")
print(f"  B^(1/4) = {B14:.1f} MeV  vs  Λ_QCD = {np.pi*m0:.1f} MeV  ({abs(B14-np.pi*m0)/(np.pi*m0)*100:.1f}%)")

# Parametric comparison
print(f"\n  Parametric comparison:")
print(f"  {'α':>6}  {'f_pk (nHz)':>11}  {'h²Ω_peak':>11}  {'note':>20}")
print(f"  {'─'*6}  {'─'*11}  {'─'*11}  {'─'*20}")
for a, note in [(0.03, "shear only"),
                (0.10, "intermediate"),
                (0.30, "strong"),
                (alpha_PT, "bag model (derived)")]:
    bH = DeltaG / T_c
    _, fp, hp = gw_spectrum(1, T_c, a, g_star, bH)
    marker = "  ◄ UPPER BOUND" if a == alpha_PT else ""
    print(f"  {a:6.3f}  {fp:11.0f}  {hp:11.2e}  {note:>20}{marker}")


# ═══════════════════════════════════════════
# PLOT
# ═══════════════════════════════════════════

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    f_nHz = np.logspace(-0.5, 3.5, 1000)
    fig, ax = plt.subplots(figsize=(10, 7))

    # Derived spectrum (upper bound)
    h2, _, _ = gw_spectrum(f_nHz, T_c, alpha_PT, g_star, beta_H)
    ax.loglog(f_nHz, h2, 'r-', lw=2.5,
              label=rf'Derived ($\alpha={alpha_PT:.2f}$, upper bound)')

    # Parametric comparisons
    for a, col, ls, lab in [(0.03, '#2980b9', '--', r'$\alpha=0.03$'),
                             (0.30, '#27ae60', '--', r'$\alpha=0.30$')]:
        h2p, _, _ = gw_spectrum(f_nHz, T_c, a, g_star, beta_H)
        ax.loglog(f_nHz, h2p, color=col, ls=ls, lw=1.5, label=lab)

    # NANOGrav 15yr approximate free spectrum
    ng_f = np.arange(1, 15) / (15 * 365.25 * 24 * 3600) * 1e9  # nHz
    ng_h2 = np.array([2e-9, 5.5e-10, 2.2e-10, 1.2e-10, 7.5e-11, 5e-11, 4e-11,
                      2.8e-11, 2.2e-11, 1.8e-11, 1.5e-11, 1.3e-11, 1.1e-11, 9e-12])
    ax.fill_between(ng_f, ng_h2/3.5, ng_h2*3.5, alpha=0.15, color='gray')
    ax.scatter(ng_f, ng_h2, color='black', s=20, zorder=5, label='NANOGrav 15yr')

    # Sensitivity bands
    ax.axvspan(1, 100, alpha=0.03, color='blue')
    ax.text(8, 5e-8, 'PTA band', fontsize=9, color='blue', alpha=0.5)
    ax.axvspan(100, 3000, alpha=0.03, color='green')
    ax.text(500, 5e-8, 'SKA projected', fontsize=9, color='green', alpha=0.5)

    ax.set_xlim(0.3, 3000)
    ax.set_ylim(1e-16, 1e-7)
    ax.set_xlabel(r'Frequency $f$ [nHz]', fontsize=13)
    ax.set_ylabel(r'$h^2 \Omega_{\mathrm{GW}}(f)$', fontsize=13)
    ax.set_title('GW spectrum from vacuum crystallisation\n'
                 r'(peak amplitude is an upper bound; $\kappa_{\rm wall} \approx 0$)',
                 fontsize=12)
    ax.legend(fontsize=10, loc='upper right', framealpha=0.9)
    ax.grid(True, which='both', alpha=0.15)
    plt.tight_layout()
    plt.savefig('/mnt/user-data/outputs/gw_spectrum_crystallisation.png', dpi=150)
    print(f"\n  Plot saved to outputs/")
except Exception as e:
    print(f"\n  (Plot skipped: {e})")
