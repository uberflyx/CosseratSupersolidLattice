#!/usr/bin/env python3
"""
gw_spectrum_crystallisation.py

Gravitational wave spectrum from the vacuum lattice crystallisation
phase transition. Computes the stochastic GW background from sound
waves, bubble collisions, and MHD turbulence using the standard
fitting formulas (Caprini et al. 2016, 2020; Hindmarsh et al. 2017;
Ellis & Lewicki 2020).

Lattice inputs:
    ell     = 2.8179 fm     (lattice spacing = classical electron radius)
    m0      = 70.025 MeV    (node mass = m_e / alpha)
    mu      = m0 c^2 / ell^3 (shear modulus)
    Theta_D = pi * m0       (Debye temperature ~ 220 MeV)
    gamma_USF = pi * m0 / ell^2 (unstable stacking fault energy)

Derived phase-transition parameters:
    T*      = nucleation temperature (~ 150-220 MeV)
    alpha   = latent heat / radiation energy density
    beta/H* = inverse duration in Hubble units
    v_w     = bubble wall velocity (Jouguet: c/sqrt(3))

Repository: https://github.com/uberflyx/CosseratSupersolidLattice
Author: Mitchell A. Cox
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

# ============================================================
# Physical constants (natural units: hbar = c = k_B = 1)
# ============================================================
hbar_c = 197.3269804  # MeV fm
c_light = 2.99792458e8  # m/s
G_Newton = 6.67430e-11  # m^3 kg^-1 s^-2
M_Pl = 1.22089e19  # GeV (Planck mass)
h_hubble = 0.674  # Planck 2018
H0_per_sec = h_hubble * 100 * 1e3 / (3.0857e22)  # H0 in s^-1

# ============================================================
# Lattice parameters
# ============================================================
alpha_em = 1.0 / 137.035999084
m_e_MeV = 0.51099895  # electron mass in MeV
m0_MeV = m_e_MeV / alpha_em  # node mass ~ 70.025 MeV
ell_fm = hbar_c / m0_MeV  # lattice spacing ~ 2.818 fm
mu_MeV_fm3 = m0_MeV / ell_fm**3  # shear modulus in MeV/fm^3
Theta_D_MeV = np.pi * m0_MeV  # Debye temperature ~ 220 MeV
gamma_USF_MeV_fm2 = np.pi * m0_MeV / ell_fm**2  # USF energy

print("=" * 60)
print("LATTICE PARAMETERS")
print("=" * 60)
print(f"  alpha_em       = 1/{1/alpha_em:.6f}")
print(f"  m_e            = {m_e_MeV:.5f} MeV")
print(f"  m_0            = {m0_MeV:.3f} MeV")
print(f"  ell            = {ell_fm:.4f} fm")
print(f"  mu             = {mu_MeV_fm3:.4f} MeV/fm^3")
print(f"  Theta_D        = {Theta_D_MeV:.2f} MeV")
print(f"  gamma_USF      = {gamma_USF_MeV_fm2:.4f} MeV/fm^2")
print()

# ============================================================
# Phase transition parameters from lattice mechanics
# ============================================================

def radiation_energy_density(T_MeV, g_star):
    """Radiation energy density in MeV/fm^3."""
    return (np.pi**2 / 30.0) * g_star * T_MeV**4 / hbar_c**3


def compute_PT_parameters(T_star_MeV, alpha_PT, g_star=62.0):
    """
    Compute the phase transition parameters from lattice inputs.

    Parameters
    ----------
    T_star_MeV : float
        Nucleation temperature in MeV.
    alpha_PT : float
        Phase transition strength (latent heat / radiation).
    g_star : float
        Effective relativistic DOF at T*.

    Returns
    -------
    dict with all GW-relevant parameters.
    """
    # CNT barrier from monograph Eq. 13.11
    DeltaG_star = (16 * np.pi / 3.0) * m0_MeV  # MeV (= 1170 MeV)

    # Inverse duration: beta/H* = DeltaG* / T*
    beta_over_H = DeltaG_star / T_star_MeV

    # Bubble wall velocity: Jouguet detonation
    v_w = 1.0 / np.sqrt(3.0)

    # Hubble rate at T* (radiation domination)
    # H* = sqrt(8 pi^3 g* / 90) * T*^2 / M_Pl
    # in natural units (GeV): M_Pl = 1.22e19 GeV
    T_star_GeV = T_star_MeV / 1000.0
    H_star = np.sqrt(8 * np.pi**3 * g_star / 90.0) * T_star_GeV**2 / M_Pl  # GeV

    # Mean bubble separation
    R_star = v_w / (beta_over_H * H_star)  # in 1/GeV

    # Efficiency factor kappa_v for Jouguet detonation
    # From Espinosa et al. (2010), Eq. 47:
    # For Jouguet, v_w = c_s = 1/sqrt(3)
    # kappa_v ~ alpha / (0.73 + 0.083*sqrt(alpha) + alpha)
    kappa_v = alpha_PT / (0.73 + 0.083 * np.sqrt(alpha_PT) + alpha_PT)

    # Sound wave lifetime suppression factor
    # Upsilon = 1 - 1/(1 + H* tau_sw)
    # tau_sw ~ R* / c_s, where c_s = 1/sqrt(3)
    c_s = 1.0 / np.sqrt(3.0)
    H_tau_sw = (H_star * R_star) / c_s  # dimensionless
    # Simpler: H*tau_sw = v_w / (beta/H* * c_s) = 1/(beta/H*)
    H_tau_sw_simple = 1.0 / beta_over_H
    Upsilon = 1.0 - 1.0 / (1.0 + H_tau_sw_simple)

    rho_rad = radiation_energy_density(T_star_MeV, g_star)

    return {
        'T_star_MeV': T_star_MeV,
        'T_star_GeV': T_star_GeV,
        'alpha': alpha_PT,
        'beta_over_H': beta_over_H,
        'v_w': v_w,
        'g_star': g_star,
        'H_star_GeV': H_star,
        'kappa_v': kappa_v,
        'Upsilon': Upsilon,
        'DeltaG_star': DeltaG_star,
        'rho_rad': rho_rad,
    }


# ============================================================
# GW spectrum: sound wave contribution
# (Hindmarsh et al. 2017; Caprini et al. 2016, 2020)
# ============================================================

def Omega_GW_sound_waves(f_Hz, params):
    """
    GW energy density spectrum from sound waves.

    Uses the broken power law fit from Hindmarsh et al. (2017)
    and the suppression factor from Ellis & Lewicki (2020).

    Parameters
    ----------
    f_Hz : array
        Frequencies in Hz.
    params : dict
        Phase transition parameters from compute_PT_parameters.

    Returns
    -------
    h2_Omega : array
        h^2 Omega_GW(f).
    """
    alpha = params['alpha']
    beta_H = params['beta_over_H']
    v_w = params['v_w']
    T_GeV = params['T_star_GeV']
    g_star = params['g_star']
    kappa_v = params['kappa_v']
    Ups = params['Upsilon']

    # Peak frequency today (Hz) - Caprini et al. 2016 Eq. 15
    f_peak = 1.9e-5 * (1.0 / v_w) * beta_H * (T_GeV / 100.0) * (g_star / 100)**0.1667

    # Peak amplitude - Caprini et al. 2016 Eq. 2.12, with Upsilon
    K = kappa_v * alpha / (1.0 + alpha)  # kinetic energy fraction
    h2_Omega_peak = (2.65e-6 * (1.0 / beta_H)**2 * K**2
                     * (100.0 / g_star)**(1.0/3.0) * v_w * Ups)

    # Spectral shape: broken power law (Caprini et al. 2016 Eq. 2.13)
    s = f_Hz / f_peak
    # C(s) = s^3 * (7 / (4 + 3 s^2))^(7/2)
    spectral_shape = s**3 * (7.0 / (4.0 + 3.0 * s**2))**3.5

    h2_Omega = h2_Omega_peak * spectral_shape

    return h2_Omega, f_peak, h2_Omega_peak


def Omega_GW_bubble_collisions(f_Hz, params):
    """
    GW energy density spectrum from bubble collisions.

    Uses the envelope approximation (Huber & Konstandin 2008).
    Subdominant for transitions with plasma coupling.
    """
    alpha = params['alpha']
    beta_H = params['beta_over_H']
    v_w = params['v_w']
    T_GeV = params['T_star_GeV']
    g_star = params['g_star']

    # Efficiency: kappa_col ~ 1 - kappa_v (uncaptured vacuum energy)
    kappa_v = params['kappa_v']
    kappa_col = 1.0 - kappa_v / (alpha / (1 + alpha))
    kappa_col = max(kappa_col * alpha / (1 + alpha), 0.01 * alpha / (1 + alpha))

    # Peak frequency
    f_peak = 1.65e-5 * (0.62 / (1.8 - 0.1 * v_w + v_w**2)) * beta_H * (T_GeV / 100.0) * (g_star / 100)**0.1667

    # Peak amplitude (envelope approximation)
    h2_Omega_peak = (1.67e-5 * (1.0 / beta_H)**2 * (kappa_col * alpha / (1 + alpha))**2
                     * (100.0 / g_star)**(1.0/3.0) * (0.11 * v_w**3 / (0.42 + v_w**2)))

    # Spectral shape
    s = f_Hz / f_peak
    a_col, b_col = 2.8, 1.0
    spectral_shape = (a_col + b_col) / (b_col * s**(-a_col/b_col) + a_col * s**(a_col))

    return h2_Omega_peak * spectral_shape


def Omega_GW_turbulence(f_Hz, params):
    """
    GW energy density spectrum from MHD turbulence.

    Uses the Caprini et al. (2009) fitting formula.
    Subdominant and uncertain.
    """
    alpha = params['alpha']
    beta_H = params['beta_over_H']
    v_w = params['v_w']
    T_GeV = params['T_star_GeV']
    g_star = params['g_star']
    kappa_v = params['kappa_v']

    # Turbulence efficiency: epsilon ~ 0.05-0.10 of kinetic energy
    epsilon_turb = 0.05
    kappa_turb = epsilon_turb * kappa_v

    # Peak frequency
    f_peak = 2.7e-5 * (1.0 / v_w) * beta_H * (T_GeV / 100.0) * (g_star / 100)**0.1667

    # Peak amplitude
    h2_Omega_peak = (3.35e-4 * (1.0 / beta_H)
                     * (kappa_turb * alpha / (1 + alpha))**(3.0/2.0)
                     * (100.0 / g_star)**(1.0/3.0) * v_w)

    # Spectral shape (Kolmogorov + Hubble cutoff)
    s = f_Hz / f_peak
    h_star = params['H_star_GeV'] * 1e9 / (2 * np.pi * 6.582e-16)  # H* in Hz (approximate)
    # Simplified: just use power law
    spectral_shape = s**3 / (1.0 + s)**(11.0/3.0) / (1.0 + 8 * np.pi * f_Hz / h_star)

    return h2_Omega_peak * np.abs(spectral_shape)


# ============================================================
# NANOGrav 15-year free spectrum data (approximate)
# From Agazie et al. (2023), Table 1
# ============================================================
nanograv_freqs_nHz = np.array([2.0, 4.0, 5.9, 7.9, 9.9, 11.9, 13.8, 15.8,
                               17.8, 19.8, 21.8, 23.7, 25.7, 27.7])
# Approximate h^2 Omega_GW from characteristic strain h_c
# h_c^2 = (2 pi^2 / 3 H0^2) * f^2 * Omega_GW
# Using NANOGrav 15yr free spectrum (approximate central values)
nanograv_h2Omega = np.array([2.2e-9, 6.3e-10, 2.5e-10, 1.4e-10, 8.0e-11,
                             5.5e-11, 4.5e-11, 3.0e-11, 2.5e-11, 2.0e-11,
                             1.8e-11, 1.5e-11, 1.2e-11, 1.0e-11])
# Upper and lower bounds (approximate 1-sigma)
nanograv_upper = nanograv_h2Omega * 3.0
nanograv_lower = nanograv_h2Omega / 3.0


# ============================================================
# PTA sensitivity curves (approximate)
# ============================================================
def pta_sensitivity(f_Hz, experiment='SKA'):
    """Approximate PTA power-law integrated sensitivity."""
    f_nHz = f_Hz * 1e9
    if experiment == 'NANOGrav':
        # Approximate NANOGrav 15yr sensitivity
        A_ref = 2.4e-15
        f_ref = 1e-8  # 10 nHz
        gamma = 13.0 / 3.0
        Omega = (2 * np.pi**2 / (3 * H0_per_sec**2)) * f_Hz**2 * A_ref**2 * (f_Hz / f_ref)**(3 - gamma)
        return h_hubble**2 * Omega
    elif experiment == 'SKA':
        # SKA projected (factor ~10 improvement)
        return pta_sensitivity(f_Hz, 'NANOGrav') / 100.0
    return None


# ============================================================
# Main computation and plotting
# ============================================================

def main():
    """Compute and plot the GW spectrum for several lattice scenarios."""

    f_Hz = np.logspace(-10, -5, 1000)  # 0.1 nHz to 10 muHz

    # Scenarios: varying alpha and T*
    scenarios = [
        {'label': r'$\alpha=0.03$, $T_*=156$ MeV (conservative)',
         'T_star': 156.0, 'alpha': 0.03, 'g_star': 47.5,
         'color': '#2980b9', 'ls': '-'},
        {'label': r'$\alpha=0.1$, $T_*=180$ MeV (bag model)',
         'T_star': 180.0, 'alpha': 0.10, 'g_star': 55.0,
         'color': '#27ae60', 'ls': '-'},
        {'label': r'$\alpha=0.3$, $T_*=200$ MeV (strong)',
         'T_star': 200.0, 'alpha': 0.30, 'g_star': 60.0,
         'color': '#c0392b', 'ls': '-'},
        {'label': r'$\alpha=1.0$, $T_*=220$ MeV (total conversion)',
         'T_star': 220.0, 'alpha': 1.00, 'g_star': 62.0,
         'color': '#8e44ad', 'ls': '--'},
    ]

    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    print("=" * 60)
    print("GW SPECTRUM PREDICTIONS")
    print("=" * 60)

    for sc in scenarios:
        params = compute_PT_parameters(sc['T_star'], sc['alpha'], sc['g_star'])

        # Sound waves (dominant)
        h2_Omega_sw, f_peak_sw, h2_peak_sw = Omega_GW_sound_waves(f_Hz, params)

        # Bubble collisions
        h2_Omega_bc = Omega_GW_bubble_collisions(f_Hz, params)

        # Turbulence
        h2_Omega_turb = Omega_GW_turbulence(f_Hz, params)

        # Total
        h2_Omega_total = h2_Omega_sw + h2_Omega_bc + h2_Omega_turb

        ax.loglog(f_Hz * 1e9, h2_Omega_total, color=sc['color'],
                  ls=sc['ls'], lw=2.0, label=sc['label'])

        print(f"\n  {sc['label']}:")
        print(f"    T*         = {params['T_star_MeV']:.1f} MeV")
        print(f"    alpha      = {params['alpha']:.3f}")
        print(f"    beta/H*    = {params['beta_over_H']:.2f}")
        print(f"    v_w        = {params['v_w']:.4f} (= 1/sqrt(3))")
        print(f"    kappa_v    = {params['kappa_v']:.4f}")
        print(f"    Upsilon    = {params['Upsilon']:.4f}")
        print(f"    DeltaG*    = {params['DeltaG_star']:.1f} MeV")
        print(f"    f_peak(sw) = {f_peak_sw*1e9:.1f} nHz")
        print(f"    h2 Omega_peak(sw) = {h2_peak_sw:.2e}")

    # NANOGrav 15-year data
    ax.fill_between(nanograv_freqs_nHz, nanograv_lower, nanograv_upper,
                    alpha=0.25, color='gray', label='NANOGrav 15yr (approx.)')
    ax.scatter(nanograv_freqs_nHz, nanograv_h2Omega, color='black',
               s=25, zorder=5, marker='o')

    # Formatting
    ax.set_xlim(0.3, 1000)
    ax.set_ylim(1e-14, 1e-7)
    ax.set_xlabel(r'Frequency $f$ [nHz]', fontsize=14)
    ax.set_ylabel(r'$h^2 \Omega_{\mathrm{GW}}(f)$', fontsize=14)
    ax.set_title('Gravitational wave spectrum from vacuum lattice crystallisation',
                 fontsize=13)
    ax.legend(fontsize=10, loc='upper right', framealpha=0.9)
    ax.grid(True, which='both', alpha=0.2)

    # Mark PTA bands
    ax.axvspan(1, 100, alpha=0.05, color='blue')
    ax.text(8, 5e-8, 'PTA band', fontsize=10, color='blue', alpha=0.6,
            ha='center')

    plt.tight_layout()
    plt.savefig('gw_spectrum_crystallisation.pdf', dpi=150)
    plt.savefig('gw_spectrum_crystallisation.png', dpi=150)
    print("\n  Saved: gw_spectrum_crystallisation.pdf/png")

    # --------------------------------------------------------
    # Second figure: decomposition for the bag-model scenario
    # --------------------------------------------------------
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 7))

    params_bag = compute_PT_parameters(180.0, 0.10, 55.0)
    h2_sw, _, _ = Omega_GW_sound_waves(f_Hz, params_bag)
    h2_bc = Omega_GW_bubble_collisions(f_Hz, params_bag)
    h2_turb = Omega_GW_turbulence(f_Hz, params_bag)
    h2_total = h2_sw + h2_bc + h2_turb

    ax2.loglog(f_Hz * 1e9, h2_total, 'k-', lw=2.5, label='Total')
    ax2.loglog(f_Hz * 1e9, h2_sw, '--', color='#2980b9', lw=1.5,
               label='Sound waves')
    ax2.loglog(f_Hz * 1e9, h2_bc, '--', color='#27ae60', lw=1.5,
               label='Bubble collisions')
    ax2.loglog(f_Hz * 1e9, h2_turb, '--', color='#c0392b', lw=1.5,
               label='MHD turbulence')

    ax2.fill_between(nanograv_freqs_nHz, nanograv_lower, nanograv_upper,
                     alpha=0.25, color='gray', label='NANOGrav 15yr (approx.)')
    ax2.scatter(nanograv_freqs_nHz, nanograv_h2Omega, color='black',
                s=25, zorder=5, marker='o')

    ax2.set_xlim(0.3, 1000)
    ax2.set_ylim(1e-14, 1e-7)
    ax2.set_xlabel(r'Frequency $f$ [nHz]', fontsize=14)
    ax2.set_ylabel(r'$h^2 \Omega_{\mathrm{GW}}(f)$', fontsize=14)
    ax2.set_title(r'Spectral decomposition: $\alpha=0.1$, $T_*=180$ MeV',
                  fontsize=13)
    ax2.legend(fontsize=11, loc='upper right', framealpha=0.9)
    ax2.grid(True, which='both', alpha=0.2)
    ax2.axvspan(1, 100, alpha=0.05, color='blue')

    plt.tight_layout()
    plt.savefig('gw_spectrum_decomposition.pdf', dpi=150)
    plt.savefig('gw_spectrum_decomposition.png', dpi=150)
    print("  Saved: gw_spectrum_decomposition.pdf/png")


if __name__ == '__main__':
    main()
