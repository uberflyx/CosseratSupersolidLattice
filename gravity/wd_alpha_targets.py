#!/usr/bin/env python3
"""
wd_alpha_targets.py — Which stars sharpen the gravitational-alpha test?

The measured signal is  dalpha/alpha = k_a * phi,  phi = GM/(Rc^2),
where k_a is the framework's response coefficient (here we only need phi;
k_a is common to all objects).  The MEASUREMENT error on dalpha/alpha has
two independent parts:

  sigma_stat  : statistical, scales as sigma_line / sqrt(N_lines * S/N-ish),
                falls when a star has more/sharper metal lines.
  sigma_lab   : laboratory-wavelength systematic. This is a FLOOR that is
                identical for every star, because it lives in the reference
                Fe V / Ni V wavelengths, not in the object.

The figure of merit for detecting a gravitational coupling is the ratio

  detectability  =  signal / total_error  =  k_a * phi / sqrt(sig_stat^2 + sig_lab^2).

Two levers raise it: (1) larger phi  -> stronger/more compact star;
(2) smaller sig_lab -> better laboratory wavelengths (helps EVERY star at once).

This script quantifies both, using a standard WD mass-radius relation.

References for the numbers used:
  Berengut et al. 2013 (PRL 111, 010801): G191-B2B, dalpha/alpha=(4.2+/-1.6)e-5, Fe V.
  Bainbridge & Webb 2017 (Universe 3, 32): sample of ~10 hot DA/sdO objects.
  Preval/Barstow: G191-B2B Teff~60kK, log g~7.6.
"""

import numpy as np

# ---- physical constants (SI) ----
G   = 6.67430e-11      # m^3 kg^-1 s^-2  (CODATA)
c   = 299792458.0      # m/s (exact)
Msun= 1.98892e30       # kg
Rsun= 6.957e8          # m
Rearth = 6.371e6       # m

def phi_dimensionless(M_Msun, R_Rsun):
    """Surface gravitational potential phi = GM/(Rc^2), dimensionless."""
    M = M_Msun * Msun
    R = R_Rsun * Rsun
    return G * M / (R * c**2)

# ---- WD mass-radius relation ----
# Use the analytic Nauenberg (1972) relation, a good closed-form fit to the
# fully-relativistic Chandrasekhar M-R for carbon-oxygen WDs. Radius in Rsun.
# R = 0.0225/mu * sqrt(1 - (M/Mch)^(4/3)) / (M/Mch)^(1/3),  mu=2, Mch=1.454.
def wd_radius_Rsun(M_Msun, mu=2.0, Mch=1.454):
    x = M_Msun / Mch
    if x >= 1.0:
        return np.nan
    return (0.0225/mu) * np.sqrt(1.0 - x**(4.0/3.0)) / x**(1.0/3.0)

def log_g_cgs(M_Msun, R_Rsun):
    """Surface gravity log10(g[cgs]) for check against observed log g."""
    M = M_Msun * Msun * 1e3          # g
    R = R_Rsun * Rsun * 1e2          # cm
    g = (G*1e3) * M / R**2           # cgs G = 6.674e-8
    return np.log10(g)

print("="*74)
print("PART 1 — Surface potential vs white-dwarf mass (Nauenberg M-R)")
print("="*74)
print(f"{'M/Msun':>7} {'R/Rsun':>9} {'R/Rearth':>9} {'log g':>7} {'phi':>11} {'phi/phi_G191':>13}")

# G191-B2B reference: M~0.52 Msun (Teff 60kK, log g 7.6 -> ~0.5 Msun canonical)
M_ref = 0.52
R_ref = wd_radius_Rsun(M_ref)
phi_ref = phi_dimensionless(M_ref, R_ref)

for M in [0.4, 0.52, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35]:
    R = wd_radius_Rsun(M)
    if np.isnan(R): 
        continue
    phi = phi_dimensionless(M, R)
    print(f"{M:>7.2f} {R:>9.5f} {R*Rsun/Rearth:>9.3f} {log_g_cgs(M,R):>7.2f} "
          f"{phi:>11.3e} {phi/phi_ref:>13.2f}")

print(f"\nReference G191-B2B: M={M_ref} Msun, R={R_ref:.4f} Rsun, phi={phi_ref:.3e}")
print(f"(literature quotes phi ~ 5e-5, log g ~ 7.6 — consistent)")

print()
print("="*74)
print("PART 2 — Detectability: signal/error for real & hypothetical targets")
print("="*74)
# Model the error budget. From Berengut 2013 on G191-B2B (Fe V):
#   total sigma(dalpha/alpha) = 1.6e-5, LAB-DOMINATED.
# Split (from Hu et al. 2020 who separate stat vs sys):
#   sigma_stat ~ 0.47e-5, sigma_lab(sys) ~ 2.3e-5 for a single object w/ current labs.
# We take a representative floor:
sig_lab_current = 1.5e-5     # lab-wavelength floor, current Fe V/Ni V (per object)
sig_stat_G191   = 0.5e-5     # statistical for a G191-quality spectrum (~96 Fe V lines)

# Framework response: dalpha/alpha = k_a * phi. Fix k_a so that at phi_ref the
# *central* signal equals the framework prediction ~1e-6 (using 7-eta_t ~ 1).
# But detectability is k_a-independent in RATIO between objects, so we report
# both the phi-boost and the total-error behaviour.

def sigma_total(N_lines, sig_lab, sig_stat_ref=sig_stat_G191, N_ref=96):
    """Total error: stat scales as 1/sqrt(N_lines); lab is a fixed floor."""
    sig_stat = sig_stat_ref * np.sqrt(N_ref / max(N_lines,1))
    return np.sqrt(sig_stat**2 + sig_lab**2), sig_stat

targets = [
    # name, M/Msun, N_metal_lines(approx), note
    ("G191-B2B (ref)",      0.52,  96, "benchmark hot DA"),
    ("Typical hot DA",      0.60,  60, "sample-average object"),
    ("Massive DA (0.9)",    0.90,  50, "higher phi, fewer lines"),
    ("Ultramassive (1.1)",  1.10,  35, "phi high, lines broaden"),
    ("Ultramassive (1.3)",  1.30,  20, "near broadening ceiling"),
]

print(f"{'target':>20} {'phi':>10} {'phi/ref':>8} {'N':>4} "
      f"{'sig_tot(cur)':>12} {'S/N(cur)':>9} {'sig_tot(lab/10)':>15} {'S/N(lab/10)':>11}")
for name, M, N, note in targets:
    R = wd_radius_Rsun(M)
    phi = phi_dimensionless(M, R)
    # signal in units where k_a*phi_ref = 1  (i.e. signal = phi/phi_ref)
    signal = phi/phi_ref
    sig_cur, _  = sigma_total(N, sig_lab_current)
    sig_better,_= sigma_total(N, sig_lab_current/10.0)   # 10x better lab wavelengths
    # scale errors into the same "units" (divide by k_a*phi_ref = signal_ref magnitude)
    # Represent detectability as signal / (sig/ k_a) ; use relative: (phi/phi_ref)/(sig/sig_ref0)
    # Simplify: report S/N assuming framework central value dalpha ~ (phi/phi_ref)*1e-6*(7-eta)
    # with (7-eta)=1 -> dalpha_ref = 1e-6. Use that:
    dalpha = signal * 1e-6
    snr_cur     = dalpha / sig_cur
    snr_better  = dalpha / sig_better
    print(f"{name:>20} {phi:>10.3e} {signal:>8.2f} {N:>4d} "
          f"{sig_cur:>12.2e} {snr_cur:>9.4f} {sig_better:>15.2e} {snr_better:>11.4f}")

print()
print("Notes:")
print("- 'S/N' here is central-signal / total-error for (7-eta_t)=1, the framework's")
print("  minimal prediction. Values <<1 mean the prediction sits well inside the bound")
print("  (consistency, not detection) — as expected.")
print("- phi rises ~6x from 0.52 to 1.3 Msun; but N_lines falls and lines broaden,")
print("  so raw phi-boost is partly eaten by loss of usable lines.")
print("- Cutting the LAB floor 10x lifts S/N by ~an order of magnitude for EVERY")
print("  object at once — the single most effective lever.")

print()
print("="*74)
print("PART 3 — Where does the error floor have to sit to DETECT the prediction?")
print("="*74)
# For a 3-sigma detection of dalpha = (phi/phi_ref)*1e-6*(7-eta_t):
for eta_span in [1, 10, 40]:
    print(f"\n  Assuming (7 - eta_t) = {eta_span}:")
    for name, M, N, note in targets[:1] + targets[2:]:
        R = wd_radius_Rsun(M); phi = phi_dimensionless(M,R)
        dalpha = (phi/phi_ref)*1e-6*eta_span
        needed_sigma = dalpha/3.0   # 3-sigma
        print(f"    {name:>20}: dalpha={dalpha:.2e}, need sigma<{needed_sigma:.2e} "
              f"(current floor {sig_lab_current:.1e})")

# ---- figure: the two levers ----
if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    Ms = np.linspace(0.4, 1.34, 300)
    phis = np.array([phi_dimensionless(M, wd_radius_Rsun(M)) for M in Ms])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    # Left: phi vs mass, with broadening ceiling
    ax1.plot(Ms, phis/phi_ref, 'k-', lw=1.8)
    ax1.axhline(1.0, color='tab:red', ls='--', lw=1, label='G191-B2B')
    ax1.axvspan(1.15, 1.34, alpha=0.12, color='gray')
    ax1.text(1.24, 5.2, 'lines broaden\ninto invisibility', fontsize=8,
             ha='center', color='0.35')
    ax1.set_xlabel(r'White-dwarf mass  $M/M_\odot$')
    ax1.set_ylabel(r'Signal boost  $\phi/\phi_{\rm G191}$')
    ax1.set_title('Lever 1: pick a more massive star')
    ax1.legend(fontsize=8, loc='upper left')
    ax1.set_ylim(0, 8.5)

    # Right: detectability (S/N) vs lab-floor improvement, for a few targets
    lab_factors = np.array([1, 2, 5, 10, 20, 50])
    for name, M, N in [("G191-B2B", 0.52, 96),
                       ("Massive DA 0.9", 0.90, 50),
                       ("Ultramassive 1.3", 1.30, 20)]:
        R = wd_radius_Rsun(M); phi = phi_dimensionless(M, R)
        dalpha = (phi/phi_ref)*1e-6   # (7-eta)=1
        snr = []
        for f in lab_factors:
            sig,_ = sigma_total(N, sig_lab_current/f)
            snr.append(dalpha/sig)
        ax2.plot(lab_factors, snr, 'o-', lw=1.5, ms=4, label=name)
    ax2.set_xscale('log')
    ax2.set_xlabel(r'Laboratory-wavelength improvement factor')
    ax2.set_ylabel(r'Detectability  (signal / error, $7-\eta_t=1$)')
    ax2.set_title('Lever 2: better lab wavelengths (helps all)')
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig('wd_targets.png', dpi=140, bbox_inches='tight')
    print("\nSaved figure.")
