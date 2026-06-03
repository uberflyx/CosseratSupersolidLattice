#!/usr/bin/env python3
r"""
Magnetic self-energy contribution to the proton-neutron mass splitting.

The neutron is electrically neutral, but its magnetic moment is large
(mu_n = -1.913 nuclear magnetons), so its magnetic self-energy is NOT
negligible.  What enters the splitting is the difference of field energies,
proportional to mu_p^2 - mu_n^2, and with the framework's derived moment ratio
this is a clean 5/9 of the proton value: the neutron removes four-ninths.

Physics
-------
The nucleon spin is a microrotation of the defect core; a rotating core sources
a strain field, and the energy held in that field is the magnetic self-energy.
For a moment mu spread with magnetisation form factor rho(q) (normalised to
rho(0)=1), the field energy is

      U_M = (mu_0 mu^2 / 6 pi^2) * INT_0^inf q^2 |rho(q)|^2 dq.                (*)

Only the transverse part of the magnetisation carries field energy; the angular
average of sin^2(theta_q) = 2/3 with the Fourier measure gives the 1/(6 pi^2).
For a uniformly magnetised sphere of radius R this integrates to the familiar
closed form U_M = mu_0 mu^2 / (4 pi R^3), used here to validate the integral.

The term that enters the splitting:

      delta_mag = U_M(n) - U_M(p) = -(1 - (mu_n/mu_p)^2) U_M(p) < 0,

proton heavier (the stronger magnet).  With the framework ratio mu_n/mu_p = -2/3
(Y-junction + Fermi statistics; appendix), 1 - (mu_n/mu_p)^2 = 5/9.

Scale and shape
---------------
R is the nucleon electromagnetic radius.  The framework derives the proton
charge radius R_p = 0.845 fm = 0.30 l from the Peierls-Nabarro core profile,
matched to the world average at 0.5%; the proton magnetic radius is measured
close to this.  The result is robust to the profile shape at the ~15% level:
the empirical dipole form factor and a uniform sphere of the same radius agree
to better than that, because the moment and the radius are fixed and only the
distribution between them is free.

Sub-lattice structure (a physical check)
----------------------------------------
The nucleon (R_p = 0.30 l) is SMALLER than the lattice spacing, so the integral
(*) draws on wavenumbers above the Brillouin-zone edge pi/l.  The magnetic
self-energy is therefore the defect's continuum near-field, the strain inside
the core, not a sum over lattice phonons.  Cutting the integral at the zone edge
removes ~9/10 of it -- the wrong thing to do, since the sub-cell structure that
holds the field is the same structure that fixes the charge radius.  The script
reports this cut as a check, not as the physical value.

Companion to nucleon_tensor_channel.py / nucleon_tensor_magnitude.py: the
magnetic self-energy is the spin-DEPENDENT part of the microrotation
(couple-stress) sector; the tensor admixture is the spin-INDEPENDENT part.
Both live in the sector that also carries the graviton and the f2(1270).
"""

import numpy as np
from scipy import integrate

# ---- constants (SI; CODATA / project files) --------------------------------
MU0_4PI = 1.0e-7              # mu_0 / 4pi  [T m / A], exact
MU_N    = 5.0507837e-27       # nuclear magneton [J/T]
FM      = 1.0e-15             # metre
MEV     = 1.602176634e-13     # joule
M_E_MEV = 0.51099895          # electron mass [MeV]

# ---- measured nucleon moments (nuclear magnetons) --------------------------
MU_P        = 2.7928473
MU_N_MOMENT = -1.9130427
RATIO_MEAS  = MU_N_MOMENT / MU_P         # ~ -0.685
FACTOR_FW   = 1.0 - (2.0 / 3.0)**2       # 5/9, framework ratio -2/3
FACTOR_MEAS = 1.0 - RATIO_MEAS**2        # measured

# ---- nucleon electromagnetic size (framework-derived) ----------------------
ELL     = 2.8179403e-15 / FM             # lattice spacing = r_e [fm]
R_P     = 0.845                          # proton charge radius [fm] = 0.30 l
R_MAG   = 0.85                           # proton magnetic radius [fm], measured
BZ_EDGE = np.pi / ELL                    # Brillouin-zone edge pi/l [1/fm]


def field_energy(mu_in_muN, form_factor, q_max):
    """Magnetic field energy U_M from Eq. (*), in MeV.

    mu_in_muN   : moment in nuclear magnetons
    form_factor : callable rho(q), q in 1/fm, rho(0)=1
    q_max       : upper integration limit [1/fm]
    """
    mu = mu_in_muN * MU_N
    prefactor = (4 * np.pi * MU0_4PI) * mu**2 / (6 * np.pi**2)   # mu_0 mu^2 / 6pi^2
    integral_invfm3, _ = integrate.quad(
        lambda q: q**2 * form_factor(q)**2, 0.0, q_max, limit=400)
    integral_SI = integral_invfm3 / FM**3                        # 1/m^3
    return prefactor * integral_SI / MEV


def ff_uniform_sphere(R):
    """Form factor of a uniformly magnetised sphere of radius R [fm]."""
    def f(q):
        x = q * R
        if x < 1e-8:
            return 1.0
        return 3.0 * (np.sin(x) - x * np.cos(x)) / x**3
    return f


def ff_dipole(R_rms):
    """Dipole form factor 1/(1+q^2/MD^2)^2 with RMS radius R_rms [fm].

    The dipole has <r^2> = 12/MD^2, the empirical nucleon shape (exponential
    magnetisation density), which the Peierls-Nabarro core profile reproduces.
    """
    M_D = np.sqrt(12.0) / R_rms          # [1/fm]
    return lambda q: 1.0 / (1.0 + (q / M_D)**2)**2


def U_sphere_analytic(mu_in_muN, R):
    """Closed form mu_0 mu^2 / (4 pi R^3) [MeV] -- validation reference."""
    mu = mu_in_muN * MU_N
    return MU0_4PI * mu**2 / (R * FM)**3 / MEV


# ---------------------------------------------------------------------------
print("=" * 74)
print("VALIDATION: numerical integral vs analytic uniformly magnetised sphere")
print("=" * 74)
for R in (0.70, 0.845, 1.00):
    num = field_energy(MU_P, ff_uniform_sphere(R), q_max=1000.0)
    ana = U_sphere_analytic(MU_P, R)
    print(f"  R = {R:.3f} fm :  numerical = {num:.4f} MeV   "
          f"analytic = {ana:.4f} MeV   (diff {abs(num - ana) / ana * 100:.2f}%)")

print()
print("=" * 74)
print("PROTON MAGNETIC SELF-ENERGY  (proper field-energy integral)")
print("=" * 74)
U_sphere = field_energy(MU_P, ff_uniform_sphere(R_P), q_max=1000.0)
U_dipole = field_energy(MU_P, ff_dipole(R_MAG),       q_max=1000.0)
print(f"  uniform sphere, R = R_p   = {R_P:.3f} fm :  U_M(p) = {U_sphere:.4f} MeV")
print(f"  dipole,        R_mag      = {R_MAG:.3f} fm :  U_M(p) = {U_dipole:.4f} MeV")
print(f"  shape spread (dipole vs sphere)            : "
      f"{abs(U_sphere - U_dipole) / U_dipole * 100:.0f}%")
U_p = 0.5 * (U_sphere + U_dipole)        # central value
print(f"  central U_M(p)                             : {U_p:.4f} MeV")

print()
print("-" * 74)
print("NEUTRON IS NOT NEGLIGIBLE")
print("-" * 74)
print(f"  measured ratio  mu_n/mu_p = {RATIO_MEAS:+.4f}   "
      f"(framework -2/3 = {-2/3:+.4f})")
print(f"  U_M(n)/U_M(p) = (mu_n/mu_p)^2 = {RATIO_MEAS**2:.3f}  "
      f"-> neutron removes {RATIO_MEAS**2 * 100:.0f}% of the proton term")
print(f"  enters as  mu_p^2 - mu_n^2 = (1 - (mu_n/mu_p)^2) mu_p^2")
print(f"    framework factor  = 5/9 = {FACTOR_FW:.4f}")
print(f"    measured factor         = {FACTOR_MEAS:.4f}")

delta_mag = -FACTOR_FW * U_p
print()
print(f"  delta_mag = -(5/9) U_M(p) = {delta_mag:+.4f} MeV "
      f"= {delta_mag / M_E_MEV:+.3f} m_e   (proton heavier)")

print()
print("-" * 74)
print("SUB-LATTICE CHECK  (the term is the defect's continuum near-field)")
print("-" * 74)
print(f"  nucleon size R_p = {R_P / ELL:.2f} l  <  lattice spacing l = {ELL:.2f} fm")
print(f"  Brillouin-zone edge  pi/l = {BZ_EDGE:.3f} 1/fm")
U_full = field_energy(MU_P, ff_dipole(R_MAG), q_max=1000.0)
U_cut  = field_energy(MU_P, ff_dipole(R_MAG), q_max=BZ_EDGE)
print(f"  U_M(p) full                 = {U_full:.4f} MeV")
print(f"  U_M(p) cut at zone edge     = {U_cut:.4f} MeV  "
      f"({U_cut / U_full * 100:.0f}% -- wrong: discards sub-cell core structure)")

print()
print("-" * 74)
print("RADIUS SENSITIVITY  (delta_mag ~ 1/R^3, dipole form factor)")
print("-" * 74)
for R in (0.80, 0.84, 0.85, 0.86, 0.90):
    U = field_energy(MU_P, ff_dipole(R), q_max=1000.0)
    d = -FACTOR_FW * U
    tag = "  <- magnetic radius" if abs(R - R_MAG) < 1e-9 else ""
    print(f"  R = {R:.2f} fm :  delta_mag = {d:+.4f} MeV "
          f"= {d / M_E_MEV:+.3f} m_e{tag}")

# ---- the full budget -------------------------------------------------------
print()
print("=" * 74)
print("BUDGET  (m_n - m_p, MeV;  + = neutron heavier)")
print("=" * 74)
strong   = 5 * M_E_MEV
coulomb  = -(5 / 3) * M_E_MEV
tensor   = -0.075                        # couple-stress, spin-independent
magnetic = delta_mag                     # couple-stress, spin-dependent
subtotal = strong + coulomb + tensor + magnetic
observed = 1.2933324
residual = observed - subtotal
rows = [("strong character  5 m_e",              strong),
        ("electric Coulomb  -5/3 m_e",           coulomb),
        ("couple-stress, spin-independent",      tensor),
        ("couple-stress, spin-dependent (mag)",  magnetic)]
for name, val in rows:
    print(f"  {name:<38s} {val:+.4f}")
print("  " + "-" * 48)
print(f"  {'derived subtotal':<38s} {subtotal:+.4f}")
print(f"  {'observed':<38s} {observed:+.4f}")
print(f"  {'residual (three-body + inelastic)':<38s} {residual:+.4f} "
      f"= {residual / M_E_MEV:+.3f} m_e")

em_total = coulomb + tensor + magnetic + residual
print()
print(f"  EM total (Coulomb + tensor + mag + residual) = {em_total:+.4f} MeV")
print(f"  dispersive EM:  WLCM 2012 -1.30(47);  needed {observed - strong:+.4f}")
print(f"  couple-stress sector (tensor + magnetic)     = {tensor + magnetic:+.4f} MeV"
      f" = {(tensor + magnetic) / M_E_MEV:+.3f} m_e")
