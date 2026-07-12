"""
Figure: the vacuum-channel budget between two conduction electrons.

Plots the interaction energy carried by each condensate channel against the
electron separation L, alongside the ordinary Coulomb energy and the Nb
pairing scale 2*Delta.  Three curves tell the story:

  1. The interaction two OPEN vortex lines would have (log-growing, hundreds
     of MeV per core length): removed exactly by the closure of the winding.
  2. The surviving dipole channel of the closed winding (1/L^3), which lands
     within an order-unity factor of...
  3. ...the magnetic dipole-dipole energy, its crystal-face twin.

The near-coincidence of 2 and 3 is the locked-faces identity; the gulf
between the dipole channel and 2*Delta is why the vacuum supplies no
pairing glue.

Output: fig_sc_channel_budget.pdf (vector, for the monograph figures/ dir).
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------- constants
alpha = 7.2973525643e-3
hbar = 1.054571817e-34            # J s
h = 6.62607015e-34                # J s
c = 2.99792458e8                  # m/s
m_e = 9.1093837139e-31            # kg
e = 1.602176634e-19               # C
mu0 = 1.25663706127e-6            # N/A^2
eps0 = 8.8541878188e-12           # F/m
eV = e

l_ = 2.8179403205e-15             # lattice spacing = classical electron radius
m0 = m_e / alpha                  # node mass
f_s = 4.0 / 5.0                   # superfluid fraction
rho_s = f_s * m0 / l_**3          # superfluid density
kappa = h / m0                    # circulation quantum
mu_B = e * hbar / (2 * m_e)       # Bohr magneton

L = np.logspace(np.log10(3e-15), np.log10(3e-7), 600)   # 3 fm to 300 nm

# channel energies [J]
E_line = (rho_s * kappa**2 / (2 * np.pi)) * l_ * np.log(np.maximum(L / l_, 1.0))
E_dip = (rho_s / (4 * np.pi)) * (kappa * np.pi * l_**2) ** 2 / L**3
E_mag = (mu0 / (4 * np.pi)) * mu_B**2 / L**3
E_coul = e**2 / (4 * np.pi * eps0 * L)
gap_Nb = 3.1e-3 * eV              # 2*Delta for niobium

fig, ax = plt.subplots(figsize=(7.2, 5.0))

ax.loglog(L * 1e9, E_line / eV, "--", color="#b22222", lw=1.8,
          label="open-line channel (removed by closure)")
ax.loglog(L * 1e9, E_coul / eV, "-", color="#666666", lw=1.4,
          label="Coulomb (crystal shear face)")
ax.loglog(L * 1e9, E_dip / eV, "-", color="#1f77b4", lw=2.2,
          label="closed-winding dipole (condensate face)")
ax.loglog(L * 1e9, E_mag / eV, ":", color="#1f77b4", lw=2.2,
          label="magnetic dipole--dipole (crystal face)")
ax.axhline(gap_Nb / eV, color="#2e8b57", lw=1.4)
ax.text(4e-6 * 1e9, 4.0e-3, r"$2\Delta$ (Nb)", color="#2e8b57", fontsize=10)

# markers: metallic spacing and Cooper-pair size
for x, tag in [(2e-10, "2 \u00c5"), (3.8e-8, r"$\xi_{\mathrm{Nb}}$")]:
    ax.axvline(x * 1e9, color="0.8", lw=0.9, zorder=0)
    ax.text(x * 1e9 * 1.1, 3e8, tag, fontsize=9, color="0.45")

ax.set_xlabel(r"electron separation $L$ [nm]")
ax.set_ylabel(r"interaction energy [eV]")
ax.set_ylim(1e-14, 3e9)
ax.set_xlim(L[0] * 1e9, L[-1] * 1e9)
ax.legend(loc="upper right", fontsize=9, framealpha=0.95)
ax.grid(alpha=0.25, which="both", lw=0.4)

fig.tight_layout()
fig.savefig("fig_sc_channel_budget.pdf")
print("wrote fig_sc_channel_budget.pdf")
