"""
Vacuum-channel budget for the superconductivity appendix (hard nodes D and H).

Computes, from framework constants only:
  1. The averted catastrophe: what an OPEN vortex-line electron-electron
     interaction would have cost (per core length), versus Coulomb.
  2. The dipole channel: condensate flow dipole-dipole energy between two
     closed-winding electrons, compared against the magnetic dipole-dipole
     energy, with the analytic ratio 4 pi^3 f_s alpha (a/l)^4.
  3. The effective ring radius a_eff forced by the locked-faces identity
     (condensate face == crystal/magnetic face of one interaction).
  4. Scale comparison against BCS quantities (gap, Debye, coherence length).
  5. Operator check: every component of S_tot annihilates the two-electron
     singlet; the difference operator (sigma1 - sigma2) does not.
  6. The channel-budget figure, fig_sc_channel_budget.pdf: interaction
     energy against electron separation for each channel of the budget.

Everything is evaluated at the adopted superfluid fraction f_s = 0.99932.
All CODATA values as in the project codata.txt (2022 adjustment).
"""

import os

import numpy as np
import sympy as sp

# ---------------------------------------------------------------- constants
alpha  = 7.2973525643e-3          # fine structure constant
hbar   = 1.054571817e-34          # J s
h      = 6.62607015e-34           # J s
c      = 2.99792458e8             # m/s
m_e    = 9.1093837139e-31         # kg
e      = 1.602176634e-19          # C
mu0    = 1.25663706127e-6         # N/A^2
eV     = e                        # J

# framework constants (bootstrap: l = r_e, m0 = m_e/alpha)
l_     = 2.8179403205e-15         # m, classical electron radius = lattice spacing
m0     = m_e / alpha              # kg, node mass (~70 MeV)
f_s    = 0.99932                  # the Leggett ceiling adopted as the working superfluid fraction
rho    = m0 / l_**3               # kg/m^3, total density
rho_s  = f_s * rho                # superfluid density
kappa  = h / m0                   # m^2/s, circulation quantum of the condensate
mu_B   = e * hbar / (2 * m_e)     # Bohr magneton

print("=== framework inputs ===")
print(f"m0 c^2         = {m0*c**2/eV/1e6:.3f} MeV")
print(f"kappa = h/m0   = {kappa:.4e} m^2/s   (check 2 pi c l = {2*np.pi*c*l_:.4e})")
print(f"rho_s          = {rho_s:.4e} kg/m^3")

# ------------------------------------------- 1. the averted line catastrophe
# Open vortex lines: interaction energy per length (rho_s kappa^2 / 2 pi) ln(R/L).
# Quote the scale at log = 1, per core length l, in MeV; compare Coulomb at l.
# With kappa = 2 pi c l and rho_s = f_s m0/l^3 this is 2 pi f_s m0 c^2, and
# Coulomb at l is alpha hbar c/l = m_e c^2, so the ratio is 2 pi f_s/alpha.
E_line_per_l = (rho_s * kappa**2 / (2*np.pi)) * l_          # J per core length
E_coul_at_l  = alpha * hbar * c / l_                          # = m_e c^2
print("\n=== 1. open-line catastrophe (what compactness prevents) ===")
print(f"line-line scale per core length: {E_line_per_l/eV/1e6:.1f} MeV (x log)")
print(f"Coulomb energy at r = l:         {E_coul_at_l/eV/1e6:.3f} MeV")
print(f"ratio (line/Coulomb at core):    {E_line_per_l/E_coul_at_l:.0f}"
      f"   (closed form 2 pi f_s/alpha = {2*np.pi*f_s/alpha:.0f})")
# consistency: the same line tension over one lattice spacing prices the
# smallest free ring, E1 = 2 pi^2 f_s m0 c^2 (the forge entry price)
E1 = 2*np.pi**2 * f_s * m0 * c**2
print(f"framework ring price E1:         {E1/eV/1e9:.3f} GeV (= pi x line scale)")

# ------------------------------------------------------- 2. dipole channel
# Closed winding of radius a: hydrodynamic dipole d = kappa * pi * a^2
# (current-loop correspondence u <-> B, kappa <-> mu0 I, d <-> m = I pi a^2).
# Cross kinetic energy of two dipoles at separation L (coefficient as in
# magnetostatics): E = (rho_s / 4 pi) * d1 d2 * f_angle / L^3, f_angle ~ 1.
def E_vac_dipole(L, a=l_):
    d = kappa * np.pi * a**2
    return (rho_s / (4*np.pi)) * d**2 / L**3

def E_mag_dipole(L):
    return (mu0 / (4*np.pi)) * mu_B**2 / L**3

print("\n=== 2. dipole channel vs magnetic dipole-dipole ===")
for L, name in [(2e-10, "2 Angstrom (metal)"), (4e-8, "40 nm (Cooper pair)")]:
    ev_, em_ = E_vac_dipole(L), E_mag_dipole(L)
    print(f"L = {name:22s}: E_vac = {ev_/eV:.3e} eV,  E_mag = {em_/eV:.3e} eV")

# analytic ratio, symbolically: E_vac/E_mag = 4 pi^3 f_s alpha (a/l)^4
a_sym, l_sym, fs_sym, al_sym, hb, m0s, cs = sp.symbols(
    'a ell f_s alpha hbar m_0 c', positive=True)
me_s   = al_sym * m0s
rho_ss = fs_sym * m0s / l_sym**3
kap_s  = 2*sp.pi*cs*l_sym                       # h/m0 with l = hbar/(m0 c)
d_s    = kap_s * sp.pi * a_sym**2
Evac_s = rho_ss * d_s**2 / (4*sp.pi)            # coefficient of 1/L^3
mu_Bs  = sp.symbols('e', positive=True) * hb / (2*me_s)
# mu0 mu_B^2 = pi alpha hbar^3 / (m_e^2 c) after e^2/(4 pi eps0) = alpha hbar c
Emag_s = sp.pi * al_sym * hb**3 / (me_s**2 * cs) / (4*sp.pi)
ratio  = sp.simplify(Evac_s / Emag_s)
ratio  = ratio.subs(hb, m0s*cs*l_sym)           # hbar = m0 c l (bootstrap)
ratio  = sp.simplify(ratio)
print(f"\nanalytic ratio E_vac/E_mag = {ratio}  [expect 4 pi^3 f_s alpha (a/l)^4]")
rnum = float(ratio.subs({a_sym: 1, l_sym: 1, fs_sym: f_s, al_sym: alpha}))
print(f"numeric at a = l:  {rnum:.4f}")

# 3. effective radius forced by the locked-faces identity (ratio == 1)
a_eff = (4*np.pi**3 * f_s * alpha)**(-0.25)
print(f"\n=== 3. locked-faces closure ===")
# The electron's own ring geometry, sqrt(6) times the spin-half column radius
# sqrt(1/(2 pi f_s)) l, gives sqrt(3/(pi f_s)) l; it does not close on a_eff.
a_ring = np.sqrt(3/(np.pi*f_s))
print(f"a_eff / l = {a_eff:.4f}")
print(f"ring geometry sqrt(3/(pi f_s)) = {a_ring:.4f} l  "
      f"({(1 - a_ring/a_eff)*100:.1f}% short of a_eff: a coincidence, not a closure)")

# --------------------------------------------- 4. scale comparison with BCS
Delta_Nb   = 1.55e-3 * eV      # Nb gap
hw_D_Nb    = 24.0e-3 * eV      # Nb Debye energy ~ 275 K
print("\n=== 4. against BCS scales (Nb) ===")
print(f"gap 2Delta = {2*Delta_Nb/eV*1e3:.2f} meV;  hbar w_D = {hw_D_Nb/eV*1e3:.1f} meV")
print(f"E_vac(2 A)/2Delta = {E_vac_dipole(2e-10)/(2*Delta_Nb):.2e}")
print(f"E_vac(xi)/2Delta  = {E_vac_dipole(4e-8)/(2*Delta_Nb):.2e}")

# ------------------------------- 5. singlet darkness: operator verification
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]])
sz = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2)
def two(op, which):
    return np.kron(op, I2) if which == 1 else np.kron(I2, op)
singlet = (np.kron([1,0],[0,1]) - np.kron([0,1],[1,0])) / np.sqrt(2)
print("\n=== 5. singlet darkness (operator level) ===")
for name, op in [("x", sx), ("y", sy), ("z", sz)]:
    tot  = two(op,1) + two(op,2)
    diff = two(op,1) - two(op,2)
    print(f"|(s1+s2)_{name} |singlet>| = {np.linalg.norm(tot @ singlet):.2e}   "
          f"|(s1-s2)_{name} |singlet>| = {np.linalg.norm(diff @ singlet):.2f}")

# compactness cross-check against the monograph's quoted suppression
print(f"\n(l / 1 cm)^2 = {(l_/1e-2)**2:.1e}  (monograph: ~1e-25 at a centimetre)")

# ---------------------------------------- 6. the channel-budget figure
# Energy against separation L for each line of the budget:
#   open-line channel   (rho_s kappa^2 / 2 pi) l ln(L/l)  (removed by closure)
#   Coulomb             alpha hbar c / L
#   closed-winding dipole (condensate face), E_vac_dipole at a = l
#   magnetic dipole-dipole (crystal face),   E_mag_dipole
# against the niobium pairing scale 2 Delta, with the metallic spacing 2 A
# and the Nb coherence length 40 nm marked.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L = np.logspace(np.log10(l_), np.log10(3e-7), 600)          # m
L_open = L[L > l_ * 1.0001]                                  # ln(L/l) > 0
nm = 1e-9
fig, ax = plt.subplots(figsize=(7.2, 5.0))
ax.loglog(L_open/nm, (rho_s*kappa**2/(2*np.pi))*l_*np.log(L_open/l_)/eV,
          "--", color="#b22222", lw=2.0, label="open-line channel (removed by closure)")
ax.loglog(L/nm, alpha*hbar*c/L/eV, "-", color="#666666", lw=1.5,
          label="Coulomb (crystal shear face)")
ax.loglog(L/nm, E_vac_dipole(L)/eV, "-", color="#1f77b4", lw=2.5,
          label="closed-winding dipole (condensate face)")
ax.loglog(L/nm, E_mag_dipole(L)/eV, ":", color="#1f77b4", lw=3.0,
          label="magnetic dipole--dipole (crystal face)")
ax.axhline(2*Delta_Nb/eV, color="#2e8b57", lw=1.5)
ax.text(3e-6, 2*Delta_Nb/eV*2.0, r"$2\Delta$ (Nb)", color="#2e8b57", fontsize=9)
for Lm in (2e-10, 4e-8):
    ax.axvline(Lm/nm, color="#bbbbbb", lw=1.0)
ax.text(2e-10/nm*1.15, 0.3, "2 \N{LATIN CAPITAL LETTER A WITH RING ABOVE}", color="#888888", fontsize=8)
ax.text(4e-8/nm*1.15, 0.3, r"$\xi_{\mathrm{Nb}}$", color="#888888", fontsize=8)
ax.set_xlim(l_/nm, 3e-7/nm)
ax.set_ylim(1e-14, 1e9)
ax.set_xlabel(r"electron separation $L$ [nm]")
ax.set_ylabel("interaction energy [eV]")
ax.grid(True, which="both", lw=0.3, alpha=0.4)
ax.legend(loc="upper right", fontsize=9)
fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig_sc_channel_budget.pdf")
fig.savefig(out, metadata={"CreationDate": None})
print(f"\nwrote {out}")
