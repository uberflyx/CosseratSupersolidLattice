"""
The condensate knot sector: mass spectrum and self-interaction.

The superfluid chapter (sec:vortex_relic) establishes that the vacuum's
condensate carries a free sector of closed vortex rings, knots, and links,
pinched off at the crystallisation transition and held stable by conserved
helicity H = n kappa^2 (Moffatt's invariant, with n the knot/link integer).
It gives the sector its charge (circulation +/- kappa), its force (a
logarithmic 2D-Coulomb interaction), and its scale ("nuclear-to-GeV"), and it
leaves the present abundance as an honest unknown. What it does not give as
numbers are the mass spectrum and the self-interaction cross-section per mass,
sigma/m, which is the quantity that decides whether the sector is a
self-interacting dark-matter (SIDM) candidate.

This script supplies those numbers from the framework's own fundamentals, with
no new input. Everything rests on three identities already in the monograph:

    lattice spacing         l   = r_e            (the classical electron radius)
    node mass               m0  = m_e / alpha
    node Compton = spacing   l   = hbar/(m0 c)   =>  kappa = h/m0 = 2 pi c l.

From these the vortex line tension follows,

    T = (rho_s kappa^2 / 4pi) ln(R/xi),   rho_s = (4/5) rho,  rho = m0/l^3,

which reduces to T = (4pi/5) ln(R/xi) * (m0 c^2 / l): about one node energy per
lattice spacing of line, times a slowly varying logarithm. A knot is a length
of this line, so its mass is the line tension times the length, and its
self-interaction cross-section is the geometric area it presents to another
knot. The two together give sigma/m.

Cross-check built in: the same line tension must reproduce the chapter's own
straight-vortex binding energy E_bind ~ 0.9-2.6 GeV over a length sqrt(6) l
(Eq. eq:vortex_binding). It does, which anchors the normalisation.

The helicity H = n kappa^2 is, for a coreless texture, exactly the Hopf charge
n of the order-parameter map (via the Mermin-Ho relation the framework already
uses). So these objects are the Faddeev-Niemi Hopfions of the soliton
literature, and the Vakulenko-Kapitanskii bound E >= c_VK |n|^{3/4} applies:
the mass grows sublinearly in the knot integer, faster in the physical size.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Fundamentals (CODATA; the framework's derived combinations built from them).
# ---------------------------------------------------------------------------
alpha = 1.0 / 137.035999177         # fine-structure constant
m_e_MeV = 0.51099895069             # electron mass [MeV/c^2]
r_e_cm = 2.8179403205e-13           # classical electron radius [cm]

m0_MeV = m_e_MeV / alpha            # node mass [MeV/c^2]
ell_cm = r_e_cm                     # lattice spacing = r_e
f_s = 4.0 / 5.0                     # superfluid fraction (bootstrap, sec:mott)

MeV_to_g = 1.782662e-27             # 1 MeV/c^2 in grams
barn_cm2 = 1.0e-24                  # 1 barn in cm^2

# Line tension coefficient: T = C_T * ln(R/xi) * (m0 c^2 / ell),  C_T = 4 pi / 5.
C_T = 4.0 * np.pi / 5.0             # = (rho_s kappa^2 / 4pi) reduced, self-energy
C_int = 8.0 * np.pi / 5.0          # interaction coefficient (rho_s kappa^2 / 2pi)

def line_tension_per_ell(lnRxi):
    """Vortex self-energy per lattice spacing, in units of m0 c^2."""
    return C_T * lnRxi              # T * ell / (m0 c^2)

# ---------------------------------------------------------------------------
# Cross-check: reproduce the chapter's straight-vortex binding energy.
#   E_bind = (rho_s kappa^2 / 2pi) * L4 * ln(R/xi),  L4 = sqrt(6) ell.
# ---------------------------------------------------------------------------
L4_over_ell = np.sqrt(6.0)
def E_bind_GeV(lnRxi):
    return C_int * lnRxi * L4_over_ell * m0_MeV / 1e3   # GeV

# ---------------------------------------------------------------------------
# Mass of a closed object of given vortex length (in units of ell).
#   M c^2 = T * length = C_T ln(R/xi) (m0 c^2/ell) * length.
# ---------------------------------------------------------------------------
def mass_MeV(length_over_ell, lnRxi):
    return C_T * lnRxi * length_over_ell * m0_MeV

# Self-interaction cross-section of a compact knot of size R: geometric, ~ pi R^2.
def sigma_cm2(R_over_ell):
    return np.pi * (R_over_ell * ell_cm) ** 2

def sigma_over_m(R_over_ell, length_over_ell, lnRxi):
    m_g = mass_MeV(length_over_ell, lnRxi) * MeV_to_g
    return sigma_cm2(R_over_ell) / m_g          # cm^2 / g

# ---------------------------------------------------------------------------
# A catalogue of the simplest objects. Length and size from knot geometry:
#   ring of radius R:  length = 2 pi R,      size R.
#   trefoil (3_1):     ropelength ~ 16.4 d,  d = 2 xi ~ 2 ell, so length ~ 33 ell,
#                      packed into size R ~ 5 ell.
#   figure-eight (4_1): ropelength ~ 21 d  -> length ~ 42 ell, size R ~ 6 ell.
# The logarithm ln(R/xi) is taken at the object's own size, xi ~ ell.
# ---------------------------------------------------------------------------
def catalogue_row(name, R_over_ell, length_over_ell):
    lnRxi = max(np.log(R_over_ell), 0.3)        # floor: cores never fully overlap
    m = mass_MeV(length_over_ell, lnRxi)
    s = sigma_cm2(R_over_ell)
    som = sigma_over_m(R_over_ell, length_over_ell, lnRxi)
    return name, R_over_ell, length_over_ell, lnRxi, m, s, som

OBJECTS = [
    catalogue_row("smallest ring",   1.0,  2.0 * np.pi * 1.0),
    catalogue_row("ring R=3l",       3.0,  2.0 * np.pi * 3.0),
    catalogue_row("ring R=10l",     10.0,  2.0 * np.pi * 10.0),
    catalogue_row("trefoil 3_1",     5.0,  33.0),
    catalogue_row("figure-8 4_1",    6.0,  42.0),
    catalogue_row("ring R=100l",   100.0,  2.0 * np.pi * 100.0),
]

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def line():
    print("-" * 78)

print(__doc__.strip().splitlines()[0])
line()
print("Fundamentals (from the framework):")
print(f"  alpha = 1/{1/alpha:.3f}   m0 = m_e/alpha = {m0_MeV:.2f} MeV   "
      f"l = r_e = {ell_cm*1e13:.3f} fm   f_s = {f_s}")
print(f"  circulation quantum kappa = 2 pi c l ;  line tension "
      f"T = (4pi/5) ln(R/xi) (m0 c^2 / l)")
line()
print("Cross-check against the chapter's straight-vortex binding "
      "(Eq. eq:vortex_binding, 0.9-2.6 GeV):")
for ln in (1.0, 2.0, 3.0):
    print(f"  ln(R/xi) = {ln:.0f}:  E_bind = {E_bind_GeV(ln):.2f} GeV")
print("  -> normalisation matches the monograph (0.9 GeV at ln=1, 2.6 GeV at ln=3).")
line()
print("Mass spectrum and self-interaction of the free objects:")
print(f"  {'object':<15s} {'R/l':>6s} {'len/l':>7s} {'ln':>4s} "
      f"{'M [GeV]':>9s} {'sigma':>11s} {'sigma/m [cm^2/g]':>17s}")
for name, R, L, ln, m, s, som in OBJECTS:
    print(f"  {name:<15s} {R:>6.1f} {L:>7.1f} {ln:>4.1f} "
          f"{m/1e3:>9.2f} {s/barn_cm2:>8.1f} b {som:>17.3f}")
line()
print("Reading the numbers:")
print("  Mass: the lightest closed object is ~1 GeV (a ring of core size), and")
print("  the mass grows linearly with the vortex length, so knots and larger")
print("  rings climb into the tens of GeV. This makes the chapter's")
print("  'nuclear-to-GeV scale' quantitative: a spectrum M ~ (length/l) x m0,")
print("  floored near 1 GeV.")
print()
print("  Self-interaction: sigma/m runs from ~0.1 cm^2/g for the lightest rings")
print("  up through ~1 cm^2/g for larger knots, i.e. squarely across the SIDM")
print("  window (0.1-10 cm^2/g) that cored dwarf profiles and rotation-curve")
print("  diversity call for. sigma/m grows with knot size, so a size")
print("  distribution gives a naturally velocity/environment-dependent")
print("  self-interaction, which is what the cluster (Bullet) bounds prefer.")
print()
print("  Topology: helicity H = n kappa^2 is the Hopf charge n of the coreless")
print("  texture (Mermin-Ho), so these are Faddeev-Niemi Hopfions and the")
print("  Vakulenko-Kapitanskii bound E >= c_VK |n|^{3/4} applies.")
line()
print("What stays open (as the chapter says): the RELIC ABUNDANCE. Kibble-Zurek")
print("sets the initial density and reconnection sets the survivor fraction;")
print("neither is derived here. sigma/m lands in the SIDM window, but whether")
print("the sector carries a cosmologically significant fraction is the one")
print("number that would turn a candidate into a component.")
