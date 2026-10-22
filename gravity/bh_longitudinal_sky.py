"""
The longitudinal sky: what the compression channel carries from black holes.

The monograph's black-holes chapter leaves open whether the condensate's
waves feel the missing sector's lapse at all.  If they do not, the hole has
no horizon on the longitudinal channel and its interior would be audible on
that channel; if they do, they freeze at r_s like every other wave.  This
script prices the first case, together with an estimate the chapter sets
aside.  Every interior quantity rests on the chapter's conditional premise
(a shear-free interior in which the condensate persists), and the interior
superfluid fraction is not derived: a shear-free interior sits between
f_s = 4/5 and 1, so interior quantities are bracketed between 4/5 and the
adopted exterior value f_s = 0.99932.

  1. The channel-horizon hierarchy of the "river" picture, in which a
     wave of speed v is trapped where an inflow reaches v, giving a
     trapping radius r_s (c/v)^2.  The framework's vacuum does not flow
     (the crystal is at rest in the static frame), so these radii are the
     river-picture numbers the chapter sets aside, kept for comparison.
  2. The standing tangle of a fed hole: the Eddington-fed Vinen steady
     state of bh_vortex_reservoir.py, with the injected line power capped
     by the accreted energy, P <= (1 - eta) Mdot c^2 (one loop of radius
     3 ell per accreted pair would cost 30 to 80 times the pair's rest
     energy).  Priced over the interior fractions 4/5 to 0.99932 and the
     line-tension logarithm cut at ten or at the inter-vortex spacing
     (about 22), the stock is 0.8 to 1.4e-19 of Mc^2 for a ten-solar-mass
     hole.  This stock bounds what the interior's spectroscopy could carry
     out; no luminosity is estimated here.
  3. The merger burst, on the rigid-array model: the remnant's Feynman
     array rotating rigidly at Omega_H = a* c / (2 r_+) holds
     J_array = (2/5) M_s r_+^2 Omega_H, with
     M_s = f_s (m0/ell^3)(4/3) pi r_+^3 the superfluid mass inside the
     horizon, and rotational energy (1/2) J_array Omega_H = M_s a*^2 c^2/20.
     For a GW150914-class remnant (62 Msun, a* = 0.67) that scale is
     1.5 to 1.9e47 J, a modelling premise rather than a derived bound.
  4. The preview timing: second sound at v_2 = sqrt((3-4N^2)/(1-N^2)) c
     = 1.592 c (N^2 = 1/pi) outruns light, so the longitudinal image of an
     event at distance D arrives (1 - c/v_2) D/c = 0.37 D/c early; first
     sound at ~1e20 c crosses a gigaparsec in about a millisecond.
  5. The spectroscopy: the interior's own lines (annihilation 2E1,
     bead-dimer bands, Kelvin hiss), quoted from the interior-chemistry
     section (interior_knot_chemistry.py computes them).
"""

import numpy as np

# ------------------------------------------------------------ constants
alpha = 1.0 / 137.035999177
m_e   = 9.1093837015e-31             # electron mass [kg]
m0    = m_e / alpha                  # node mass [kg]
ell   = 2.8179403205e-15             # lattice spacing r_e [m]
c     = 2.99792458e8
G     = 6.67430e-11
Msun  = 1.98892e30                   # [kg]
m_p   = 1.67262192e-27               # proton mass [kg]
yr    = 3.156e7                      # [s]
pc    = 3.0857e16                    # [m]
kappa = 2.0 * np.pi * ell * c        # circulation quantum h/m0 [m^2/s]
rho_n = m0 / ell**3                  # node mass density [kg/m^3]

F_S_ADOPTED = 0.99932                # exterior working superfluid fraction
F_S_INT_LOW = 4 / 5                  # interior lower bound (not derived)
FS_BRACKET  = (F_S_INT_LOW, F_S_ADOPTED)

N2 = 1.0 / np.pi                     # Cosserat coupling number
v2 = np.sqrt((3.0 - 4.0 * N2) / (1.0 - N2)) * c   # second sound, 1.592 c
vp = 1e20 * c                        # first sound / pilot wave (order only)

# ---------------------------------------------- 1. channel-horizon sizes
print("=== 1. river-picture trapping radii (10 Msun); set aside by the chapter ===")
M = 10 * Msun
rs = 2 * G * M / c**2
for name, v in [("shear (light, matter)", c), ("second sound", v2),
                ("first sound", vp)]:
    rh = rs * (c / v) ** 2
    note = "" if rh > ell else "  < lattice spacing: NO horizon"
    print(f"  {name:22s}: r_h = {rh:9.3e} m = r_s/{rs/rh:9.3g}{note}")

# ----------------------------------------- 2. the standing tangle (stock)
# Vinen steady state under Eddington feed, the model of bh_vortex_reservoir.py,
# with the feed capped by energy.  A loop of radius 3 ell per accreted charge
# pair would cost 6 pi^2 f_s ln m0 c^2, 30 to 80 times the pair's rest energy,
# so the injected line power is capped at P = (1 - eta) Mdot c^2.  Then
#   dL/dt = P/(V T) - chi kappa L^2  ->  L* = sqrt(P / (V T chi kappa)),
#   E = T L* V = sqrt(P T V / (chi kappa)),
# with T = pi f_s ln m0 c^2 / ell the line tension and chi ~ 0.1 the
# superfluid-turbulence Vinen constant.  MODEL INPUTS: eta, chi, the cap.
print("\n=== 2. the standing tangle of a fed hole (10 Msun at Eddington) ===")
eta      = 0.1                       # radiative efficiency of accretion
L_Edd    = 1.26e31 * (M / Msun)      # Eddington luminosity [W]
Mdot     = L_Edd / (eta * c**2)      # accretion rate [kg/s]
P_cap    = (1 - eta) * Mdot * c**2   # injected line power at the cap [W]
V_int    = (4 / 3) * np.pi * rs**3   # interior volume [m^3]
chi      = 0.1
T_line   = lambda fs, ln: np.pi * fs * ln * m0 * c**2 / ell   # [J/m]


def stock(fs, ln):
    """Capped Vinen stock as a share of Mc^2, with its line density."""
    L_star = np.sqrt(P_cap / (V_int * T_line(fs, ln) * chi * kappa))
    return T_line(fs, ln) * L_star * V_int / (M * c**2), L_star


def ln_self_consistent(fs):
    """Iterate the tension logarithm to its value at the inter-vortex spacing."""
    ln = 20.0
    for _ in range(50):
        _, L_star = stock(fs, ln)
        ln = np.log(1 / np.sqrt(L_star) / ell)
    return ln


shares = []
for fs in FS_BRACKET:
    for ln in (10.0, ln_self_consistent(fs)):
        s, L_star = stock(fs, ln)
        shares.append(s)
        print(f"  f_s = {fs:7.5f}, log cut {ln:5.2f}: L* = {L_star:.2e} m^-2,"
              f" tangle share of Mc^2 = {s:.2e}")
print(f"  -> standing stock {min(shares):.1e} to {max(shares):.1e} of Mc^2,"
      f" about 1e-19,")
print("     the bound on what the interior's spectroscopy could carry out.")

# ------------------------------------------------------ 3. merger burst
print("\n=== 3. merger-burst scale on the rigid-array model (GW150914-class) ===")
Mf, af = 62 * Msun, 0.67
r_plus  = (1 + np.sqrt(1 - af**2)) * G * Mf / c**2      # horizon radius [m]
Omega_H = af * c / (2 * r_plus)                         # horizon angular velocity
E_rot   = Mf * c**2 * (1 - np.sqrt(r_plus * c**2 / (2 * G * Mf)))  # above M_irr
E_gw    = 3.0 * Msun * c**2                             # radiated in GWs (input)
M_s     = lambda fs: fs * rho_n * (4 / 3) * np.pi * r_plus**3
E_array = lambda fs: M_s(fs) * af**2 * c**2 / 20        # (1/2) J_array Omega_H
# consistency: (1/2)(2/5) M_s r_+^2 Omega_H^2 must equal M_s a*^2 c^2 / 20
assert np.isclose(0.5 * 0.4 * M_s(1.0) * r_plus**2 * Omega_H**2, E_array(1.0))
print(f"  r_+ = {r_plus/1e3:.1f} km, Omega_H = {Omega_H:.0f} rad/s")
print(f"  extractable rotational energy: {E_rot:.2e} J "
      f"({E_rot/(Mf*c**2)*100:.1f}% of Mf c^2)")
for fs in FS_BRACKET:
    print(f"  f_s = {fs:7.5f}: M_s/M = {M_s(fs)/Mf:.2f}, "
          f"E = M_s a*^2 c^2/20 = {E_array(fs):.2e} J = "
          f"{E_array(fs)/E_rot:.2f} of E_rot = {E_array(fs)/E_gw:.2f} of E_GW")
print(f"  (E_GW = 3.0 Msun c^2 = {E_gw:.2e} J)")

# --------------------------------------------------- 4. preview timing
print("\n=== 4. the preview channel ===")
lead_frac = 1 - c / v2
print(f"  v_2 = sqrt((3-4N^2)/(1-N^2)) c = {v2/c:.4f} c  ->  lead = "
      f"(1 - c/v_2) D/c = {lead_frac:.4f} D/c")
for name, D in [("GW150914 (410 Mpc)", 410e6 * pc),
                ("GW170817 ( 40 Mpc)", 40e6 * pc),
                ("galactic centre (8 kpc)", 8e3 * pc)]:
    print(f"  {name:24s}: light time = {D/c/yr:10.3e} yr; "
          f"second-sound lead = {lead_frac*D/c/yr:10.3e} yr; "
          f"first-sound transit = {D/vp*1e3:7.3f} ms")
print(f"  first sound at 1e20 c crosses 1 Gpc in {1e9*pc/vp*1e3:.2f} ms")

# ------------------------------------------------------ 5. the colours
print("\n=== 5. what the interior sounds like (from interior_knot_chemistry) ===")
print("  annihilation line: 2 E1 = 2.21 GeV (bead-antibead, f_s = 4/5)")
print("  dimer bands: bond 419 MeV at contact, vibrational quantum 49 MeV,")
print("               ~9 levels: a sparse nuclear-style spectrum")
print("  Kelvin hiss: soft quanta ~0.4 MeV, quadratically dispersive")
print("  silence map: quiescent Schwarzschild interiors empty in ~Myr;")
print("               the longitudinal sky maps FED holes and spin events")
