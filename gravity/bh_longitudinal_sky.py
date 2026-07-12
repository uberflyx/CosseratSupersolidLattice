"""
The longitudinal sky: what the compression channel carries from black holes.

The monograph establishes (sec:interior_knot_chemistry) that the event
horizon is a shear-sector horizon only: first sound rides the condensate
bulk stiffness, its impedance is continuous across the melt boundary, and
interior compression waves transmit outward with an order-unity
coefficient. This script prices the consequences:

  1. The channel-horizon hierarchy: trapping radius scales as (c/v)^2,
     so each propagation channel sees a different sized hole.
  2. The steady hum of a fed hole: standing-tangle stock over a
     residence time. The stock is 1e-19 of Mc^2 (tex and the
     bh_vortex_reservoir.py computation agree at 5e-19; the earlier
     1e-12 in that script's summary string was a stale figure, since
     corrected).
  3. The merger burst: an upper bound from the Feynman array's share
     of the extractable rotational energy released at spin
     reconfiguration.
  4. The preview timing: second sound at 3.65c outruns light, so the
     longitudinal image of any event arrives 0.726 D/c early; first
     sound at ~1e20 c makes the whole observable universe one
     millisecond-latency 'now'.
  5. The spectroscopy: the interior's own lines (annihilation 2E1,
     bead-dimer bands, Kelvin hiss), quoted from the chemistry section.
"""

import numpy as np

# ------------------------------------------------------------ constants
c = 2.99792458e8
G = 6.67430e-11
Msun = 1.98892e30
yr = 3.156e7
pc = 3.0857e16

v2 = 3.65 * c                    # second sound group velocity (derived)
vp = 1e20 * c                    # first sound / pilot-wave speed (order)

# ---------------------------------------------- 1. channel-horizon sizes
print("=== 1. one hole, three horizons (10 Msun) ===")
M = 10 * Msun
rs = 2 * G * M / c**2
for name, v in [("shear (light, matter)", c), ("second sound", v2),
                ("first sound", vp)]:
    rh = rs * (c / v) ** 2
    note = "" if rh > 2.8e-15 else "  < lattice spacing: NO horizon"
    print(f"  {name:22s}: r_h = {rh:9.3e} m = r_s/{rs/rh:9.3g}{note}")

# ------------------------------------------------- 2. the steady-state hum
print("\n=== 2. the hum of a fed hole (10 Msun at Eddington) ===")
Mc2 = M * c**2
t_salpeter = 4.5e7 * yr          # Eddington e-folding (residence) time
L_edd = 1.26e31 * 10             # W
share = 1e-19                    # resolved stock: tex == reservoir-script computation
L_hum = share * Mc2 / t_salpeter
F_1kpc = L_hum / (4 * np.pi * (1e3 * pc) ** 2)
print(f"  tangle share {share:.0e}: L_hum ~ {L_hum:8.2e} W "
      f"= {L_hum/L_edd:.1e} L_Edd; flux at 1 kpc {F_1kpc:.1e} W/m^2")
print("  -> a fed hole hums faintly in the longitudinal channel; only")
print("     the monopole transducer could hear it.")

# ------------------------------------------------------ 3. merger burst
print("\n=== 3. merger burst upper bound (GW150914-like) ===")
Mf, af = 62 * Msun, 0.67
E_rot = Mf * c**2 * (1 - np.sqrt(0.5 * (1 + np.sqrt(1 - af**2))))
array_share = 6e-4               # array's share of extractable rot. energy
E_burst = array_share * E_rot    # released at spin reconfiguration (upper)
E_gw = 3 * Msun * c**2
print(f"  extractable rotational energy: {E_rot:.2e} J "
      f"({E_rot/(Mf*c**2)*100:.1f}% of Mf c^2)")
print(f"  array share {array_share:.0e} -> E_burst <= {E_burst:.1e} J "
      f"= {E_burst/E_gw:.1e} of the GW energy")

# --------------------------------------------------- 4. preview timing
print("\n=== 4. the preview channel ===")
for name, D in [("GW150914 (410 Mpc)", 410e6 * pc),
                ("GW170817 ( 40 Mpc)", 40e6 * pc),
                ("galactic centre (8 kpc)", 8e3 * pc)]:
    lead2 = (D / c) * (1 - 1 / 3.65)
    t1 = D / vp
    print(f"  {name:24s}: second-sound lead = {lead2/yr:10.3e} yr; "
          f"first-sound transit = {t1*1e3:7.3f} ms")
print("  -> today's second-sound sky previews events whose light arrives")
print("     up to 0.73 D/c in the future; the first-sound sky is live.")

# ------------------------------------------------------ 5. the colours
print("\n=== 5. what the hum sounds like (from interior_knot_chemistry) ===")
print("  annihilation line: 2 E1 = 2.21 GeV (bead-antibead)")
print("  dimer bands: bond 419 MeV at contact, vibrational quantum 49 MeV,")
print("               ~9 levels: a sparse nuclear-style spectrum")
print("  Kelvin hiss: soft quanta ~0.4 MeV, quadratically dispersive")
print("  silence map: quiescent Schwarzschild interiors empty in ~Myr;")
print("               the longitudinal sky maps FED holes and spin events")
