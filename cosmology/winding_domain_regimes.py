#!/usr/bin/env python3
"""
winding_domain_regimes.py -- what sets the winding-domain size at the
deconfinement crossover, now that bubbles are gone.

The baryon asymmetry in this framework is a per-domain choice: each causally
coherent region of the stacking order picks a winding chirality, the
propagation chirality biases the pick by delta = 0.34, and the observed
n_B/s = 8.7e-11 then demands a domain size

    xi_req = (delta / (n_B/s * s))^(1/3) = 663 fm = 235 ell        (at T_c)

This script asks which physical mechanism can deliver that number, by pricing
every candidate regime with the framework's own scales.  Four regimes:

  A. Standard Kibble-Zurek at a critical point.  Requires critical slowing
     down.  The measured transition is a crossover (potts_fcc_field_endpoint),
     so the correlation length and relaxation time saturate at xi_max, tau_max
     instead of diverging, and the first question is whether the quench is
     even fast enough to fall out of equilibrium at all.

  B. Adiabatic passage.  If the system equilibrates through the crossover, the
     domain structure at the transition is the equilibrium one, xi ~ xi_max,
     of order 1/T_c.  Overproduces catastrophically.

  C. Free coarsening after the crossover.  Domain walls under tension in a
     low-viscosity medium accelerate toward c; the network coarsens toward the
     causal bound.  Underproduces catastrophically by the time the chirality
     window closes at T_geom.

  D. Coarsening arrested by Peierls pinning.  A stacking-chirality wall is a
     crystallographic object and sees the lattice corrugation.  Curvature
     pressure falls as sigma/L, so growth self-arrests at L_pin where the
     pressure equals the pinning stress.  The arrest length is a pure lattice
     number, and the question becomes what wall width the required 235 ell
     implies through the Peierls-Nabarro exponential.

The point of the exercise is that A, B and C fail by four to sixty orders in
three different directions, which forces D or something like it: the domain
size is not a relic of the quench at all, but a static property of the lattice.
The baryon-to-photon ratio then stops being a cosmological accident and
becomes a crystallographic one.

All numbers per the corpus: T_c = 156.1 MeV, T_geom = 28.6 MeV, delta = 0.34,
ell = 2.818 fm, g* = 61.75 above the crossover.
"""
import numpy as np

hbar_c = 197.3269804          # MeV fm
hbar_s = 6.582119569e-22      # MeV s
c_fm_s = 2.99792458e23        # fm / s
M_Pl = 1.22089e22             # MeV
T_c, T_geom = 156.1, 28.6     # MeV
g_star = 61.75
ell = 2.8179403262            # fm
delta = 0.34
nB_s_obs = 8.7e-11

s_ent = (2 * np.pi**2 / 45) * g_star * T_c**3 / hbar_c**3     # fm^-3
xi_req = (delta / (nB_s_obs * s_ent))**(1 / 3)                # fm


def hubble(T):
    return np.sqrt(8 * np.pi**3 * g_star / 90) * T**2 / M_Pl  # MeV


def t_of(T):
    return hbar_s / (2 * hubble(T))                           # s (rad era)


if __name__ == "__main__":
    print("What the observed asymmetry demands")
    print(f"  s(T_c) = {s_ent:.2f} fm^-3   ->   xi_req = {xi_req:.0f} fm "
          f"= {xi_req/ell:.0f} ell\n")

    H = hubble(T_c)
    tau_Q = hbar_s / H
    tau_0 = hbar_s / T_c

    # ---------------- A: standard KZ, and whether it even applies ----------
    print("A. Standard Kibble-Zurek (needs a critical point)")
    for nu, z, tag in ((0.5, 2.0, "mean field"), (0.63, 2.0, "3D Ising")):
        e = nu / (1 + nu * z)
        xi = ell * (tau_Q / tau_0)**e
        print(f"   {tag:10s}: xi_KZ = {xi:9.3g} fm   "
              f"overshoots xi_req by {xi/xi_req:8.1f}x  "
              f"-> n_B/s low by {(xi/xi_req)**3:.1e}")
    # crossover rounding: does the system fall out of equilibrium at all?
    # the field h ~ 0.3 (an order past the endpoint h_c ~ 0.03) rounds the
    # transition at eps_h ~ h^(1/(beta delta)); take the QCD crossover width
    # ~15 MeV as the empirical rounding.
    eps_h = 15.0 / T_c
    xi_max = ell * eps_h**(-0.63)
    tau_max = tau_0 * (xi_max / ell)**2
    rate = 0.63 / eps_h * H / T_c * T_c        # |dln xi/dt| peak ~ nu/eps_h * H
    print(f"\n   Crossover rounding: eps_h ~ {eps_h:.2f}, xi_max ~ "
          f"{xi_max:.1f} fm, tau_max ~ {tau_max:.1e} s")
    print(f"   adiabaticity parameter tau_max * (nu/eps_h) * H = "
          f"{tau_max*0.63/eps_h*H/hbar_s:.1e}")
    print("   << 1 by seventeen orders: the system NEVER falls out of"
          " equilibrium.\n   Standard KZ does not apply at this crossover"
          " at all.\n")

    # ---------------- B: adiabatic passage ---------------------------------
    print("B. Adiabatic passage (equilibrium domains at the crossover)")
    nB_B = delta / xi_max**3 / s_ent
    print(f"   xi ~ xi_max = {xi_max:.1f} fm  ->  n_B/s = {nB_B:.1e}"
          f"   overproduces by {nB_B/nB_s_obs:.1e}\n")

    # ---------------- C: free coarsening to the window edge ----------------
    print("C. Free coarsening until the chirality window closes at T_geom")
    t_window = t_of(T_geom) - t_of(T_c)
    L_causal = c_fm_s * t_window
    nB_C = delta / L_causal**3 / s_ent
    print(f"   window duration = {t_window*1e6:.0f} us,  causal bound "
          f"L = {L_causal:.2e} fm")
    print(f"   n_B/s = {nB_C:.1e}   underproduces by {nB_s_obs/nB_C:.1e}\n")

    # ---------------- D: Peierls-arrested coarsening ------------------------
    print("D. Coarsening arrested by lattice pinning")
    print("   Growth stops when curvature pressure sigma/L falls below the")
    print("   Peierls pinning stress gamma_P/ell:  L_pin = ell (sigma/gamma_P).")
    ratio = xi_req / ell
    print(f"   Landing L_pin = xi_req needs   gamma_P/sigma = 1/{ratio:.0f}"
          f" = {1/ratio:.1e}")
    w = ell * np.log(ratio) / (2 * np.pi)
    print(f"   Through the Peierls-Nabarro form gamma_P/sigma ~ exp(-2 pi w/ell),")
    print(f"   that pinning ratio corresponds to a wall width")
    print(f"       w = ell ln({ratio:.0f})/(2 pi) = {w/ell:.2f} ell")
    print("   against the corpus's dislocation core width w = 0.45 ell:")
    print(f"   a factor {w/(0.45*ell):.1f}, i.e. the right regime, not yet the"
          " right number.\n")

    # timing checks for D
    print("   Timing checks for regime D:")
    # walls must be removable by the bias before wall domination / today.
    # bias pressure is a volume energy difference ~ delta * T^4-scale acting
    # on a pinned wall; it exceeds curvature pressure once L > L_pin, so the
    # network is swept by the bias unless the pinning also beats the bias.
    print("   - the same pinning must NOT beat the chirality bias pressure,")
    print("     or the wall network survives the window and over-closes the")
    print("     universe; the bias must depin and sweep the minority domains")
    print("     between T_c and T_geom.  Sweep time at speed v across L_pin:")
    for v in (1e-3, 1e-6):
        t_sweep = xi_req / (v * c_fm_s)
        print(f"       v = {v:g} c  ->  {t_sweep:.1e} s  "
              f"(window is {t_window:.1e} s)")
    print("     even v = 1e-6 c clears the network with twelve orders to"
          " spare,\n     so the ordering gamma_P: curvature < gamma_P < bias"
          " is the whole\n     requirement, and it spans many orders of"
          " admissible gamma_P.")
