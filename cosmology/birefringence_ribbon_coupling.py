#!/usr/bin/env python3
"""
birefringence_ribbon_coupling.py
================================
The accumulating birefringence mechanism, after the sign audit
(Confrontation chapter, "Cosmic birefringence").

The sign audit (birefringence_sign_audit.py) shows the twist handedness of
the walls alternates along the line of sight, so any rotation ODD in the wall
twist telescopes. This script does two things:

1. Structure theorem. Split the per-wall rotation into parts odd and even in
   the chirality jump s_w = +/-1. Over N strictly alternating walls the odd
   part telescopes to an endpoint bounded by theta_ch; only the EVEN part
   accumulates as N * (even part). An even part exists only if mirror
   symmetry is broken by something other than the wall itself.

2. The even coupling from the fault ribbons. The lattice supplies the broken
   mirror: the global A -> B -> C stacking handedness, uniform across every
   domain. Each wall screw dissociates on {111} into Shockley partials
   bounding an intrinsic stacking-fault ribbon whose type is fixed by that
   global stacking, for EITHER sign of the parent Burgers vector (standard
   FCC metallurgy). The ribbons are the even, globally handed material. Their
   areal fraction per crossing sets the accumulating rotation, and it carries
   the theta_ch scaling automatically because the wall's dislocation line
   density is proportional to theta_ch.

Result: beta = c1 * theta_ch * N with c1 = (8 w / b) * g, where w is the PN
core half-width, b the screw Burgers magnitude, and g the polarisation
rotation per fault-core traversal (the one remaining O(1) unknown, a
Maurel-Lund scattering problem off a KNOWN core profile). The sign is the
global stacking sense, single across all walls, tying beta to the baryon
asymmetry.

Units: theta_ch dimensionless; lengths in the lattice spacing ell.
"""
import numpy as np

alpha = 1 / 137.035999177
theta_ch = alpha**2 / (2 * np.pi)          # rad
theta_ch_deg = np.degrees(theta_ch)


def structure_theorem(N=940, seed=1, n_trials=3000):
    """Odd part telescopes, even part accumulates. Demonstrate numerically."""
    rng = np.random.default_rng(seed)
    s = np.array([(-1) ** i for i in range(N)])            # alternating jump sign
    a, b = 1.0, 1.0                                        # arbitrary odd, even coeffs
    # per-wall rotation = a * s_w (odd) + b (even), in units theta_ch
    odd_sum = (a * s).sum()
    even_sum = (b * np.ones(N)).sum()
    print("structure theorem (rotation = a*s_w + b, units theta_ch):")
    print(f"  odd part  sum over N={N}: {odd_sum:+.0f}  (telescopes; |.| <= 1)")
    print(f"  even part sum over N={N}: {even_sum:+.0f}  (accumulates as N*b)")
    print("  => only the even part survives; it needs a mirror broken off-wall.\n")


def ribbon_coupling(N=940):
    """The even coupling from globally handed Shockley-partial fault ribbons."""
    b_B = 1.0                                  # screw Burgers magnitude [ell]
    D_grid = b_B / (2 * theta_ch)              # twist-grid line spacing [ell]
    d_p = 1 / np.sqrt(3)                       # partial hop distance d = ell/sqrt3
    w = (np.pi / 4) * d_p                      # PN core half-width, w/d = pi/4

    # chiral-material areal fraction per crossing: two line families, each of
    # line density 1/D_grid, each line dressed by a core/ribbon of width ~2w.
    f_core = 2 * (2 * w) / D_grid
    c1_prefactor = 8 * w / b_B                 # so that f_core = c1_prefactor * theta_ch

    print("ribbon coupling (even, globally handed):")
    print(f"  twist-grid spacing D = b/(2 theta_ch) = {D_grid:.3e} ell (dilute grid)")
    print(f"  PN core half-width  w = (pi/4) ell/sqrt3 = {w:.4f} ell")
    print(f"  chiral areal fraction per wall f = {f_core:.3e} = {f_core/theta_ch:.2f} theta_ch")
    print(f"  => c1 = (8 w / b) * g = {c1_prefactor:.2f} * g\n")

    print("  per-core rotation g        beta = c1 theta_ch N   implied by data")
    for g, lab in [(1/np.pi, "g = 1/pi"), (0.5, "g = 1/2"), (1.0, "g = 1")]:
        c1 = c1_prefactor * g
        beta = c1 * theta_ch_deg * N
        print(f"    {lab:9s}: c1 = {c1:.2f}, beta = {beta:.2f} deg")
    # invert: what g does the measured band want?
    for name, bmeas in [("Eskilt-Komatsu 0.342", 0.342),
                        ("Planck PR4 0.30", 0.30),
                        ("ACT DR6 0.215", 0.215)]:
        c1_need = bmeas / (theta_ch_deg * N)
        g_need = c1_need / c1_prefactor
        print(f"  {name:22s} deg -> c1 = {c1_need:.2f}, g = {g_need:.3f}")
    print("\n  theta_ch scaling is automatic; sign = global A->B->C stacking;")
    print("  the one open O(1) is g, a Maurel-Lund scattering integral off a known core.")


if __name__ == "__main__":
    structure_theorem()
    ribbon_coupling()
