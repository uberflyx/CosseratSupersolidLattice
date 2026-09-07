#!/usr/bin/env python3
"""Transverse energy partition of the anti-plane Cosserat medium.

Purpose
-------
Establish, symbolically and then numerically, which combination of the Cosserat
moduli is the fraction of a transverse wave's elastic energy that the
microrotation channel carries away, and confirm that this combination is the
coupling number

    N^2 = kappa_c / [2 (mu + kappa_c)] .

Two dimensionless ratios can both be called "the rotational fraction", and they
differ by 47 per cent at the operating point.  They are distinguished only by
the reference energy in the denominator:

    N^2 = (relaxation energy) / (energy of the same wave with phi frozen)
    g_c = (relaxation energy) / (energy of the same wave with phi relaxed)
        = N^2 / (1 - N^2) .

The script shows that the first is the coupling number and the second is not.

Method
------
The anti-plane problem carries a transverse displacement u_3(x_1) and a
microrotation phi_2(x_1).  Its strain energy density is

    W = (mu/2) u'^2 + (kappa_c/2) [ (u' + phi)^2 + phi^2 ] + (gamma/2) phi'^2

where u' = du_3/dx_1 is the shear strain [-], phi = phi_2 is the microrotation
angle [rad], mu is the Cosserat shear modulus [Pa], kappa_c is the coupling
modulus that resists a mismatch between a node's own spin and the twist of the
medium around it [Pa], and gamma is the curvature modulus [Pa m^2].  The two
squared terms inside the kappa_c bracket are the two independent antisymmetric
strain components that involve phi_2, which is why the local restoring torque
in the field equations carries 2 kappa_c rather than kappa_c.

Varying W reproduces the printed field equations, so the energy and the
equations agree before any partition is attempted.

Author: M. A. Cox
"""

from __future__ import annotations

import numpy as np
import sympy as sp

# ----------------------------------------------------------------------------
# Symbols.  All positive: moduli, wavenumber and amplitude are physical.
# ----------------------------------------------------------------------------
mu, kc, gam, k, U = sp.symbols("mu kappa_c gamma k U", positive=True)


def strain_energy_density():
    """Return W(u', phi, phi') and the two Euler-Lagrange field equations.

    The field equations are returned in the form (LHS_u, LHS_phi) where each is
    the elastic force (or torque) density that balances the inertia term.  They
    should read

        rho u_3,tt   = (mu + kappa_c) u_3,11 + kappa_c phi_2,1
        J   phi_2,tt = gamma phi_2,11 - kappa_c u_3,1 - 2 kappa_c phi_2 .
    """
    x = sp.symbols("x", real=True)
    u = sp.Function("u")(x)
    ph = sp.Function("phi")(x)

    W = (
        mu / 2 * sp.diff(u, x) ** 2
        + kc / 2 * ((sp.diff(u, x) + ph) ** 2 + ph ** 2)
        + gam / 2 * sp.diff(ph, x) ** 2
    )

    # Euler-Lagrange: force density = d/dx (dW/du') - dW/du, with a sign so that
    # the result is the right-hand side of the equation of motion.
    force_u = sp.expand(sp.diff(sp.diff(W, sp.diff(u, x)), x))
    force_phi = sp.expand(sp.diff(sp.diff(W, sp.diff(ph, x)), x) - sp.diff(W, ph))
    return W, force_u, force_phi


def cycle_average(A, B):
    """Cycle average of Re(A e^{ikx}) Re(B e^{ikx}) for complex amplitudes A, B."""
    return sp.re(sp.expand(A * sp.conjugate(B))) / 2


def effective_moduli():
    """Return (mu_eff(k), mu_frozen, f_phi(k), g_c(k)).

    A plane wave u_3 = U e^{ikx} slaves the microrotation to the displacement
    through the second field equation in the static limit,

        Phi(k) = -i kappa_c k U / (gamma k^2 + 2 kappa_c),

    which is the elimination that turns the coupled pair into the single
    Cosserat kernel.  Substituting it back into W and cycle-averaging gives an
    effective shear modulus mu_eff(k) through W_avg = mu_eff k^2 U^2 / 4.
    """
    D = gam * k ** 2 + 2 * kc                       # rotational denominator [Pa]
    Phi = -sp.I * kc * k * U / D                    # slaved microrotation [rad]
    up = sp.I * k * U                               # shear strain amplitude [-]

    W_relaxed = (
        mu / 2 * cycle_average(up, up)
        + kc / 2 * (cycle_average(up + Phi, up + Phi) + cycle_average(Phi, Phi))
        + gam / 2 * cycle_average(sp.I * k * Phi, sp.I * k * Phi)
    )
    # Same wave, microrotation held at zero: the reference the coupling is
    # switched on from.
    W_frozen = mu / 2 * cycle_average(up, up) + kc / 2 * cycle_average(up, up)

    mu_eff = sp.simplify(4 * W_relaxed / (k ** 2 * U ** 2))
    mu_frozen = sp.simplify(4 * W_frozen / (k ** 2 * U ** 2))

    f_phi = sp.simplify((mu_frozen - mu_eff) / mu_frozen)   # against frozen
    g_c = sp.simplify((mu_frozen - mu_eff) / mu_eff)        # against relaxed
    return mu_eff, mu_frozen, f_phi, g_c


def report():
    W, force_u, force_phi = strain_energy_density()
    print("=" * 74)
    print("1. Strain energy density and the field equations it generates")
    print("=" * 74)
    print("  W        =", sp.expand(W))
    print("  force_u  =", force_u, "   [expect (mu+kappa_c) u'' + kappa_c phi']")
    print("  force_phi=", force_phi, "  [expect gamma phi'' - kappa_c u' - 2 kappa_c phi]")

    mu_eff, mu_frozen, f_phi, g_c = effective_moduli()
    mubar = mu + kc / 2                              # relaxed (low-k) modulus
    mutot = mu + kc                                  # frozen (high-k) modulus
    N2 = kc / (2 * (mu + kc))                        # coupling number

    print()
    print("=" * 74)
    print("2. Effective transverse modulus and its two limits")
    print("=" * 74)
    print("  mu_eff(k)   =", sp.factor(mu_eff))
    print("  mu_eff(0)   =", sp.limit(mu_eff, k, 0), " ; mubar - mu_eff(0) =",
          sp.simplify(mubar - sp.limit(mu_eff, k, 0)))
    print("  mu_eff(inf) =", sp.limit(mu_eff, k, sp.oo), " ; mutot - mu_eff(inf) =",
          sp.simplify(mutot - sp.limit(mu_eff, k, sp.oo)))
    print("  mu_frozen   =", mu_frozen, " (equals mu_tot, as it must)")

    print()
    print("=" * 74)
    print("3. The rotational fraction, against each reference energy")
    print("=" * 74)
    lorentzian = N2 / (1 + gam * k ** 2 / (2 * kc))
    print("  f_phi(k)                        =", sp.simplify(f_phi))
    print("  f_phi(k) - N^2/[1+gamma k^2/2kc]=", sp.simplify(f_phi - lorentzian),
          "  <- zero means the printed k-dependence is exact")
    print("  f_phi(0)                        =", sp.simplify(sp.limit(f_phi, k, 0)))
    print("  f_phi(0) - N^2                  =", sp.simplify(sp.limit(f_phi, k, 0) - N2))
    print("  g_c(0)                          =", sp.simplify(sp.limit(g_c, k, 0)))
    print("  g_c(0) - N^2/(1-N^2)            =", sp.simplify(sp.limit(g_c, k, 0) - N2 / (1 - N2)))

    # Two equivalent readings of N^2 worth recording.
    print()
    print("  N^2 = 1 - mubar/mutot ?  residual =",
          sp.simplify(N2 - (1 - mubar / mutot)))
    print("  g_c = mutot/mubar - 1 ?  residual =",
          sp.simplify(kc / (2 * mu + kc) - (mutot / mubar - 1)))

    print()
    print("=" * 74)
    print("4. Numbers at the rolling-contact operating point")
    print("=" * 74)
    N2v = 1.0 / np.pi
    kcv = 2 * N2v / (1 - 2 * N2v)                    # kappa_c in units of mu
    mubarv, mutotv = 1 + kcv / 2, 1 + kcv
    print(f"  N^2 (target)              = {N2v:.9f}   = 1/pi")
    print(f"  kappa_c / mu              = {kcv:.9f}   = 2/(pi-2)")
    print(f"  kappa_c / [2(mu+kappa_c)] = {kcv/(2*(1+kcv)):.9f}   = 1/pi        <- coupling number")
    print(f"  kappa_c / (2 mu+kappa_c)  = {kcv/(2+kcv):.9f}   = 1/(pi-1)    <- g_c, not N^2")
    print(f"  ratio of the two          = {(kcv/(2+kcv))/(kcv/(2*(1+kcv))):.9f}   = pi/(pi-1)")
    print(f"  1 - mubar/mutot           = {1-mubarv/mutotv:.9f}")
    print(f"  mutot/mubar - 1           = {mutotv/mubarv-1:.9f}")
    print(f"  discrepancy if confused   = {(kcv/(2+kcv))/N2v - 1:+.2%}")


if __name__ == "__main__":
    report()
