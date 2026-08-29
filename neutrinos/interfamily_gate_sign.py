#!/usr/bin/env python3
"""The sign of the inter-family gate, derived analytically. It is bonding.

THE QUESTION
interfamily_hopping.py bracketed |t| between 7e-6 and 1.5e-2 meV and left the
sign open, noting that an antibonding gate (t > 0) would close most of the
solar splitting's -0.66 sigma residual while a bonding gate (t < 0) would
deepen it. The sign was called the sector's cheapest live derivation target.
It is cheaper than that: it needs no computation at all. Three independent
arguments fix it, and they agree.

ARGUMENT 1: SUPEREXCHANGE POSITIVITY (the sign of the vertex is irrelevant)
t is not a bare matrix element. The inter-family hop is a virtual passage
through the constricted intermediate at gap Delta_f = 28.6 MeV above the
glide manifold, so degenerate second-order perturbation theory gives

    t = sum_c <A|V|c><c|V|B> / (E_0 - E_c) = -(1/Delta_f) sum_c <A|V|c><c|V|B>.

The denominator is negative because the intermediate lies ABOVE: that is what
"virtual passage through a barrier state" means. If the intermediate is even
under the A <-> B exchange then <A|V|c> = <B|V|c> = g_x and

    t = -g_x^2 / Delta_f  <  0    for either sign of g_x.

The unknown gate vertex enters squared, so its own sign never mattered. What
mattered was the parity of the intermediate, and that is fixed: the
constricted state is the perfect, undissociated dislocation sitting on the
<110> line where the two {111} planes intersect. That line and its Burgers
vector both lie in the mirror that exchanges the two planes, so the
intermediate is even, and t < 0.

ARGUMENT 2: ISOTROPY REQUIRES IT (the strongest of the three)
Four families coupled pairwise form K4, the tetrahedron. Its spectrum is
exactly (verified in sympy below)

    A_1 singlet at 3t,    T_2 triplet at -t.

The framework selects the A_1 singlet on the grounds that the propagating
neutrino is isotropic, equal on all four families. A_1 is the ground state
iff 3t < -t, i.e. iff t < 0. Were t > 0, the lowest four-family state would
be the threefold T_2 multiplet, which is NOT equal on all four families: the
neutrino would acquire a preferred direction in the crystal. That is a hard
Lorentz violation the framework cannot tolerate. So the isotropy the
framework already assumes is equivalent to the bonding sign, and the
projection step was never sign-neutral.

ARGUMENT 3: VARIATIONAL FLOOR
Any second-order admixture of states above the manifold lowers the manifold's
energy (the standard variational statement). The A_1 combination is the one
that admits the constriction coherently from all four families, so it must be
the one lowered most, which again places it at the bottom: t < 0.

THE CONSEQUENCE, AND IT IS UNFAVOURABLE
delta H_eff = 3t I shifts every mass by 3t, so with t < 0 all three masses
fall, and the solar splitting falls with them:

    delta(Dm2_21) = 6 t (m_2 - m_1) < 0.

The prediction already sits 0.66 sigma BELOW the measured solar splitting, so
the derived sign deepens the residual rather than closing it. The earlier
statement that an antibonding gate might supply +0.45 sigma is withdrawn: no
antibonding gate is available.

WHAT THAT BUYS ANYWAY
A one-sided constraint, which is worth more than a free parameter. Requiring
the degradation to stay inside 0.1 sigma caps |t| at 3.4e-3 meV, which
excludes the top quarter-decade of the mechanics bracket on data grounds and
tells the gate-vertex computation what it must come out under. The mass sum
moves the other way, 9t < 0 lowering sum m_nu slightly and easing the
cosmological comparison by a hair.
"""

import numpy as np
import sympy as sp

# ---- sector numbers [meV, meV^2]
M1, M2 = 5.03, 10.00
DM21_PRED, DM21_EXP, DM21_ERR_LO = 74.711, 75.37, 1.00   # meV^2, low-side sigma
T_LO, T_HI = 7e-6, 1.5e-2                                # |t| bracket, meV
SUM_PRED = 65.515                                        # meV


def k4_spectrum():
    """Exact spectrum of four pairwise-coupled families."""
    t = sp.symbols('t', real=True)
    H = t * (sp.ones(4, 4) - sp.eye(4))
    return t, H.eigenvals()


def solar_shift(t):
    """delta(Dm2_21) = 6 t (m2 - m1), from m_k -> m_k + 3t."""
    return 6 * t * (M2 - M1)


def pull(dm21):
    return (dm21 - DM21_EXP) / DM21_ERR_LO


def main():
    t, spec = k4_spectrum()
    print("=" * 68)
    print("The inter-family gate is bonding: t < 0")
    print("=" * 68)
    print(f"  K4 spectrum: {spec}")
    print("  A_1 (isotropic, equal on all four families) sits at 3t;")
    print("  T_2 (threefold, anisotropic) sits at -t. A_1 is the ground")
    print("  state iff t < 0, so the framework's isotropy assumption and")
    print("  the bonding sign are the same statement.\n")

    print("  superexchange: t = -g_x^2 / Delta_f for an even intermediate,")
    print("  negative for either sign of the gate vertex g_x.\n")

    print("consequence for the solar splitting (prediction already -0.66 sigma):")
    print(f"  {'|t| [meV]':>12s} {'d(Dm2_21)':>11s} {'Dm2_21':>9s} {'pull':>7s}")
    for tm in [T_LO, 1e-4, 1e-3, 3.4e-3, T_HI]:
        d = solar_shift(-tm)
        print(f"  {tm:12.1e} {d:11.4f} {DM21_PRED + d:9.3f} {pull(DM21_PRED + d):+7.2f}")

    cap = 0.1 * DM21_ERR_LO / (6 * (M2 - M1))
    print(f"\n  keeping the degradation inside 0.1 sigma caps |t| < {cap:.1e} meV,")
    print(f"  which excludes the top of the mechanics bracket ({T_HI:.1e} meV)")
    print(f"  on data grounds and sets a target for the gate-vertex computation.")
    print(f"\n  mass sum moves the other way: 9t = {9*-T_HI:+.3f} meV at the cap's")
    print(f"  worst, easing sum m_nu = {SUM_PRED:.1f} meV by a hair.")


if __name__ == "__main__":
    main()
