#!/usr/bin/env python3
"""
second_order_index_bond_shape.py -- the post-Newtonian beta of the lattice's
weak field, from the bond's shape, the pre-stress law and the source.

Speed law (Born bond-stretch sum with the Compton identity applied locally):
    c/c0 = [V''(s l)/V''(l)] s^3                       (exact, any bond)
Bond to second order, x = s - 1:
    V''(s l)/V''(l) = 1 + xi x + zeta x^2/2,  xi = l V'''/V'',  zeta = l^2 V''''/V''
With y = ln s:
    ln(c/c0) = (xi + 3) y + D y^2/2,           D = zeta - xi^2 + xi
Pre-stress map and source:
    y = Phi/(3 kappa) + m Phi^2,  Phi = eps + sigma eps^2,  kappa = K_cr/mubar
With mu'_n = 2 (xi + 3 = -3 kappa) and g00 = n^-2:
    beta = 1 - sigma - 3 kappa m + D/(18 kappa^2)
Pre-stress linear in the potential, bulk modulus following the bonds:
    m = (1 - xi)/(18 kappa^2)
Barometric pre-stress (ln P linear in Phi) with a power-law contact:
    m = 0,  and n = exp(eps) to all orders.
"""
import numpy as np
import sympy as sp
from scipy.optimize import brentq

xi, zeta, kap, sig, eps, y, m, r, l = sp.symbols("xi zeta kappa sigma epsilon y m r ell")
KAPPA = float((5 - 8 / sp.pi) / (3 * (1 - 1 / sp.pi)))
XI_G = -3 * (1 + KAPPA)


def derive():
    s = sp.exp(y)
    lnc = sp.series(sp.log((1 + xi * (s - 1) + zeta * (s - 1) ** 2 / 2) * s ** 3), y, 0, 3).removeO()
    print("ln(c/c0) =", sp.collect(sp.expand(lnc), y))
    D = sp.Symbol("D")
    Phi = eps + sig * eps ** 2
    lnn = -((xi + 3) * y + D * y ** 2 / 2).subs(y, Phi / (3 * kap) + m * Phi ** 2)
    lnn = sp.expand(sp.series(lnn, eps, 0, 3).removeO().subs(xi, -3 - 3 * kap))
    g00 = sp.expand(sp.series(sp.exp(-2 * lnn), eps, 0, 3).removeO())
    beta = sp.simplify(g00.coeff(eps, 2) / 2)
    print("ln n     =", sp.collect(lnn, eps))
    print("beta     =", beta)
    # map for a pre-stress linear in the potential, K(s) = K V''(s l)/(V''(l) s)
    P = 3 * (y + (xi - 1) * y ** 2 / 2)            # |P|/K_cr, integrated
    w = sp.Symbol("w")                               # w = Phi/(3 kappa)
    y2 = sp.solve(sp.Eq(sp.series(P.subs(y, w + sp.Symbol("c2") * w ** 2), w, 0, 3).removeO(), 3 * w),
                  sp.Symbol("c2"))
    print("linear pre-stress: y = w + c2 w^2 with c2 =", y2, " -> m = (1 - xi)/(18 kappa^2)")
    return beta


def bond_families(xv):
    print(f"\nbond families at xi = {xv:.4f}:")
    fams = {"power law r^-p": xv ** 2 - xv, "Morse at minimum": 7 * xv ** 2 / 9,
            "exponential e^(-r/rho)": xv ** 2, "Koide flat quartic": -(18 * xv + 9) / 3}
    for name, z in fams.items():
        print(f"   {name:24s} zeta = {z:7.3f}  D = {z - xv ** 2 + xv:+8.3f}")
    return fams


def stiffer_than_power(xv):
    print("\nbonds stiffer than any power, tuned to the same xi:")
    def shape(V):
        V2 = sp.diff(V, r, 2)
        return (float((r * sp.diff(V, r, 3) / V2).subs(r, 1)),
                float((r ** 2 * sp.diff(V, r, 4) / V2).subs(r, 1)))
    for rc in (0.1, 0.2, 0.3):
        q = brentq(lambda qv: shape((r - rc) ** (-qv))[0] - xv, 0.05, 40)
        x_, z_ = shape((r - rc) ** (-q))
        print(f"   hard-core soft sphere r_c = {rc}: D = {z_ - x_ ** 2 + x_:+.2f}")
    b = brentq(lambda bv: shape(sp.exp(bv / r))[0] - xv, 0.1, 20)
    x_, z_ = shape(sp.exp(b / r))
    print(f"   exp(b/r), b = {b:.2f}: D = {z_ - x_ ** 2 + x_:+.2f}")


def main():
    beta = derive()
    fams = bond_families(XI_G)
    mlin = (1 - XI_G) / (18 * KAPPA ** 2)
    for label, mv in (("pre-stress linear in potential", mlin), ("barometric pre-stress", 0.0)):
        print(f"\n{label}: m = {mv:.4f}")
        for sv in (0, 0.5):
            need = float(sp.solve(sp.Eq(beta.subs({kap: KAPPA, m: mv, sig: sv}), 1), sp.Symbol("D"))[0])
            row = "  ".join(f"{k.split()[0]}: {float(beta.subs({kap: KAPPA, m: mv, sig: sv, sp.Symbol('D'): z - XI_G ** 2 + XI_G})):+.3f}"
                            for k, z in fams.items())
            print(f"   sigma = {sv}: beta = 1 needs D = {need:.2f};  {row}")
    stiffer_than_power(XI_G)
    s_ = sp.symbols("s", positive=True)
    n_exact = s_ ** (-(XI_G + 3))
    print(f"\npower law, barometric: n = s^{-(XI_G + 3):.4f} with 3 kappa ln s = eps -> "
          f"n = exp({float(-(XI_G + 3) / (3 * KAPPA)):.6f} eps)")
    print(f"p = 3 kappa + 1 = {3 * KAPPA + 1:.4f} = 3(2pi-3)/(pi-1) = {3 * (2 * np.pi - 3) / (np.pi - 1):.4f}")


if __name__ == "__main__":
    main()
