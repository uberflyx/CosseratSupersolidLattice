#!/usr/bin/env python3
"""
second_order_index_bond_shape.py -- the second-order refractive index of the
weak field against the shape of the internode bond.

The first-order index rests on the speed law c ~ V''(l) l^3 (Born bond-stretch
sum with the Compton identity applied locally). Its second order needs the
bond's quartic shape parameter zeta = l^2 V''''/V'' alongside the anharmonicity
xi = l V'''/V'', and two conventions the first order never chose:

  map   : how pre-stress sets spacing beyond first order
          'linear'        3 kappa ln s = Phi
          'recrystallised' dP = -3 K(s) dy with K(s) ~ V''(s l)/s
  source: Phi = eps + sigma eps^2, sigma = 1/2 with the field-energy iteration,
          sigma = 0 without it.

With mu'_n = 2 imposed (xi + 3 = -3 kappa, kappa = K_cr/mubar), the index is
n = 1 + eps + a2 eps^2 and g00 = -n^-2 gives beta = 3/2 - a2. The script prints
beta as a function of zeta for each convention, evaluates it for Morse,
exponential and power-law bonds at the xi that G fixes, and checks that a
power-law bond with the linear map gives n = exp(eps) to all orders.
"""
import sympy as sp

xi, zeta, kap, sig, eps = sp.symbols("xi zeta kappa sigma epsilon")
KAPPA = (5 - 8 / sp.pi) / (3 * (1 - 1 / sp.pi))       # K_cr/mubar at N^2 = 1/pi
XI_G = -3 * (1 + KAPPA)                               # xi fixed by G (mu'_n = 2)


def beta_expr(recrystallised):
    Phi = eps + sig * eps ** 2
    y = Phi / (3 * kap)
    if recrystallised:
        y -= (xi - 1) / 2 * (Phi / (3 * kap)) ** 2
    y = sp.series(y, eps, 0, 3).removeO()
    x = y + y ** 2 / 2                                  # s - 1
    c_ratio = (1 + xi * x + zeta * x ** 2 / 2) * (1 + 3 * y + sp.Rational(9, 2) * y ** 2)
    n = sp.series(1 / c_ratio, eps, 0, 3).removeO()
    a2 = sp.expand(n.coeff(eps, 2)).subs(xi, -3 * (1 + kap))
    return sp.simplify(sp.Rational(3, 2) - a2)


def main():
    xv = float(XI_G)
    families = {"Morse (7 xi^2/9)": 7 * xv ** 2 / 9,
                "exponential / Toda (xi^2)": xv ** 2,
                "power law r^-p (xi^2 - xi)": xv ** 2 - xv,
                "Koide flat quartic -(18 xi + 9)/3": -(18 * xv + 9) / 3}
    print(f"xi = {xv:.4f}, kappa = {float(KAPPA):.4f}\n")
    for recr in (False, True):
        b = beta_expr(recr)
        print(("recrystallised K(s)" if recr else "linear map") + ":  beta =",
              sp.collect(sp.expand(b), [zeta, sig]))
        for s_ in (0, sp.Rational(1, 2)):
            zreq = float(sp.solve(sp.Eq(b.subs(sig, s_), 1), zeta)[0].subs(kap, KAPPA))
            print(f"   sigma = {s_}: beta = 1 needs zeta = {zreq:.2f}")
            for name, z in families.items():
                bv = float(b.subs({kap: KAPPA, sig: s_, zeta: z}))
                print(f"      {name:36s} zeta = {z:6.2f}   beta = {bv:+.3f}")
        print(f"   d beta/d zeta = {float(sp.diff(b, zeta).subs(kap, KAPPA)):.4f}; "
              f"LLR |beta-1| < 5e-4 pins zeta to +/- {5e-4 / float(sp.diff(b, zeta).subs(kap, KAPPA)):.3f}\n")
    # all-orders check for the power-law bond with the linear map and sigma = 0
    s = sp.symbols("s", positive=True)
    n_exact = s ** (-(XI_G + 3))                        # c/c0 = s^(xi+3) exactly
    print("power law, linear map: n =", sp.simplify(n_exact.subs(s, sp.exp(eps / (3 * KAPPA)))),
          f";  p = -xi - 2 = {float(-XI_G - 2):.4f} = 3(2pi-3)/(pi-1) = {float(3*(2*sp.pi-3)/(sp.pi-1)):.4f}")


if __name__ == "__main__":
    main()
