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
    (sigma = -s/4: the field-energy correction, +1/4 for negative field energy,
     -1/4 for positive, from Poisson's equation with the field energy as source)
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


def deflection_2pn():
    """Second-order light deflection through the exponential and the isotropic
    Schwarzschild metrics, by quadrature of Fermat's principle in N = sqrt(g_ij/g00)."""
    from scipy.integrate import quad
    N_exp = lambda rr, mm: np.exp(2 * mm / rr)
    N_sch = lambda rr, mm: (1 + mm / (2 * rr)) ** 3 / (1 - mm / (2 * rr))
    print("\nsecond-order deflection coefficient (d - 4 m/b)/(m/b)^2:")
    for name, N in (("exponential", N_exp), ("Schwarzschild", N_sch)):
        vals = []
        for r0 in (4e3, 1.6e4):
            L = N(r0, 1.0) * r0
            f = lambda u: (L / r0) / np.sqrt(max(N(r0 / u, 1.0) ** 2 - (L * u / r0) ** 2, 1e-300))
            d = 2 * quad(f, 0, 1, limit=400, epsabs=1e-15, epsrel=1e-13)[0] - np.pi
            vals.append((d - 4 / L) * L ** 2)
        print(f"   {name:14s} {vals[0]:.4f} -> {vals[1]:.4f}")
    x = 1476.625 / 6.957e8
    uas = 180 / np.pi * 3600e6
    print(f"   4 pi = {4 * np.pi:.4f}, 15 pi/4 = {15 * np.pi / 4:.4f}; solar limb: GR term "
          f"{15 * np.pi / 4 * x ** 2 * uas:.2f} uas, difference {np.pi / 4 * x ** 2 * uas:.3f} uas")


def relativistic_prestress():
    """Inertia of a stressed medium by a Lorentz boost, relativistic hydrostatics
    dP/dPhi = -(eps + P), and the resulting beta for rest-energy-dominated and
    pressure-proportional media."""
    e_, P, v, c = sp.symbols("varepsilon P v c", positive=True)
    g = 1 / sp.sqrt(1 - v ** 2 / c ** 2)
    Lam = sp.Matrix([[g, g * v / c, 0, 0], [g * v / c, g, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    T = Lam * sp.diag(e_, P, P, P) * Lam.T
    print("\nboosted T^01/c =", sp.simplify(T[0, 1] / c))
    Phi = sp.Symbol("Phi")
    Pf = sp.Function("P")
    sol = sp.dsolve(sp.Eq(Pf(Phi).diff(Phi), -(e_ + Pf(Phi))), ics={Pf(0): 0}).rhs
    print("case 1: P(Phi) =", sp.simplify(sol), "; |P|/eps =", sp.series(-sol / e_, Phi, 0, 3).removeO())
    D = sp.Symbol("D")
    beta_rel = 1 - sig - 2 / (3 * kap) + D / (18 * kap ** 2)
    print("case 1 beta = 1 - sigma - 2/(3 kappa) + D/(18 kappa^2); beta = 1 needs D = 12 kappa =",
          f"{12 * KAPPA:.3f}")
    xv = float(XI_G)
    for name, z in (("power law", xv ** 2 - xv), ("exponential", xv ** 2), ("Morse", 7 * xv ** 2 / 9),
                    ("Koide flat quartic", -(18 * xv + 9) / 3)):
        print(f"   {name:20s} beta = {float(beta_rel.subs({kap: KAPPA, sig: 0, D: z - xv ** 2 + xv})):+.3f}")
    for sv in (0.25, -0.25):
        need = float(sp.solve(sp.Eq(beta_rel.subs({kap: KAPPA, sig: sv}), 1), D)[0])
        print(f"   sigma = {sv:+.2f}: beta = 1 needs D = {need:.2f}")
    print(f"case 2 (eps = c_e P, barometric): exact exponential; mu'_n = 2 needs c_e = 4/(3 kappa) = "
          f"{4 / (3 * KAPPA):.3f}; power-law virial gives 3/p = {3 / (3 * KAPPA + 1):.3f}")


def deflection_orbit_check():
    """The orbit u = C + A cos(w phi) solves the ray equation exactly for N^2 quadratic
    in u, and its swept angle expands to pi + 4 e + pi (2 + a2) e^2."""
    u, ph, mm, L, a2, e = sp.symbols("u phi m L a_2 e", positive=True)
    N2 = 1 + 4 * mm * u + (4 + 2 * a2) * mm ** 2 * u ** 2
    w2 = 1 - (4 + 2 * a2) * mm ** 2 / L ** 2
    C = (2 * mm / L ** 2) / w2
    A = sp.sqrt(1 / (L ** 2 * w2) + C ** 2)
    uu = C + A * sp.cos(sp.sqrt(w2) * ph)
    ok = sp.simplify(sp.expand(sp.diff(uu, ph) ** 2 - (N2.subs(u, uu) / L ** 2 - uu ** 2))) == 0
    defl = sp.series(((2 / sp.sqrt(w2)) * sp.acos(-C / A) - sp.pi).subs(mm, e * L), e, 0, 3).removeO()
    print(f"\norbit satisfies the first integral: {ok}; deflection = {sp.simplify(defl)}")


def field_energy_potential():
    """Potential of the gravitating field energy rho_f = s G M^2/(8 pi c^2 r^4):
    enclosed shells plus exterior shells, source radius R, then the source's own
    field energy absorbed into M."""
    G, M, c, r, R, rp, s_ = sp.symbols("G M c r R r' s", positive=True)
    rho = s_ * G * M ** 2 / (8 * sp.pi * c ** 2 * rp ** 4)
    m_enc = sp.integrate(4 * sp.pi * rp ** 2 * rho, (rp, R, r))
    ext = -G * sp.integrate(4 * sp.pi * rp * rho, (rp, r, sp.oo))
    Phi_f = sp.expand(-G * m_enc / r + ext)
    renorm = sp.expand(Phi_f + G / r * s_ * G * M ** 2 / (2 * c ** 2 * R))
    print("\nfield-energy potential after mass renormalisation:", sp.simplify(renorm),
          "; Poisson check:", sp.simplify(sp.diff(r ** 2 * sp.diff(renorm, r), r) / r ** 2 - 4 * sp.pi * G * rho.subs(rp, r)) == 0)


def main():
    import warnings
    warnings.filterwarnings("ignore")
    beta = derive()
    fams = bond_families(XI_G)
    mlin = (1 - XI_G) / (18 * KAPPA ** 2)
    for label, mv in (("pre-stress linear in potential", mlin), ("barometric pre-stress", 0.0)):
        print(f"\n{label}: m = {mv:.4f}")
        for sv in (0, 0.25, -0.25):
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
    field_energy_potential()
    deflection_orbit_check()
    deflection_2pn()
    relativistic_prestress()


if __name__ == "__main__":
    main()
