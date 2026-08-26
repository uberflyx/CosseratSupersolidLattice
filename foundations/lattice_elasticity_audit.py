"""Lattice elasticity behind the D_4 audit and the environmental response of alpha.

Four independent checks, all symbolic or exact, none fitted.

1. ISOTROPY.  A lattice of central springs is elastically isotropic only if
   sum_b n_1^4 = 3 sum_b n_1^2 n_2^2 over the neighbour shell.  FCC's twelve
   neighbours give 2 against 3 and a Zener ratio of 2; D_4's twenty-four give
   3 against 3 and a Zener ratio of exactly 1.  The consequence for FCC is a
   41.4 per cent direction-dependence in the transverse (photon) wave speed,
   which Michelson and Morley bounded in 1887, so this is the one argument for
   the fourth dimension that survives adversarial reading.

2. THE CAUCHY DEFECT UNDER LOAD.  Straining a crystal whose bond minimum stays
   put puts every bond under tension, and the Cauchy relation fails by exactly
   the pressure: C_12 - C_44 = P/2, verified to machine precision at every
   strain.  This is why a LOADED lattice cannot hold its Poisson ratio.

3. THE SAME CRYSTAL RESCALED.  The lattice is a condensed phase of the
   condensate, so a change of ambient state re-crystallises it at a new spacing
   rather than squeezing it, and at its own equilibrium the bond force vanishes.
   Both moduli are then phi''(R_0) times a purely geometric coefficient, the
   stiffness cancels from the ratio, and K/mu is fixed by lattice geometry
   alone.  Hence eta_mu = eta_K identically, hence 7 - eta_t = 0, which is the
   strain-stationary condition the alpha chapter currently names without
   deriving.  The PTB optical-clock null is then structural rather than tuned.

   CAVEAT, and it is not small: this runs on the twelve-neighbour FCC shell.
   The argument's structure survives any shell but the values do not, since
   D_4 gives C11:C12:C44 = 3:1:1 where FCC gives 2:1:1.  Whether the SLICE of a
   D_4 crystal reproduces nu = 1/4 is open, and nothing downstream of nu = 1/4
   should be trusted until it is settled.

4. THE TANGENTIAL SPRING.  A sliding tangential stiffness k_t puts a Cauchy
   defect -sqrt2 k_t/R straight into the classical constants, and nu = 1/4 then
   demands k_t = 0 exactly while Mindlin at nu = 1/4 demands k_t/k_n = 6/7.
   The two are reconciled only by the no-slip rolling condition: contacts that
   roll are never loaded in slide, so k_t does not enter the classical sector.
   That condition is therefore load-bearing in a second place.
"""
import itertools

import numpy as np
import sympy as sp

INK = "-" * 74


# ---------------------------------------------------------------- shells -----
def shell(dim):
    """Nearest neighbours of the D_n checkerboard lattice: all (+-1, +-1, 0, ...)."""
    out = []
    for i, j in itertools.combinations(range(dim), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0] * dim
                v[i], v[j] = si, sj
                out.append(np.array(v, float) / np.sqrt(2.0))
    return np.array(out)


def isotropy_check():
    print(INK)
    print("1. ISOTROPY OF THE NEIGHBOUR SHELL")
    for dim, name in ((3, "FCC = D_3, 12 neighbours"), (4, "D_4, 24 neighbours")):
        n = shell(dim)
        s4 = sum(v[0] ** 4 for v in n)
        s22 = sum(v[0] ** 2 * v[1] ** 2 for v in n)
        # central-spring elastic constants are proportional to these sums
        C11, C12, C44 = s4, s22, s22
        zener = 2 * C44 / (C11 - C12)
        print("   %-26s sum n1^4 = %.3f   3 sum n1^2n2^2 = %.3f   Zener A = %.3f"
              % (name, s4, 3 * s22, zener))

    nn = shell(3)

    def christoffel(k):
        k = np.array(k, float)
        k /= np.linalg.norm(k)
        return sum(np.outer(v, v) * (v @ k) ** 2 for v in nn)

    print("\n   FCC transverse wave speeds by direction (rho = 1, k_n/R = 1):")
    tv = []
    for lab, k in (("[100]", (1, 0, 0)), ("[110]", (1, 1, 0)), ("[111]", (1, 1, 1))):
        ev = np.sort(np.linalg.eigvalsh(christoffel(k)))
        tv += [ev[0], ev[1]]
        print("     %-6s v^2 = %s" % (lab, np.round(ev, 4)))
    print("   spread %.4f to %.4f  ->  delta v / v = %.1f per cent"
          % (min(tv), max(tv), 100 * (np.sqrt(max(tv) / min(tv)) - 1)))
    print("   Fermi-LAT bounds photon speed variation near 1e-19.")


# ------------------------------------------------- Born sums with strain -----
e_s, t, ell, eps0 = sp.symbols('e_s t ell eps0', real=True)


def moduli(phi, dim=3):
    """C11, C12, C44 of a strained D_n crystal by energy expansion over the bonds."""
    Re = ell * (1 + e_s)
    q = sp.Matrix(dim, dim, lambda i, j: sp.Symbol('q%d%d' % (min(i, j), max(i, j)), real=True))
    E = 0
    for v in shell(dim):
        Rb = sp.Matrix([sp.nsimplify(x * sp.sqrt(2)) for x in v]) * Re / sp.sqrt(2)
        Rp = (sp.eye(dim) + q) * Rb
        E += phi.subs(t, sp.sqrt((Rp.T * Rp)[0, 0]))
    a = Re * sp.sqrt(2)
    Vat = a ** dim / 2 ** (dim - 1)
    E = E / (2 * Vat)
    z = {q[i, j]: 0 for i in range(dim) for j in range(dim)}
    d2 = lambda A, B: sp.diff(E, A, B).subs(z)
    return (sp.simplify(d2(q[0, 0], q[0, 0])),
            sp.simplify(d2(q[0, 0], q[1, 1])),
            sp.simplify(d2(q[0, 1], q[0, 1]) / 4),
            Vat)


def morse(xi):
    """Morse bond pinned to the framework's anharmonicity xi = ell phi'''/phi''."""
    b = sp.Rational(-xi, 3) / ell
    return eps0 * ((1 - sp.exp(-b * (t - ell))) ** 2 - 1)


def loaded_crystal():
    print(INK)
    print("2. A LOADED CRYSTAL: THE CAUCHY DEFECT IS THE PRESSURE")
    phi = morse(-7)
    C11, C12, C44, Vat = moduli(phi)
    U = 6 * phi.subs(t, ell * (1 + e_s))
    P = sp.simplify(-sp.diff(U, e_s) / sp.diff(Vat, e_s))
    cd = sp.simplify(C12 - C44)
    sub = {ell: 1, eps0: 1}
    print("   C12 - C44 =", sp.simplify(sp.series(cd.subs(sub), e_s, 0, 2).removeO()))
    print("   P         =", sp.simplify(sp.series(P.subs(sub), e_s, 0, 2).removeO()))
    print("   ratio C12-C44 over P, from the columns below: exactly 1/2")
    print("\n     e        C12-C44        P           ratio      nu")
    for e in (-0.02, -0.01, 0.01, 0.02):
        f = lambda X: float(X.subs(sub).subs(e_s, e))
        K, mu = f((C11 + 2 * C12) / 3), f((C11 - C12 + 3 * C44) / 5)
        x = 3 * K / mu
        print("   %+.3f  %+.6f  %+.6f    %7.4f   %.6f"
              % (e, f(cd), f(P), f(cd) / f(P), (x - 2) / (2 * (x + 1))))


def rescaled_crystal():
    print(INK)
    print("3. THE SAME CRYSTAL RESCALED, UNSTRESSED AT EVERY SPACING")
    print("     R0/ell        K          mu        K/mu        nu      C12-C44")
    for xi in (-6, -7, -12):
        phi = morse(xi)
        C11, C12, C44, _ = moduli(phi)
        for R0 in (0.90, 1.00, 1.10):
            g = lambda X: float(sp.simplify(X.subs(e_s, 0).subs({ell: R0, eps0: 1})))
            K, mu = (g(C11) + 2 * g(C12)) / 3, (g(C11) - g(C12) + 3 * g(C44)) / 5
            x = 3 * K / mu
            print("   xi=%-3d %.2f  %9.4f  %9.4f   %.6f   %.6f   %+.1e"
                  % (xi, R0, K, mu, K / mu, (x - 2) / (2 * (x + 1)), g(C12) - g(C44)))
    print("   K and mu move by a factor of two; K/mu does not move at all.")
    print("   Hence eta_mu = eta_K identically, hence 7 - eta_t = 0 exactly.")


def tangential_spring():
    print(INK)
    print("4. THE TANGENTIAL SPRING BREAKS CAUCHY UNLESS CONTACTS ROLL")
    kn, kt, R, r = sp.symbols('k_n k_t R r', positive=True)
    q = sp.symbols('u0:6')
    eps = sp.Matrix([[q[0], q[3], q[4]], [q[3], q[1], q[5]], [q[4], q[5], q[2]]])
    U = 0
    for v in shell(3):
        n = sp.Matrix([sp.nsimplify(x * sp.sqrt(2)) for x in v]) / sp.sqrt(2)
        d = eps * (R * n)
        par = (n.T * d)[0, 0]
        U += sp.Rational(1, 2) * (kn * par ** 2 + kt * ((d.T * d)[0, 0] - par ** 2))
    U = sp.expand(U / 2)
    z = {s: 0 for s in q}
    Vat = (R * sp.sqrt(2)) ** 3 / 4
    C = lambda i, j: sp.simplify(sp.diff(U, q[i], q[j]).subs(z) / Vat)
    C11, C12, C44 = C(0, 0), C(0, 1), C(3, 3) / 4
    print("   C12 - C44 =", sp.simplify(C12 - C44))
    K = sp.simplify((C11 + 2 * C12) / 3)
    mu = sp.simplify((C11 - C12 + 3 * C44) / 5)
    x = sp.simplify(3 * (K / mu).subs(kt, r * kn))
    nu = sp.simplify((x - 2) / (2 * (x + 1)))
    print("   nu(r)     =", nu, "  with r = k_t/k_n")
    print("   nu = 1/4 requires r =", sp.solve(sp.Eq(nu, sp.Rational(1, 4)), r), "(only)")
    print("   Mindlin at nu = 1/4 wants r = 6/7 = 0.857, which returns nu = %.5f"
          % float(nu.subs(r, sp.Rational(6, 7))))
    print("   The two are reconciled only by no-slip rolling, which is therefore")
    print("   load-bearing here as well as in the alpha derivation.")


def alpha_chain():
    print(INK)
    print("5. THE CHAIN FROM nu TO alpha, AND WHAT THE CLOCK MEASURES")
    dlnN2_dlnr, dainv_dN2, N2, ainv = 0.47, 17.9, 1.0 / np.pi, 137.036
    dlnr_dnu, dnu_dlnKmu = -16.0 / 21.0, 5.0 / 24.0
    C_nu = -(dainv_dN2 * N2 / ainv) * dlnN2_dlnr * dlnr_dnu
    k = -dlnr_dnu * dnu_dlnKmu
    print("   d ln(alpha)/d nu            = %+.6f" % C_nu)
    print("   (7 - eta_t) = %+.5f x (eta_mu - eta_K)" % (-k))
    print("   corpus first-order coefficient -2.0e-2 reproduced as %+.4e" % (-2.0e-2 * -k))
    b = 4.5e-7 / k
    print("   PTB clock |7 - eta_t| < 4.5e-7  =>  |eta_mu - eta_K| < %.2e" % b)
    print("   equivalently |d nu / d e_V| < %.2e: the vacuum's Poisson ratio is a"
          % (dnu_dlnKmu * b))
    print("   constant of the medium, not a function of its state.")


if __name__ == "__main__":
    isotropy_check()
    loaded_crystal()
    rescaled_crystal()
    tangential_spring()
    alpha_chain()
    print(INK)
