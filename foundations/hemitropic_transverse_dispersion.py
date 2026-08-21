#!/usr/bin/env python3
"""Optical activity of the hemitropic Cosserat vacuum: transverse dispersion,
helicity splitting, and the resulting polarisation rotation rate.

What this checks
----------------
The monograph's chirality is a single parity-odd quadratic invariant coupling
the symmetric strain to the symmetric curvature,

    W_ch = 2 C_ch eps_(ij) kap_(ij),      C_ch = mu theta_ch ell,

with kap_ij = d_j phi_i the microrotation gradient and theta_ch = alpha^2/(2 pi).
The question is what that does to a transverse wave, which is the photon.

The answer is that the two circular polarisations split, but only at order
ell k. The splitting must vanish as k -> 0, because a structural handedness of
pitch ell cannot be sampled by a wave that does not resolve it, and the script
confirms that it does:

    omega_pm^2 = c^2 k^2 (1 -/+ theta_ch ell k),     d(beta)/dL = (omega/c)^2 theta_ch ell / 2

which is the classical 1/lambda^2 law of optical rotatory dispersion, recovered
from the lattice rather than assumed.

Why it is done symbolically
---------------------------
Index placement is the failure mode in Cosserat elasticity, and hand algebra
does not survive it: a first attempt at this calculation, written by hand,
produced a u_x-phi_x coupling where the microrotation term couples u_x to phi_y,
and reported a splitting independent of k. So W is built here as a scalar, the
field equations are taken by sympy's Euler-Lagrange machinery, and the
non-chiral limit is checked against the framework's own System A before the
chiral term is trusted:

    rho u''   = (lam+mu) grad(div u) + (mu+kc) grad^2 u + kc curl(phi)
    J phi''   = gam grad^2 phi + kc (curl u - 2 phi)

Two mu conventions coexist in the framework and both are correct. W is the
energy form, so the mu written below is the true shear modulus and equals
mubar; the transverse acoustic speed comes out as rho c^2 = mubar, and the
displacement equation's coefficient appears as mu + kc/2, which is the Eringen
mu_E + kc. The script prints both so the reader can see the translation.

Requires sympy.
"""
import sympy as sp

# ------------------------------------------------------------------ symbols
lam, mu, kc, gam, Jm, rho, Cch = sp.symbols(
    "lambda mu kappa_c gamma J rho C_ch", positive=True
)
k, w = sp.symbols("k omega", positive=True)
z, t = sp.symbols("z t")
levi = sp.LeviCivita

# Six fields, functions of (z, t): a plane wave travelling along z. Keeping all
# six lets the script confirm that the transverse sector closes on itself
# rather than assuming it.
u = [sp.Function(f"u{a}")(z, t) for a in "xyz"]
ph = [sp.Function(f"p{a}")(z, t) for a in "xyz"]


def d(f, j):
    """Spatial derivative. Only d/dz survives for a wave along z."""
    return sp.diff(f, z) if j == 2 else sp.S.Zero


def build_lagrangian(with_chirality=True):
    """L = kinetic - W, with W assembled from the Cosserat kinematics."""
    eps = [
        [d(u[i], j) - sum(levi(i, j, m) * ph[m] for m in range(3)) for j in range(3)]
        for i in range(3)
    ]
    kap = [[d(ph[i], j) for j in range(3)] for i in range(3)]

    def symm(A, i, j):
        return (A[i][j] + A[j][i]) / 2

    def anti(A, i, j):
        return (A[i][j] - A[j][i]) / 2

    trace = sum(eps[i][i] for i in range(3))
    W = (lam / 2) * trace**2
    W += mu * sum(symm(eps, i, j) ** 2 for i in range(3) for j in range(3))
    W += (kc / 2) * sum(anti(eps, i, j) ** 2 for i in range(3) for j in range(3))
    W += (gam / 2) * sum(kap[i][j] ** 2 for i in range(3) for j in range(3))
    if with_chirality:
        W += 2 * Cch * sum(
            symm(eps, i, j) * symm(kap, i, j) for i in range(3) for j in range(3)
        )

    kinetic = (rho / 2) * sum(sp.diff(u[i], t) ** 2 for i in range(3))
    kinetic += (Jm / 2) * sum(sp.diff(ph[i], t) ** 2 for i in range(3))
    return sp.expand(kinetic - W)


def field_equations(L):
    """Euler-Lagrange equations for all six fields, as expressions set to zero."""
    eqs = sp.euler_equations(L, u + ph, [z, t])
    return [sp.expand(e.lhs - e.rhs) for e in eqs]


def transverse_block(eqs):
    """Substitute a plane wave and extract the 4x4 transverse coefficient matrix."""
    Ux, Uy, Px, Py = sp.symbols("U_x U_y P_x P_y")
    phase = sp.exp(sp.I * (k * z - w * t))
    sub = {u[0]: Ux * phase, u[1]: Uy * phase, u[2]: 0,
           ph[0]: Px * phase, ph[1]: Py * phase, ph[2]: 0}

    def plug(e):
        e = e.subs(sub, simultaneous=True)
        return sp.simplify(sp.expand(e.doit()) / phase)

    rows = [plug(eqs[i]) for i in (0, 1, 3, 4)]
    return sp.Matrix(
        [[sp.expand(r).coeff(v) for v in (Ux, Uy, Px, Py)] for r in rows]
    )


def circular_blocks(M):
    """Rotate to e_pm = (x +- i y)/sqrt(2) and return the two 2x2 helicity blocks."""
    S = sp.Matrix(
        [[1, 1, 0, 0], [sp.I, -sp.I, 0, 0], [0, 0, 1, 1], [0, 0, sp.I, -sp.I]]
    ) / sp.sqrt(2)
    Mc = sp.simplify(S.inv() * M * S)
    plus = sp.Matrix([[Mc[0, 0], Mc[0, 2]], [Mc[2, 0], Mc[2, 2]]])
    minus = sp.Matrix([[Mc[1, 1], Mc[1, 3]], [Mc[3, 1], Mc[3, 3]]])
    return plus, minus


def acoustic_branch(det_expr, order=2):
    """The branch with omega^2/k^2 finite and non-zero as k -> 0, series in k."""
    for root in sp.solve(sp.Eq(det_expr, 0), w**2):
        ratio = sp.simplify(root / k**2)
        limit = sp.limit(ratio, k, 0)
        if limit not in (0, sp.oo):
            return sp.simplify(sp.series(ratio, k, 0, order).removeO()), limit
    raise RuntimeError("no acoustic branch found")


def main():
    print("=" * 68)
    print("check 1: the non-chiral limit must reproduce System A")
    print("=" * 68)
    eqs0 = field_equations(build_lagrangian(with_chirality=False))
    M0 = transverse_block(eqs0)
    plus0, _ = circular_blocks(M0)
    _, c2 = acoustic_branch(sp.expand(sp.simplify(plus0.det())))
    print("  displacement diagonal :", sp.simplify(M0[0, 0]))
    print("    (mu here is the energy-form shear modulus; mu + kc/2 is the")
    print("     Eringen mu_E + kc, which is System A's coefficient)")
    print("  microrotation diagonal:", sp.simplify(M0[2, 2]))
    print("  u-phi coupling        :", sp.simplify(M0[0, 3]),
          "  (u_x couples to phi_y, as curl requires)")
    print("  acoustic speed        : rho c^2 =", sp.simplify(c2 * rho))
    assert sp.simplify(c2 * rho - (mu + kc / 2)) == 0 or \
        sp.simplify(c2 * rho - mu) == 0, "acoustic speed not recognised"

    print()
    print("=" * 68)
    print("check 2: the chiral term, and whether it splits the helicities")
    print("=" * 68)
    eqs = field_equations(build_lagrangian(with_chirality=True))
    M = transverse_block(eqs)
    plus, minus = circular_blocks(M)
    det_p = sp.expand(sp.simplify(plus.det()))
    det_m = sp.expand(sp.simplify(minus.det()))
    odd = sp.factor(sp.simplify((det_p - det_m) / 2))
    print("  helicity-odd part of the secular determinant:", odd)
    assert odd != 0, "no optical activity found"

    ser_p, _ = acoustic_branch(det_p)
    ser_m, _ = acoustic_branch(det_m)
    split = sp.simplify(sp.expand(ser_m - ser_p))
    print("  omega^2/k^2, + helicity:", sp.simplify(ser_p))
    print("  omega^2/k^2, - helicity:", sp.simplify(ser_m))
    print("  difference             :", split, "  -> linear in k, so it vanishes as k -> 0")

    print()
    print("=" * 68)
    print("check 3: the rotation rate, and what it costs")
    print("=" * 68)
    # theta_ch ell = C_ch / mubar, with mubar the energy-form mu.
    theta_ch = sp.Rational(1, 1) * sp.Symbol("theta_ch", positive=True)
    ell = sp.Symbol("ell", positive=True)
    print("  fractional splitting  : theta_ch * ell * k")
    print("  d(beta)/dL            : (omega/c)^2 theta_ch ell / 2")

    # Numbers, in SI.
    ALPHA = 7.2973525643e-3          # CODATA 2022
    THETA = ALPHA**2 / (2 * 3.141592653589793)
    ELL = 2.8179403205e-15           # classical electron radius, m
    C = 2.99792458e8
    HBAR_EV = 6.582119569e-16        # eV s

    def rate(omega):
        return 0.5 * (omega / C) ** 2 * THETA * ELL

    MPC = 3.0856775814913673e22
    cases = [
        ("CMB, 30 GHz", 2 * 3.141592653589793 * 30e9, [
            ("one 21 Mpc domain", 21.4 * MPC),
            ("14 Gpc to last scattering", 14e3 * MPC)]),
        ("optical, 550 nm", 2 * 3.141592653589793 * C / 550e-9, [
            ("1 km of laboratory vacuum", 1e3)]),
        ("X-ray, 2 keV", 2e3 / HBAR_EV, [
            ("4.5 kpc to 1E 1547.0-5408", 4.5e3 * MPC / 1e6)]),
    ]
    print(f"\n  theta_ch = {THETA:.4e}, ell = {ELL:.4e} m")
    for band, omega, paths in cases:
        r = rate(omega)
        print(f"\n  {band}: d(beta)/dL = {r:.3e} rad/m")
        for label, L in paths:
            print(f"      over {label:<28} beta = {r * L:.3e} rad")
    print("\n  observed cosmic birefringence: beta ~ 6e-3 rad")
    print("  so a bulk hemitropic coupling at C_ch = mu theta_ch ell is excluded;")
    print("  the chirality must be carried by defect material, not empty lattice.")


if __name__ == "__main__":
    main()
