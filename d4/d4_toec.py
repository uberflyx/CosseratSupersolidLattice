#!/usr/bin/env python3
"""
d4_toec.py -- third-order elastic constants of the D_4 shell once the contact
carries a tangential spring, and the fate of the C_456 selection rule.

The monograph advertises C_456 = 0 as a selection rule of the shell geometry,
derived from the central-force Born sum where the third-order constant carries
S_222 = sum n1^2 n2^2 n3^2 and every D_4 bond has at least two vanishing spatial
components.  That derivation assumes central forces.  The rolling constraint
puts a tangential spring on every contact, which is what generates kappa_c and
hence the fine structure constant, so the selection rule has to be re-derived
with k_t on before it can be quoted.

Until this week the calculation was not well posed, because a tangential spring
at finite strain needs a law for how k_t moves with bond length, and that law
was a free function.  It is now fixed: eta_t = 7 means k_t tracks the normal
curvature, k_t(L) = r V''(L), at every spacing.  With that, the bond energy is
determined and the expansion is unambiguous.

Model, per bond with reference spatial part R_s and compact part R_4:

    L^2      = |F R_s|^2 + R_4^2 = L_0^2 + 2 R_s.eta.R_s        (eta Lagrangian)
    Delta    = (F - I) R_s                    relative centre displacement
    E_bond   = V(L) + (1/2) r V''(L) [ |Delta|^2 - (Delta.n_hat)^2 ]

with F = (I + 2 eta)^(1/2) the symmetric square root, which makes every
expression a function of eta alone and so manifestly rotation invariant:

    |Delta|^2      = 2 R_s.(I + eta - F).R_s
    Delta.n_hat    = R_s.(I + 2 eta - F).R_s / L

Strain energy per unit reference volume, W = (1/2 Omega_0) sum_bonds E_bond,
expands as

    W = (1/2) C_ij eta_i eta_j + (1/6) C_ijk eta_i eta_j eta_k     (Voigt)

so C_456 = d^3 W / d(eta_4) d(eta_5) d(eta_6), with eta_4 = 2 eta_23 etc.

The second-order block doubles as a validation gate and carries a result of its
own: with k_t on, the lattice violates the Cauchy relation by exactly half the
Cosserat coupling,

    C_12 - C_44 = -kappa_c/2 = -N^2 mu_tot,

on the 24-bond D_4 shell and the 12-bond slice alike, at every r.  Since
N^2 = 1/pi is what fixes alpha, the Cauchy discrepancy measured against the
lattice-scale shear modulus IS the fine structure constant's coupling number.
Note that mu_tot = mubar + kappa_c/2 is the lattice-scale modulus and C44 is
the relaxed one; conflating them is a factor (1+5r)/(1+2r) on D_4.

Method: exact numerical energy as a function of the three shear components,
differentiated by a mixed central-difference stencil.  Richardson extrapolation
in the step size separates a true zero from a small nonzero value, which is the
whole question here.
"""
import itertools

import numpy as np

S = 1.0 / np.sqrt(2.0)


def shell(kind="d4"):
    """Reference bonds: spatial 3-vectors and compact components."""
    sp, cp = [], []
    for i, j in itertools.combinations(range(3), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = np.zeros(3)
                v[i], v[j] = si * S, sj * S
                sp.append(v)
                cp.append(0.0)
    if kind == "d4":
        for i in range(3):
            for si in (1, -1):
                for s4 in (1, -1):
                    v = np.zeros(3)
                    v[i] = si * S
                    sp.append(v)
                    cp.append(s4 * S)
    return np.array(sp), np.array(cp)


def morse_derivs(a_m):
    """V, V', V'', V''' for a Morse contact normalised to V''(1) = 1."""
    D = 1.0 / (2.0 * a_m**2)

    def V(L):
        e = np.exp(-a_m * (L - 1.0))
        return D * ((1.0 - e) ** 2 - 1.0)

    def V2(L):
        e = np.exp(-a_m * (L - 1.0))
        return 2.0 * D * a_m**2 * (2.0 * e * e - e)

    return V, V2


def sym_sqrt(M):
    """Symmetric square root of a symmetric positive-definite matrix."""
    w, Q = np.linalg.eigh(M)
    return (Q * np.sqrt(w)) @ Q.T


def energy(eta, sp, cp, r, a_m):
    """Total bond energy per node for Lagrangian strain eta (3x3 symmetric)."""
    V, V2 = morse_derivs(a_m)
    I = np.eye(3)
    F = sym_sqrt(I + 2.0 * eta)
    A = I + eta - F                       # for |Delta|^2 = 2 R.A.R
    B = I + 2.0 * eta - F                 # for Delta.n = R.B.R / L
    tot = 0.0
    for k in range(len(cp)):
        R = sp[k]
        L2 = R @ R + cp[k] ** 2 + 2.0 * (R @ eta @ R)
        L = np.sqrt(L2)
        d2 = 2.0 * (R @ A @ R)            # |Delta|^2
        dn = (R @ B @ R) / L              # Delta . n_hat
        perp2 = d2 - dn * dn
        tot += V(L) + 0.5 * r * V2(L) * perp2
    return 0.5 * tot                      # bonds shared between two nodes


def eta_shear(a, b, c):
    """Pure-shear Lagrangian strain: eta_23 = a, eta_13 = b, eta_12 = c."""
    return np.array([[0.0, c, b], [c, 0.0, a], [b, a, 0.0]])


def mixed_third(f, h):
    """d^3 f / da db dc by the 8-point mixed central stencil."""
    tot = 0.0
    for sa in (1, -1):
        for sb in (1, -1):
            for sc in (1, -1):
                tot += sa * sb * sc * f(sa * h, sb * h, sc * h)
    return tot / (8.0 * h**3)


def C456(sp, cp, r, a_m, vol0, h=2e-3):
    """C_456 in Voigt convention, Richardson-extrapolated in the step size."""
    def f(a, b, c):
        return energy(eta_shear(a, b, c), sp, cp, r, a_m) / vol0

    d3_h = mixed_third(f, h)
    d3_2h = mixed_third(f, 2 * h)
    d3 = (4.0 * d3_h - d3_2h) / 3.0       # Richardson, kills the O(h^2) term
    return d3 / 8.0                       # eta_4 = 2 eta_23, three factors of 2


def second_order(sp, cp, r, a_m, vol0, h=1e-4):
    """C_11, C_12, C_44 by second differences, as a validation gate."""
    def W(eta):
        return energy(eta, sp, cp, r, a_m) / vol0

    E0 = W(np.zeros((3, 3)))

    def e_diag(i, x):
        m = np.zeros((3, 3))
        m[i, i] = x
        return m

    def e_off(i, j, x):
        m = np.zeros((3, 3))
        m[i, j] = m[j, i] = x
        return m

    C11 = (W(e_diag(0, h)) - 2 * E0 + W(e_diag(0, -h))) / h**2
    Wpp = W(e_diag(0, h) + e_diag(1, h))
    Wpm = W(e_diag(0, h) + e_diag(1, -h))
    Wmp = W(e_diag(0, -h) + e_diag(1, h))
    Wmm = W(e_diag(0, -h) + e_diag(1, -h))
    C12 = (Wpp - Wpm - Wmp + Wmm) / (4 * h**2)
    C44 = (W(e_off(1, 2, h)) - 2 * E0 + W(e_off(1, 2, -h))) / (4 * h**2)
    return C11, C12, C44


if __name__ == "__main__":
    VOL0 = 1.0 / np.sqrt(2.0)
    A_M = 7.0 / 3.0
    r_star = 1.0 / (3.0 * np.pi - 5.0)
    LAM = 1.0 / (2.0 * VOL0)

    for tag, kind in (("D_4 (24 bonds)", "d4"), ("FCC slice (12)", "fcc")):
        sp, cp = shell(kind)
        print(f"\n=== {tag} ===")
        # Closed forms, in units of Lambda, with the microrotation relaxed:
        #   D_4 : C11 = 3(1+r), C12 = 1-r, C44 = mubar = 1+2r, kappa_c = 6r
        #   slice: C11 = 2(1+r), C12 = 1-r, C44 = mubar = 1+r,  kappa_c = 4r
        # The lattice-scale modulus is mu_tot = mubar + kappa_c/2, which is
        # 1+5r on D_4 and 1+3r on the slice; it is NOT C44.  The Cauchy defect
        # is C12 - C44 = -kappa_c/2 on both shells and at every r.
        kap = 6.0 * (kind == "d4") + 4.0 * (kind == "fcc")
        c11 = 3.0 if kind == "d4" else 2.0
        c44 = 2.0 if kind == "d4" else 1.0
        print("\nValidation, second order (expect "
              f"C11={c11:.0f}+{c11:.0f}r, C12=1-r, C44=1+{c44:.0f}r, "
              "in units Lambda)\n")
        print(f"{'r':>8}{'C11/L':>10}{'C12/L':>10}{'C44/L':>10}"
              f"{'C12-C44':>11}{'-kap_c/2':>11}{'mu_tot':>10}{'N^2':>9}")
        for r in (0.0, 0.1, r_star, 0.3):
            C11, C12, C44 = second_order(sp, cp, r, A_M, VOL0)
            kc = kap * r * LAM
            mu_tot = C44 + kc / 2.0
            print(f"{r:8.4f}{C11/LAM:10.5f}{C12/LAM:10.5f}{C44/LAM:10.5f}"
                  f"{(C12-C44)/LAM:11.5f}{-kc/(2*LAM):11.5f}"
                  f"{mu_tot/LAM:10.5f}{kc/(2*mu_tot):9.5f}")

        print("\nThird order: the C_456 selection rule\n")
        print(f"{'r':>8}{'C_456':>16}{'relative to C44':>18}")
        for r in (0.0, 0.05, 0.1, r_star, 0.3, 0.6):
            c456 = C456(sp, cp, r, A_M, VOL0)
            _, _, C44 = second_order(sp, cp, r, A_M, VOL0)
            print(f"{r:8.4f}{c456:16.3e}{c456/C44:18.3e}")
