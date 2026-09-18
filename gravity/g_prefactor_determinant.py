#!/usr/bin/env python3
"""
g_prefactor_determinant.py -- the order-one prefactor of Newton's constant,
as the rigid-kink instanton of the nineteen-node Born cluster defines it.

The rigid kink fixes the exponent of G = C_G alpha^19 hbar c / m0^2 because its
action is nineteen single-node actions. Its prefactor comes from the eighteen
internal fluctuation operators

    H_n = -M d^2/dtau^2 + V''_PN(u(tau)) + k1 lambda_n ,

each the reflectionless Poschl-Teller kink operator shifted by k1 lambda_n.
Their Gel'fand-Yaglom ratios to the kink-free operators have the closed form
(kappa_n - omega0)/(kappa_n + omega0), kappa_n^2 = omega0^2 + k1 lambda_n / M, so
the internal modes contribute

    P_int = prod_n [ (kappa_n + omega0) / (kappa_n - omega0) ]^(1/2),

a function of the single stiffness ratio s = k1 / (M omega0^2).

The script
  1. checks the closed form against the Gel'fand-Yaglom initial-value problem;
  2. builds the Born cluster's spring network (centre, 12 nearest and 6
     next-nearest neighbours of FCC) and evaluates P_int against s;
  3. tries to fix s from the framework's own numbers, reading the single-node
     factor either as a WKB tunnelling amplitude (S = pi^2/2) or as a Peierls
     barrier of the usual size (V0 ~ alpha k1 b^2), and reports that neither
     reading reproduces both.

Units for step 3: m0 = c = ell = hbar = 1 (the Compton identity hbar = m0 c ell),
Burgers vector b = ell, bond stiffness k1 = 2 m0 c^2 / ell^2 from the Born
bond-stretch sum c^2 = S2 V'' ell^2 / (2 m0) at S2 = 1.
"""
import itertools

import numpy as np
from scipy.integrate import solve_ivp

ALPHA = 1 / 137.035999177


def gy_ratio(kappa, T=12.0):
    """det(-d2 + kappa^2 - 2 sech^2)/det(-d2 + kappa^2) by Gel'fand-Yaglom (omega0 = 1)."""
    def rhs(t, y, with_kink):
        V = 2 / np.cosh(t) ** 2 if with_kink else 0.0
        return [y[1], (kappa ** 2 - V) * y[0]]
    ends = [solve_ivp(rhs, [-T, T], [0.0, 1.0], args=(w,), rtol=1e-12, atol=1e-14).y[0, -1]
            for w in (True, False)]
    return ends[0] / ends[1]


def born_cluster_laplacian(k2_over_k1):
    """Nonzero graph-Laplacian eigenvalues of the 19-node Born cluster, units of k1."""
    shell1 = {p for v in ([1, 1, 0], [1, -1, 0], [-1, -1, 0]) for p in itertools.permutations(v)}
    shell2 = {(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2)}
    pts = [np.zeros(3)] + [np.array(p, float) for p in sorted(shell1)] \
        + [np.array(p, float) for p in sorted(shell2)]
    assert len(pts) == 19
    L = np.zeros((19, 19))
    for i, j in itertools.combinations(range(19), 2):
        d2 = np.sum((pts[i] - pts[j]) ** 2)       # nearest neighbours at d2 = 2
        k = 1.0 if d2 == 2 else (k2_over_k1 if d2 == 4 else 0.0)
        L[i, j] -= k; L[j, i] -= k; L[i, i] += k; L[j, j] += k
    return np.linalg.eigvalsh(L)[1:]              # drop the rigid mode


def internal_prefactor(s, lam):
    kappa = np.sqrt(1.0 + s * lam)
    return np.prod(np.sqrt((kappa + 1) / (kappa - 1)))


def main():
    print("1. Gel'fand-Yaglom closed form (kappa-1)/(kappa+1):")
    for kappa in (1.2, 2.0, 5.0):
        print(f"   kappa = {kappa}: numerical {gy_ratio(kappa):.12f}, "
              f"closed {(kappa - 1) / (kappa + 1):.12f}")

    print("\n2. Internal-mode prefactor on the Born cluster:")
    for r in (0.0, 1.0):
        lam = born_cluster_laplacian(r)
        row = ", ".join(f"s={s:g}: {internal_prefactor(s, lam):.3g}"
                        for s in (1, 10, 137, 1e4))
        print(f"   k2/k1 = {r}: lambda in [{lam.min():.2f}, {lam.max():.2f}];  {row}")

    print("\n3. Fixing s from the framework's numbers (m0 = c = ell = hbar = 1):")
    k1 = 2.0
    S = lambda V0: (2 / np.pi) * np.sqrt(2 * V0)          # WKB action, one period
    w0sq = lambda V0: 2 * np.pi ** 2 * V0                  # well curvature
    lam0, lam1 = born_cluster_laplacian(0.0), born_cluster_laplacian(1.0)
    V0A = (np.pi ** 3 / 4) ** 2 / 2                        # S = pi^2/2
    sA = k1 / w0sq(V0A)
    print(f"   A: S = pi^2/2 needs V0 = {V0A:.2f} m0c^2 = {V0A / k1:.1f} bond energies; "
          f"s = {sA:.4f}; P_int = {internal_prefactor(sA, lam1):.2g}"
          f"-{internal_prefactor(sA, lam0):.2g}")
    V0B = ALPHA * k1                                       # Peierls barrier ~ alpha bond
    sB = k1 / w0sq(V0B)
    print(f"   B: V0 = alpha k1 b^2 gives S = {S(V0B):.3f} (e^-S = {np.exp(-S(V0B)):.2f}); "
          f"s = {sB:.2f}; P_int = {internal_prefactor(sB, lam1):.3g}"
          f"-{internal_prefactor(sB, lam0):.3g}")
    print("   Neither reading gives a node that both tunnels with probability alpha and "
          "sits behind a Peierls-sized barrier.")


if __name__ == "__main__":
    main()
