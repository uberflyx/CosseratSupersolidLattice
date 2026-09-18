#!/usr/bin/env python3
"""
g_prefactor_dyson.py -- the coherent nineteen-node hop as a Dyson process.

Chapter 4 defines alpha as the ratio of the Peierls locking energy to an elastic
reference energy E_ref, and iterates it as the dimensionless zero-momentum
Peierls-Nabarro amplitude of a Born (Dyson) series. In that reading the coherent
compression event of the Born cluster is a nineteenth-order process: nineteen
single-node hops of amplitude alpha*E_ref, separated by eighteen intermediate
states in which some nodes have hopped and the rest have not. Such a state costs
the spring energy of the bonds strained between the two groups,

    E(S) = E_cut * cut(S),   E_cut = k1 b^2 / 2,
    cut(S) = sum_{i in S, j not in S} w_ij,
    w_ij = [(n_ij . bhat)^2 + r (1 - (n_ij . bhat)^2)] * (1 or k2/k1),

with n_ij the bond direction, bhat the hop direction, r the tangential contact
ratio and k2/k1 the second-shell spring ratio. The effective amplitude between
the two degenerate rigid states is

    A / E_ref = alpha^19 * (E_ref / E_cut)^18 * P,
    P = sum over the 19! orderings of prod_{k=1}^{18} 1 / cut(S_k),

and P is evaluated exactly by dynamic programming over the 2^19 subsets.
The Peierls-Nabarro functional pins E_ref: Poisson resummation of the discrete
misfit sum gives every harmonic of the Peierls potential as W_m/W_0 = alpha^m with
W_0 = E_mis = Gamma d^2 w/(2 pi), so per node E_ref = E_mis * ell and, with the hop
over the misfit period d, E_ref/E_cut = w/(2 pi d111) = 1/(8 sqrt 2). The script
checks the harmonics against a direct row sum and evaluates the Dyson prefactor
with that ratio; it falls short of the 1.31 the measured G needs by 1e16 or more.
"""
import itertools
import math

import numpy as np

ALPHA = 1 / 137.035999177
N = 19


def born_cluster_weights(bhat, r, k2_over_k1):
    """Bond weights w_ij of the Born cluster for a hop along bhat."""
    shell1 = sorted({p for v in ([1, 1, 0], [1, -1, 0], [-1, -1, 0])
                     for p in itertools.permutations(v)})
    shell2 = [(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2)]
    pts = [np.zeros(3)] + [np.array(p, float) for p in shell1] \
        + [np.array(p, float) for p in shell2]
    W = np.zeros((N, N))
    for i, j in itertools.combinations(range(N), 2):
        d = pts[i] - pts[j]
        d2 = d @ d
        if d2 not in (2.0, 4.0):
            continue
        c2 = ((d / np.sqrt(d2)) @ bhat) ** 2
        W[i, j] = W[j, i] = (1.0 if d2 == 2 else k2_over_k1) * (c2 + r * (1 - c2))
    return W


def ordering_sum(W):
    """P = sum over orderings of prod 1/cut(S_k), by DP over subsets by size."""
    idx = np.arange(1 << N)
    cut = np.zeros(1 << N)
    for i in range(N):
        for j in range(i + 1, N):
            if W[i, j] == 0:
                continue
            m = ((idx >> i) & 1) ^ ((idx >> j) & 1)
            cut[m == 1] += W[i, j]
    pop = np.array([bin(s).count("1") for s in range(1 << N)])
    g = np.zeros(1 << N)
    g[0] = 1.0
    for k in range(1, N):                      # g(S) = sum_{i in S} g(S - i) / cut(S)
        Ss = np.nonzero(pop == k)[0]
        acc = np.zeros(len(Ss))
        for i in range(N):
            has = ((Ss >> i) & 1) == 1
            acc[has] += g[Ss[has] ^ (1 << i)]
        g[Ss] = acc / cut[Ss]
    full = (1 << N) - 1                        # last hop lands on the degenerate state
    return sum(g[full ^ (1 << i)] for i in range(N)), cut


def peierls_harmonics(w_over_d=np.pi / 4, n_rows=200000, n_u=512):
    """Harmonic ratios W_m/W_0 of the discrete Lorentzian row sum (d = 1)."""
    n = np.arange(-n_rows, n_rows + 1)
    u = np.linspace(0.0, 1.0, n_u, endpoint=False)
    W = np.array([np.sum(1.0 / ((n - ui) ** 2 + w_over_d ** 2)) for ui in u])
    F = np.abs(np.fft.rfft(W))
    return [F[m] / F[0] for m in (1, 2, 3)], np.exp(-2 * np.pi * w_over_d)


def vertex_and_profile():
    """Vertex weight at zero wavevector against the misfit wavevector, the
    anharmonicity that mu'_n = 2 then requires, and rigid versus independent
    nineteen-node form factors for three core shapes (width set so each shape's
    single-node form factor at G equals alpha)."""
    from scipy.optimize import brentq
    N2, kc = 1 / np.pi, 2 / (np.pi - 2)           # coupling number, kappa_c/mu
    q2 = 2 * kc                                   # q^2 = 2 kappa_c / gamma, gamma = mu l^2
    Gw = 2 * np.pi * np.sqrt(3)                   # misfit wavevector 2 pi/d, d = l/sqrt3
    fG = N2 / (1 + Gw ** 2 / q2)
    print(f"C0 at k -> 0: {1 + N2:.4f};  C0 at G = 2 pi/d: {1 + fG:.4f} "
          f"(ratio {(1 + fG) / (1 + N2):.3f})")
    G = 2 * np.pi
    shapes = {"Lorentzian": lambda k, w: np.exp(-abs(k) * w),
              "sech (sine-Gordon)": lambda k, w: 1 / np.cosh(np.pi * k * w / 2),
              "Gaussian": lambda k, w: np.exp(-(k * w) ** 2 / 2)}
    for name, f in shapes.items():
        w = brentq(lambda w: f(G, w) - ALPHA, 1e-3, 5)
        print(f"   {name:20s} rigid/independent 19-node form factor = "
              f"{f(19 * G, w) / f(G, w) ** 19:.3e}")


def compton_inertia_and_budget():
    """Index coefficient with each node's mass set locally by the Compton identity
    m v l = hbar, jointly with the Born sum v^2 = S2 V''(l) l^2/(2m); the xi that
    mu'_n = 2 requires; walk-count radius of the cluster series; precision budget."""
    import sympy as sp
    s_, xi = sp.symbols("s xi", positive=True)
    V2 = sp.Function("V2")
    v = s_ ** 3 * V2(s_)                           # joint solution, up to constants
    rho = 1 / (v * s_ ** 4)                        # m/(s l)^3 with m = hbar/(v s l)
    sub = {sp.Derivative(V2(s_), s_): xi * V2(s_) / s_}
    dv = sp.simplify((s_ * sp.diff(sp.log(v), s_)).subs(sub))
    dr = sp.simplify((s_ * sp.diff(sp.log(rho), s_)).subs(sub))
    print(f"d ln v/d ln s = {dv},  d ln rho_n/d ln s = {sp.expand(dr)}")
    N2 = 1 / np.pi
    K = (5 - 8 * N2) / (3 * (1 - N2))              # K_cr/mubar
    print(f"mu'_n = -2(xi+3)/(3K/mu): {-2 * (-7 + 3) / 3 / K:.4f} at xi = -7; "
          f"mu'_n = 2 at xi = -(8pi-11)/(pi-1) = {-(8 * np.pi - 11) / (np.pi - 1):.5f}")
    W = born_cluster_weights(np.array([1.0, 0, 0]), 1.0, 1.0)   # r = 1: every bond present
    A = (W > 0).astype(float)
    lam = np.linalg.eigvalsh(A).max()
    print(f"cluster adjacency lambda_max = {lam:.3f}; radius {1 / lam:.3f}; "
          f"alpha inside by {1 / lam / ALPHA:.1f}x")
    print(f"budget: 17a/18 = {17 * ALPHA / 18 * 1e6:.0f} ppm; a^3 = {ALPHA ** 3 * 1e6:.2f} ppm; "
          f"inputs ~ {21 * 0.15 + 2 * 0.3:.1f} ppb; xi from 22 ppm: +/- {22e-6 * 3.5992:.1e}")


def main():
    compton_inertia_and_budget()
    print()
    vertex_and_profile()
    print()
    ratios, a0 = peierls_harmonics()
    print("Peierls harmonics W_m/W_0 against alpha_0^m:")
    for m, r in enumerate(ratios, 1):
        print(f"   m={m}: {r:.8e}  {a0 ** m:.8e}")
    d, d111 = 1 / np.sqrt(3), np.sqrt(2 / 3)
    eratio = (np.pi / 4 * d) / (2 * np.pi * d111)
    print(f"E_ref/E_cut = w/(2 pi d111) = {eratio:.6f} (1/(8 sqrt2) = {1 / (8 * np.sqrt(2)):.6f})\n")

    cases = [("central springs, hop along <110>", [1, 1, 0], 0.0, 0.0),
             ("central springs, hop along <112>", [1, 1, -2], 0.0, 0.0),
             ("rolling contact r = 0.2264, k2 = k1, <112>", [1, 1, -2], 0.2264, 1.0),
             ("rolling contact, hop along <110>", [1, 1, 0], 0.2264, 1.0)]
    print(f"{'spring model':44s} {'P':>10s} {'Dyson C':>12s} {'E_ref/E_cut for 1.31':>22s}")
    for label, b, r, k2 in cases:
        bhat = np.array(b, float)
        bhat /= np.linalg.norm(bhat)
        P, cut = ordering_sum(born_cluster_weights(bhat, r, k2))
        need = (1.3092 / P) ** (1 / 18)
        print(f"{label:44s} {P:10.3e} {eratio ** 18 * P:12.3e} {need:22.3f}")
    print(f"\ntarget alpha_G = G m0^2/(hbar c) = 3.3e-41; 19! = {math.factorial(N):.3e}")


if __name__ == "__main__":
    main()
