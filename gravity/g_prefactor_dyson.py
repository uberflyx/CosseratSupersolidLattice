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
The result is a pure number for each spring model; what the monograph does not
yet supply is E_ref, which enters to the eighteenth power.
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


def main():
    cases = [("central springs, nearest neighbours only", [1, 1, 0], 0.0, 0.0),
             ("rolling contact r = 0.2264, k2 = k1", [1, 1, 0], 0.2264, 1.0),
             ("rolling contact, hop along [100]", [1, 0, 0], 0.2264, 1.0)]
    print(f"{'spring model':44s} {'P':>10s} {'alpha^19 P':>12s} {'E_ref/E_cut for 1.31':>22s}")
    for label, b, r, k2 in cases:
        bhat = np.array(b, float)
        bhat /= np.linalg.norm(bhat)
        P, cut = ordering_sum(born_cluster_weights(bhat, r, k2))
        need = (1.3092 / P) ** (1 / 18)
        print(f"{label:44s} {P:10.3e} {ALPHA**19 * P:12.3e} {need:22.3f}")
    print(f"\ntarget alpha_G = G m0^2/(hbar c) = 3.3e-41; 19! = {math.factorial(N):.3e}")


if __name__ == "__main__":
    main()
