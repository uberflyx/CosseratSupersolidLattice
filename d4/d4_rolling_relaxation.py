#!/usr/bin/env python3
"""
d4_rolling_relaxation.py — the rolling-relaxation residual eta_0 of a
{111} glide wall in the D4 lattice
=====================================================================

A glide of the upper half-crystal excites the tangential (k_t) springs
of every contact crossing the wall.  The spheres can relax by rolling:
layer-uniform rotations Omega_m lower the slip energy, at the cost of
slipping against their own-side contacts.  This script shows:

  1. SYMMETRY: the relaxation is pure rolling — only the t^w bivector
     (glide direction ^ wall normal) activates, one angle per layer,
     mirror-symmetric about the wall.

  2. THE RESIDUAL: the minimised tangential energy is a geometric
     fraction of the frozen value,

         eta_0 = 0.688340330958...,

     independent of the stiffness ratio r = k_t/k_n (the relaxation
     problem is k_t-homogeneous).  Verified two ways: converged
     least-squares on slabs of increasing depth, and an exact symbolic
     solve of the mirror-reduced rolling chain (the truncated value is
     exactly rational at every depth; the limit is algebraic, defined
     by the chain recurrence, with oscillatory tails).

  3. THE HARMONIC BAND: the relaxed harmonic misfit curvature is
     (2 + 10 eta_0 r)/(3(1+2r)) = 0.816 of the Frenkel value at the
     rolling point, against the frozen-contact 2pi/(3(pi-1)) = 0.978.
     The band [0.82, 0.98] is spanned by how the rotational response
     is partitioned between the elastic kernel and the misfit
     potential — the partition rests on the tangential contact-
     handover law, the open axiom identified in the monograph.

Usage: python3 d4_rolling_relaxation.py     Requires numpy, sympy.

Author: Mitchell A. Cox
Date:   June 2026
"""

import numpy as np
import sympy as sp

s = 1/np.sqrt(2)
w_hat = np.array([1, 1, 1, 0])/np.sqrt(3)
t_hat = np.array([1, 1, -2, 0])/np.sqrt(6)
R = 0.5


def star():
    """24 D4 bonds tagged by layer step dm along the wall normal."""
    out = []
    for i in range(4):
        for j in range(i+1, 4):
            for si in (-1, 1):
                for sj in (-1, 1):
                    v = np.zeros(4); v[i] = si/np.sqrt(2); v[j] = sj/np.sqrt(2)
                    out.append((int(round((v@w_hat)/(1/np.sqrt(6)))), v))
    return out


def eta0_numeric(L):
    """Least-squares rolling relaxation on a slab of L layers per side."""
    S = star()
    gens = []
    for i in range(4):
        for j in range(i+1, 4):
            G = np.zeros((4, 4)); G[i, j] = 1; G[j, i] = -1; gens.append(G)
    ms = list(range(-L+1, L+1)); idx = {m: k for k, m in enumerate(ms)}
    nv = 6*len(ms); rows = []; rhs = []

    def add(m, dm, n, wgt):
        crossing = (m <= 0) and (m+dm >= 1)
        dv = t_hat if crossing else np.zeros(4)
        dvp = dv - (dv@n)*n
        row = np.zeros((4, nv))
        for layer in (m, m+dm):
            if layer in idx:
                base = 6*idx[layer]
                for k, G in enumerate(gens):
                    row[:, base+k] += -R*(G@n)
        sw = np.sqrt(wgt); rows.append(sw*row); rhs.append(sw*dvp)

    for m in range(ms[0]-2, ms[-1]+1):
        for dm, n in S:
            if dm < 0:
                continue
            if dm == 0:
                if m in idx:
                    add(m, 0, n, 0.5)
            elif (m in idx) or ((m+dm) in idx):
                add(m, dm, n, 1.0)
    A = np.vstack(rows); b = np.hstack(rhs)
    Ef = 0.5*(b@b)
    sol, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    rv = A@sol - b
    return 0.5*(rv@rv)/Ef, sol, idx, gens


def eta0_exact(M):
    """Exact mirror-reduced rolling chain: theta_{1-m} = theta_m."""
    s6 = sp.sqrt(6)
    th = sp.symbols(f'a0:{M}', real=True)

    def theta(m):
        k = (m-1) if m >= 1 else (-m)
        return th[k] if k < M else sp.Integer(0)

    Rs = sp.Rational(1, 2)
    E = sp.Integer(0)
    p_st, p_sh = 2/s6, 1/s6
    q2_st = [sp.Rational(4, 12), sp.Rational(1, 12), sp.Rational(1, 12)]
    q2_sh = [sp.Rational(1, 12)]*2 + [sp.Rational(4, 12)]
    for m in range(-M-2, M+3):
        Th = theta(m) + theta(m+2)
        dc = 1 if (m <= 0 and m+2 >= 1) else 0
        for q2 in q2_st:
            E += sp.Rational(1, 2)*(dc*(1-q2) - 2*dc*Rs*Th*p_st + (Rs*Th)**2*(p_st**2+q2))
        Th = theta(m) + theta(m+1)
        dc = 1 if (m <= 0 and m+1 >= 1) else 0
        for q2 in q2_sh:
            E += 2*sp.Rational(1, 2)*(dc*(1-q2) - 2*dc*Rs*Th*p_sh + (Rs*Th)**2*(p_sh**2+q2))
        E += sp.Rational(1, 4)*((Rs*2*theta(m))**2*3)
    eqs = [sp.diff(sp.expand(E), v) for v in th]
    vals = list(sp.linsolve(eqs, list(th)))[0]
    return sp.simplify(E.subs(dict(zip(th, vals)))/5), vals


if __name__ == '__main__':
    print("numeric convergence:")
    for L in (4, 6, 8, 10):
        e, sol, idx, gens = eta0_numeric(L)
        print(f"  L={L:2d}: eta_0 = {e:.9f}")
    # symmetry: dominant components
    co = sol[6*idx[0]:6*idx[0]+6]
    print("active components at the wall layer (t^w is index of the",
          "(t,w) generator in the i<j basis): ",
          [f"{c:+.4f}" for c in co])
    eta, vals = eta0_exact(12)
    print(f"\nexact (M=12 truncation): eta_0 = {sp.N(eta, 12)}  (exactly rational: {eta.is_Rational})")
    print("first rolling angles:", [float(sp.N(v, 6)) for v in vals[:4]])
    r = 1/(3*np.pi - 5)
    e0 = float(sp.N(eta, 12))
    print(f"\nharmonic band at the rolling point:")
    print(f"  frozen : 2pi/(3(pi-1))          = {2*np.pi/(3*(np.pi-1)):.4f}")
    print(f"  relaxed: (2+10 eta0 r)/(3(1+2r)) = {(2+10*e0*r)/(3*(1+2*r)):.4f}")
