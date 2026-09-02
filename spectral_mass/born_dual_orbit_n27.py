#!/usr/bin/env python3
"""
born_dual_orbit_n27.py
======================
The sixth fully symmetric host below the third-shell closure: the Born
cluster with both tetrahedral void orbits, N = 27, point group O_h.

The host-admissibility ladder (host_admissibility.py) closes the list of
O_h and T_d hosts below N = 43 at {13, 17, 19, 21, 23, 27}; five are in use
and 27 is not.  Its mass band, 1836 to 2070 MeV, holds the Delta(1920) 3/2+
and the Delta(1910) 1/2+, so this script applies the framework's adiabatic
continuation rule there and reports what it returns.

Procedure (the one sec-delta1600-spectral uses): linear interpolation of the
Cosserat dynamical matrix between the parent cluster with the new atoms
decoupled and the full cluster, forty-one tau steps, global Hungarian
assignment of eigenvectors between steps.  The parent mass mode's degenerate
subspace is followed as a whole.

Validation first: the Delta(1232) T_1 root at 9.0515 on N = 17 and its
descendant on the dual-orbit N = 21, which lands at 10.1085 in T_2u (the
label subduction requires, since T_2u of O_h restricts to T_1 of T_d).

Result: on N = 27 the Delta(1600)'s root reaches T_2u at 10.924 (1986 MeV)
when the second shell is switched on; the Delta(1232)'s root reaches the
T_1g partner at 11.400 (1993 MeV) when routed through Born plus one
tetrahedron; the Sigma's isovector A_2 root reaches A_2g at 9.498 (1967 MeV).
All three sit about 3.5 % above the Breit-Wigner masses of the Delta(1920)
and Delta(1910), outside the framework's alpha m band, and both J = 3/2
descendants keep under three quarters of their weight on the original
shell.  The host exists; the continuation rule does not read those states
on it.  A failed attempt on one host, not a no-go.

Author: M. Cox, with Claude (Anthropic).  License: MIT.
"""

import sys, numpy as np
from scipy.linalg import eigh
from scipy.optimize import linear_sum_assignment
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spectral_classifier import cluster_coord_shell, A_LAT
from hadron_spectral_mass import mass_from_lambda
from delta_first_principles import (build_cosserat_matrix_two_d, generate_Td,
                                    build_td_projector, TD_CHARS, cluster_delta)
from spectral_classifier import generate_Oh
from delta_first_principles import build_cosserat_rep

s = A_LAT / 4.0
VP = np.array([[+s, +s, +s], [+s, -s, -s], [-s, +s, -s], [-s, -s, +s]])
VM = -VP
SHELL, _ = cluster_coord_shell()
SECOND = np.array([[sg * A_LAT if ax == k else 0.0 for k in range(3)]
                   for ax in range(3) for sg in (+1, -1)])


def H_interp(coords, tau, n_old):
    N = len(coords)
    Hf = build_cosserat_matrix_two_d(coords)
    Hi = build_cosserat_matrix_two_d(coords[:n_old])
    Hp = np.zeros((6 * N, 6 * N))
    Hp[:3*n_old, :3*n_old] = Hi[:3*n_old, :3*n_old]
    Hp[3*N:3*N+3*n_old, 3*N:3*N+3*n_old] = Hi[3*n_old:, 3*n_old:]
    Hp[:3*n_old, 3*N:3*N+3*n_old] = Hi[:3*n_old, 3*n_old:]
    Hp[3*N:3*N+3*n_old, :3*n_old] = Hi[3*n_old:, :3*n_old]
    return (1 - tau) * Hp + tau * Hf


def embed(v_old, n_old, N):
    v = np.zeros(6 * N)
    v[:3*n_old] = v_old[:3*n_old]
    v[3*N:3*N+3*n_old] = v_old[3*n_old:]
    return v


def track(coords, n_old, seed_vecs, n_steps=41):
    """Global Hungarian tracking; returns final lambdas and vectors of the seeds."""
    N = len(coords)
    taus = np.linspace(0, 1, n_steps)
    vals, vecs = eigh(H_interp(coords, 0.0, n_old))
    # locate the seed vectors inside the tau=0 eigenbasis
    ov = (np.abs(vecs.T @ seed_vecs) ** 2).sum(axis=1)   # weight of each tau=0 vector in the seed subspace
    k = seed_vecs.shape[1]
    idx = list(np.argsort(ov)[-k:])
    assert ov[idx].sum() > 0.95 * k, ov[idx]
    labels = np.arange(6 * N)
    for tau in taus[1:]:
        v2, V2 = eigh(H_interp(coords, tau, n_old))
        O = np.abs(vecs.T @ V2) ** 2
        r, c = linear_sum_assignment(-O)
        perm = np.zeros(6 * N, dtype=int); perm[r] = c
        vals, vecs = v2[perm], V2[:, perm]
    return vals[idx], vecs[:, idx]


def chars(v, N, inner=13, new=()):
    phi = np.sum(v[3*N:]**2) / np.sum(v**2)
    per = np.array([np.sum(v[3*i:3*i+3]**2) + np.sum(v[3*N+3*i:3*N+3*i+3]**2) for i in range(N)])
    per /= per.sum()
    return phi, per[:inner].sum(), per[list(new)].sum() if new else 0.0


OH_CHARS = {  # classes: E 8C3 3C2 6C4 6C2' i 8S6 3sigma_h 6S4 6sigma_d
    'A_1g': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 'A_2g': [1, 1, 1, -1, -1, 1, 1, 1, -1, -1],
    'E_g': [2, -1, 2, 0, 0, 2, -1, 2, 0, 0], 'T_1g': [3, 0, -1, 1, -1, 3, 0, -1, 1, -1],
    'T_2g': [3, 0, -1, -1, 1, 3, 0, -1, -1, 1], 'A_1u': [1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
    'A_2u': [1, 1, 1, -1, -1, -1, -1, -1, 1, 1], 'E_u': [2, -1, 2, 0, 0, -2, 1, -2, 0, 0],
    'T_1u': [3, 0, -1, 1, -1, -3, 0, 1, -1, 1], 'T_2u': [3, 0, -1, -1, 1, -3, 0, 1, 1, -1]}


def oh_class(R):
    """Conjugacy class of a signed permutation matrix, by determinant, trace and diagonality."""
    d, t = round(np.linalg.det(R)), round(np.trace(R))
    diag = np.allclose(R, np.diag(np.diag(R)))
    if d == 1:
        return {3: 0, 0: 1, 1: 3}[t] if t != -1 else (2 if diag else 4)
    return {-3: 5, 0: 6, -1: 8}[t] if t != 1 else (7 if diag else 9)


def oh_content(coords, V):
    """Mean O_h irrep content of the columns of V, on a cluster with full O_h symmetry."""
    n = len(coords)
    out = {}
    for nm, ch in OH_CHARS.items():
        P = np.zeros((6 * n, 6 * n))
        for R in generate_Oh():
            D = build_cosserat_rep(R, coords)
            P += ch[oh_class(R)] * D
        P *= ch[0] / 48
        c = float(np.mean([v @ (P @ v) / (v @ v) for v in V.T]))
        if c > 0.01:
            out[nm] = round(c, 3)
    return out


# ---- parent: Delta(1232) T_1 triplet at 9.0515 on N=17
c17 = cluster_delta()
w17, V17 = eigh(build_cosserat_matrix_two_d(c17))
Ptd = build_td_projector(c17, generate_Td(), TD_CHARS['T_1'])
k17 = [k for k in range(len(w17)) if abs(w17[k] - 9.0515) < 1e-3 and V17[:, k] @ (Ptd @ V17[:, k]) > 0.5]
print("parent T_1 triplet:", [round(w17[k], 4) for k in k17])

# ---- validation: N=21 dual orbit
c21 = np.vstack([SHELL, VP, VM])
seeds = np.column_stack([embed(V17[:, k], 17, 21) for k in k17])
lam21, sub21 = track(c21, 17, seeds)
ph, sh, nw = chars(sub21[:, 0], 21, new=range(17, 21))
print(f"Delta(1600) N=21: lambda {lam21.mean():.4f} (spread {np.ptp(lam21):.1e}), "
      f"m {mass_from_lambda(21, lam21.mean()):.2f}; phi {ph:.2f} shell {sh:.2f}; ",
      oh_content(c21, sub21))
print("   expected: 10.1085 / 1536.08, T_2u")

# ---- N=27 route (a): N=21 -> + second shell
c27a = np.vstack([SHELL, VP, VM, SECOND])
seeds = np.column_stack([embed(sub21[:, k], 21, 27) for k in range(3)])
lam_a, sub_a = track(c27a, 21, seeds)
ph, sh, nw = chars(sub_a[:, 0], 27, new=range(21, 27))
print(f"\nN=27 route (a) Delta(1600) + Born orbit: lambda {lam_a.mean():.4f} (spread {np.ptp(lam_a):.1e}), "
      f"m {mass_from_lambda(27, lam_a.mean()):.2f}; phi {ph:.2f} inner {sh:.2f} 2nd {nw:.2f}; ",
      oh_content(c27a, sub_a))

# ---- N=27 route (b): N=17 -> + Born orbit (N=23) -> + second tetrahedron
c23 = np.vstack([SHELL, VP, SECOND])
seeds = np.column_stack([embed(V17[:, k], 17, 23) for k in k17])
lam23, sub23 = track(c23, 17, seeds)
ph, sh, nw = chars(sub23[:, 0], 23, new=range(17, 23))
print(f"\nroute (b) step 1, N=23 Born + one tetrahedron: lambda {lam23.mean():.4f} (spread {np.ptp(lam23):.1e}), "
      f"m {mass_from_lambda(23, lam23.mean()):.2f}; phi {ph:.2f} inner {sh:.2f} 2nd {nw:.2f}")
c27b = np.vstack([SHELL, VP, SECOND, VM])
seeds = np.column_stack([embed(sub23[:, k], 23, 27) for k in range(3)])
lam_b, sub_b = track(c27b, 23, seeds)
ph, sh, nw = chars(sub_b[:, 0], 27, new=range(23, 27))
print(f"route (b) step 2, N=27: lambda {lam_b.mean():.4f} (spread {np.ptp(lam_b):.1e}), "
      f"m {mass_from_lambda(27, lam_b.mean()):.2f}; phi {ph:.2f} inner {sh:.2f}; ",
      oh_content(c27b, sub_b))

# ---- the isovector J=1/2 (Sigma-type A_2) parent on N=17 at 4.624, same two routes
Pa2 = build_td_projector(c17, generate_Td(), TD_CHARS['A_2'])
ka = [k for k in range(len(w17)) if abs(w17[k] - 4.6236) < 1e-3 and V17[:, k] @ (Pa2 @ V17[:, k]) > 0.5]
print("\nisovector J=1/2 parent A_2 at", [round(w17[k], 4) for k in ka])
seeds = np.column_stack([embed(V17[:, k], 17, 27) for k in ka])
# straight from N=17 to N=27 (both additions together) and via each route
for tag, coords, n_old, sd in (("N=17 -> N=27 direct", c27b, 17, seeds),):
    lam, sub = track(coords, n_old, sd)
    ph, sh, nw = chars(sub[:, 0], 27)
    print(f"  {tag}: lambda {lam.mean():.4f}, m {mass_from_lambda(27, lam.mean()):.2f}; phi {ph:.2f} inner {sh:.2f}; ",
          oh_content(coords, sub))


# ---------------------------------------------------------------------------
# Born with one tetrahedron (N = 23): the isovector J = 1/2 channel followed
# along both extension orders, symmetry-restricted so that only same-channel
# avoided crossings can occur.  Both orders return the same three roots.
# ---------------------------------------------------------------------------
def channel_curves(coords, n_old, irrep, n_steps=201):
    """Eigenvalue curves of one T_d channel as the new atoms (index >= n_old)
    are switched on, computed in the channel's own subspace."""
    P = build_td_projector(coords, generate_Td(), TD_CHARS[irrep])
    w, V = np.linalg.eigh(P)
    B = V[:, w > 0.5]
    taus = np.linspace(0, 1, n_steps)
    curves = np.array([np.linalg.eigvalsh(B.T @ H_interp(coords, t, n_old) @ B) for t in taus])
    return taus, curves


print("\nN = 23 isovector A_2 channel along both extension orders (tau = 0, 0.5, 1):")
for tag, coords, n_old in (("Born, then voids", np.vstack([SHELL, SECOND, VP]), 19),
                           ("void host, then second shell", np.vstack([SHELL, VP, SECOND]), 17)):
    taus, cv = channel_curves(coords, n_old, 'A_2')
    gap = np.diff(cv, axis=1).min()
    print(f"  {tag}: smallest same-channel gap {gap:.3f}")
    for i in range(cv.shape[1]):
        print("     " + "  ".join(f"{cv[j, i]:7.4f}" for j in (0, 100, 200)))
print("  Sigma_b reads 9.031 by either order; 7.545 is the Sigma(1660)'s root by either order.")
