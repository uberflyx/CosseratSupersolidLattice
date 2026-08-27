#!/usr/bin/env python3
"""
potts_fcc_field_endpoint.py -- where the deconfinement jump dies.

Baryogenesis in this framework needs the stacking transition to be genuinely
first order, because a crossover has no bubble walls for the propagation
chirality to bias across.  The Polyakov derivation of T_c says in the same
breath that dynamical quarks soften the pure-gauge transition into the observed
crossover, and those two statements cannot both stand unexamined.

The question is decidable, and it is a classical one wearing a new hat.  Quarks
in the fundamental representation break the Z_3 centre symmetry explicitly,
which in the Svetitsky-Yaffe mapping is an external field h coupled to one
Potts state.  A field on a first-order line weakens it and terminates it at a
critical endpoint h_c; beyond that there is only a crossover.  This is the
same structure as the Columbia plot, where heavy quarks leave deconfinement
first order and light quarks wash it out.

So the framework's baryogenesis survives if and only if the effective quark
field sits below h_c, and h_c for the FCC lattice is not in the literature.
This script measures it.

Model:  E = -J sum_<ij> delta(s_i,s_j) - H sum_i delta(s_i,0)
        K = J/kT,  h = H/kT,  three states, twelve neighbours, FCC.

Method: for each h, reweight the bond-count histogram to the coupling of equal
peak heights and measure the latent heat as the peak separation.  The
first-order line ends where that separation collapses to the width of a single
peak.  Reported alongside is the transition strength alpha_PT the measured
latent heat implies, since a jump that survives but shrinks is no more use to
baryogenesis than one that vanishes.
"""
import time

import numpy as np
from numba import njit

hbar_c = 197.3269804
K_C0 = 0.2578
T_C = 156.1
ELL = hbar_c / (0.51099895 / (1 / 137.035999177))


def build_fcc(L):
    N = 4 * L**3
    basis = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    nn = []
    for s1 in (-.5, .5):
        for s2 in (-.5, .5):
            nn += [[s1, s2, 0], [s1, 0, s2], [0, s1, s2]]
    nn = np.array(nn)
    pos = np.zeros((N, 3))
    sub = np.zeros(N, dtype=np.int32)
    for iz in range(L):
        for iy in range(L):
            for ix in range(L):
                for s in range(4):
                    idx = ((iz * L + iy) * L + ix) * 4 + s
                    pos[idx] = [ix + basis[s][0], iy + basis[s][1],
                                iz + basis[s][2]]
                    sub[idx] = s
    nbrs = np.zeros((N, 12), dtype=np.int32)
    for i in range(N):
        np_ = (pos[i] + nn) % L
        k = 0
        for dv in range(12):
            frac = np_[dv] % 1.0
            cell = np.floor(np_[dv]).astype(np.int32) % L
            for s2 in range(4):
                if np.allclose(frac, basis[s2], atol=1e-6):
                    nbrs[i, k] = ((cell[2] * L + cell[1]) * L + cell[0]) * 4 + s2
                    k += 1
                    break
        assert k == 12
    masks = [np.where(sub == s)[0].astype(np.int32) for s in range(4)]
    return nbrs, N, masks


@njit(cache=True)
def sweep(spins, nbrs, K, h, mask, r):
    """Heat bath with a Z_3-breaking field on state 0."""
    for ii in range(mask.shape[0]):
        i = mask[ii]
        c0 = 0; c1 = 0; c2 = 0
        for k in range(12):
            s = spins[nbrs[i, k]]
            if s == 0:
                c0 += 1
            elif s == 1:
                c1 += 1
            else:
                c2 += 1
        w0 = np.exp(K * c0 + h); w1 = np.exp(K * c1); w2 = np.exp(K * c2)
        x = r[ii] * (w0 + w1 + w2)
        if x < w0:
            spins[i] = 0
        elif x < w0 + w1:
            spins[i] = 1
        else:
            spins[i] = 2


@njit(cache=True)
def bonds(spins, nbrs, N):
    t = 0
    for i in range(N):
        si = spins[i]
        for k in range(12):
            if spins[nbrs[i, k]] == si:
                t += 1
    return t // 2


def series(L, K, h, n_therm, n_meas, seed, hot):
    nbrs, N, masks = build_fcc(L)
    rng = np.random.default_rng(seed)
    spins = (rng.integers(0, 3, N) if hot else np.zeros(N)).astype(np.int8)
    for _ in range(n_therm):
        for m in masks:
            sweep(spins, nbrs, K, h, m, rng.random(m.shape[0]))
    out = np.empty(n_meas, dtype=np.int64)
    for t in range(n_meas):
        for m in masks:
            sweep(spins, nbrs, K, h, m, rng.random(m.shape[0]))
        out[t] = bonds(spins, nbrs, N)
    return out, N


def latent_heat(L, K_run, h, n_meas=40000, seed=11):
    """Peak separation of the reweighted bond-count histogram, per site."""
    a, N = series(L, K_run, h, 3000, n_meas, seed, True)
    b, _ = series(L, K_run, h, 3000, n_meas, seed + 1, False)
    M = np.concatenate([a, b])
    lo, hi = M.min(), M.max()
    Mx = np.arange(lo, hi + 1)
    H, _ = np.histogram(M, bins=np.arange(lo - .5, hi + 1.5))
    logH = np.where(H > 0, np.log(np.maximum(H, 1)), -np.inf)

    def w_of(K):
        lw = logH + (K - K_run) * Mx
        lw -= lw.max()
        return np.exp(lw)

    def split(w):
        mid = int(np.sum(w * Mx) / np.sum(w))
        im = np.searchsorted(Mx, mid)
        if im < 4 or im > len(Mx) - 4:
            return None
        return im, w[:im].max(), w[im:].max()

    klo, khi = K_run - 0.03, K_run + 0.03
    for _ in range(60):
        km = .5 * (klo + khi)
        sp = split(w_of(km))
        if sp is None:
            break
        _, plo, phi = sp
        if phi > plo:
            khi = km
        else:
            klo = km
    Ks = .5 * (klo + khi)
    w = w_of(Ks)
    sp = split(w)
    if sp is None:
        return Ks, 0.0, 0.0
    im, _, _ = sp
    ilo = int(np.argmax(w[:im])); ihi = im + int(np.argmax(w[im:]))
    imin = ilo + int(np.argmin(w[ilo:ihi + 1]))
    de = (Mx[ihi] - Mx[ilo]) / N
    depth = np.log(np.sqrt(w[ilo] * w[ihi]) / max(w[imin], 1e-300))
    return Ks, de, depth


def alpha_PT(de):
    """Transition strength implied by a latent bond count de (J per site)."""
    J = K_C0 * T_C
    n_site = np.sqrt(2) / ELL**3
    L_v = de * J * n_site
    g = 61.75
    rho = (np.pi**2 / 30) * g * T_C**4 / hbar_c**3
    return L_v / rho


if __name__ == "__main__":
    print("FCC three-state Potts under a Z_3-breaking field: where the jump dies\n")
    print("  The field h is what dynamical quarks look like in the")
    print("  Svetitsky-Yaffe mapping. h = 0 is pure gauge.\n")
    print(f"{'h':>6}{'K*':>9}{'latent de [J/site]':>21}{'valley depth':>14}"
          f"{'alpha_PT':>11}")
    L = 8
    rows = []
    for h in (0.0, 0.02, 0.05, 0.10, 0.20, 0.40):
        t0 = time.time()
        Ks, de, depth = latent_heat(L, 0.2585 - 0.35 * h, h)
        rows.append((h, Ks, de, depth))
        print(f"{h:6.2f}{Ks:9.5f}{de:21.4f}{depth:14.2f}{alpha_PT(de):11.2e}"
              f"   ({time.time()-t0:.0f}s)")
    print("\n  A valley depth below about 1 means the two peaks have merged:")
    print("  no coexistence, no bubbles, crossover.")
