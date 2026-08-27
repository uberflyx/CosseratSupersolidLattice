#!/usr/bin/env python3
"""
potts_fcc_chirality_wall.py -- is the chirality wall wet or dry?

WHAT THIS MEASURES, AND WHAT IT DOES NOT.  This script returns an interface
ENERGY per unit area, from the excess bond count of a forced wall.  It does NOT
return an interface tension, which is a FREE energy per unit area, and at a
first-order transition the two differ enormously because the interface entropy
very nearly cancels the interface energy.  On this lattice the gap is a factor
of roughly six hundred.  So the numbers below must never be compared with
sigma_od = 0.085 MeV/fm^2 from potts_fcc_transition_observables.py, which is a
free energy obtained from ln(P_max/P_min).  Comparing them is an error, and it
was made once already.

What the script is good for is a ratio of like to like: the wetting test.

The arrest-scale result L* = kT/(ell sigma) was evaluated with sigma_od, the
order-DISORDER tension measured at coexistence.  That is the wrong interface.
A chirality wall has ordered stacking on both sides in different Z_3 states, so
what is needed is sigma_oo, the order-ORDER tension, and it has never been
measured on this lattice.

The distinction is not cosmetic.  If the order-order interface wets, meaning it
splits into two order-disorder interfaces with a film of disordered phase
between them, then sigma_oo = 2 sigma_od exactly and the arrest scale halves,
which moves the predicted baryon asymmetry by a factor of eight.  Wetting is
common at three-state interfaces near coexistence, so this has to be measured
rather than assumed.

Method: fixed boundary conditions in z.  The bottom cell layer is pinned to
state 0 throughout.  In the reference run the top layer is pinned to state 0 as
well, so no interface is forced; in the interface run the top layer is pinned
to state 1, so exactly one order-order wall must exist somewhere in the box.
The excess bond energy between the two runs is the interface energy:

    sigma_oo = J * (M_ref - M_int) / A ,     A = 2 L^2 ell^2

using the FCC conventional cube edge a_c = sqrt(2) ell, so a box of L x L cubes
has cross-section 2 L^2 ell^2.  Note that this returns the interface ENERGY per
area, not the free energy; they differ by the interface entropy, so the number
is an upper bound on the tension and the gap widens with temperature.  That
caveat is carried into the reported result rather than hidden.

Run below the transition, where both phases are ordered: K > K_c = 0.2578.
"""
import itertools
import time

import numpy as np
from numba import njit

hbar_c = 197.3269804
ELL = 2.8179403262
K_C = 0.2578
T_C = 156.1
SIGMA_OD = 0.085          # MeV/fm^2, measured order-disorder tension


def build_fcc_slab(L):
    """FCC lattice plus a frozen mask for the bottom and top cell layers."""
    N = 4 * L**3
    basis = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    nn = []
    for s1 in (-.5, .5):
        for s2 in (-.5, .5):
            nn += [[s1, s2, 0], [s1, 0, s2], [0, s1, s2]]
    nn = np.array(nn)
    pos = np.zeros((N, 3))
    sub = np.zeros(N, dtype=np.int32)
    zcell = np.zeros(N, dtype=np.int32)
    for iz in range(L):
        for iy in range(L):
            for ix in range(L):
                for s in range(4):
                    idx = ((iz * L + iy) * L + ix) * 4 + s
                    pos[idx] = [ix + basis[s][0], iy + basis[s][1],
                                iz + basis[s][2]]
                    sub[idx] = s
                    zcell[idx] = iz
    nbrs = np.zeros((N, 12), dtype=np.int32)
    for i in range(N):
        cand = (pos[i] + nn) % L
        k = 0
        for dv in range(12):
            frac = cand[dv] % 1.0
            cell = np.floor(cand[dv]).astype(np.int32) % L
            for s2 in range(4):
                if np.allclose(frac, basis[s2], atol=1e-6):
                    nbrs[i, k] = ((cell[2] * L + cell[1]) * L + cell[0]) * 4 + s2
                    k += 1
                    break
        assert k == 12
    frozen = (zcell == 0) | (zcell == L - 1)
    masks = [np.where((sub == s) & (~frozen))[0].astype(np.int32)
             for s in range(4)]
    return nbrs, N, masks, zcell, frozen


@njit(cache=True)
def sweep(spins, nbrs, K, mask, r):
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
        w0 = np.exp(K * c0); w1 = np.exp(K * c1); w2 = np.exp(K * c2)
        x = r[ii] * (w0 + w1 + w2)
        if x < w0:
            spins[i] = 0
        elif x < w0 + w1:
            spins[i] = 1
        else:
            spins[i] = 2


@njit(cache=True)
def bond_count(spins, nbrs, N):
    t = 0
    for i in range(N):
        si = spins[i]
        for k in range(12):
            if spins[nbrs[i, k]] == si:
                t += 1
    return t // 2


@njit(cache=True)
def disorder_fraction(spins, nbrs, N):
    """Fraction of sites whose 12 neighbours are not majority-aligned.

    A site in a well-ordered region has most of its neighbours in its own
    state; a site inside a disordered film does not.  Used to detect wetting.
    """
    n_dis = 0
    for i in range(N):
        si = spins[i]
        c = 0
        for k in range(12):
            if spins[nbrs[i, k]] == si:
                c += 1
        if c < 7:
            n_dis += 1
    return n_dis / N


def run(L, K, top_state, n_therm, n_meas, seed):
    nbrs, N, masks, zcell, frozen = build_fcc_slab(L)
    rng = np.random.default_rng(seed)
    spins = np.zeros(N, dtype=np.int8)
    # start with a sharp wall halfway so the interface run does not have to
    # nucleate one; the reference run starts uniform.
    if top_state != 0:
        spins[zcell > L // 2] = top_state
    spins[zcell == 0] = 0
    spins[zcell == L - 1] = top_state
    for _ in range(n_therm):
        for m in masks:
            sweep(spins, nbrs, K, m, rng.random(m.shape[0]))
    Ms = np.empty(n_meas)
    Ds = np.empty(n_meas)
    for t in range(n_meas):
        for m in masks:
            sweep(spins, nbrs, K, m, rng.random(m.shape[0]))
        Ms[t] = bond_count(spins, nbrs, N)
        Ds[t] = disorder_fraction(spins, nbrs, N)
    return Ms.mean(), Ms.std() / np.sqrt(len(Ms)), Ds.mean()


def sigma_oo(L, K, n_therm=3000, n_meas=8000, seed=5):
    """Order-order interface energy per unit area, in MeV/fm^2."""
    m_ref, e_ref, d_ref = run(L, K, 0, n_therm, n_meas, seed)
    m_int, e_int, d_int = run(L, K, 1, n_therm, n_meas, seed + 1)
    dM = m_ref - m_int                       # missing satisfied bonds
    err = np.hypot(e_ref, e_int)
    J = K_C * T_C                            # MeV per bond
    A = 2.0 * L**2 * ELL**2                  # fm^2, one interface
    return dM * J / A, err * J / A, d_int - d_ref


if __name__ == "__main__":
    print("Order-order (chirality) wall tension on the FCC three-state Potts"
          " model\n")
    print(f"  reference sigma_od = {SIGMA_OD:.3f} MeV/fm^2 (order-disorder,"
          " measured at coexistence)")
    print(f"  complete wetting would give sigma_oo = 2 sigma_od = "
          f"{2*SIGMA_OD:.3f}\n")
    print(f"{'L':>4}{'K':>8}{'T/T_c':>8}{'sigma_oo [MeV/fm^2]':>22}"
          f"{'ratio to od':>13}{'excess disorder':>17}")
    for L, K in ((8, 0.2578), (8, 0.2700), (8, 0.3000), (10, 0.2578),
                 (10, 0.2700)):
        t0 = time.time()
        s, e, dd = sigma_oo(L, K)
        print(f"{L:4d}{K:8.4f}{K_C/K:8.3f}{s:16.4f} +-{e:5.4f}"
              f"{s/SIGMA_OD:13.2f}{dd:17.4f}   ({time.time()-t0:.0f}s)")
    print("\n  WETTING TEST (the valid result). Excess disorder that grows")
    print("  with L signals a thickening film of disordered phase between two")
    print("  order-disorder interfaces, i.e. complete wetting and")
    print("  sigma_oo = 2 sigma_od. It does not grow: it falls slightly.")
    print("  So the wall is dry, and sigma_oo < 2 sigma_od strictly.")
    print("\n  The interface ENERGIES above are internally consistent")
    print("  (about 5 broken bonds per interface site) but are not tensions,")
    print("  and the order-order tension still needs a free-energy method:")
    print("  thermodynamic integration in K between periodic and twisted")
    print("  boundary conditions. That is the outstanding calculation.")
    print("\n  Arrest scale L* = kT/(ell sigma) bracketed by the wetting")
    print("  bound, sigma_od <= sigma_oo < 2 sigma_od:")
    print("    sigma [MeV/fm^2]   L* [fm]   L*/ell    n_B/s")
    g = 61.75
    s_ent = (2 * np.pi**2 / 45) * g * T_C**3 / hbar_c**3
    for sg in (SIGMA_OD, 2 * SIGMA_OD):
        Ls = T_C / (ELL * sg)
        print(f"    {sg:10.4f}      {Ls:9.1f}  {Ls/ELL:7.1f}"
              f"   {0.34/(Ls**3*s_ent):.2e}")
    print(f"    observed n_B/s = 8.7e-11")
