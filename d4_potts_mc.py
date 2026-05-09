#!/usr/bin/env python3
"""
d4_potts_mc.py — 4-state Potts model on the D4 lattice.

Determines the critical temperature T_c for compact-direction ordering
in the D4 vacuum lattice. The result is used in the Cosserat Supersolid
Lattice framework to assess whether thermal fluctuations of the compact
direction contribute to the chirality domain-wall energy.

The D4 lattice consists of integer points (x1,x2,x3,x4) with even
coordinate sum. Each node has 24 nearest neighbours (all permutations
of (±1, ±1, 0, 0)). Each node carries a Potts variable σ ∈ {0,1,2,3}
representing one of four equivalent compact directions.

Energy: E = -J Σ_{<ij>} δ(σ_i, σ_j)

Result: first-order transition at kT_c/J ≈ 7.1.
The crystallisation temperature is kT/J = π (exact: Θ_D = π m₀c²).
Since π ≈ 3.14 << 7.1, the compact direction is deeply ordered
at crystallisation. No thermal misorientation contributes to the
domain-wall energy.

Requires: numpy, numba
Usage:    python3 d4_potts_mc.py
"""

import numpy as np
from numba import njit
import time


# ─────────────────────────────────────────────────────────────
# Lattice construction
# ─────────────────────────────────────────────────────────────

def build_d4_lattice(L):
    """
    Build a D4 lattice in a periodic box of side L.

    Returns
    -------
    N : int
        Number of D4 sites.
    nbr : ndarray of shape (N, 24), dtype int32
        Neighbour table: nbr[i, k] is the index of the k-th
        nearest neighbour of site i.
    """
    idx_map = -np.ones((L, L, L, L), dtype=np.int32)
    sites = []
    for x1 in range(L):
        for x2 in range(L):
            for x3 in range(L):
                for x4 in range(L):
                    if (x1 + x2 + x3 + x4) % 2 == 0:
                        idx_map[x1, x2, x3, x4] = len(sites)
                        sites.append((x1, x2, x3, x4))
    N = len(sites)

    # 24 nearest-neighbour offsets: all permutations of (±1, ±1, 0, 0)
    offsets = []
    for i in range(4):
        for j in range(i + 1, 4):
            for si in [1, -1]:
                for sj in [1, -1]:
                    v = [0, 0, 0, 0]
                    v[i] = si
                    v[j] = sj
                    offsets.append(v)
    assert len(offsets) == 24

    nbr = np.zeros((N, 24), dtype=np.int32)
    for idx, (x1, x2, x3, x4) in enumerate(sites):
        for k, (d1, d2, d3, d4) in enumerate(offsets):
            n1 = (x1 + d1) % L
            n2 = (x2 + d2) % L
            n3 = (x3 + d3) % L
            n4 = (x4 + d4) % L
            nbr[idx, k] = idx_map[n1, n2, n3, n4]

    return N, nbr


# ─────────────────────────────────────────────────────────────
# Numba-accelerated Monte Carlo
# ─────────────────────────────────────────────────────────────

@njit
def mc_sweep(states, nbr, N, beta, site_rng, new_rng, acc_rng):
    """
    One Metropolis sweep: N single-site update attempts.

    Parameters
    ----------
    states : ndarray of int8, shape (N,)
        Current Potts states (0, 1, 2, or 3).
    nbr : ndarray of int32, shape (N, 24)
        Precomputed neighbour table.
    N : int
        Number of sites.
    beta : float
        Inverse temperature J/kT.
    site_rng, new_rng, acc_rng : ndarray of float64, shape (N,)
        Pre-generated uniform random numbers in [0, 1).
    """
    for step in range(N):
        site = int(site_rng[step] * N)
        if site >= N:
            site = N - 1
        old = states[site]
        new = int(new_rng[step] * 3)
        if new >= old:
            new += 1
        # Energy change from flipping site: old → new
        nn_old = 0
        nn_new = 0
        for k in range(24):
            s = states[nbr[site, k]]
            if s == old:
                nn_old += 1
            if s == new:
                nn_new += 1
        dE = -(nn_new - nn_old)  # in units of J
        if dE <= 0:
            states[site] = np.int8(new)
        elif acc_rng[step] < np.exp(-beta * dE):
            states[site] = np.int8(new)


@njit
def measure(states, nbr, N):
    """
    Measure energy per site and Potts order parameter.

    Returns
    -------
    E_per_site : float
        Energy per site in units of J (each bond counted once).
    m : float
        Order parameter: (q × max_fraction - 1) / (q - 1),
        where q = 4. Ranges from 0 (disordered) to 1 (ordered).
    """
    E = 0.0
    counts = np.zeros(4, dtype=np.int64)
    for i in range(N):
        s = states[i]
        counts[s] += 1
        for k in range(24):
            if states[nbr[i, k]] == s:
                E -= 0.5  # each bond counted from both ends
    mx = 0
    for j in range(4):
        if counts[j] > mx:
            mx = counts[j]
    m = (4.0 * mx / N - 1.0) / 3.0
    return E / N, m


# ─────────────────────────────────────────────────────────────
# Simulation driver
# ─────────────────────────────────────────────────────────────

def run_scan(L, temps, n_eq=200, n_meas=200, cold_start=True):
    """
    Scan a range of temperatures and measure E, m at each.

    Parameters
    ----------
    L : int
        Box side length (4D periodic box).
    temps : array-like
        Values of kT/J to scan.
    n_eq : int
        Number of equilibration sweeps.
    n_meas : int
        Number of measurement sweeps.
    cold_start : bool
        If True, start from fully ordered state at each temperature.
    """
    N, nbr = build_d4_lattice(L)
    print(f"  L={L}, N={N} D4 nodes, {n_eq} eq + {n_meas} meas sweeps")
    print(f"  Start: {'cold (ordered)' if cold_start else 'hot (random)'}")
    print(f"  {'kT/J':>7s}  {'E/N':>9s}  {'m':>7s}  {'χ':>9s}  {'C_v':>9s}")
    print(f"  " + "─" * 46)

    rng = np.random.RandomState(42)
    results = []
    for kT in temps:
        beta = 1.0 / kT
        if cold_start:
            states = np.zeros(N, dtype=np.int8)
        else:
            states = rng.randint(0, 4, size=N).astype(np.int8)

        for _ in range(n_eq):
            mc_sweep(states, nbr, N, beta,
                     rng.random(N), rng.random(N), rng.random(N))

        Es = np.zeros(n_meas)
        Ms = np.zeros(n_meas)
        for i in range(n_meas):
            mc_sweep(states, nbr, N, beta,
                     rng.random(N), rng.random(N), rng.random(N))
            Es[i], Ms[i] = measure(states, nbr, N)

        E_avg = np.mean(Es)
        m_avg = np.mean(Ms)
        chi = N * np.var(Ms) / kT
        Cv = N * np.var(Es) / kT**2
        results.append((kT, E_avg, m_avg, chi, Cv))
        print(f"  {kT:7.2f}  {E_avg:9.4f}  {m_avg:7.4f}  {chi:9.1f}  {Cv:9.1f}")

    return results


# ─────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("4-STATE POTTS MODEL ON THE D4 LATTICE")
    print("Finding T_c for compact-direction ordering")
    print("=" * 60)

    # Compile numba functions with a tiny lattice
    print("\n  Compiling numba kernels...", end=" ", flush=True)
    t0 = time.time()
    N_t, nbr_t = build_d4_lattice(4)
    st = np.zeros(N_t, dtype=np.int8)
    r = np.random.random(N_t)
    mc_sweep(st, nbr_t, N_t, 1.0, r, r, r)
    measure(st, nbr_t, N_t)
    print(f"done ({time.time() - t0:.1f}s)")

    # Coarse scan: L=8
    print("\n--- COARSE SCAN (L=8, cold start) ---")
    t0 = time.time()
    temps_coarse = np.arange(4.0, 10.0, 0.5)
    res_c = run_scan(8, temps_coarse, n_eq=200, n_meas=200, cold_start=True)
    print(f"  Wall time: {time.time() - t0:.0f}s")

    # Find the transition
    for i in range(1, len(res_c)):
        if res_c[i][2] < 0.1 and res_c[i-1][2] > 0.5:
            T_lo = res_c[i-1][0]
            T_hi = res_c[i][0]
            print(f"\n  Transition between kT/J = {T_lo:.1f} and {T_hi:.1f}")
            break

    # Fine scan around transition
    print(f"\n--- FINE SCAN (L=8, cold start) ---")
    t0 = time.time()
    temps_fine = np.arange(T_lo - 0.2, T_hi + 0.2, 0.1)
    res_f = run_scan(8, temps_fine, n_eq=300, n_meas=300, cold_start=True)
    print(f"  Wall time: {time.time() - t0:.0f}s")

    # Identify T_c as midpoint of the energy jump
    for i in range(1, len(res_f)):
        if res_f[i][2] < 0.1 and res_f[i-1][2] > 0.3:
            Tc = 0.5 * (res_f[i-1][0] + res_f[i][0])
            break

    # Summary
    print(f"\n{'=' * 60}")
    print(f"RESULT")
    print(f"{'=' * 60}")
    print(f"  First-order transition at kT_c/J ≈ {Tc:.2f}")
    print(f"  Crystallisation temperature:  kT/J = π ≈ {np.pi:.4f}")
    print(f"  Ratio: T_c / T_cryst = {Tc / np.pi:.2f}")
    print(f"")
    print(f"  The compact direction is ORDERED at crystallisation")
    print(f"  (T_cryst is {Tc/np.pi:.1f}× below T_c).")
    print(f"  No thermal misorientation of the compact direction")
    print(f"  survives the crystallisation transition.")
    print(f"{'=' * 60}")
