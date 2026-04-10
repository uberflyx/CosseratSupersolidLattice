#!/usr/bin/env python3
"""
potts_fcc_numba.py — Production MC for 3D 3-state Potts on FCC.
Numba JIT for the inner loop; sublattice-parallel updates.
"""
import numpy as np
from numba import njit, prange
from scipy.optimize import brentq
import time

# ============================================================
# Lattice construction (numpy, runs once)
# ============================================================
def build_fcc(L):
    N = 4 * L**3
    basis = np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]])
    nn_vecs = []
    for s1 in [-.5,.5]:
        for s2 in [-.5,.5]:
            nn_vecs += [[s1,s2,0],[s1,0,s2],[0,s1,s2]]
    nn_vecs = np.array(nn_vecs)
    
    pos = np.zeros((N, 3))
    sub = np.zeros(N, dtype=np.int32)
    for iz in range(L):
        for iy in range(L):
            for ix in range(L):
                for s in range(4):
                    idx = ((iz*L + iy)*L + ix)*4 + s
                    pos[idx] = [ix+basis[s][0], iy+basis[s][1], iz+basis[s][2]]
                    sub[idx] = s
    
    nbrs = np.zeros((N, 12), dtype=np.int32)
    for i in range(N):
        nn_pos = (pos[i] + nn_vecs) % L
        k = 0
        for dv_idx in range(12):
            npos = nn_pos[dv_idx]
            frac = npos % 1.0
            cell = np.floor(npos).astype(np.int32) % L
            for s2 in range(4):
                if np.allclose(frac, basis[s2], atol=1e-6):
                    nbrs[i, k] = ((cell[2]*L + cell[1])*L + cell[0])*4 + s2
                    k += 1
                    break
        assert k == 12
    
    masks = [np.where(sub == s)[0].astype(np.int32) for s in range(4)]
    return nbrs, N, masks


# ============================================================
# Numba-accelerated sweep
# ============================================================
@njit(cache=True)
def mc_sweep_numba(spins, nbrs, K, mask, rand_vals):
    """Heat-bath update for one sublattice. Fully compiled."""
    n = mask.shape[0]
    for ii in range(n):
        i = mask[ii]
        # Count neighbours in each state
        c0 = 0; c1 = 0; c2 = 0
        for k in range(12):
            s = spins[nbrs[i, k]]
            if s == 0: c0 += 1
            elif s == 1: c1 += 1
            else: c2 += 1
        
        # Heat-bath weights
        w0 = np.exp(K * c0)
        w1 = np.exp(K * c1)
        w2 = np.exp(K * c2)
        total = w0 + w1 + w2
        
        r = rand_vals[ii] * total
        if r < w0:
            spins[i] = 0
        elif r < w0 + w1:
            spins[i] = 1
        else:
            spins[i] = 2


@njit(cache=True)
def compute_energy_numba(spins, nbrs, N):
    """Total energy per site (each bond counted once)."""
    total = 0
    for i in range(N):
        si = spins[i]
        for k in range(12):
            if spins[nbrs[i, k]] == si:
                total += 1
    return -total / (2 * N)  # each bond counted twice


# ============================================================
# Simulation driver
# ============================================================
def simulate(L, K_values, n_therm=3000, n_meas=8000, seed=42):
    print(f"  Building FCC L={L} ({4*L**3} sites)...", end=" ", flush=True)
    nbrs, N, masks = build_fcc(L)
    rng = np.random.default_rng(seed)
    print("done.")
    
    # Warm up numba
    dummy_spins = np.zeros(N, dtype=np.int32)
    dummy_rand = rng.random(masks[0].shape[0])
    mc_sweep_numba(dummy_spins, nbrs, 0.1, masks[0], dummy_rand)
    compute_energy_numba(dummy_spins, nbrs, N)
    
    results = {'K': [], 'E': [], 'C': [], 'E_err': [], 'L': L, 'N': N}
    
    for K in K_values:
        spins = rng.integers(0, 3, size=N).astype(np.int32)
        
        # Thermalise
        for _ in range(n_therm):
            for m_idx in range(4):
                rv = rng.random(masks[m_idx].shape[0])
                mc_sweep_numba(spins, nbrs, K, masks[m_idx], rv)
        
        # Measure
        E_arr = np.zeros(n_meas)
        for meas in range(n_meas):
            for m_idx in range(4):
                rv = rng.random(masks[m_idx].shape[0])
                mc_sweep_numba(spins, nbrs, K, masks[m_idx], rv)
            E_arr[meas] = compute_energy_numba(spins, nbrs, N)
        
        E_mean = E_arr.mean()
        E_var = E_arr.var()
        C = K**2 * N * E_var
        # Jackknife error on C
        n_jack = 20
        C_jack = np.zeros(n_jack)
        block = n_meas // n_jack
        for j in range(n_jack):
            mask_j = np.concatenate([E_arr[:j*block], E_arr[(j+1)*block:]])
            C_jack[j] = K**2 * N * mask_j.var()
        C_err = np.sqrt((n_jack-1) * np.var(C_jack))
        
        results['K'].append(K)
        results['E'].append(E_mean)
        results['C'].append(C)
        results['E_err'].append(C_err)
        print(f"    K={K:.4f}: <E>/N={E_mean:.5f}, C/N={C:.3f} ± {C_err:.3f}")
    
    return results


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 60)
    print("3D 3-STATE POTTS ON FCC — NUMBA PRODUCTION RUN")
    print("=" * 60)
    
    m0 = 0.51099895 / (1/137.035999084)
    eps_SF = np.pi * np.sqrt(3)/2 * m0
    N_tau = 3
    
    # Phase 1: L=6 coarse scan
    print(f"\n--- L=6 ({4*216} sites), coarse ---")
    K_coarse = np.arange(0.20, 0.32, 0.010)
    t0 = time.time()
    r6c = simulate(6, K_coarse, n_therm=2000, n_meas=5000, seed=42)
    K_peak6 = r6c['K'][np.argmax(r6c['C'])]
    print(f"  K_peak(L=6) ≈ {K_peak6:.3f} [{time.time()-t0:.1f}s]")
    
    # Phase 2: L=8 fine scan
    print(f"\n--- L=8 ({4*512} sites), fine ---")
    K_fine = np.arange(K_peak6 - 0.020, K_peak6 + 0.025, 0.003)
    t0 = time.time()
    r8 = simulate(8, K_fine, n_therm=3000, n_meas=8000, seed=123)
    K_peak8 = r8['K'][np.argmax(r8['C'])]
    print(f"  K_peak(L=8) ≈ {K_peak8:.4f} [{time.time()-t0:.1f}s]")
    
    # Phase 3: L=10 fine scan
    print(f"\n--- L=10 ({4*1000} sites), fine ---")
    K_fine10 = np.arange(K_peak8 - 0.012, K_peak8 + 0.015, 0.003)
    t0 = time.time()
    r10 = simulate(10, K_fine10, n_therm=3000, n_meas=8000, seed=456)
    K_peak10 = r10['K'][np.argmax(r10['C'])]
    print(f"  K_peak(L=10) ≈ {K_peak10:.4f} [{time.time()-t0:.1f}s]")
    
    # Phase 4: L=12 confirmation
    print(f"\n--- L=12 ({4*1728} sites), fine ---")
    K_fine12 = np.arange(K_peak10 - 0.009, K_peak10 + 0.012, 0.003)
    t0 = time.time()
    r12 = simulate(12, K_fine12, n_therm=4000, n_meas=10000, seed=789)
    K_peak12 = r12['K'][np.argmax(r12['C'])]
    print(f"  K_peak(L=12) ≈ {K_peak12:.4f} [{time.time()-t0:.1f}s]")
    
    # Best estimate
    K_c = K_peak12
    K_c_err = max(abs(K_peak12 - K_peak10), 0.003)
    
    # Convert to T_c
    def J_eff(T):
        x = np.exp(-eps_SF / T)
        return N_tau * np.log((1 + 2*x) / (1 - x))
    
    Z_eff = 3 + np.sqrt(3)
    J_bond = eps_SF / Z_eff
    T_c = J_bond / K_c
    T_c_lo = J_bond / (K_c + K_c_err)
    T_c_hi = J_bond / (K_c - K_c_err)
    
    print(f"\n{'='*60}")
    print(f"PRODUCTION RESULT")
    print(f"{'='*60}")
    print(f"  K_c(FCC) = {K_c:.4f} ± {K_c_err:.4f}")
    print(f"  K_c(L=6)  = {K_peak6:.4f}")
    print(f"  K_c(L=8)  = {K_peak8:.4f}")
    print(f"  K_c(L=10) = {K_peak10:.4f}")
    print(f"  K_c(L=12) = {K_peak12:.4f}")
    print(f"  Finite-size drift: {abs(K_peak12-K_peak8):.4f}")
    print(f"\n  Z_eff = 3 + sqrt(3) = {Z_eff:.4f}")
    print(f"  T_c = {T_c:.1f} MeV ({T_c_lo:.0f}--{T_c_hi:.0f})")
    print(f"  T_c(QCD) = 156 ± 3 MeV")
    print(f"  Discrepancy: {abs(T_c-156)/156*100:.1f}%")
    print(f"  T_c/Lambda_QCD = {T_c/(np.pi*m0):.4f}")
    
    # Check first-order: C_max should scale as N
    C_max = [max(r6c['C']), max(r8['C']), max(r10['C']), max(r12['C'])]
    N_vals = [4*6**3, 4*8**3, 4*10**3, 4*12**3]
    print(f"\n  First-order check (C_max should scale ~ N):")
    for n, c in zip(N_vals, C_max):
        print(f"    N={n:5d}: C_max/N = {c:.3f}, C_max = {c*n:.0f}")
    
    np.savez('potts_fcc_production.npz',
             K_c=K_c, K_c_err=K_c_err, T_c=T_c, Z_eff=Z_eff,
             K8=r8['K'], C8=r8['C'], E8=r8['E'],
             K10=r10['K'], C10=r10['C'], E10=r10['E'],
             K12=r12['K'], C12=r12['C'], E12=r12['E'])
    print(f"\n  Saved: potts_fcc_production.npz")


if __name__ == '__main__':
    main()
