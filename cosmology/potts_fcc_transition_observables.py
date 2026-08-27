#!/usr/bin/env python3
"""
potts_fcc_transition_observables.py -- the latent heat and interface tension
of the deconfinement transition, measured, and the two observables they price.

The framework fixes the deconfinement transition with no freedom: a three-state
Potts model on the FCC lattice at K_c = 0.258, mapped to T_c = 156.1 MeV by the
stacking-fault energy (monograph, sec:polyakov_Tc).  The transition is first
order, so as the universe cools through it, bubbles nucleate, collide, and
radiate gravitational waves.  gw_spectrum_crystallisation.py already carries
the spectrum, but with two estimated inputs: the strength alpha_PT from a bag
count, and beta/H = 7.5 from a placeholder.  Both are set by two numbers this
script measures from the Potts model itself:

    L_v      the latent heat per unit volume, from the separation of the
             ordered and disordered peaks of the energy histogram at
             coexistence (single-histogram reweighting in K);
    sigma_od the order-disorder interface tension, from the depth of the
             valley between the peaks (Binder's histogram method,
             2 sigma A = k_B T ln[P_max / P_min] with two interfaces in a
             periodic box).

The same two numbers price the Gibbs-Thomson coefficient of a curved melt
boundary, which is the Hawking-temperature check of the droplet paper:
T_GT/T_H = 8 pi sigma k_B T_c / (L_v hbar c) at T_m = T_c.

Units: Potts coupling J per bond; K = J/(k_B T); J = K_c k_B T_c = 40.3 MeV.
FCC conventional cube edge a_c = sqrt(2) ell, four sites per cube, so the site
density is sqrt(2)/ell^3 and a box of side L cubes has cross-section
A = 2 L^2 ell^2.

Caveat stated up front: at a strongly first-order transition a local heat-bath
algorithm tunnels between phases rarely, and the tunnelling rate falls
exponentially with cross-section.  The script therefore reports the tunnel
count per run, uses modest sizes (L = 8, 10, 12), and treats the L-trend of
sigma as the honest error bar.  A multicanonical run is the upgrade path.
"""
import time

import numpy as np
from numba import njit

# ---------------------------------------------------------------- lattice ---
# (build_fcc, mc_sweep_numba, compute_energy_numba reproduced from
#  potts_fcc_mc.py so this script runs standalone)


def build_fcc(L):
    N = 4 * L**3
    basis = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    nn_vecs = []
    for s1 in (-.5, .5):
        for s2 in (-.5, .5):
            nn_vecs += [[s1, s2, 0], [s1, 0, s2], [0, s1, s2]]
    nn_vecs = np.array(nn_vecs)
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
        nn_pos = (pos[i] + nn_vecs) % L
        k = 0
        for dv in range(12):
            npos = nn_pos[dv]
            frac = npos % 1.0
            cell = np.floor(npos).astype(np.int32) % L
            for s2 in range(4):
                if np.allclose(frac, basis[s2], atol=1e-6):
                    nbrs[i, k] = ((cell[2] * L + cell[1]) * L + cell[0]) * 4 + s2
                    k += 1
                    break
        assert k == 12
    masks = [np.where(sub == s)[0].astype(np.int32) for s in range(4)]
    return nbrs, N, masks


@njit(cache=True)
def mc_sweep_numba(spins, nbrs, K, mask, rand_vals):
    n = mask.shape[0]
    for ii in range(n):
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
        r = rand_vals[ii] * (w0 + w1 + w2)
        if r < w0:
            spins[i] = 0
        elif r < w0 + w1:
            spins[i] = 1
        else:
            spins[i] = 2


@njit(cache=True)
def bond_count(spins, nbrs, N):
    """Number of satisfied bonds M (each once); E = -J M."""
    total = 0
    for i in range(N):
        si = spins[i]
        for k in range(12):
            if spins[nbrs[i, k]] == si:
                total += 1
    return total // 2


def run_histogram(L, K_run, n_therm, n_meas, seed, start):
    """Time series of the bond count M at coupling K_run."""
    nbrs, N, masks = build_fcc(L)
    rng = np.random.default_rng(seed)
    if start == "cold":
        spins = np.zeros(N, dtype=np.int8)
    elif start == "hot":
        spins = rng.integers(0, 3, N).astype(np.int8)
    else:                                    # slab: half ordered, half hot
        spins = rng.integers(0, 3, N).astype(np.int8)
        spins[: N // 2] = 0
    for _ in range(n_therm):
        for m in masks:
            mc_sweep_numba(spins, nbrs, K_run, m, rng.random(m.shape[0]))
    Ms = np.empty(n_meas, dtype=np.int64)
    for t in range(n_meas):
        for m in masks:
            mc_sweep_numba(spins, nbrs, K_run, m, rng.random(m.shape[0]))
        Ms[t] = bond_count(spins, nbrs, N)
    return Ms, N


def reweight_to_coexistence(Ms, K_run, N):
    """
    Reweight P(M) from K_run to the K* of equal peak heights.
    P_K(M) ~ H(M) exp[(K - K_run) M].
    Returns K*, the reweighted histogram, its bin centres, and the peak split.
    """
    lo, hi = Ms.min(), Ms.max()
    bins = np.arange(lo - 0.5, hi + 1.5)
    H, _ = np.histogram(Ms, bins=bins)
    M_ax = np.arange(lo, hi + 1)
    logH = np.where(H > 0, np.log(np.maximum(H, 1)), -np.inf)

    def peaks(K):
        lw = logH + (K - K_run) * M_ax
        lw -= lw.max()
        w = np.exp(lw)
        # split at the minimum between the two lobes of the running mean
        mid = int(np.sum(w * M_ax) / np.sum(w))
        i_mid = np.searchsorted(M_ax, mid)
        if i_mid < 3 or i_mid > len(M_ax) - 3:
            return None
        p_lo = w[:i_mid].max()
        p_hi = w[i_mid:].max()
        return w, p_lo, p_hi

    # bisection on equal peak heights
    Klo, Khi = K_run - 0.02, K_run + 0.02
    for _ in range(60):
        Km = 0.5 * (Klo + Khi)
        out = peaks(Km)
        if out is None:
            break
        _, p_lo, p_hi = out
        # higher K favours order (large M): if ordered peak too big, lower K
        if p_hi > p_lo:
            Khi = Km
        else:
            Klo = Km
    Kstar = 0.5 * (Klo + Khi)
    w, _, _ = peaks(Kstar)
    return Kstar, w, M_ax


def analyse(L, K_run, n_therm=4000, n_meas=60000, seed=7):
    t0 = time.time()
    Ms_h, N = run_histogram(L, K_run, n_therm, n_meas, seed, "hot")
    Ms_c, _ = run_histogram(L, K_run, n_therm, n_meas, seed + 1, "cold")
    Ms = np.concatenate([Ms_h, Ms_c])
    # tunnelling diagnostic: sign changes of (M - median)
    med = np.median(Ms_h)
    tunnels = int(np.sum(np.diff(np.sign(Ms_h - med)) != 0))
    Kstar, w, M_ax = reweight_to_coexistence(Ms, K_run, N)
    # peak positions and valley
    mid = int(np.sum(w * M_ax) / np.sum(w))
    i_mid = np.searchsorted(M_ax, mid)
    i_lo = int(np.argmax(w[:i_mid]))
    i_hi = i_mid + int(np.argmax(w[i_mid:]))
    sl = slice(i_lo, i_hi + 1)
    i_min = i_lo + int(np.argmin(w[sl]))
    Pmax = np.sqrt(w[i_lo] * w[i_hi])        # geometric mean of the two peaks
    Pmin = w[i_min]
    de = (M_ax[i_hi] - M_ax[i_lo]) / N       # latent bond count per site [J]
    # interface free energy: two interfaces of area A = 2 L^2 ell^2 each
    two_sigma_A = np.log(Pmax / max(Pmin, 1e-300))       # in units k_B T
    sigma_lat = two_sigma_A / (2 * 2 * L**2)             # k_B T per ell^2
    dt = time.time() - t0
    return dict(L=L, Kstar=Kstar, de=de, sigma=sigma_lat,
                tunnels=tunnels, secs=dt)


# ------------------------------------------------------- physical pricing ---
def price(de_J, sigma_lat, verbose=True):
    """Convert lattice measurements to alpha_PT, beta/H, f_peak, and Hawking."""
    hbar_c = 197.3269804                      # MeV fm
    M_Pl = 1.22089e22                         # MeV
    alpha_em = 1 / 137.035999177
    m0 = 0.51099895 / alpha_em                # MeV
    ell = hbar_c / m0                         # fm
    K_c = 0.258
    T_c = 156.1                               # MeV
    J = K_c * T_c                             # MeV per bond
    n_site = np.sqrt(2) / ell**3              # sites per fm^3

    L_v = de_J * J * n_site                   # MeV / fm^3
    sigma = sigma_lat * T_c / ell**2          # MeV / fm^2

    g_hot = 61.75
    rho_rad = (np.pi**2 / 30) * g_hot * T_c**4 / hbar_c**3
    alpha_PT = L_v / rho_rad

    # thin-wall nucleation: S3 = 16 pi sigma^3 / (3 (L_v eta)^2), eta = 1-T/Tc
    def S3_over_T(eta):
        T = T_c * (1 - eta)
        return 16 * np.pi * sigma**3 / (3 * (L_v * eta)**2) / T

    def log_Gamma_over_H4(eta):
        T = T_c * (1 - eta)
        s = S3_over_T(eta)
        H = np.sqrt(8 * np.pi**3 * g_hot / 90) * T**2 / M_Pl   # MeV
        return 4 * np.log(T / H) + 1.5 * np.log(max(s, 1e-12) / (2 * np.pi)) - s

    etas = np.linspace(1e-4, 0.5, 20000)
    vals = np.array([log_Gamma_over_H4(e) for e in etas])
    i_star = int(np.argmin(np.abs(vals)))
    eta_s = etas[i_star]
    # beta/H = -d(S3/T)/dln T = (1-eta) d(S3/T)/d eta
    h = 1e-5
    beta_H = (1 - eta_s) * (S3_over_T(eta_s + h) - S3_over_T(eta_s - h)) / (2 * h)
    beta_H = abs(beta_H)
    T_star = T_c * (1 - eta_s)

    # sound-wave GW peak (formulas as in gw_spectrum_crystallisation.py)
    v_w = 1 / np.sqrt(3)
    kv = alpha_PT / (0.73 + 0.083 * np.sqrt(alpha_PT) + alpha_PT)
    Ups = 1 - 1 / np.sqrt(1 + 2 / beta_H) if beta_H > 0 else 1.0
    f_pk = 1.9e-5 * (1 / v_w) * beta_H * (T_star * 1e-3 / 100) \
        * (g_hot / 100)**(1 / 6) * 1e9        # nHz
    Kfac = kv * alpha_PT / (1 + alpha_PT)
    h2 = 2.65e-6 * beta_H**(-1) * Kfac**2 * (100 / g_hot)**(1 / 3) * v_w * Ups

    # Gibbs-Thomson against Hawking, at T_m = T_c
    R_hawking = 8 * np.pi * sigma * T_c / (L_v * hbar_c)

    if verbose:
        print(f"\n  latent heat        L_v      = {L_v:10.4g} MeV/fm^3"
              f"   ({de_J:.3f} J per site)")
        print(f"  interface tension  sigma    = {sigma:10.4g} MeV/fm^2"
              f"   ({sigma_lat:.4f} k_B T_c per ell^2)")
        print(f"  dimensionless      sigma/(L_v ell) = {sigma/(L_v*ell):.4f}"
              f"   (Turnbull ~ 0.45)")
        print(f"  transition strength alpha_PT = {alpha_PT:.3f}"
              f"      (script previously assumed 0.94)")
        print(f"  undercooling       eta*     = {eta_s:.4f}"
              f"    ->  T* = {T_star:.1f} MeV")
        print(f"  inverse duration   beta/H*  = {beta_H:.1f}"
              f"      (script previously assumed 7.5)")
        print(f"  GW peak            f_peak   = {f_pk:.1f} nHz"
              f"     h^2 Omega_peak = {h2:.2e}")
        print(f"  Hawking check      T_GT/T_H = {R_hawking:.3f}"
              f"   (1 means Gibbs-Thomson IS Hawking)")
    return dict(L_v=L_v, sigma=sigma, alpha_PT=alpha_PT, eta=eta_s,
                beta_H=beta_H, f_pk=f_pk, h2=h2, R=R_hawking)


if __name__ == "__main__":
    print("Three-state Potts on FCC at coexistence: latent heat and interface"
          " tension\n")
    results = []
    for L, K_run, nm in ((8, 0.2585, 80000), (10, 0.2582, 80000),
                         (12, 0.2580, 60000)):
        r = analyse(L, K_run, n_meas=nm)
        results.append(r)
        print(f"  L={L:2d}: K* = {r['Kstar']:.5f}   de = {r['de']:.4f} J/site"
              f"   sigma = {r['sigma']:.4f} kT/ell^2   tunnels = "
              f"{r['tunnels']:4d}   ({r['secs']:.0f}s)")

    de = results[-1]["de"]
    sig = results[-1]["sigma"]
    print("\nPricing the observables with the largest-lattice values:")
    price(de, sig)
