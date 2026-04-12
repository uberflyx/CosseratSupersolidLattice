#!/usr/bin/env python3
"""
d4_transfer_matrix_20x20.py
============================
D4 Cosserat PN transfer matrix with k₄-dependent coupling.

The 20×20 system decomposes into 4 INDEPENDENT blocks at k₂=k₃=0:
  Block CG: (u₃, φ₂, u₃', φ₂', φ₃₄, φ₃₄')  6×6 — GIVES α
  Block DF: (u₂, φ₃, u₂', φ₃', φ₂₄, φ₂₄')  6×6 — 2nd EM polarisation
  Block AE: (u₁, u₁', u₄, φ₁₄, u₄', φ₁₄')  6×6 — compression + compact
  Block Bt: (φ₁, φ₁')                         2×2 — torsion

At k₄=0 these reduce to the existing 4 sectors (A,D,C,Bt) ⊕ System B.
At k₄≠0 the compact-direction fields couple in, the eigenvalues flow.

M. A. Cox, University of the Witwatersrand (2026)
"""

import numpy as np
from scipy.integrate import solve_ivp

# ═══════════════════════════════════════════
# MODULI (units: μ = 1, ℓ = 1)
# ═══════════════════════════════════════════

N2 = 1.0 / np.pi
mu = 1.0
ell = 1.0
kappa_c = 2 * N2 * mu / (1.0 - N2)
mu_tot = mu + kappa_c
gamma_c = mu * ell**2
K_over_mu = 1e6
K = K_over_mu * mu
lam = K - 2*mu/3
M_long = lam + 2*mu + kappa_c
q2 = 2 * kappa_c / gamma_c
alpha_target = 1.0 / 137.035999177

V0 = 1.0
def V_PN(x):
    return V0 * np.sin(np.pi * x / ell)**2


# ═══════════════════════════════════════════
# BLOCK MATRICES
# ═══════════════════════════════════════════

def block_CG(x, k4):
    """Anti-plane shear + anti-plane mixed rotation.
    State: (u₃, φ₂, u₃', φ₂', φ₃₄, φ₃₄')
    At k₄=0: first 4 components = existing Sector C (gives α).
    At k₄≠0: u₃ ↔ φ₃₄ coupling through iκ_c k₄."""
    V = V_PN(x)
    k4sq = k4**2
    M = np.zeros((6, 6), dtype=complex)
    # u₃' = u₃'
    M[0, 2] = 1.0
    # φ₂' = φ₂'
    M[1, 3] = 1.0
    # u₃'' = [V + μ_tot k₄²] u₃/μ_tot - κ_c φ₂'/μ_tot - iκ_c k₄ φ₃₄/μ_tot
    M[2, 0] = (V + mu_tot * k4sq) / mu_tot
    M[2, 3] = -kappa_c / mu_tot
    M[2, 4] = -1j * kappa_c * k4 / mu_tot
    # φ₂'' = [2κ_c + γk₄²] φ₂/γ + κ_c u₃'/γ
    M[3, 1] = (2 * kappa_c + gamma_c * k4sq) / gamma_c
    M[3, 2] = kappa_c / gamma_c
    # φ₃₄' = φ₃₄'
    M[4, 5] = 1.0
    # φ₃₄'' = [2κ_c + γk₄²] φ₃₄/γ + iκ_c k₄ u₃/γ
    M[5, 4] = (2 * kappa_c + gamma_c * k4sq) / gamma_c
    M[5, 0] = 1j * kappa_c * k4 / gamma_c
    return M


def block_DF(x, k4):
    """In-plane shear + in-plane mixed rotation.
    State: (u₂, φ₃, u₂', φ₃', φ₂₄, φ₂₄')
    At k₄=0: first 4 = existing Sector D (2nd EM polarisation).
    At k₄≠0: u₂ ↔ φ₂₄ coupling."""
    V = V_PN(x)
    k4sq = k4**2
    M = np.zeros((6, 6), dtype=complex)
    M[0, 2] = 1.0
    M[1, 3] = 1.0
    # u₂'' = [V + μ_tot k₄²] u₂/μ_tot + κ_c φ₃'/μ_tot + iκ_c k₄ φ₂₄/μ_tot
    M[2, 0] = (V + mu_tot * k4sq) / mu_tot
    M[2, 3] = kappa_c / mu_tot
    M[2, 4] = 1j * kappa_c * k4 / mu_tot
    # φ₃'' = [2κ_c + γk₄²] φ₃/γ - κ_c u₂'/γ
    M[3, 1] = (2 * kappa_c + gamma_c * k4sq) / gamma_c
    M[3, 2] = -kappa_c / gamma_c
    M[4, 5] = 1.0
    # φ₂₄'' = [2κ_c + γk₄²] φ₂₄/γ + iκ_c k₄ u₂/γ
    M[5, 4] = (2 * kappa_c + gamma_c * k4sq) / gamma_c
    M[5, 0] = 1j * kappa_c * k4 / gamma_c
    return M


def block_AE(x, k4):
    """Compression + compact displacement + propagation-dir mixed rotation.
    State: (u₁, u₁', u₄, φ₁₄, u₄', φ₁₄')
    At k₄=0: (u₁,u₁') = Sector A; (u₄,φ₁₄,u₄',φ₁₄') = Sector E."""
    V = V_PN(x)
    k4sq = k4**2
    M = np.zeros((6, 6), dtype=complex)
    # u₁' = u₁'
    M[0, 1] = 1.0
    # u₁'' = [V + μ_tot k₄²] u₁/M_long - (λ+μ)ik₄ u₄'/M_long + iκ_c k₄ φ₁₄/M_long
    M[1, 0] = (V + mu_tot * k4sq) / M_long
    M[1, 4] = -(lam + mu) * 1j * k4 / M_long
    M[1, 3] = 1j * kappa_c * k4 / M_long
    # u₄' = u₄'
    M[2, 4] = 1.0
    # φ₁₄' = φ₁₄'
    M[3, 5] = 1.0
    # u₄'' = [V + M_long k₄²] u₄/μ_tot + κ_c φ₁₄'/μ_tot - (λ+μ)ik₄ u₁'/μ_tot
    M[4, 2] = (V + M_long * k4sq) / mu_tot
    M[4, 5] = kappa_c / mu_tot
    M[4, 1] = -(lam + mu) * 1j * k4 / mu_tot
    # φ₁₄'' = [2κ_c + γk₄²] φ₁₄/γ - κ_c u₄'/γ + iκ_c k₄ u₁/γ
    M[5, 3] = (2 * kappa_c + gamma_c * k4sq) / gamma_c
    M[5, 4] = -kappa_c / gamma_c
    M[5, 0] = 1j * kappa_c * k4 / gamma_c
    return M


def block_Bt(x, k4):
    """Axial torsion. State: (φ₁, φ₁'). No new coupling at k₄≠0."""
    M = np.zeros((2, 2), dtype=complex)
    M[0, 1] = 1.0
    M[1, 0] = (2 * kappa_c + gamma_c * k4**2) / gamma_c
    return M


# ═══════════════════════════════════════════
# INTEGRATOR
# ═══════════════════════════════════════════

def integrate_block(matrix_func, dim, k4, n_steps=2000):
    """Integrate dΨ/dx = M(x,k₄)Ψ through one period."""
    def rhs(x, Psi_flat):
        Psi = Psi_flat.reshape(dim, dim)
        return (matrix_func(x, k4) @ Psi).flatten()

    Psi0 = np.eye(dim, dtype=complex).flatten()
    sol = solve_ivp(rhs, [0, ell], Psi0, method='DOP853',
                    rtol=1e-12, atol=1e-14, max_step=ell/n_steps,
                    t_eval=[ell])
    if not sol.success:
        raise RuntimeError(f"Integration failed: {sol.message}")
    return sol.y[:, -1].reshape(dim, dim)


def get_evanescent(T):
    """Return sorted list of |λ| < 1 eigenvalues."""
    evals = np.linalg.eigvals(T)
    return sorted([abs(e) for e in evals if abs(e) < 0.999])


# ═══════════════════════════════════════════
# STEP 1: CALIBRATE V₀
# ═══════════════════════════════════════════

print("=" * 72)
print("  D4 TRANSFER MATRIX (BLOCK-DIAGONAL, k₄-DEPENDENT)")
print("=" * 72)
print(f"\n  N² = 1/π = {N2:.6f},  κ_c = {kappa_c:.6f},  μ_tot = {mu_tot:.6f}")
print(f"  γ = {gamma_c:.6f},  K/μ = {K_over_mu:.0e},  √(2κ_c/γ) = {np.sqrt(q2):.6f}")

print("\n  Calibrating V₀ from block CG at k₄=0...")

def get_lam_CG(V0_val, k4_val=0.0):
    global V0
    V0 = V0_val
    T = integrate_block(block_CG, 6, k4_val)
    return min(abs(e) for e in np.linalg.eigvals(T))

# Coarse scan
best_V0, best_err = 1.0, 1.0
for lv in np.linspace(1, 6, 40):
    try:
        lam_val = get_lam_CG(np.exp(lv))
        err = abs(lam_val - alpha_target)
        if err < best_err:
            best_err = err
            best_V0 = np.exp(lv)
    except Exception:
        pass

# Binary search
lo, hi = best_V0 * 0.7, best_V0 * 1.3
for _ in range(50):
    mid = (lo + hi) / 2
    lam_val = get_lam_CG(mid)
    if lam_val < alpha_target:
        hi = mid
    else:
        lo = mid

V0 = (lo + hi) / 2
lam_cal = get_lam_CG(V0)
print(f"  V₀ = {V0:.6f},  |λ_CG| = {lam_cal:.10f},  α = {alpha_target:.10f}")
print(f"  Residual: {abs(lam_cal - alpha_target):.2e}")


# ═══════════════════════════════════════════
# STEP 2: k₄ = 0 VERIFICATION
# ═══════════════════════════════════════════

print(f"\n{'=' * 72}")
print("  STEP 2: k₄ = 0 — verify block structure")
print(f"{'=' * 72}")

for name, func, dim in [("CG (anti-plane → α)", block_CG, 6),
                          ("DF (in-plane → α)", block_DF, 6),
                          ("AE (compression)", block_AE, 6),
                          ("Bt (torsion)", block_Bt, 2)]:
    T = integrate_block(func, dim, k4=0.0)
    evals = np.linalg.eigvals(T)
    evals_sorted = sorted(evals, key=lambda x: abs(x))
    evan = [abs(e) for e in evals_sorted if abs(e) < 0.999]
    print(f"\n  Block {name}:")
    for i, e in enumerate(evals_sorted):
        tag = ""
        if abs(e) < 0.999:
            tag = f"  ← tunnelling"
        elif abs(abs(e) - 1.0) < 0.01:
            tag = f"  (propagating)"
        print(f"    λ_{i} = {abs(e):14.8f}  (phase {np.angle(e)/np.pi:+.4f}π){tag}")


# ═══════════════════════════════════════════
# STEP 3: EIGENVALUE FLOW WITH k₄
# ═══════════════════════════════════════════

print(f"\n{'=' * 72}")
print("  STEP 3: Eigenvalue flow — CG block (anti-plane, gives α)")
print(f"{'=' * 72}")

k4_max = 2 * np.pi / (3 * ell)
n_k4 = 15
k4_values = np.linspace(0, k4_max, n_k4)

cg_flow = []
df_flow = []

print(f"\n  {'k₄ℓ':>7}  {'|λ₁|_CG':>14}  {'|λ₂|_CG':>14}  {'|λ₁|/α':>11}  {'Δα/α %':>9}  {'|λ₁|_DF':>14}  {'CG-DF':>10}")
print(f"  {'─'*7}  {'─'*14}  {'─'*14}  {'─'*11}  {'─'*9}  {'─'*14}  {'─'*10}")

for k4 in k4_values:
    try:
        # CG block
        T_cg = integrate_block(block_CG, 6, k4)
        ev_cg = sorted([abs(e) for e in np.linalg.eigvals(T_cg)])
        evan_cg = [e for e in ev_cg if e < 0.999]

        # DF block
        T_df = integrate_block(block_DF, 6, k4)
        ev_df = sorted([abs(e) for e in np.linalg.eigvals(T_df)])
        evan_df = [e for e in ev_df if e < 0.999]

        cg_flow.append((k4, evan_cg))
        df_flow.append((k4, evan_df))

        if evan_cg and evan_df:
            l1_cg = evan_cg[0]
            l2_cg = evan_cg[1] if len(evan_cg) > 1 else float('nan')
            l1_df = evan_df[0]
            shift = (l1_cg - alpha_target) / alpha_target * 100
            split = abs(l1_cg - l1_df) / alpha_target
            print(f"  {k4*ell:7.4f}  {l1_cg:14.10f}  {l2_cg:14.10f}  "
                  f"{l1_cg/alpha_target:11.8f}  {shift:+9.5f}  "
                  f"{l1_df:14.10f}  {split:10.2e}")
        else:
            print(f"  {k4*ell:7.4f}  (insufficient evanescent modes)")
    except Exception as ex:
        print(f"  {k4*ell:7.4f}  ERROR: {ex}")
        cg_flow.append((k4, []))
        df_flow.append((k4, []))


# ═══════════════════════════════════════════
# STEP 4: COMPRESSION SECTOR FLOW
# ═══════════════════════════════════════════

print(f"\n{'=' * 72}")
print("  STEP 4: Compression block (AE) flow")
print(f"{'=' * 72}")
print(f"\n  {'k₄ℓ':>7}  {'|λ₁|':>14}  {'|λ₂|':>14}  {'n_evan':>7}")
print(f"  {'─'*7}  {'─'*14}  {'─'*14}  {'─'*7}")

for k4 in k4_values[::3]:  # every 3rd point
    try:
        T_ae = integrate_block(block_AE, 6, k4)
        ev_ae = sorted([abs(e) for e in np.linalg.eigvals(T_ae)])
        evan_ae = [e for e in ev_ae if e < 0.999]
        if evan_ae:
            print(f"  {k4*ell:7.4f}  {evan_ae[0]:14.8f}  "
                  f"{evan_ae[1] if len(evan_ae)>1 else float('nan'):14.8f}  "
                  f"{len(evan_ae):7d}")
        else:
            print(f"  {k4*ell:7.4f}  (barrier too high — no tunnelling)")
    except Exception as ex:
        print(f"  {k4*ell:7.4f}  ERROR: {ex}")


# ═══════════════════════════════════════════
# STEP 5: SUMMARY
# ═══════════════════════════════════════════

print(f"\n{'=' * 72}")
print("  SUMMARY")
print(f"{'=' * 72}")

if cg_flow and cg_flow[0][1] and cg_flow[-1][1]:
    l0 = cg_flow[0][1][0]
    lk = cg_flow[-1][1][0]
    shift = (lk - l0) / l0 * 100

    print(f"""
  EM tunnelling amplitude (block CG):
    k₄ = 0:          |λ| = {l0:.10f}  (calibrated to α)
    k₄ = 2π/(3ℓ):    |λ| = {lk:.10f}
    Shift:            {shift:+.6f}%
""")

    if abs(shift) < 0.001:
        print("  → α is INSENSITIVE to the compact direction at this level.")
    elif abs(shift) < 1.0:
        print(f"  → α shifts by {abs(shift):.4f}% at the KK scale.")
        print(f"     This is a {abs(shift)/100:.1e} KK correction to the fine structure constant.")
    else:
        print(f"  → α shifts by {abs(shift):.2f}% — significant KK correction!")

if cg_flow and df_flow and cg_flow[-1][1] and df_flow[-1][1]:
    lcg = cg_flow[-1][1][0]
    ldf = df_flow[-1][1][0]
    split = abs(lcg - ldf) / alpha_target * 100
    print(f"\n  CG vs DF splitting at k₄ = 2π/(3ℓ):")
    print(f"    |λ_CG| = {lcg:.10f}")
    print(f"    |λ_DF| = {ldf:.10f}")
    print(f"    Splitting: {split:.6f}% of α")

    if split < 0.001:
        print("    → The two EM polarisations remain DEGENERATE.")
        print("       No chirality breaking at this level.")
    else:
        print(f"    → CHIRALITY BREAKING: the two polarisations split by {split:.4f}%.")
        print(f"       This breaks the SU(2)_L × SU(2)_R symmetry.")
