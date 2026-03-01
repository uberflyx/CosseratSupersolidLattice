#!/usr/bin/env python3
"""
pn_variational.py — Peierls–Nabarro variational calculation for α and G
=========================================================================

Reproduces the calculation in:
  M. A. Cox, "The fine-structure constant and Newton's constant from
  Peierls–Nabarro tunnelling in an FCC Cosserat lattice" (PRL, 2026).

Physics overview
----------------
A screw dislocation in an FCC Cosserat lattice tunnels through the
Peierls–Nabarro (PN) barrier.  The tunnelling amplitude

    α_PN = exp(−2π w/ℓ)

where w is the equilibrium core half-width, plays the role of the QED
coupling constant.

The core width is set by a variational competition:
  - Elastic self-energy (wants a wide, diffuse core)
  - Lattice misfit potential (wants a narrow, localised core)

In a Cosserat solid the elastic kernel has both shear (μ) and micropolar
(κ) contributions.  The coupling number N² = κ/(2μ+κ) controls the
rotational channel.  Rolling without slip fixes N² = 1/π.

Sections
--------
  1. FCC Born–Huang homogenisation → Cosserat moduli
  2. Anti-plane Cosserat kernel in k-space (exact: Si, Ci functions)
  3. PN equilibrium: solve for core width w(N²)
  4. Self-consistent perturbation series for α⁻¹ (Eq. 2 of paper)
  5. Convergence analysis
  6. Skyrmion numbers (gauge-mode correction)
  7. N² variational scan & sensitivity
  8. Compression channel → Newton's constant G

Usage
-----
    python3 pn_variational.py

Requirements: numpy, scipy (standard Anaconda / pip install).

Author:  Mitchell A. Cox
Date:    February 2026
"""

import numpy as np
from scipy.special import sici
from scipy.optimize import brentq
from scipy.integrate import quad


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  Constants
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# CODATA 2022 — Rev. Mod. Phys. 97, 025002 (2025)
ALPHA_INV_EXP = 137.035_999_177       # α⁻¹
ALPHA_UNC     = 0.000_000_021         # 1σ uncertainty
G_CODATA      = 6.67430e-11           # m³ kg⁻¹ s⁻² (22 ppm)
G_UNC         = 0.00015e-11

# SI constants (2019 exact + CODATA)
HBAR = 1.054_571_817e-34   # J·s
C    = 299_792_458.0        # m/s
M_E  = 9.109_383_7015e-31   # kg

# FCC {111} slip geometry
D_OVER_ELL = 1.0 / np.sqrt(3.0)   # inter-planar spacing / lattice spacing


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  SECTION 1 — FCC Born–Huang homogenisation
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
#  The 12 FCC nearest-neighbour bond vectors are permutations of
#  (1/√2)(±1, ±1, 0).  Each contact has stiffness matrix
#      K_ij = k_n n̂_i n̂_j + k_t (δ_ij − n̂_i n̂_j)
#
#  Born–Huang homogenisation (Born & Huang 1954; Suiker et al. 2001)
#  yields the continuum Cosserat moduli:
#      μ     = √2 k_n (1 − r)        shear modulus
#      κ     = 4√2 k_n r              micropolar modulus
#      N²    = 2r / (1 + 3r)          coupling number
#      γ     = μ ℓ²                   curvature modulus  (ℓ_c = ℓ/2)
#
#  where r = k_t/k_n.  Stability requires r < 1 (i.e. μ > 0).
#
#  The Cosserat characteristic length ℓ_c = ℓ/2 (sphere radius) gives
#  γ = 4μ (ℓ/2)² = μ ℓ².

def cosserat_params(N2, mu=1.0, ell_c=0.5):
    """
    Compute all anti-plane Cosserat kernel parameters from N².

    The anti-plane screw dislocation problem couples displacement u_z
    and microrotation φ_y.  Eliminating φ from the coupled equilibrium
    equations (Eringen 1999, Ch. 5) gives a k-dependent effective kernel

        K̂(k) = (μ_tot/2) · k · (k² + p²) / (k² + q²)

    where μ_tot = μ + κ and p², q² are characteristic wavenumbers
    encoding the Cosserat dispersion.

    Parameters
    ----------
    N2 : float
        Cosserat coupling number κ/(2μ+κ).
    mu : float
        Shear modulus (default 1, dimensionless units).
    ell_c : float
        Cosserat characteristic length in units of ℓ (default 0.5).

    Returns
    -------
    dict with: kappa, gamma, mu_bar, mu_tot, p2, q2, p, q
    """
    if N2 < 1e-15:
        return {"kappa": 0, "gamma": 4*mu*ell_c**2,
                "mu_bar": mu, "mu_tot": mu,
                "p2": 0, "q2": 0, "p": 0, "q": 0}

    kappa  = 2 * N2 * mu / (1 - 2*N2)
    gamma  = 4 * mu * ell_c**2
    mu_bar = mu + kappa/2
    mu_tot = mu + kappa
    p2     = kappa * (2*mu + kappa) / (mu_tot * gamma)
    q2     = mu_tot * p2 / mu_bar
    return {"kappa": kappa, "gamma": gamma,
            "mu_bar": mu_bar, "mu_tot": mu_tot,
            "p2": p2, "q2": q2,
            "p": np.sqrt(p2), "q": np.sqrt(q2)}


def N2_from_r(r):
    """Coupling number from stiffness ratio: N² = 2r/(1+3r)."""
    return 2*r / (1 + 3*r)


def r_from_N2(N2):
    """Stiffness ratio from coupling number: r = N²/(2−3N²)."""
    return N2 / (2 - 3*N2)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  SECTION 2 — Anti-plane Cosserat kernel (exact integral via Si, Ci)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
#  The PN equilibrium condition equates the elastic restoring stress at
#  the core to the maximum slope of the sinusoidal lattice potential.
#
#  For an arctan core u(x) = (b/2π) arctan(x/w), the strain profile is
#  Lorentzian: ε(x) = (b/2π) · w/(x²+w²), with Fourier transform
#  ε̂(k) = (b/2) exp(−|k|w).
#
#  The elastic stress in Fourier space is σ̂(k) = K̂(k) ε̂(k), and the
#  equilibrium integral becomes
#
#      LHS(w) = ∫₀^∞ K̂(k) exp(−2kw) dk
#
#  where K̂(k) = (μ_tot/2) k(k²+p²)/(k²+q²).
#
#  This splits as:
#      LHS = (μ_tot/2) [ I₁ + (p²−q²) I₂ ]
#
#  with I₁ = ∫k exp(−2kw) dk = 1/(4w²)  and
#       I₂ = ∫k/(k²+q²) exp(−2kw) dk   (evaluated via Si, Ci).

def _I2_exact(q, w):
    """
    Evaluate ∫₀^∞ k/(k²+q²) exp(−2kw) dk via Si, Ci special functions.

    This is a standard Laplace-transform result.  For z = 2qw:
        I₂ = −cos(z) Ci(z) + sin(z)(π/2 − Si(z))
    """
    if q < 1e-15:
        return np.inf
    z = 2.0 * q * w
    Si_z, Ci_z = sici(z)
    return -np.cos(z) * Ci_z + np.sin(z) * (np.pi/2 - Si_z)


def kernel_integral(w, N2, mu=1.0, ell_c=0.5):
    """
    Analytical evaluation of the PN equilibrium LHS integral:
        LHS(w) = ∫₀^∞ K̂(k) exp(−2kw) dk
    using the Si/Ci closed form.
    """
    if N2 < 1e-15:
        return mu / (8 * w**2)

    P     = cosserat_params(N2, mu, ell_c)
    Delta = P["p2"] - P["q2"]
    I1    = 1.0 / (4 * w**2)
    I2    = _I2_exact(P["q"], w)
    return (P["mu_tot"] / 2) * (I1 + Delta * I2)


def kernel_integral_numerical(w, N2, mu=1.0, ell_c=0.5):
    """
    Numerical verification of the kernel integral via direct quadrature.
    Use to cross-check the analytical Si/Ci result.
    """
    if N2 < 1e-15:
        return mu / (8 * w**2)

    P = cosserat_params(N2, mu, ell_c)

    def integrand(k):
        if k < 1e-15:
            return 0.0
        return (P["mu_tot"]/2) * k * (k**2 + P["p2"]) / (k**2 + P["q2"]) \
               * np.exp(-2*k*w)

    result, _ = quad(integrand, 0, max(50/w, 500), limit=500, epsrel=1e-12)
    return result


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  SECTION 3 — Solve PN equilibrium for the core width w
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
#  The PN equilibrium equates the elastic kernel integral to the maximum
#  restoring force of the sinusoidal lattice potential:
#
#      LHS(w) = μ b / (4d)
#
#  where b = ℓ (Burgers vector) and d = ℓ/√3 (FCC {111} spacing).

def find_w(N2, mu=1.0, ell=1.0, ell_c=0.5, d_over_ell=D_OVER_ELL):
    """
    Solve the PN equilibrium for the core half-width w.
    Returns w/ℓ (dimensionless).
    """
    b   = ell
    d   = d_over_ell * ell
    RHS = mu * b / (4 * d)

    def residual(w):
        return kernel_integral(w, N2, mu, ell_c) - RHS

    try:
        w_sol = brentq(residual, 0.01*ell, 20*ell, xtol=1e-15, rtol=1e-15)
        return w_sol / ell
    except ValueError:
        return np.nan


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  SECTION 4 — Self-consistent perturbation series for α⁻¹
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
#  Eq. (2) of the paper:
#
#      α⁻¹ = α₀⁻¹ − 2 − α − α/[π(1−α)] − 6α³/π²
#
#  Each term:
#      α₀⁻¹            bare PN tunnelling
#      −2               gauge-mode (skyrmion number, topological)
#      −α               translational self-energy
#      −α/[π(1−α)]     rotational self-energy (Dyson-resummed)
#      −6α³/π²          inter-valley scattering (6 NN × N⁴)

def alpha_inv_rhs(alpha, bare_inv):
    """Right-hand side of the self-consistent equation."""
    gauge   = 2.0
    trans   = alpha
    rot     = alpha / (np.pi * (1.0 - alpha))
    ivalley = 6.0 * alpha**3 / np.pi**2
    return bare_inv - gauge - trans - rot - ivalley


def solve_alpha(bare_inv):
    """Solve α⁻¹ = f(α) by Brent root-finding on g(α) = 1/α − f(α) = 0."""
    def residual(alpha):
        return 1.0/alpha - alpha_inv_rhs(alpha, bare_inv)

    return brentq(residual, 1e-8, 0.5, xtol=1e-15, rtol=1e-15)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  SECTION 5 — Convergence analysis
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def convergence_table(alpha, bare_inv):
    """
    Show how each correction reduces the residual.
    Returns list of (label, α⁻¹ partial sum, residual in ppb).
    """
    a = alpha
    stages = [
        ("Bare: exp(2πw/ℓ)",                    bare_inv),
        ("− Gauge (−2)",                         bare_inv - 2),
        ("− Translational (−α)",                 bare_inv - 2 - a),
        ("− Rotational leading (−α/π)",          bare_inv - 2 - a - a/np.pi),
        ("− Rotational resummed (−α/[π(1−α)])", bare_inv - 2 - a
                                                  - a/(np.pi*(1-a))),
        ("− Inter-valley (−6α³/π²)",             bare_inv - 2 - a
                                                  - a/(np.pi*(1-a))
                                                  - 6*a**3/np.pi**2),
    ]
    rows = []
    for label, val in stages:
        ppb = (val - ALPHA_INV_EXP) / ALPHA_INV_EXP * 1e9
        rows.append((label, val, ppb))
    return rows


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  SECTION 6 — Skyrmion number (gauge-mode correction)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def skyrmion_number(with_pn_pinning=True):
    """Topological charge per polarisation channel (analytic formula)."""
    m = 1
    cos_0   = np.cos(np.pi)
    cos_inf = np.cos(0.0) if with_pn_pinning else np.cos(np.pi/2)
    return -m/2 * (cos_inf - cos_0)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  SECTION 7 — Compression channel: Newton's constant G
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_G(alpha):
    """
    Newton's constant from the 19-node FCC cluster tunnelling amplitude.

    Cluster: 1 central + 12 NN + 6 NNN = 19 nodes.
    Purely-elastic-mass: ℓ = r_e, m₀ = m_e/α.
    G = (1+1/π)(1−17α/18) × (ℏc/m₀²) × α^19
    """
    m_0       = M_E / alpha
    cosserat  = 1.0 + 1.0/np.pi
    trans_cor = 1.0 - 17.0*alpha/18.0
    prefactor = cosserat * trans_cor * HBAR * C / m_0**2
    return prefactor * alpha**19


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  MAIN
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    sep  = "=" * 72
    sep2 = "-" * 72

    print(sep)
    print("  Peierls-Nabarro variational calculation for alpha and G")
    print("  PRL: Cox (2026)")
    print(sep)

    mu, ell, ell_c = 1.0, 1.0, 0.5

    # ── STAGE 1: Verify analytical kernel vs quadrature ──────────────

    print(f"\nSTAGE 1: Kernel verification (analytical Si,Ci vs quadrature)")
    print(sep2)

    w_test = np.pi / 4
    print(f"  Test point: w/ell = pi/4 = {w_test:.8f}\n")
    print(f"  {'N^2':>10s}  {'I_analytical':>16s}  {'I_numerical':>16s}  {'rel error':>12s}")
    print(f"  {'-'*58}")

    for N2 in [0.0, 0.1, 0.2, 1/np.pi, 1/3, 0.4, 0.49]:
        I_ana = kernel_integral(w_test, N2, mu, ell_c)
        I_num = kernel_integral_numerical(w_test, N2, mu, ell_c)
        rel   = abs(I_ana - I_num) / abs(I_num) if abs(I_num) > 0 else 0
        tag   = " <-- 1/pi" if abs(N2 - 1/np.pi) < 0.002 else ""
        print(f"  {N2:10.6f}  {I_ana:16.10f}  {I_num:16.10f}  {rel:12.2e}{tag}")

    # ── STAGE 2: Equilibrium core width w(N^2) ───────────────────────

    print(f"\n\nSTAGE 2: Equilibrium core width w/ell vs N^2")
    print(sep2)

    RHS = mu * ell / (4 * D_OVER_ELL * ell)
    print(f"  RHS = mu*b/(4*d) = {RHS:.8f}")
    print(f"  d/ell = 1/sqrt(3) = {D_OVER_ELL:.8f}\n")
    print(f"  {'N^2':>10s}  {'w/ell':>12s}  {'alpha_0^-1':>14s}  {'alpha_0^-1 - 2':>16s}")
    print(f"  {'-'*56}")

    for N2 in [0.001, 0.10, 0.20, 0.25, 0.30, 1/np.pi, 1/3, 0.35, 0.40, 0.45, 0.49]:
        w = find_w(N2, mu, ell, ell_c)
        if np.isnan(w):
            continue
        a0_inv = np.exp(2*np.pi*w)
        tag = " <-- 1/pi" if abs(N2 - 1/np.pi) < 0.002 else ""
        print(f"  {N2:10.6f}  {w:12.8f}  {a0_inv:14.4f}  {a0_inv-2:16.4f}{tag}")

    # ── STAGE 3: N^2 that gives w/ell = pi/4 exactly ────────────────

    print(f"\n\nSTAGE 3: Variational optimum N^2* such that w/ell = pi/4")
    print(sep2)

    target_w = np.pi / 4

    def _res_N2(N2):
        return find_w(N2, mu, ell, ell_c) - target_w

    N2_star = brentq(_res_N2, 0.1, 0.49, xtol=1e-14, rtol=1e-14)
    w_star  = find_w(N2_star, mu, ell, ell_c)
    r_star  = r_from_N2(N2_star)

    print(f"  N^2*          = {N2_star:.10f}")
    print(f"  1/pi          = {1/np.pi:.10f}")
    print(f"  Offset        = {(N2_star - 1/np.pi)/(1/np.pi)*100:+.4f}%")
    print(f"  w*/ell        = {w_star:.10f}")
    print(f"  pi/4          = {np.pi/4:.10f}")
    print(f"  r* = k_t/k_n  = {r_star:.8f}")

    # At N^2 = 1/pi exactly
    w_rolling    = find_w(1/np.pi, mu, ell, ell_c)
    bare_rolling = np.exp(2*np.pi*w_rolling)
    bare_analytic = np.exp(np.pi**2/2)

    print(f"\n  At N^2 = 1/pi exactly:")
    print(f"    w/ell       = {w_rolling:.10f}")
    print(f"    pi/4        = {np.pi/4:.10f}")
    print(f"    deviation   = {(w_rolling - np.pi/4)/(np.pi/4)*100:+.4f}%")
    print(f"    alpha_0^-1  = {bare_rolling:.8f}")
    print(f"    exp(pi^2/2) = {bare_analytic:.8f}")
    print(f"    deviation   = {(bare_rolling - bare_analytic)/bare_analytic*100:+.6f}%")

    # ── STAGE 4: Self-consistent solution ────────────────────────────

    print(f"\n\nSTAGE 4: Self-consistent perturbation series")
    print(sep2)

    # (a) Variational bare value
    alpha_var   = solve_alpha(bare_rolling)
    inv_var     = 1.0 / alpha_var
    res_var_ppb = (inv_var - ALPHA_INV_EXP) / ALPHA_INV_EXP * 1e9
    res_var_sig = (inv_var - ALPHA_INV_EXP) / ALPHA_UNC

    # (b) Analytic bare
    alpha_an    = solve_alpha(bare_analytic)
    inv_an      = 1.0 / alpha_an
    res_an_ppb  = (inv_an - ALPHA_INV_EXP) / ALPHA_INV_EXP * 1e9
    res_an_sig  = (inv_an - ALPHA_INV_EXP) / ALPHA_UNC

    print(f"  (a) Using variational bare (N^2=1/pi -> w/ell={w_rolling:.8f}):")
    print(f"      alpha_0^-1  = {bare_rolling:.8f}")
    print(f"      alpha_PN^-1 = {inv_var:.12f}")
    print(f"      CODATA      = {ALPHA_INV_EXP:.9f} +/- {ALPHA_UNC}")
    print(f"      Residual    = {res_var_ppb:+.3f} ppb  ({res_var_sig:+.3f} sigma)")

    print(f"\n  (b) Using analytic bare exp(pi^2/2):")
    print(f"      alpha_0^-1  = {bare_analytic:.8f}")
    print(f"      alpha_PN^-1 = {inv_an:.12f}")
    print(f"      CODATA      = {ALPHA_INV_EXP:.9f} +/- {ALPHA_UNC}")
    print(f"      Residual    = {res_an_ppb:+.3f} ppb  ({res_an_sig:+.3f} sigma)")

    # ── STAGE 5: Convergence ─────────────────────────────────────────

    print(f"\n\nSTAGE 5: Convergence of perturbation series")
    print(sep2)

    print(f"  {'Term':<42s} {'alpha^-1':>18s} {'Residual':>14s}")
    print(f"  {'-'*76}")
    for label, val, ppb in convergence_table(alpha_an, bare_analytic):
        print(f"  {label:<42s} {val:>18.9f} {ppb:>+13.3f} ppb")

    # ── STAGE 6: Skyrmion numbers ────────────────────────────────────

    print(f"\n\nSTAGE 6: Topological gauge correction (skyrmion numbers)")
    print(sep2)

    Q_sky = skyrmion_number(with_pn_pinning=True)
    Q_mer = skyrmion_number(with_pn_pinning=False)

    print(f"  With PN pinning (skyrmion):   Q = {Q_sky:.1f}  per channel")
    print(f"  Without pinning (meron):      Q = {Q_mer:.1f}  per channel")
    print(f"  Gauge correction = 2 x |Q|    = {2*abs(Q_sky):.1f}")

    # ── STAGE 7: N^2 sensitivity scan ────────────────────────────────

    print(f"\n\nSTAGE 7: N^2 sensitivity scan")
    print(sep2)

    print(f"  {'N^2':>10s}  {'w/ell':>10s}  {'bare^-1':>12s}"
          f"  {'alpha_PN^-1':>18s}  {'Residual':>14s}")
    print(f"  {'-'*68}")

    for N2_trial in np.linspace(0.20, 0.45, 26):
        r_t = r_from_N2(N2_trial)
        if r_t <= 0 or r_t >= 1:
            continue
        w_t = find_w(N2_trial, mu, ell, ell_c)
        if np.isnan(w_t):
            continue
        bare_t = np.exp(2*np.pi*w_t)
        try:
            a_t   = solve_alpha(bare_t)
            inv_t = 1.0/a_t
        except ValueError:
            continue
        ppb = (inv_t - ALPHA_INV_EXP) / ALPHA_INV_EXP * 1e9
        tag = " <-- 1/pi" if abs(N2_trial - 1/np.pi) < 0.006 else ""
        print(f"  {N2_trial:>10.5f}  {w_t:>10.6f}  {bare_t:>12.4f}"
              f"  {inv_t:>18.6f}  {ppb:>+13.1f} ppb{tag}")

    # Sensitivity derivative
    dN2 = 1e-7
    w_p = find_w(N2_star + dN2, mu, ell, ell_c)
    w_m = find_w(N2_star - dN2, mu, ell, ell_c)
    dw_dN2 = (w_p - w_m) / (2*dN2)
    dalpha_inv_dN2 = 2*np.pi * np.exp(2*np.pi*target_w) * dw_dN2

    print(f"\n  Sensitivity at N^2*:")
    print(f"    dw/dN^2            = {dw_dN2:.6f}")
    print(f"    d(alpha^-1)/d(N^2) = {dalpha_inv_dN2:.2f}")
    print(f"    1% error in N^2 -> delta(alpha^-1) = "
          f"{0.01*N2_star*dalpha_inv_dN2:.2f}"
          f"  ({0.01*N2_star*dalpha_inv_dN2/ALPHA_INV_EXP*100:.2f}%)")

    # ── STAGE 8: Compression → G ────────────────────────────────────

    print(f"\n\nSTAGE 8: Compression channel -> Newton's constant")
    print(sep2)

    G_calc = compute_G(alpha_an)
    G_ppm  = (G_calc - G_CODATA) / G_CODATA * 1e6

    print(f"  Cluster: 1 + 12 + 6 = 19 nodes")
    print(f"  alpha^19        = {alpha_an**19:.6e}")
    print(f"  Node mass m_0   = m_e/alpha = {M_E/alpha_an:.4e} kg"
          f"  ({M_E/alpha_an * C**2 / 1.602e-13:.1f} MeV/c^2)")
    print(f"  Cosserat factor = 1 + 1/pi = {1+1/np.pi:.6f}")
    print(f"  Trans. corr.    = 1 - 17alpha/18 = {1-17*alpha_an/18:.6f}")
    print(f"\n  G_PN    = {G_calc:.6e} m^3 kg^-1 s^-2")
    print(f"  G_COD   = {G_CODATA:.6e} +/- {G_UNC:.2e}")
    print(f"  Offset  = {G_ppm:+.1f} ppm  (expt. unc. = 22 ppm)")

    # Hierarchy
    alpha_G   = G_calc * M_E**2 / (HBAR * C)
    hierarchy = alpha_an / alpha_G

    print(f"\n  alpha_PN  = {alpha_an:.8e}")
    print(f"  alpha_G   = {alpha_G:.8e}")
    print(f"  Ratio     = {hierarchy:.3e}")
    print(f"  alpha^18  = {alpha_an**18:.3e}")

    # ── SUMMARY ──────────────────────────────────────────────────────

    print(f"\n\n{sep}")
    print(f"{'SUMMARY':^72s}")
    print(sep)
    print(f"  Rolling constraint:  N^2 = 1/pi = {1/np.pi:.8f}")
    print(f"  Variational optimum: N^2*       = {N2_star:.8f}"
          f"  (offset {(N2_star-1/np.pi)/(1/np.pi)*100:+.3f}%)")
    print(f"  Core width at 1/pi: w/ell       = {w_rolling:.8f}"
          f"  (pi/4 = {np.pi/4:.8f})")
    print(f"  Bare amplitude:     alpha_0^-1  = {bare_rolling:.6f}"
          f"  (exp(pi^2/2) = {bare_analytic:.6f})")
    print(f"\n  alpha_PN^-1  = {inv_an:.12f}")
    print(f"  CODATA       = {ALPHA_INV_EXP:.9f} +/- {ALPHA_UNC}")
    print(f"  Residual     = {res_an_ppb:+.3f} ppb  ({res_an_sig:+.3f} sigma)")
    print(f"\n  G_PN         = {G_calc:.4e} m^3 kg^-1 s^-2")
    print(f"  G_CODATA     = {G_CODATA:.4e} +/- {G_UNC:.2e}")
    print(f"  Offset       = {G_ppm:+.1f} ppm")
    print(sep)


if __name__ == "__main__":
    main()
