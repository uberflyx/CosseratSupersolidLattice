#!/usr/bin/env python3
"""
gamma4_final_v2.py  —  γ₄ from quantum dressing: γ₃ = -7
===================================================================

Two independent methods:
  A) Analytic moment expansion (valid for σ/d ≲ 0.3)
  B) Numerical quadrature with polynomial derivative extraction

Mitchell A. Cox, 2026
"""
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from numpy.polynomial.laguerre import laggauss

# ═══════════════════════════════════════════════════════════════════
# CORRECTED MORSE PARAMETERS
# ═══════════════════════════════════════════════════════════════════
a = 7.0/3.0       # alpha (with d = 1)
D = 1.0/(2*a**2)  # = 9/98 ≈ 0.09184

print(f"Morse parameters: a*d = {a:.4f}, D = {D:.6f} = 9/98")
print(f"  gamma_3 = {-3*a:.4f}, gamma_4 = {7*a**2:.4f} = 343/9")
print(f"  gamma_4/gamma_3^2 = {7*a**2/(-3*a)**2:.6f} = 7/9")
print()

# ═══════════════════════════════════════════════════════════════════
# A. ANALYTIC MOMENT EXPANSION
# ═══════════════════════════════════════════════════════════════════

print("="*72)
print("  A. ANALYTIC MOMENT EXPANSION for Morse (a=7/3, D=9/98)")
print("="*72)
print()

def morse_deriv(n):
    """n-th derivative of Morse at r = d = 1."""
    return D * a**n * ((-2)**n - 2*(-1)**n)

V0 = morse_deriv(0)
V1 = morse_deriv(1)
V2 = morse_deriv(2)
V3 = morse_deriv(3)
V4 = morse_deriv(4)
V5 = morse_deriv(5)
V6 = morse_deriv(6)

print(f"  V(1) = {V0:.4f},  V'(1) = {V1:.1f}")
print(f"  V''  = {V2:.4f},  V''' = {V3:.4f},  V'''' = {V4:.4f}")
print(f"  V''''' = {V5:.2f},  V'''''' = {V6:.2f}")
print(f"  gamma_3 = V'''/V'' = {V3/V2:.4f}")
print(f"  gamma_4 = V''''/V'' = {V4/V2:.4f}")
print()

# Updated eta_4 decomposition:
# W4^(V''') = -|gamma_3| * C_3 = -7 * 0.03967 = -0.2777
C_3 = 0.03967  # geometric lattice sum ratio (unchanged)
C_4 = 0.0142   # geometric lattice sum ratio (unchanged)
W4_Vpp = 0.668
W4_Vppp = -7 * C_3  # = -0.2777
eta4_roll = -0.019
eta4_base = W4_Vpp + W4_Vppp + eta4_roll  # = 0.3713

print(f"  eta_4 decomposition (corrected):")
print(f"    W4^(V'')  = {W4_Vpp:+.4f}")
print(f"    W4^(V''') = {W4_Vppp:+.4f}  (was -0.238 with gamma_3=-6)")
print(f"    eta_4^(roll) = {eta4_roll:+.4f}")
print(f"    Subtotal  = {eta4_base:+.4f}  (was +0.411)")
print(f"    C_4       = {C_4}")
print(f"    Koide target: eta_4 = 0.62 => gamma_4_eff = {(0.62-eta4_base)/C_4:.1f}")
print()

print(f"  {'σ/d':>6s} {'V″_eff':>9s} {'V″″_eff':>10s} {'γ₃_eff':>8s}"
      f" {'γ₄_eff':>8s} {'γ₄/γ₃²':>8s} {'η₄':>7s}")
print("  "+"-"*60)

for sig100 in range(0, 36, 2):
    s = sig100/100.0
    s2 = s**2; s4 = s**4
    
    V2_eff = V2 + s2*(0.5*V4 + V2) + s4*(0.125*V6 + 0.5*V4 + 0.75*V2)
    V4_eff = V4 + s2*(0.5*V6 + 6*V4) + s4*(0.125*V6*28 + 3*V6)
    V3_eff = V3 + s2*(0.5*V5 + 3*V3)
    
    g3_eff = V3_eff/V2_eff if abs(V2_eff) > 0.01 else 0
    g4_eff = V4_eff/V2_eff if abs(V2_eff) > 0.01 else 0
    rat = g4_eff/g3_eff**2 if abs(g3_eff) > 0.01 else 0
    eta4 = eta4_base + g4_eff*C_4
    
    if sig100 % 4 == 0:
        mark = " ←" if abs(eta4 - 0.62) < 0.02 else ""
        print(f"  {s:6.3f} {V2_eff:9.4f} {V4_eff:10.3f} {g3_eff:8.3f}"
              f" {g4_eff:8.2f} {rat:8.4f} {eta4:7.3f}{mark}")

# ═══════════════════════════════════════════════════════════════════
# B. NUMERICAL QUADRATURE
# ═══════════════════════════════════════════════════════════════════

def morse_func(r):
    r = np.asarray(r, dtype=float)
    return D*(np.exp(np.clip(-2*a*(r-1), -80, 80))
              - 2*np.exp(np.clip(-a*(r-1), -80, 80)))

def V_dressed_num(R, sigma, nq=50):
    xi, wh = hermegauss(nq)
    tl, wl = laggauss(nq)
    uz = sigma*xi; s = 2*sigma**2*tl
    wz = wh/np.sqrt(2*np.pi)
    R = np.atleast_1d(np.asarray(R, dtype=float))
    Rz = R[:,None,None] + uz[None,:,None]
    r = np.sqrt(Rz**2 + s[None,None,:])
    r = np.maximum(r, 0.02)
    W = wz[None,:,None]*wl[None,None,:]
    return np.sum(morse_func(r)*W, axis=(1,2))

print(f"\n{'='*72}")
print("  B. NUMERICAL QUADRATURE (50-pt Gauss-Hermite×Laguerre)")
print("="*72)
print()

Rg = np.linspace(0.3, 2.5, 601)

print(f"  {'σ/d':>6s} {'d_eff':>7s} {'V″_eff':>9s} {'γ₃_eff':>8s}"
      f" {'γ₄_eff':>8s} {'γ₄/γ₃²':>8s} {'η₄':>7s}")
print("  "+"-"*60)

num_data = []
for sig100 in range(5, 38):
    sig = sig100/100.0
    Ve = V_dressed_num(Rg, sig, nq=50)
    imin = 10+np.argmin(Ve[10:-10])
    hw = min(25, imin, len(Rg)-1-imin)
    if hw < 8: continue
    x = Rg[imin-hw:imin+hw+1] - Rg[imin]
    c = np.polyfit(x, Ve[imin-hw:imin+hw+1], 8)
    dm = Rg[imin]; V2n = 2*c[-3]; V3n = 6*c[-4]; V4n = 24*c[-5]
    if abs(V2n) < 1e-10: continue
    g3n = dm*V3n/V2n; g4n = dm**2*V4n/V2n
    ratn = g4n/g3n**2 if abs(g3n) > 0.01 else 0
    eta4n = eta4_base + g4n*C_4
    num_data.append((sig, dm, V2n, g3n, g4n, ratn, eta4n))
    if sig100 % 5 == 0:
        mark = " ←" if abs(eta4n - 0.62) < 0.02 else ""
        print(f"  {sig:6.3f} {dm:7.4f} {V2n:9.5f} {g3n:8.3f}"
              f" {g4n:8.2f} {ratn:8.4f} {eta4n:7.3f}{mark}")

num_data = np.array(num_data)

# ═══════════════════════════════════════════════════════════════════
# C. COMBINED ANALYSIS
# ═══════════════════════════════════════════════════════════════════

print(f"\n{'='*72}")
print("  C. COMBINED ANALYSIS")
print("="*72)
print()

etas = num_data[:, 6]
g4s  = num_data[:, 4]
sigs = num_data[:, 0]

print("  γ₄_eff trend (numerical, method B):")
for i in range(0, len(num_data), 5):
    print(f"    σ/d = {sigs[i]:.3f}: γ₄ = {g4s[i]:.2f}, η₄ = {etas[i]:.3f}")

if len(num_data) > 10:
    mask = (sigs >= 0.25) & (sigs <= 0.35)
    if np.sum(mask) >= 3:
        s_fit = sigs[mask]
        e_fit = etas[mask]
        p = np.polyfit(s_fit, e_fit, 2)
        p_solve = p.copy(); p_solve[2] -= 0.62
        roots = np.roots(p_solve)
        real_roots = roots[np.isreal(roots)].real
        valid = real_roots[(real_roots > 0.2) & (real_roots < 0.6)]
        if len(valid) > 0:
            sig_extrap = valid[0]
            pg4 = np.polyfit(s_fit, g4s[mask], 2)
            g4_extrap = np.polyval(pg4, sig_extrap)
            print(f"\n  Quadratic extrapolation of η₄(σ) → 0.62:")
            print(f"    σ_Koide / d = {sig_extrap:.3f}")
            print(f"    γ₄_eff (extrapolated) = {g4_extrap:.1f}")
            sig_E = 2**(-7/8)
            print(f"    Jastrow: 1-(σ/σ_E)² = {1-(sig_extrap/sig_E)**2:.3f}")
        else:
            print("\n  No valid crossing found")

# ═══════════════════════════════════════════════════════════════════
# D. MONOGRAPH TABLE (LaTeX)
# ═══════════════════════════════════════════════════════════════════

print(f"\n{'='*72}")
print("  D. MONOGRAPH TABLE (LaTeX)")
print("="*72)
print()

table_sigmas = [0.10, 0.20, 0.25, 0.30, 0.35]
print("\\begin{center}")
print("\\begin{tabular}{cccccc}")
print("\\toprule")
print("$\\sigma/d$ & $d_{\\text{eff}}/d$ & $V''_{\\text{eff}}$ & $\\gamma_{3,\\text{eff}}$ & $\\gamma_{4,\\text{eff}}$ & $\\gamma_4/\\gamma_3^2$ \\\\")
print("\\midrule")
print(f"$0.00$ & $1.000$ & $1.000$ & ${V3/V2:.2f}$ & ${V4/V2:.1f}$ & ${(V4/V2)/(V3/V2)**2:.3f}$ \\\\")

for ts in table_sigmas:
    idx = np.argmin(np.abs(num_data[:,0] - ts))
    sig, dm, v2, g3, g4, rat, eta = num_data[idx]
    print(f"${sig:.2f}$ & ${dm:.3f}$ & ${v2:.3f}$ & ${g3:.2f}$ & ${g4:.1f}$ & ${rat:.3f}$ \\\\")

print("\\bottomrule")
print("\\end{tabular}")
print("\\end{center}")

# ═══════════════════════════════════════════════════════════════════
# E. FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════

print(f"\n{'='*72}")
print("  FINAL SUMMARY")
print("="*72)
print(f"""
  ┌──────────────────────────────────────────────────────────────┐
  │  CORRECTED PARAMETERS (μ' = (3-ξ)/5, ξ = -7):              │
  │  BARE POTENTIAL:   γ₃ = -7,  γ₄ = 343/9 ≈ 38.1            │
  │                    γ₄/γ₃² = 7/9 ≈ 0.778                    │
  │  D = 9/98 ≈ 0.0918,  a·d = 7/3 ≈ 2.333                    │
  ├──────────────────────────────────────────────────────────────┤
  │  η₄ decomposition (corrected):                              │
  │    W₄(V'')  = +0.668   (harmonic, unchanged)                │
  │    W₄(V''') = -0.278   (was -0.238 with γ₃ = -6)           │
  │    η₄(roll) = -0.019   (Cosserat, unchanged)                │
  │    Subtotal = +0.371   (was +0.411)                          │
  │    Koide target η₄ = 0.62 → γ₄_eff = {(0.62-eta4_base)/C_4:.1f}               │
  └──────────────────────────────────────────────────────────────┘
""")
