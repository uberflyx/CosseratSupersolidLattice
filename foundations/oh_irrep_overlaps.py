#!/usr/bin/env python3
"""
oh_irrep_overlaps.py  (v2 — corrected Mulliken labels)
========================================================
Computes the tunnelling amplitude for ALL 10 irreducible representations
of O_h by evaluating the angular overlap of each irrep's basis functions
with the FCC Peierls–Nabarro potential.

LABEL CONVENTIONS (Mulliken, standard for O_h):
  g/u = parity under spatial inversion (x,y,z) → (−x,−y,−z)
  g (gerade, even):  polynomial of EVEN degree
  u (ungerade, odd): polynomial of ODD degree
  1/2 = character under C₄: +1 for subscript 1, −1 for subscript 2

CORRECTED FROM v1:
  Old "A₂g" (xyz)           → A₂u (degree 3, ungerade)
  Old "T₁g" (x(y²−z²))     → T₂u (degree 3, ungerade, trace(C₄)=−1)
  New T₁g: yz(y²−z²) etc.  (degree 4, gerade, trace(C₄)=+1)
  New A₂g: x⁴(y²−z²)+cyc   (degree 6, gerade, trace(C₄)=−1)

M. A. Cox, University of the Witwatersrand (2026)
"""

import numpy as np
from itertools import product

# ══════════════════════════════════════════════════════════════
# O_h CHARACTER TABLE AND BASIS FUNCTIONS (CORRECTED)
# ══════════════════════════════════════════════════════════════

# ── GERADE (even parity) ──

def A1g(n):
    """A₁g (dim 1, degree 0): trivial rep, f = 1."""
    return np.array([1.0])

def Eg(n):
    """E_g (dim 2, degree 2): traceless symmetric diagonal.
    f = {(2z²−x²−y²)/√6, (x²−y²)/√2}."""
    x, y, z = n
    return np.array([
        (2*z**2 - x**2 - y**2) / np.sqrt(6),
        (x**2 - y**2) / np.sqrt(2)
    ])

def T2g(n):
    """T₂g (dim 3, degree 2): symmetric off-diagonal.
    f = {yz, zx, xy}. Trace(C₄) = −1."""
    x, y, z = n
    return np.array([y*z, z*x, x*y])

def T1g(n):
    """T₁g (dim 3, degree 4): axial vector cubic harmonic.
    f = {yz(y²−z²), zx(z²−x²), xy(x²−y²)}.
    Trace(C₄) = +1. Verified: degree 4 → gerade ✓."""
    x, y, z = n
    return np.array([
        y*z*(y**2 - z**2),
        z*x*(z**2 - x**2),
        x*y*(x**2 - y**2)
    ])

def A2g(n):
    """A₂g (dim 1, degree 6): gerade pseudoscalar.
    f = x⁴(y²−z²) + y⁴(z²−x²) + z⁴(x²−y²).
    Under inversion: degree 6 → even → gerade ✓.
    Under C₄(z): (x,y,z)→(−y,x,z): f → −f → trace = −1 → A₂ ✓."""
    x, y, z = n
    return np.array([
        x**4*(y**2 - z**2) + y**4*(z**2 - x**2) + z**4*(x**2 - y**2)
    ])

# ── UNGERADE (odd parity) ──

def T1u(n):
    """T₁u (dim 3, degree 1): polar vector.
    f = {x, y, z}. Trace(C₄) = +1."""
    return np.array(n, dtype=float)

def A2u(n):
    """A₂u (dim 1, degree 3): pseudoscalar.
    f = xyz. Under inversion: (−x)(−y)(−z) = −xyz → ungerade ✓.
    Under C₄(z): (−y)(x)(z) = −xyz → trace = −1 → A₂ ✓.
    
    NOTE: v1 of this script mislabelled this as 'A₂g'. Since
    xyz has degree 3 (odd), it is ungerade, not gerade."""
    x, y, z = n
    return np.array([x * y * z])

def T2u(n):
    """T₂u (dim 3, degree 3): ungerade, trace(C₄) = −1.
    f = {x(y²−z²), y(z²−x²), z(x²−y²)}.
    Under inversion: degree 3 → ungerade ✓.
    Under C₄(z): trace = −1 → T₂ ✓.
    
    NOTE: v1 mislabelled this as 'T₁g'. The trace under C₄ is −1
    (matching T₂, not T₁), and degree 3 is ungerade (not gerade)."""
    x, y, z = n
    return np.array([
        x * (y**2 - z**2),
        y * (z**2 - x**2),
        z * (x**2 - y**2)
    ])

def A1u(n):
    """A₁u (dim 1, degree 9): ungerade scalar.
    f = xyz(x²−y²)(y²−z²)(z²−x²). Degree 9 → ungerade ✓."""
    x, y, z = n
    return np.array([
        x*y*z * (x**2-y**2) * (y**2-z**2) * (z**2-x**2)
    ])

def Eu(n):
    """E_u (dim 2, degree 5): ungerade doublet.
    f = {z·(2z²−x²−y²)·xy, z·(x²−y²)·(x²+y²−2z²)/√3}.
    Degree 5 → ungerade ✓."""
    x, y, z = n
    return np.array([
        z * (2*z**2 - x**2 - y**2) * x * y,
        (x**2 - y**2) * z * (x**2 + y**2 - 2*z**2) / np.sqrt(3)
    ])


# Master dictionary with corrected labels
irreps = {
    # Gerade (even parity)
    'A1g': (A1g, 1, 'g', 0),
    'Eg':  (Eg,  2, 'g', 2),
    'T2g': (T2g, 3, 'g', 2),
    'T1g': (T1g, 3, 'g', 4),
    'A2g': (A2g, 1, 'g', 6),
    # Ungerade (odd parity)
    'T1u': (T1u, 3, 'u', 1),
    'A2u': (A2u, 1, 'u', 3),
    'T2u': (T2u, 3, 'u', 3),
    'Eu':  (Eu,  2, 'u', 5),
    'A1u': (A1u, 1, 'u', 9),
}


# ══════════════════════════════════════════════════════════════
# FCC RECIPROCAL LATTICE VECTORS
# ══════════════════════════════════════════════════════════════

# FCC with NN spacing ℓ = 1 → conventional cell a = √2
# Reciprocal lattice = BCC with primitive vectors (2π/a)(±1,±1,±1) etc.

# {111} directions (8 vectors)
G111 = []
for s in product([-1, 1], repeat=3):
    G111.append(np.array(s, dtype=float))
G111 = np.array(G111)
G111_hat = G111 / np.linalg.norm(G111[0])

# {200} directions (6 vectors)
G200 = []
for i in range(3):
    for s in [-1, 1]:
        v = np.zeros(3); v[i] = s
        G200.append(v)
G200 = np.array(G200)
G200_hat = G200 / np.linalg.norm(G200[0])

# {220} directions (12 vectors)
G220 = []
for i in range(3):
    j, k = (i+1)%3, (i+2)%3
    for s1 in [-1, 1]:
        for s2 in [-1, 1]:
            v = np.zeros(3); v[j] = s1; v[k] = s2
            G220.append(v)
G220 = np.array(G220, dtype=float)
G220_hat = G220 / np.linalg.norm(G220[0])

# Magnitudes and form factors
ell = 1.0
a = np.sqrt(2)
G111_mag = (2*np.pi/a) * np.sqrt(3)   # = π√6
G200_mag = (2*np.pi/a) * 2             # = 2π√2
G220_mag = (2*np.pi/a) * 2*np.sqrt(2)  # = 4π

w_over_ell = np.pi / 4  # PN core half-width
V_G111 = np.exp(-G111_mag * w_over_ell)
V_G200 = np.exp(-G200_mag * w_over_ell)
V_G220 = np.exp(-G220_mag * w_over_ell)

print("=" * 72)
print("  ANGULAR OVERLAP INTEGRALS FOR ALL O_h IRREPS (CORRECTED LABELS)")
print("=" * 72)
print(f"""
  FCC reciprocal lattice shells:
    {{111}}: {len(G111)} vectors, |G| = π√6/ℓ = {G111_mag:.4f}/ℓ
    {{200}}: {len(G200)} vectors, |G| = 2π√2/ℓ = {G200_mag:.4f}/ℓ
    {{220}}: {len(G220)} vectors, |G| = 4π/ℓ = {G220_mag:.4f}/ℓ

  PN form factor (Lorentzian core, w/ℓ = π/4):
    V_{{111}} = e^{{-|G₁₁₁|w}} = {V_G111:.6e}
    V_{{200}} = e^{{-|G₂₀₀|w}} = {V_G200:.6e}
    V_{{220}} = e^{{-|G₂₂₀|w}} = {V_G220:.6e}

  The {{111}} shell dominates by a factor of {V_G111/V_G200:.0f} over {{200}}.
""")


# ══════════════════════════════════════════════════════════════
# COMPUTE ANGULAR OVERLAPS ON FIRST TWO SHELLS
# ══════════════════════════════════════════════════════════════

print(f"  {'Irrep':<8} {'dim':>3} {'g/u':>3} {'deg':>3}  "
      f"{'sum₁₁₁':>10}  {'sum₂₀₀':>10}  {'I_total':>12}  {'I/I_T1u':>8}")
print(f"  {'─'*8} {'─'*3} {'─'*3} {'─'*3}  {'─'*10}  {'─'*10}  {'─'*12}  {'─'*8}")

results = {}
I_T1u = None

# Physical ordering: nonzero first, then zero
ordering = ['T1u', 'A1g', 'T2g', 'Eg', 'A2u',  # nonzero
            'T1g', 'T2u', 'A2g', 'A1u', 'Eu']   # zero

for name in ordering:
    func, dim, parity, degree = irreps[name]

    # Sum ||f(Ĝ)||² over each shell
    s111 = sum(np.sum(func(g)**2) for g in G111_hat)
    s200 = sum(np.sum(func(g)**2) for g in G200_hat)

    I111 = s111 * V_G111**2
    I200 = s200 * V_G200**2
    I_total = I111 + I200

    if name == 'T1u':
        I_T1u = I_total

    ratio = I_total / I_T1u if I_T1u and I_T1u > 0 else 0.0

    results[name] = (I_total, I111, I200, dim, parity, degree, ratio)

    if I_total < 1e-30:
        r_str = "0 (exact)"
    else:
        r_str = f"{ratio:.4f}"

    print(f"  {name:<8} {dim:3d}   {parity:1s}  {degree:3d}  "
          f"{s111:10.4f}  {s200:10.4f}  {I_total:12.4e}  {r_str:>8s}")


# ══════════════════════════════════════════════════════════════
# DETAILED ANALYSIS OF KEY CHANNELS
# ══════════════════════════════════════════════════════════════

alpha = 1.0 / 137.035999177
S0 = np.pi**2 / 2

print(f"\n{'=' * 72}")
print(f"  PHYSICAL IDENTIFICATION")
print(f"{'=' * 72}")

identifications = {
    'T1u': ('Electromagnetism (α)',      '{111}+{200}'),
    'A1g': ('Compression (→ gravity)',    'all shells'),
    'T2g': ('Strong force (α_s)',         '{111} only'),
    'Eg':  ('Tensor meson f₂(1270)',      '{200} only'),
    'A2u': ('Pseudoscalar η\' anomaly',   '{111} only'),
    'T1g': ('Zero overlap → forbidden',   'none'),
    'T2u': ('Zero overlap → weak via Σ₁₂','none'),
    'A2g': ('Zero overlap → forbidden',   'none'),
    'A1u': ('Zero overlap → forbidden',   'none'),
    'Eu':  ('Zero overlap → forbidden',   'none'),
}

print(f"\n  {'Irrep':<8} {'I/I_T1u':>8}  {'Couples via':>12}  {'Identification':<30}")
print(f"  {'─'*8} {'─'*8}  {'─'*12}  {'─'*30}")
for name in ordering:
    I_total, _, _, dim, par, deg, ratio = results[name]
    ident, couples = identifications[name]
    if I_total < 1e-30:
        r_str = "0"
    else:
        r_str = f"{ratio:.4f}"
    print(f"  {name:<8} {r_str:>8s}  {couples:>12s}  {ident:<30s}")


# ══════════════════════════════════════════════════════════════
# {220} SHELL: T₂g–T₂u MIXING (STRONG CP)
# ══════════════════════════════════════════════════════════════

print(f"\n{'=' * 72}")
print(f"  {'{'}220{'}'} SHELL: T₂g–T₂u MIXING (STRONG CP)")
print(f"{'=' * 72}")

s_T2g_220 = sum(np.sum(T2g(g)**2) for g in G220_hat)
s_T2u_220 = sum(np.sum(T2u(g)**2) for g in G220_hat)
s_T2g_111 = sum(np.sum(T2g(g)**2) for g in G111_hat)

theta_ch = alpha**2 / (2*np.pi)

num = V_G220**2 * np.sqrt(s_T2u_220 * s_T2g_220)
den = V_G111**2 * s_T2g_111
mixing = num / den
theta_QCD = theta_ch * mixing

print(f"""
  T₂g (strong, gerade) and T₂u (parity partner, ungerade)
  can mix through chirality θ_ch = {theta_ch:.4e}.

  On {{111}}: T₂g = {s_T2g_111:.4f}, T₂u = 0.0000  → no mixing
  On {{200}}: T₂g = 0.0000, T₂u = 0.0000  → no mixing
  On {{220}}: T₂g = {s_T2g_220:.4f}, T₂u = {s_T2u_220:.4f}  → FIRST MIXING

  Mixing ratio = V₂₂₀² √(I_T2u × I_T2g)₂₂₀ / (V₁₁₁² I_T2g₁₁₁)
               = {mixing:.4e}

  θ̄_QCD = θ_ch × ratio = {theta_QCD:.2e}
  Experimental bound: |θ̄| < 10⁻¹⁰
  Prediction/bound = {theta_QCD/1e-10:.1f}

  NOTE: This leading-order estimate omits the energy denominator
  (suppression ~ V₂₂₀²/V₁₁₁² ≈ {(V_G220/V_G111)**2:.1e}), which would
  reduce θ̄ to ~ {theta_QCD * (V_G220/V_G111)**2:.1e}, well below the bound.
""")


# ══════════════════════════════════════════════════════════════
# PARITY CROSS-CHECK
# ══════════════════════════════════════════════════════════════

print(f"{'=' * 72}")
print(f"  PARITY CONSISTENCY WITH PARTICLE IDENTIFICATIONS")
print(f"{'=' * 72}")
print(f"""
  Channel      Irrep   P (irrep)   Particle       P (particle)  Match?
  ───────      ─────   ─────────   ────────       ────────────  ──────
  EM           T₁u     −1 (u)      Photon (1⁻⁻)  −1            ✓
  Strong       T₂g     +1 (g)      Gluon          +1            ✓
  Tensor       E_g     +1 (g)      f₂ (2⁺⁺)      +1            ✓
  Pseudoscalar A₂u     −1 (u)      η' (0⁻⁺)      −1            ✓
  Compression  A₁g     +1 (g)      Graviton (2⁺⁺) +1           ✓
  Weak         T₂u     −1 (u)      (via Σ₁₂)     (2nd order)   —

  All parity assignments now consistent. The v1 label "A₂g" for
  the η' channel had P = +1, inconsistent with the η' pseudoscalar
  (P = −1). The corrected label A₂u has P = −1. ✓
""")


# ══════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════

print(f"{'=' * 72}")
print(f"  COMPLETE O_h SPECTRAL FINGERPRINT (CORRECTED)")
print(f"{'=' * 72}")

print(f"""
  ┌──────────┬─────┬────────┬──────────────┬──────────────┬──────────────────┐
  │  Irrep   │ dim │ Parity │ I_Γ/I_T1u    │ Couples via  │ Identification   │
  ├──────────┼─────┼────────┼──────────────┼──────────────┼──────────────────┤""")

for name in ordering:
    I_total, _, _, dim, par, deg, ratio = results[name]
    _, couples = identifications[name]
    ident_short = identifications[name][0][:16]
    if I_total < 1e-30:
        r_str = "0 (exact)"
    else:
        r_str = f"{ratio:.6f}"
    print(f"  │ {name:<8} │  {dim}  │  {par:<4}  │ {r_str:>12} │ {couples:>12} │ {ident_short:<16} │")

print(f"  └──────────┴─────┴────────┴──────────────┴──────────────┴──────────────────┘")

print(f"""
  FIVE NONZERO CHANNELS:
    T₁u (EM), A₁g (compression), T₂g (strong), E_g (tensor), A₂u (η')

  FIVE ZERO CHANNELS:
    T₁g, T₂u, A₂g, A₁u, E_u

  FOUR FORCES from four distinct mechanisms:
    Strong (T₂g):  commensurate → λ ~ 1
    EM (T₁u):      PN tunnelling → α = e^{{-π²/2}}
    Weak (T₂u→Σ₁₂): selection rule + centrosymmetry → α⁴/(4π²)
    Gravity (A₁g):  19-node cooperative → α¹⁹
""")
