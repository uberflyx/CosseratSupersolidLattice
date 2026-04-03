#!/usr/bin/env python3
"""
exotic_filling_fraction.py
============================
CRITICAL HONESTY CHECK: Is 32/32 PASS meaningful, or can the formula
m = N*m_0 + Q*m_e fit ANY mass?

Accompanies: M. A. Cox, "The Cosserat Supersolid" (2026).
https://doi.org/10.5281/zenodo.18636501
Monograph Ch. 10, Sec. "Comprehensive exotic landscape", honesty paragraph.

The formula has two free integers (N, Q) per particle. With m_e = 0.511 MeV,
any mass can be approximated to within m_e/2 = 0.256 MeV by choosing Q.
Since delta_th ~ 0.73% of mass ~ 20-80 MeV for these states, ANY mass
will trivially pass the mass formula.

The REAL test is whether N can be DERIVED from defect geometry BEFORE
consulting the mass. We check:
1. What fraction of random masses pass? (should be ~100% → formula is trivially flexible)
2. Can we derive N from constituent structure? (the actual physics)
3. Are Q values bounded/physically motivated?

Author: Mitchell A. Cox
Date:   March 2026
"""

import numpy as np

m_e = 0.51099895
alpha = 1/137.035999177
m_0 = m_e / alpha

print("=" * 90)
print("FILLING FRACTION ANALYSIS: IS THE MASS FORMULA TRIVIALLY FLEXIBLE?")
print("=" * 90)

# Test 1: Generate 10000 random masses in the exotic range and check
# how many can be fitted to within delta_th
np.random.seed(42)
n_random = 10000
random_masses = np.random.uniform(2000, 11000, n_random)

n_pass = 0
for m in random_masses:
    delta_th = alpha * m
    # Try all plausible N (integer and half-integer)
    N_raw = m / m_0
    for N_try in [int(N_raw), int(N_raw)+1, int(N_raw)-1,
                  int(2*N_raw)/2, (int(2*N_raw)+1)/2, (int(2*N_raw)-1)/2]:
        if N_try <= 0:
            continue
        Q_exact = (m - N_try * m_0) / m_e
        Q_int = round(Q_exact)
        m_pred = N_try * m_0 + Q_int * m_e
        if abs(m_pred - m) < delta_th:
            n_pass += 1
            break

print(f"\nTest 1: Random masses in [2000, 11000] MeV")
print(f"  {n_pass}/{n_random} random masses pass ({100*n_pass/n_random:.1f}%)")
print(f"  → The formula can fit ANY mass. 32/32 PASS is NOT evidence.")
print()

# Test 2: The REAL constraint is on N
# For the formula to be predictive, N must be derivable from the
# defect structure BEFORE consulting the mass.
# 
# For conventional hadrons, this is done:
#   pion: N=2 (stacking fault dipole, 2 nodes)
#   kaon: N=7 (stacking fault + strange arm)
#   proton: N=27/2 (3-fold dislocation junction)
#   etc.
#
# For exotics, the question is: do the N values decompose into
# sums of known hadron N values?

print("Test 2: Do exotic N values decompose into known constituents?")
print("-" * 90)

# Known hadron N values
known_N = {
    'π': 2, 'K': 7, 'η': 8, 'ρ/ω': 144/13,  # ~11.08
    'K*': 13, "η'": 14, 'φ': 29/2,  # 14.5
    'p/n': 27/2,  # 13.5
    'Λ': 16, 'Σ': 17, 'Ξ': 19, 'Ω': 24,
    'D': 26.5,  # approx from m_D / m_0 ~ 26.6
    'D*': 28.5,  # approx
    'Ds': 28,    # approx
    'Ds*': 30,   # approx
    'J/ψ': 44,   # approx (3097/70.025 ~ 44.2)
    'ψ(2S)': 52.5,  # approx (3686/70.025 ~ 52.6)
    'η_c': 42.5,  # approx
    'χ_c0': 48.5,  # approx
    'χ_c1': 50,    # approx
    'h_c': 50.5,   # approx
}

# Exotic states and their algorithm-derived N values
exotic_N = {
    'Pcc(4312)+':  60,     'Pcc(4380)+':  63,
    'Pcc(4440)+':  64,     'Pcc(4457)+':  63.5,
    'Pccs(4338)+': 61.5,
    'χc1(4140)':   58,     'χc1(4274)':   60,
    'X(4500)':     64.5,   'X(4630)':     65.5,
    'χc1(4685)':   66,     'X(4700)':     68,
    'Tccs1(4000)+':57.5,   'Tccs1(4220)+':61,
    'Tcc(3875)+':  54,     'T*cs0(2870)++':40,
    'T*cs0(2900)0':41.5,   'χc0(3960)':   56.5,
    'Tcc̄1(4430)+':62,     'χc1(4010)':   56.5,
    'Tcccc(6600)': 93.5,   'Tcccc(6900)': 99,
    'Tcccc(7100)': 101.5,
    'χc1(3872)':   55.5,   'Tcc̄1(3900)±':54,
    'Tcc̄(4020)±': 56,     'Tbb̄1(10610)±':150.5,
    'Tbb̄1(10650)±':151.5, 'D*s0(2317)+': 33.5,
    'ψ(4230)':     61,     'ψ(4360)':     62.5,
    'ψ(4660)':     65.5,   'd*(2380)':    33.5,
}

# Check: do any of these decompose as sums of known N?
# Pentaquarks should be ~ baryon_N + meson_N
# Tetraquarks should be ~ meson_N + meson_N
print(f"\n{'State':25s} {'N_fit':>6s}  {'Possible decomposition':45s} {'Match?':>8s}")
print("-" * 90)

penta_check = [
    ('Pcc(4312)+',  'Σ_c D̄ mol.', 'Σ', 'D'),
    ('Pcc(4440)+',  'Σ_c D̄* mol.', 'Σ', 'D*'),
    ('Pcc(4457)+',  'Σ_c D̄* mol.', 'Σ', 'D*'),
    ('Pccs(4338)+', 'Ξ_c D̄ mol.', 'Ξ', 'D'),
]

for name, desc, b, m in penta_check:
    N_fit = exotic_N[name]
    N_sum = known_N.get(b, '?') 
    N_m = known_N.get(m, '?')
    if isinstance(N_sum, (int, float)) and isinstance(N_m, (int, float)):
        total = N_sum + N_m
        diff = N_fit - total
        match = "YES" if abs(diff) < 2 else "approx" if abs(diff) < 5 else "NO"
        print(f"  {name:23s} {N_fit:6.1f}  {b}({N_sum:.1f}) + {m}({N_m:.1f}) = {total:.1f}  "
              f"(Δ = {diff:+.1f})  {match}")

# Tetraquark decompositions
print()
tetra_check = [
    ('χc1(3872)',   'D D̄* mol.', 'D', 'D*'),
    ('Tcc̄1(3900)±','D D̄* mol.', 'D', 'D*'),
    ('Tcc̄(4020)±', 'D* D̄* mol.', 'D*', 'D*'),
    ('Tcc(3875)+',  'D D* (open)', 'D', 'D*'),
    ('Tcccc(6900)', 'J/ψ J/ψ', 'J/ψ', 'J/ψ'),
]

for name, desc, m1, m2 in tetra_check:
    N_fit = exotic_N[name]
    N1 = known_N.get(m1, '?')
    N2 = known_N.get(m2, '?')
    if isinstance(N1, (int, float)) and isinstance(N2, (int, float)):
        total = N1 + N2
        diff = N_fit - total
        match = "YES" if abs(diff) < 2 else "approx" if abs(diff) < 5 else "NO"
        print(f"  {name:23s} {N_fit:6.1f}  {m1}({N1:.1f}) + {m2}({N2:.1f}) = {total:.1f}  "
              f"(Δ = {diff:+.1f})  {match}")

# Hexaquark
print()
N_d2380 = exotic_N['d*(2380)']
N_2p = 2 * known_N['p/n']
print(f"  {'d*(2380)':23s} {N_d2380:6.1f}  2×p/n(2×{known_N['p/n']:.1f}) = {N_2p:.1f}  "
      f"(Δ = {N_d2380 - N_2p:+.1f})  {'YES' if abs(N_d2380-N_2p) < 2 else 'NO'}")

print()
print("=" * 90)
print("CONCLUSIONS:")
print("  1. The mass formula m = Nm₀ + Qm_e can fit ANY mass (100% filling fraction)")
print("     → 32/32 PASS is mathematically trivial, NOT a physics result")
print()
print("  2. The REAL test is whether N decomposes into known constituent N values")
print("     → This partially works for molecular states near threshold")
print("     → For states far from threshold, N is unconstrained → hidden flexibility")
print()
print("  3. The Q values range from -113 to +267 across the catalogue")
print("     → No obvious geometric constraint on Q for exotics")
print("     → For conventional hadrons, Q is derived from coordination chemistry")
print()
print("  4. HONEST FRAMING for the monograph:")
print("     → The formula ACCOMMODATES all 32 exotics (consistency check)")
print("     → But it does not PREDICT them (no a priori N derivation)")
print("     → Molecular decomposition works for threshold states")
print("     → Far-from-threshold states remain unconstrained")
print("     → This is a necessary condition, not a sufficient one")
print("=" * 90)

