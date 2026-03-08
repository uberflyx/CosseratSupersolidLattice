#!/usr/bin/env python3
"""
exotic_catalogue.py
====================
Comprehensive exotic hadron catalogue: blind Cosserat lattice algorithm.

Accompanies: M. A. Cox, "The Cosserat Supersolid" (2026).
https://doi.org/10.5281/zenodo.18636501
Monograph Ch. 10, Sec. "Comprehensive exotic landscape".

For each exotic state, the algorithm:
1. Computes raw node count N_raw = m_obs / m_0
2. Identifies the nearest allowed FCC defect (integer or half-integer N)
3. Computes Q = round((m_obs - N*m_0) / m_e)
4. Checks if mass lies within delta_th = alpha * m_pred of an S-wave
   threshold pair (two conventional hadrons whose masses sum to ~m_obs)
5. Verdict: PASS (within delta_th), MARGINAL, FAIL, or ABOVE-THRESHOLD

Author: Mitchell A. Cox
Date:   March 2026
"""

import numpy as np

# ── Physical constants ──────────────────────────────────────────────
m_e   = 0.51099895       # MeV, electron mass
alpha = 1/137.035999177  # lattice-derived fine structure constant
m_0   = m_e / alpha      # ~70.025 MeV, node mass

# Theoretical uncertainty
delta_th_frac = alpha     # ~0.73%

# ── Conventional hadron masses (PDG 2024) for threshold pairs ───────
# These are the "building blocks" for molecular thresholds
had = {
    'pi+':    139.57039,   'pi0':    134.9768,
    'K+':     493.677,     'K0':     497.611,
    'eta':    547.862,     'eta_p':  957.78,
    'rho':    775.26,      'omega':  782.66,
    'phi':    1019.461,
    'D0':     1864.84,     'D+':     1869.66,
    'D*0':    2006.85,     'D*+':    2010.26,
    'Ds+':    1968.35,     'Ds*+':   2112.2,
    'D0_bar': 1864.84,     'D-':     1869.66,
    'D*0_bar':2006.85,     'D*-':    2010.26,
    'Ds-':    1968.35,     'Ds*-':   2112.2,
    'J/psi':  3096.900,    'psi2S':  3686.097,
    'eta_c':  2983.9,      'chi_c0': 3414.71,
    'chi_c1': 3510.67,     'h_c':    3525.37,
    'p':      938.272,     'n':      939.565,
    'Lambda':  1115.683,   'Sigma+': 1189.37,
    'Sigma0':  1192.642,   'Sigma-': 1197.449,
    'Sigma_c++': 2453.97,  'Sigma_c0': 2453.75,
    'Lambda_c+': 2286.46,
    'Xi_c+':   2467.71,    'Xi_c0':   2470.44,
    'Upsilon': 9460.30,
    'B0':      5279.66,    'B+':      5279.34,
    'B*0':     5324.71,    'B*+':     5324.71,
    'Bs':      5366.92,    'Bs*':     5415.4,
    'D1_2420': 2420.8,     'D2_2460': 2461.1,
    'Ds1_2536':2535.11,    'Ds2_2573':2569.1,
    'Ds1_2460':2459.5,
    'f0_980':  990.0,      'a0_980':  980.0,
    'K*_892':  891.67,
}

# ── S-wave threshold pairs ──────────────────────────────────────────
# For each exotic, we list candidate two-hadron thresholds
# Only S-wave (orbital L=0) pairs are "natural" molecular candidates
threshold_pairs = {}
# Build all meson-meson and meson-baryon pairs relevant for exotics
charm_mesons = ['D0','D+','D*0','D*+','Ds+','Ds*+','D0_bar','D-','D*0_bar','D*-','Ds-','Ds*-']
light_mesons = ['pi+','pi0','K+','K0','eta','rho','omega','phi','eta_p']
charmonia = ['J/psi','psi2S','eta_c','chi_c0','chi_c1','h_c']
baryons = ['p','n','Lambda','Sigma+','Sigma0','Sigma-','Sigma_c++','Sigma_c0','Lambda_c+','Xi_c+','Xi_c0']

# Generate all relevant pairs
all_pairs = []
# Charm-anticharm meson pairs (for hidden-charm exotics near threshold)
for m1 in ['D0','D+','D*0','D*+','Ds+','Ds*+']:
    for m2 in ['D0_bar','D-','D*0_bar','D*-','Ds-','Ds*-']:
        thr = had[m1] + had[m2]
        all_pairs.append((m1, m2, thr))

# Charmonium + light meson (for some states)
for cc in charmonia:
    for lm in light_mesons:
        thr = had[cc] + had[lm]
        all_pairs.append((cc, lm, thr))

# Charmonium + baryon (for pentaquarks)
for cc in charmonia:
    for b in baryons:
        thr = had[cc] + had[b]
        all_pairs.append((cc, b, thr))

# Charm meson + baryon (for pentaquarks)
for cm in charm_mesons[:6]:
    for b in baryons:
        thr = had[cm] + had[b]
        all_pairs.append((cm, b, thr))

# Charm meson + light meson (for open-charm exotics)
for cm in charm_mesons[:6]:
    for lm in light_mesons:
        thr = had[cm] + had[lm]
        all_pairs.append((cm, lm, thr))

# Di-charmonium (for fully heavy tetraquarks)
for c1 in charmonia:
    for c2 in charmonia:
        thr = had[c1] + had[c2]
        all_pairs.append((c1, c2, thr))


def find_nearest_thresholds(mass, n_results=3):
    """Find the n closest S-wave thresholds to a given mass."""
    diffs = [(abs(mass - thr), m1, m2, thr) for m1, m2, thr in all_pairs]
    diffs.sort()
    return diffs[:n_results]


def run_algorithm(name, mass_obs, mass_unc, quark_content, jpc="?",
                  discovered_by="", year="", notes=""):
    """
    Run the blind Cosserat lattice algorithm on one exotic state.
    
    Returns a dict with all results.
    """
    # Step 1: raw node count
    N_raw = mass_obs / m_0
    
    # Step 2: try integer and half-integer N
    candidates = []
    for N_try in [int(N_raw), int(N_raw)+1, 
                  int(N_raw)-1,
                  int(2*N_raw)/2, (int(2*N_raw)+1)/2,
                  (int(2*N_raw)-1)/2]:
        if N_try <= 0:
            continue
        Q_exact = (mass_obs - N_try * m_0) / m_e
        Q_int = round(Q_exact)
        m_pred = N_try * m_0 + Q_int * m_e
        resid_ppm = (m_pred / mass_obs - 1) * 1e6
        resid_pct = (m_pred / mass_obs - 1) * 100
        delta_th = alpha * m_pred
        within = abs(m_pred - mass_obs) < delta_th
        candidates.append({
            'N': N_try, 'Q': Q_int, 'Q_exact': Q_exact,
            'm_pred': m_pred, 'resid_ppm': resid_ppm,
            'resid_pct': resid_pct, 'delta_th': delta_th,
            'within': within,
            'Q_frac': abs(Q_exact - Q_int)
        })
    
    # Sort by |residual|
    candidates.sort(key=lambda x: abs(x['resid_ppm']))
    best = candidates[0]
    
    # Step 3: find nearest thresholds
    nearest_thr = find_nearest_thresholds(mass_obs, 3)
    
    # Step 4: classify
    closest_thr_dist = nearest_thr[0][0] if nearest_thr else 9999
    delta_th_val = alpha * mass_obs
    
    if best['within']:
        if closest_thr_dist < delta_th_val:
            verdict = "PASS (at threshold)"
        else:
            verdict = "PASS (no nearby threshold)"
    elif abs(best['resid_pct']) < 1.0:
        if closest_thr_dist < delta_th_val:
            verdict = "MARGINAL (at threshold)"
        else:
            verdict = "MARGINAL"
    elif closest_thr_dist < 10:  # within 10 MeV of a threshold
        verdict = "ABOVE-THRESHOLD resonance"
    else:
        verdict = "FAIL / OPEN"
    
    return {
        'name': name,
        'mass_obs': mass_obs,
        'mass_unc': mass_unc,
        'quark_content': quark_content,
        'jpc': jpc,
        'discovered_by': discovered_by,
        'year': year,
        'notes': notes,
        'N_raw': N_raw,
        'best': best,
        'all_candidates': candidates[:4],
        'nearest_thresholds': nearest_thr,
        'closest_thr_dist': closest_thr_dist,
        'delta_th_mass': delta_th_val,
        'verdict': verdict
    }


# ══════════════════════════════════════════════════════════════════════
# COMPLETE EXOTIC HADRON CATALOGUE
# Sources: CERN Courier Nov 2024, PDG 2024, LHCb/Belle/BESIII papers
# ══════════════════════════════════════════════════════════════════════

exotics = []

# ── LHC PENTAQUARKS (5) ──────────────────────────────────────────────
exotics.append(run_algorithm("Pcc(4312)+",  4311.9, 0.7, "cc̄uud", "1/2-",
    "LHCb", "2019", "Σc D̄ threshold"))
exotics.append(run_algorithm("Pcc(4380)+",  4380, 30, "cc̄uud", "3/2-?",
    "LHCb", "2015", "Broad; superseded by 4440/4457 split"))
exotics.append(run_algorithm("Pcc(4440)+",  4440.3, 1.3, "cc̄uud", "1/2-",
    "LHCb", "2019", "Σc D̄* threshold"))
exotics.append(run_algorithm("Pcc(4457)+",  4457.3, 0.6, "cc̄uud", "3/2-",
    "LHCb", "2019", "Σc D̄* threshold"))
exotics.append(run_algorithm("Pccs(4338)+", 4338.2, 0.7, "cc̄uds", "?",
    "LHCb", "2022", "First strange pentaquark; Ξc D̄ threshold"))

# ── LHC TETRAQUARKS: hidden-charm hidden-strange (J/ψ φ) (6) ───────
exotics.append(run_algorithm("χc1(4140)",   4146.8, 2.4, "cc̄ss̄", "1++",
    "CMS/CDF", "2013", "J/ψ φ threshold cusp?"))
exotics.append(run_algorithm("χc1(4274)",   4274, 8, "cc̄ss̄", "1++",
    "LHCb", "2016", ""))
exotics.append(run_algorithm("X(4500)",     4506, 11, "cc̄ss̄", "0++",
    "LHCb", "2016", ""))
exotics.append(run_algorithm("X(4630)",     4626, 16, "cc̄ss̄?", "1-?",
    "LHCb", "2021", "Near Ds D̄s1(2536) threshold"))
exotics.append(run_algorithm("χc1(4685)",   4684, 7, "cc̄ss̄?", "1++?",
    "LHCb", "2021", ""))
exotics.append(run_algorithm("X(4700)",     4704, 10, "cc̄ss̄", "0++",
    "LHCb", "2016", ""))

# ── LHC TETRAQUARKS: hidden-charm with strangeness (2) ──────────────
exotics.append(run_algorithm("Tccs1(4000)+", 4003, 6, "cc̄us̄", "1+",
    "LHCb", "2021", "Argand diagram confirmed; J/ψ K+"))
exotics.append(run_algorithm("Tccs1(4220)+", 4220, 15, "cc̄us̄", "1+",
    "LHCb", "2021", "J/ψ K+"))

# ── LHC TETRAQUARKS: open-charm (4) ─────────────────────────────────
exotics.append(run_algorithm("Tcc(3875)+",  3874.83, 0.11, "ccūd̄", "1+",
    "LHCb", "2021", "D0D*+ threshold; double open charm"))
exotics.append(run_algorithm("T*cs0(2870)++", 2870, 7, "cs̄ūd̄", "0+?",
    "LHCb", "2022", "Doubly charged; Ds π"))
exotics.append(run_algorithm("T*cs0(2900)0",  2900, 7, "cs̄ud̄?", "0+?",
    "LHCb", "2022", "Neutral partner of 2870++"))
exotics.append(run_algorithm("χc0(3960)",   3956, 8, "cc̄?", "0++",
    "LHCb", "2022", "Ds+ Ds- threshold; could be conventional χc0(2P)"))

# ── LHC TETRAQUARKS: hidden-charm no strangeness (3) ────────────────
exotics.append(run_algorithm("Tcc̄1(4430)+", 4478, 17, "cc̄ud̄", "1+-",
    "Belle/LHCb", "2007/2014", "Pathfinder exotic; ψ(2S)π+"))
exotics.append(run_algorithm("χc1(4010)",   4010, 5, "cc̄?", "?",
    "LHCb", "2022", "D*+ D-"))

# ── LHC TETRAQUARKS: fully charmed (3) ──────────────────────────────
exotics.append(run_algorithm("Tcccc(6600)", 6552, 20, "cc̄cc̄", "?",
    "LHCb/CMS/ATLAS", "2020", "J/ψ J/ψ; fully charmed"))
exotics.append(run_algorithm("Tcccc(6900)", 6886, 11, "cc̄cc̄", "0++/2++",
    "LHCb/CMS/ATLAS", "2020", "J/ψ J/ψ; CMS Nature 2025 spin-parity"))
exotics.append(run_algorithm("Tcccc(7100)", 7100, 50, "cc̄cc̄?", "?",
    "CMS", "2023", "Broad; needs confirmation"))

# ── PRE-LHC EXOTICS ─────────────────────────────────────────────────
exotics.append(run_algorithm("χc1(3872)",   3871.65, 0.06, "cc̄(uū/dd̄)", "1++",
    "Belle", "2003", "X(3872); D0 D̄*0 threshold"))
exotics.append(run_algorithm("Tcc̄1(3900)±", 3887.1, 2.6, "cc̄ud̄", "1+-",
    "BESIII/Belle", "2013", "Zc(3900); DD̄* threshold"))
exotics.append(run_algorithm("Tcc̄(4020)±",  4024.1, 1.9, "cc̄ud̄", "?",
    "BESIII", "2013", "Zc(4020); D*D̄* threshold"))
exotics.append(run_algorithm("Tbb̄1(10610)±",10607.2, 2.0, "bb̄ud̄", "1+-",
    "Belle", "2011", "Zb(10610); BB̄* threshold"))
exotics.append(run_algorithm("Tbb̄1(10650)±",10652.2, 1.5, "bb̄ud̄", "1+-",
    "Belle", "2011", "Zb(10650); B*B̄* threshold"))
exotics.append(run_algorithm("D*s0(2317)+", 2317.8, 0.5, "cs̄(+qq̄?)", "0+",
    "BaBar", "2003", "Near DK threshold; exotic or shifted cs"))
exotics.append(run_algorithm("ψ(4230)",     4222.5, 2.4, "cc̄?", "1--",
    "BaBar/BESIII", "2005", "Y(4260); D1(2420)D̄ molecular?"))
exotics.append(run_algorithm("ψ(4360)",     4368, 13, "cc̄?", "1--",
    "BaBar/Belle", "2007", "Y(4360)"))
exotics.append(run_algorithm("ψ(4660)",     4630, 10, "cc̄?", "1--",
    "Belle", "2007", "Y(4660); ψ(2S) f0(980)?"))

# ── NON-COLLIDER EXOTIC ─────────────────────────────────────────────
exotics.append(run_algorithm("d*(2380)",    2380, 10, "uuuddd", "3+",
    "WASA-at-COSY", "2011", "Dibaryon / hexaquark"))


# ══════════════════════════════════════════════════════════════════════
# OUTPUT
# ══════════════════════════════════════════════════════════════════════

print("=" * 120)
print("COMPREHENSIVE EXOTIC HADRON CATALOGUE — COSSERAT LATTICE BLIND ALGORITHM")
print(f"m_0 = {m_0:.3f} MeV,  m_e = {m_e:.5f} MeV,  α = 1/{1/alpha:.6f}")
print(f"δ_th = α × m_pred ≈ 0.73% of mass")
print("=" * 120)
print()

# Summary counters
n_pass = 0
n_marginal = 0
n_above = 0
n_fail = 0

for ex in exotics:
    b = ex['best']
    thr = ex['nearest_thresholds']
    
    verdict_short = ex['verdict']
    if 'PASS' in verdict_short: n_pass += 1
    elif 'MARGINAL' in verdict_short: n_marginal += 1
    elif 'ABOVE' in verdict_short: n_above += 1
    else: n_fail += 1
    
    print(f"{'─'*100}")
    print(f"  {ex['name']:25s}  m = {ex['mass_obs']:.1f} ± {ex['mass_unc']:.1f} MeV")
    print(f"  Quarks: {ex['quark_content']:12s}  J^PC = {ex['jpc']:8s}  "
          f"Disc: {ex['discovered_by']} ({ex['year']})")
    print(f"  Notes: {ex['notes']}")
    print(f"  N_raw = {ex['N_raw']:.3f}")
    print(f"  Best fit:  N = {b['N']:.1f},  Q = {b['Q']:+d}  "
          f"(Q_exact = {b['Q_exact']:+.3f}, fractional = {b['Q_frac']:.3f})")
    print(f"  m_pred = {b['m_pred']:.2f} MeV,  residual = {b['resid_pct']:+.4f}%  "
          f"({b['resid_ppm']:+.0f} ppm)")
    print(f"  δ_th = {b['delta_th']:.2f} MeV,  |Δm| = {abs(b['m_pred']-ex['mass_obs']):.2f} MeV,  "
          f"within δ_th: {b['within']}")
    print(f"  Nearest thresholds:")
    for d, m1, m2, thr_val in thr[:3]:
        print(f"    {m1} + {m2} = {thr_val:.2f} MeV  (Δ = {d:+.1f} MeV)")
    print(f"  ▶ VERDICT: {verdict_short}")
    print()

print("=" * 120)
print("SUMMARY")
print(f"  Total states analysed:  {len(exotics)}")
print(f"  PASS:                   {n_pass}")
print(f"  MARGINAL:               {n_marginal}")
print(f"  ABOVE-THRESHOLD:        {n_above}")
print(f"  FAIL/OPEN:              {n_fail}")
print()

# Compact table for LaTeX
print("\n" + "=" * 120)
print("COMPACT TABLE (for LaTeX integration)")
print(f"{'State':25s} {'m_obs':>10s} {'N':>6s} {'Q':>5s} {'m_pred':>10s} {'Δ%':>8s} {'Verdict':>30s}")
print("-" * 100)
for ex in exotics:
    b = ex['best']
    N_str = f"{b['N']:.1f}" if b['N'] != int(b['N']) else f"{int(b['N'])}"
    print(f"{ex['name']:25s} {ex['mass_obs']:10.1f} {N_str:>6s} {b['Q']:+5d} "
          f"{b['m_pred']:10.2f} {b['resid_pct']:+8.4f} {ex['verdict']:>30s}")

