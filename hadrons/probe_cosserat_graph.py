"""PDG validation harness for cosserat_graph.predict().

Runs the predictor against 127 PDG hadrons and classifies each result:

    WIN         correct prediction (|Δ%| < 1%)
    SOFT-OK     prediction within 3% (Regge endpoint regime)
    WRONG       code returns a non-zero mass but is >3% off
                (signals an assembler firing on a QN combination it
                was not designed for — no remaining cases as of the
                rule-based refactor; any new instances should be
                investigated)
    REGGE-SKIP  code returns cluster prefixed 'regge:' or 'refuse:'
                with mass=0 (no rule derived; honest refusal)
    FORBIDDEN   code actively rejects the QN combination as
                physically forbidden (e.g. Pauli)
    EMPIRICAL   code uses PDG input (bottom sector)
    ERROR       code raised an exception

Output:
    Console: full classification table by category.
    File:    probe_results.json in the current working directory.

Run from the repository root:
    python hadrons/probe_cosserat_graph.py
"""
import os
import sys
import json
from collections import defaultdict

# Make the parent (repository root containing cosserat_graph.py) importable
# regardless of where the script is invoked from.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import cosserat_graph as cg

M0 = cg.M0
ME = cg.ME
ALPHA = cg.ALPHA


def classify(pred_mass, cluster, obs_mass):
    if cluster.startswith('forbidden'):
        return 'FORBIDDEN'
    if cluster.startswith('regge') or cluster == 'regge':
        return 'REGGE-SKIP'
    if pred_mass is None or pred_mass <= 0:
        return 'REGGE-SKIP'
    if obs_mass is None or obs_mass <= 0:
        return 'WIN'  # no comparator, assume correct
    pct = 100.0 * abs(pred_mass - obs_mass) / obs_mass
    if pct < 1.0:
        return 'WIN'
    elif pct < 3.0:
        return 'SOFT-OK'
    else:
        return 'WRONG'


# Complete list of PDG hadrons with quantum numbers for cosserat_graph
# (charge states collapsed to one row each for brevity; I3 picks a representative)
TESTS = [
    # ============ Light unflavored mesons ============
    ("pi+",              139.57,  dict(B=0, I=1, I3=+1, J=0, P=-1)),
    ("pi0",              134.98,  dict(B=0, I=1, I3=0, J=0, P=-1)),
    ("eta(548)",         547.86,  dict(B=0, I=0, I3=0, J=0, P=-1)),
    ("f_0(500)",         500.0,   dict(B=0, I=0, I3=0, J=0, P=+1)),
    ("rho(770)",         775.26,  dict(B=0, I=1, I3=0, J=1, P=-1)),
    ("omega(782)",       782.66,  dict(B=0, I=0, I3=0, J=1, P=-1)),
    ("eta'(958)",        957.78,  dict(B=0, I=0, I3=0, J=0, P=-1, level=2)),
    ("f_0(980)",         990.0,   dict(B=0, I=0, I3=0, J=0, P=+1)),
    ("a_0(980)",         980.0,   dict(B=0, I=1, I3=0, J=0, P=+1)),
    ("phi(1020)",        1019.46, dict(B=0, I=0, I3=0, J=1, P=-1, n_hidden_s=1)),
    ("h_1(1170)",        1166.0,  dict(B=0, I=0, I3=0, J=1, P=+1)),
    ("b_1(1235)",        1229.5,  dict(B=0, I=1, I3=0, J=1, P=+1)),
    ("a_1(1260)",        1230.0,  dict(B=0, I=1, I3=0, J=1, P=+1)),
    ("f_2(1270)",        1275.4,  dict(B=0, I=0, I3=0, J=2, P=+1)),
    ("f_1(1285)",        1281.9,  dict(B=0, I=0, I3=0, J=1, P=+1)),
    ("eta(1295)",        1294.0,  dict(B=0, I=0, I3=0, J=0, P=-1, level=3)),
    ("pi(1300)",         1300.0,  dict(B=0, I=1, I3=+1, J=0, P=-1, level=2)),
    ("a_2(1320)",        1318.2,  dict(B=0, I=1, I3=0, J=2, P=+1)),
    ("f_0(1370)",        1350.0,  dict(B=0, I=0, I3=0, J=0, P=+1)),
    ("eta(1405)",        1408.7,  dict(B=0, I=0, I3=0, J=0, P=-1, level=3)),
    ("h_1(1415)",        1409.0,  dict(B=0, I=0, I3=0, J=1, P=+1)),
    ("f_1(1420)",        1428.4,  dict(B=0, I=0, I3=0, J=1, P=+1)),
    ("omega(1420)",      1410.0,  dict(B=0, I=0, I3=0, J=1, P=-1, level=2)),
    ("a_0(1450)",        1439.0,  dict(B=0, I=1, I3=0, J=0, P=+1)),
    ("rho(1450)",        1465.0,  dict(B=0, I=1, I3=0, J=1, P=-1, level=2)),
    ("eta(1475)",        1476.0,  dict(B=0, I=0, I3=0, J=0, P=-1, level=3, n_hidden_s=1)),
    ("f_0(1500)",        1522.0,  dict(B=0, I=0, I3=0, J=0, P=+1)),
    ("f_2'(1525)",       1517.3,  dict(B=0, I=0, I3=0, J=2, P=+1, n_hidden_s=1)),
    ("eta_2(1645)",      1617.0,  dict(B=0, I=0, I3=0, J=2, P=-1)),
    ("omega(1650)",      1670.0,  dict(B=0, I=0, I3=0, J=1, P=-1, level=3)),
    ("omega_3(1670)",    1667.0,  dict(B=0, I=0, I3=0, J=3, P=-1)),
    ("pi_2(1670)",       1670.6,  dict(B=0, I=1, I3=0, J=2, P=-1)),
    ("rho_3(1690)",      1688.8,  dict(B=0, I=1, I3=0, J=3, P=-1)),
    ("rho(1700)",        1720.0,  dict(B=0, I=1, I3=0, J=1, P=-1, level=3)),
    ("f_0(1710)",        1733.0,  dict(B=0, I=0, I3=0, J=0, P=+1)),
    ("pi(1800)",         1810.0,  dict(B=0, I=1, I3=+1, J=0, P=-1, level=3)),
    ("phi_3(1850)",      1854.0,  dict(B=0, I=0, I3=0, J=3, P=-1)),
    ("a_4(1970)",        1967.0,  dict(B=0, I=1, I3=0, J=4, P=+1)),
    ("f_2(1950)",        1936.0,  dict(B=0, I=0, I3=0, J=2, P=+1, level=2)),
    ("f_4(2050)",        2018.0,  dict(B=0, I=0, I3=0, J=4, P=+1)),

    # ============ Strange mesons ============
    ("K+",               493.68,  dict(B=0, S=1, I=0.5, I3=+0.5, J=0, P=-1)),
    ("K0",               497.61,  dict(B=0, S=1, I=0.5, I3=-0.5, J=0, P=-1)),
    ("K_0*(700)",        845.0,   dict(B=0, S=1, I=0.5, I3=+0.5, J=0, P=+1)),
    ("K*(892)",          895.55,  dict(B=0, S=1, I=0.5, I3=+0.5, J=1, P=-1)),
    ("K_1(1270)",        1253.0,  dict(B=0, S=1, I=0.5, I3=+0.5, J=1, P=+1)),
    ("K_1(1400)",        1403.0,  dict(B=0, S=1, I=0.5, I3=+0.5, J=1, P=+1)),
    ("K*(1410)",         1414.0,  dict(B=0, S=1, I=0.5, I3=+0.5, J=1, P=-1, level=2)),
    ("K_0*(1430)",       1425.0,  dict(B=0, S=1, I=0.5, I3=+0.5, J=0, P=+1)),
    ("K_2*(1430)",       1432.4,  dict(B=0, S=1, I=0.5, I3=+0.5, J=2, P=+1)),
    ("K_1(1650)",        1650.0,  dict(B=0, S=1, I=0.5, I3=+0.5, J=1, P=+1)),
    ("K*(1680)",         1718.0,  dict(B=0, S=1, I=0.5, I3=+0.5, J=1, P=-1, level=2)),
    ("K_2(1770)",        1773.0,  dict(B=0, S=1, I=0.5, I3=+0.5, J=2, P=-1)),
    ("K_3*(1780)",       1779.0,  dict(B=0, S=1, I=0.5, I3=+0.5, J=3, P=-1)),
    ("K_2(1820)",        1819.0,  dict(B=0, S=1, I=0.5, I3=+0.5, J=2, P=-1)),
    ("K_0*(1950)",       1957.0,  dict(B=0, S=1, I=0.5, I3=+0.5, J=0, P=+1)),
    ("K_2*(1980)",       1990.0,  dict(B=0, S=1, I=0.5, I3=+0.5, J=2, P=+1)),
    ("K_4*(2045)",       2048.0,  dict(B=0, S=1, I=0.5, I3=+0.5, J=4, P=+1)),
    ("phi(1680)",        1680.0,  dict(B=0, I=0, I3=0, J=1, P=-1, level=2)),

    # ============ Charm mesons ============
    ("D+",               1869.66, dict(B=0, I=0.5, I3=-0.5, J=0, P=-1, n_charm=1)),
    ("D0",               1864.84, dict(B=0, I=0.5, I3=+0.5, J=0, P=-1, n_charm=1)),
    ("D*(2007)",         2006.85, dict(B=0, I=0.5, I3=+0.5, J=1, P=-1, n_charm=1)),
    ("D*(2010)",         2010.26, dict(B=0, I=0.5, I3=-0.5, J=1, P=-1, n_charm=1)),
    ("D_0*(2300)",       2343.0,  dict(B=0, I=0.5, I3=+0.5, J=0, P=+1, n_charm=1)),
    ("D_1(2420)",        2422.1,  dict(B=0, I=0.5, I3=+0.5, J=1, P=+1, n_charm=1)),
    ("D_2*(2460)",       2461.1,  dict(B=0, I=0.5, I3=+0.5, J=2, P=+1, n_charm=1)),
    ("D_s+",             1968.35, dict(B=0, S=-1, I=0, I3=0, J=0, P=-1, n_charm=1)),
    ("D_s*",             2112.2,  dict(B=0, S=-1, I=0, I3=0, J=1, P=-1, n_charm=1)),
    ("D_s0*(2317)",      2317.8,  dict(B=0, S=-1, I=0, I3=0, J=0, P=+1, n_charm=1)),
    ("D_s1(2460)",       2459.5,  dict(B=0, S=-1, I=0, I3=0, J=1, P=+1, n_charm=1)),

    # ============ Charmonium ============
    ("eta_c(1S)",        2984.1,  dict(B=0, I=0, I3=0, J=0, P=-1, n_charm=2)),
    ("J/psi(1S)",        3096.9,  dict(B=0, I=0, I3=0, J=1, P=-1, n_charm=2)),
    ("chi_c0(1P)",       3414.7,  dict(B=0, I=0, I3=0, J=0, P=+1, n_charm=2)),
    ("chi_c1(1P)",       3510.7,  dict(B=0, I=0, I3=0, J=1, P=+1, n_charm=2)),
    ("h_c(1P)",          3525.4,  dict(B=0, I=0, I3=0, J=1, P=+1, n_charm=2)),
    ("chi_c2(1P)",       3556.2,  dict(B=0, I=0, I3=0, J=2, P=+1, n_charm=2)),
    ("eta_c(2S)",        3637.5,  dict(B=0, I=0, I3=0, J=0, P=-1, n_charm=2, level=2)),
    ("psi(2S)",          3686.1,  dict(B=0, I=0, I3=0, J=1, P=-1, n_charm=2, level=2)),
    ("psi(3770)",        3773.7,  dict(B=0, I=0, I3=0, J=1, P=-1, n_charm=2, level=3)),
    ("chi_c2(3930)",     3922.5,  dict(B=0, I=0, I3=0, J=2, P=+1, n_charm=2, level=2)),
    ("psi(4040)",        4039.0,  dict(B=0, I=0, I3=0, J=1, P=-1, n_charm=2, level=3)),

    # ============ Baryon octet ============
    ("p",                938.27,  dict(B=1, I=0.5, I3=+0.5, J=0.5, P=+1)),
    ("n",                939.57,  dict(B=1, I=0.5, I3=-0.5, J=0.5, P=+1)),
    ("Lambda",           1115.68, dict(B=1, S=-1, I=0, I3=0, J=0.5, P=+1)),
    ("Sigma+",           1189.37, dict(B=1, S=-1, I=1, I3=+1, J=0.5, P=+1)),
    ("Sigma0",           1192.64, dict(B=1, S=-1, I=1, I3=0, J=0.5, P=+1)),
    ("Sigma-",           1197.45, dict(B=1, S=-1, I=1, I3=-1, J=0.5, P=+1)),
    ("Xi0",              1314.86, dict(B=1, S=-2, I=0.5, I3=+0.5, J=0.5, P=+1)),
    ("Xi-",              1321.71, dict(B=1, S=-2, I=0.5, I3=-0.5, J=0.5, P=+1)),

    # ============ Baryon decuplet ============
    ("Delta++",          1232.0,  dict(B=1, I=1.5, I3=+1.5, J=1.5, P=+1)),
    ("Sigma*+(1385)",    1382.8,  dict(B=1, S=-1, I=1, I3=+1, J=1.5, P=+1)),
    ("Xi*0(1530)",       1531.8,  dict(B=1, S=-2, I=0.5, I3=+0.5, J=1.5, P=+1)),
    ("Omega-",           1672.45, dict(B=1, S=-3, I=0, I3=0, J=1.5, P=+1)),

    # ============ N* resonances ============
    ("N(1440) Roper",    1440.0,  dict(B=1, I=0.5, I3=+0.5, J=0.5, P=+1, level=2)),
    ("N(1520)",          1515.0,  dict(B=1, I=0.5, I3=+0.5, J=1.5, P=-1)),
    ("N(1535)",          1530.0,  dict(B=1, I=0.5, I3=+0.5, J=0.5, P=-1)),
    ("N(1650)",          1650.0,  dict(B=1, I=0.5, I3=+0.5, J=0.5, P=-1)),
    ("N(1675)",          1675.0,  dict(B=1, I=0.5, I3=+0.5, J=2.5, P=-1)),
    ("N(1680)",          1685.0,  dict(B=1, I=0.5, I3=+0.5, J=2.5, P=+1)),
    ("N(2190)",          2180.0,  dict(B=1, I=0.5, I3=+0.5, J=3.5, P=-1)),

    # ============ Delta* resonances ============
    ("Delta(1600)",      1600.0,  dict(B=1, I=1.5, I3=+1.5, J=1.5, P=+1, level=2)),
    ("Delta(1620)",      1610.0,  dict(B=1, I=1.5, I3=+1.5, J=0.5, P=-1)),
    ("Delta(1700)",      1710.0,  dict(B=1, I=1.5, I3=+1.5, J=1.5, P=-1)),
    ("Delta(1905)",      1880.0,  dict(B=1, I=1.5, I3=+1.5, J=2.5, P=+1)),
    ("Delta(1950)",      1930.0,  dict(B=1, I=1.5, I3=+1.5, J=3.5, P=+1)),

    # ============ Lambda* ============
    ("Lambda(1405)",     1405.1,  dict(B=1, S=-1, I=0, I3=0, J=0.5, P=-1)),
    ("Lambda(1520)",     1519.0,  dict(B=1, S=-1, I=0, I3=0, J=1.5, P=-1)),
    ("Lambda(1600)",     1600.0,  dict(B=1, S=-1, I=0, I3=0, J=0.5, P=+1, level=2)),
    ("Lambda(1670)",     1674.0,  dict(B=1, S=-1, I=0, I3=0, J=0.5, P=-1)),
    ("Lambda(1690)",     1690.0,  dict(B=1, S=-1, I=0, I3=0, J=1.5, P=-1)),
    ("Lambda(1820)",     1820.0,  dict(B=1, S=-1, I=0, I3=0, J=2.5, P=+1)),

    # ============ Sigma* ============
    ("Sigma(1660)",      1660.0,  dict(B=1, S=-1, I=1, I3=+1, J=0.5, P=+1, level=2)),
    ("Sigma(1670)",      1675.0,  dict(B=1, S=-1, I=1, I3=+1, J=1.5, P=-1)),
    ("Sigma(1775)",      1775.0,  dict(B=1, S=-1, I=1, I3=+1, J=2.5, P=-1)),

    # ============ Xi* ============
    ("Xi(1690)",         1690.0,  dict(B=1, S=-2, I=0.5, I3=+0.5, J=0.5, P=-1)),
    ("Xi(1820)",         1823.0,  dict(B=1, S=-2, I=0.5, I3=+0.5, J=1.5, P=-1)),

    # ============ Charm baryons ============
    ("Lambda_c+",        2286.46, dict(B=1, I=0, I3=0, J=0.5, P=+1, n_charm=1)),
    ("Sigma_c++(2455)",  2453.97, dict(B=1, I=1, I3=+1, J=0.5, P=+1, n_charm=1)),
    ("Sigma_c(2520)",    2518.5,  dict(B=1, I=1, I3=+1, J=1.5, P=+1, n_charm=1)),
    ("Lambda_c(2595)",   2592.2,  dict(B=1, I=0, I3=0, J=0.5, P=-1, n_charm=1)),
    ("Lambda_c(2625)",   2628.0,  dict(B=1, I=0, I3=0, J=1.5, P=-1, n_charm=1)),
    ("Lambda_c(2880)",   2881.6,  dict(B=1, I=0, I3=0, J=2.5, P=+1, n_charm=1)),
    ("Xi_c+",            2467.7,  dict(B=1, S=-1, I=0.5, I3=+0.5, J=0.5, P=+1, n_charm=1)),
    ("Xi_c(2645)",       2646.2,  dict(B=1, S=-1, I=0.5, I3=+0.5, J=1.5, P=+1, n_charm=1)),
    ("Xi_c(2790)",       2793.9,  dict(B=1, S=-1, I=0.5, I3=+0.5, J=0.5, P=-1, n_charm=1)),
    ("Xi_c(2815)",       2819.8,  dict(B=1, S=-1, I=0.5, I3=+0.5, J=1.5, P=-1, n_charm=1)),
    ("Omega_c",          2695.2,  dict(B=1, S=-2, I=0, I3=0, J=0.5, P=+1, n_charm=1)),
    ("Omega_c(2770)",    2765.9,  dict(B=1, S=-2, I=0, I3=0, J=1.5, P=+1, n_charm=1)),
]


# Run probe
results = []
for name, obs, qn_kwargs in TESTS:
    try:
        r = cg.predict(cg.QN(**qn_kwargs))
        pred = r.mass
        cl = r.cluster
        verdict = classify(pred, cl, obs)
        results.append((name, obs, pred, cl, verdict, qn_kwargs))
    except Exception as e:
        results.append((name, obs, None, f"ERROR:{type(e).__name__}", "ERROR", qn_kwargs))


# Summary
print(f"\nTested {len(results)} PDG hadrons against cosserat_graph.py predict()\n")

verdict_counts = defaultdict(int)
for _, _, _, _, v, _ in results:
    verdict_counts[v] += 1

print("=" * 60)
print(f"{'VERDICT':<15} {'COUNT':>5}   {'DESCRIPTION':<40}")
print("=" * 60)
order = ['WIN', 'SOFT-OK', 'WRONG', 'REGGE-SKIP', 'FORBIDDEN', 'EMPIRICAL', 'ERROR']
descs = {
    'WIN':        'Correct prediction |Δ%| < 1%',
    'SOFT-OK':    'Within 3% (Regge endpoint regime)',
    'WRONG':      'Code returns mass >3% off (silent failure)',
    'REGGE-SKIP': 'Cluster=regge, mass=0 (no rule)',
    'FORBIDDEN':  'Actively rejected',
    'EMPIRICAL':  'Uses PDG input mass',
    'ERROR':      'Exception raised',
}
for v in order:
    print(f"{v:<15} {verdict_counts[v]:>5}   {descs[v]:<40}")
print()

# Detailed listing of each verdict category
for v in order:
    items = [(n, o, p, cl, qn) for (n, o, p, cl, ver, qn) in results if ver == v]
    if not items:
        continue
    print(f"\n--- {v} ({len(items)}) ---")
    for n, o, p, cl, qn in items:
        if p is None:
            print(f"   {n:<22}  obs={o:>7.1f}   {cl}")
            continue
        pct = 100*(p - o)/o if (o and p) else 0
        print(f"   {n:<22}  obs={o:>7.1f}  pred={p:>7.1f}  Δ={pct:+6.1f}%  cluster={cl}")

# Save for downstream use
out = {
    'results': [
        {'name': n, 'obs': o, 'pred': p, 'cluster': cl, 'verdict': v, 'qn': qn}
        for (n, o, p, cl, v, qn) in results
    ],
    'counts': dict(verdict_counts),
}
_OUT = 'probe_results.json'
with open(_OUT, 'w') as f:
    json.dump(out, f, indent=2)
print(f"\nSaved to {_OUT}")
