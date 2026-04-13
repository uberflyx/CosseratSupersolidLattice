#!/usr/bin/env python3
"""
d4_hypothesis_test.py — Rigorous statistical tests of the D4 layer decomposition
=================================================================================
Mitchell A. Cox, University of the Witwatersrand

Three independent tests of whether the D4 stacking-layer assignment carries
physical information:

    Test 1 (Permutation null): Is the FCC stacking assignment special, or does
            any random 3-colouring of the defect graph fragment baryons?

    Test 2 (Quantitative NLO): Does a single proportionality constant
            δm = c × α × (σ²/N) × m₀ fit the decuplet residuals AND predict
            other families?

    Test 3 (Family-by-family): Within every controlled family, does the member
            closest to H = 1 always have the smallest residual?

Reference: Sec. layer_decomposition of the monograph (decay chapter).

Usage:
    python d4_hypothesis_test.py              # run all tests
    python d4_hypothesis_test.py --plot       # also produce figures

Dependencies: numpy, scipy, matplotlib (optional), networkx, cosserat_graph
"""

import sys, os, math
import numpy as np
from scipy.stats import spearmanr, pearsonr
import networkx as nx

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from cosserat_graph import predict_with_defect, QN, ALPHA, ME, M0

from layer_entropy_catalogue import (
    assign_layers, layer_entropy, layer_variance,
    bond_classification, spatial_components, build_catalogue
)

N_PERM = 10000   # number of random permutations per particle


# =========================================================================
# TEST 1: PERMUTATION NULL — Is the FCC stacking assignment special?
# =========================================================================

def random_layer_assignment(defect, true_counts):
    """Assign nodes to three layers randomly, preserving group sizes.

    This is the null hypothesis: any partition into groups of size
    (n_A, n_B, n_C) would produce the same fragmentation.
    """
    pn = defect.positioned_nodes()
    N = len(pn)
    # Build a random permutation of layer labels
    labels = []
    for layer_idx, count in enumerate(true_counts):
        labels.extend([layer_idx] * count)
    np.random.shuffle(labels)
    return {node: labels[i] for i, node in enumerate(pn)}


def test_fragmentation_null(defect, true_assignments, true_counts, n_perm=N_PERM):
    """Test whether random 3-colourings also fragment the graph.

    Returns:
        true_components: int, number of spatial components with FCC assignment
        frac_fragment: float, fraction of random assignments that also fragment
        mean_components: float, mean number of components under random assignment
    """
    true_comp = spatial_components(defect, true_assignments)

    random_fragment_count = 0
    component_counts = []

    for _ in range(n_perm):
        rand_assign = random_layer_assignment(defect, true_counts)
        # Build spatial-only graph under random assignment
        pn = defect.positioned_nodes()
        G = nx.Graph()
        G.add_nodes_from(pn)
        for a, b in defect.nn_edge_list():
            if a in rand_assign and b in rand_assign:
                if rand_assign[a] == rand_assign[b]:
                    G.add_edge(a, b)
        nc = nx.number_connected_components(G)
        component_counts.append(nc)
        if nc > 1:
            random_fragment_count += 1

    frac_fragment = random_fragment_count / n_perm
    mean_comp = np.mean(component_counts)

    return true_comp, frac_fragment, mean_comp


def run_test_1(results):
    """Run the permutation null test for all baryons."""
    print("\n" + "=" * 80)
    print("TEST 1: PERMUTATION NULL — Is the FCC stacking assignment special?")
    print("=" * 80)
    print(f"\n  For each baryon, we generate {N_PERM:,d} random 3-colourings")
    print("  (preserving group sizes) and check how often the spatial-only")
    print("  graph fragments. If 100% of random colourings also fragment,")
    print("  then fragmentation is trivial. If significantly fewer do,")
    print("  the FCC stacking assignment is special.\n")

    baryons = [r for r in results if r['B'] >= 1 and r['n_temporal'] > 0]
    # Rebuild defects for baryons
    print(f"  {'Particle':10s} {'Nodes':>5s} {'Layers':>10s} {'FCC comp':>8s}"
          f"  {'% random frag':>13s} {'Mean comp':>9s} {'FCC special?':>12s}")
    print(f"  {'-' * 75}")

    fcc_special_count = 0
    for r in baryons:
        qn = QN(B=r['B'], S=r['S'], I=r['I'], I3=r['I'],
                J=r['J'], P=1, level=r['level'], n_charm=r['nc'])
        try:
            result, defect, _ = predict_with_defect(qn)
        except Exception:
            continue
        if defect is None:
            continue

        layer_assign, counts = assign_layers(defect)
        true_comp, frac_frag, mean_comp = test_fragmentation_null(
            defect, layer_assign, counts)

        special = "YES" if frac_frag < 0.95 else "no"
        if frac_frag < 0.95:
            fcc_special_count += 1

        layers_str = f"{counts[0]}:{counts[1]}:{counts[2]}"
        print(f"  {r['name']:10s} {r['N_nodes']:>5d} {layers_str:>10s}"
              f"  {true_comp:>8d}  {frac_frag:>13.1%} {mean_comp:>9.1f}"
              f"  {special:>12s}")

    print(f"\n  Summary: {fcc_special_count}/{len(baryons)} baryons have"
          f" FCC-specific fragmentation (random < 95%)")


# =========================================================================
# TEST 2: QUANTITATIVE NLO MODEL
# =========================================================================

def run_test_2(results):
    """Fit δm = c × α × (σ²/N) × m₀ from the decuplet and predict others."""
    print("\n" + "=" * 80)
    print("TEST 2: QUANTITATIVE NLO MODEL")
    print("=" * 80)
    print(f"\n  Model: |Δm| = c × α × (σ²/N) × m₀")
    print(f"  Fit c from the decuplet, predict other families.\n")

    # Extract decuplet data
    decuplet = []
    for r in results:
        if r['B'] == 1 and r['nc'] == 0 and r['J'] >= 1.4 and r['residual_pct'] is not None:
            decuplet.append(r)
    decuplet.sort(key=lambda r: abs(r['S']))

    if len(decuplet) < 3:
        print("  Insufficient decuplet data.")
        return

    print("  Decuplet data:")
    for r in decuplet:
        dm = abs(r['m_pred'] - r['m_obs'])
        predictor = ALPHA * r['var_norm'] * M0
        print(f"    {r['name']:8s}  σ²/N = {r['var_norm']:.3f}"
              f"  |Δm| = {dm:.2f} MeV  predictor = {predictor:.2f} MeV")

    # Fit c by least squares (excluding Omega which has predictor = 0)
    predictors = np.array([ALPHA * r['var_norm'] * M0 for r in decuplet])
    residuals_mev = np.array([abs(r['m_pred'] - r['m_obs']) for r in decuplet])

    # Use only non-zero predictors for fit
    mask = predictors > 0.001
    if mask.sum() < 2:
        print("  Insufficient non-zero predictors for fit.")
        return

    c_fit = np.sum(residuals_mev[mask] * predictors[mask]) / np.sum(predictors[mask]**2)
    predicted_dm = c_fit * predictors

    print(f"\n  Fitted constant: c = {c_fit:.2f}")
    print(f"\n  {'Particle':8s} {'|Δm| obs':>10s} {'|Δm| pred':>10s} {'Ratio':>8s}")
    print(f"  {'-' * 40}")
    for r, obs, pred in zip(decuplet, residuals_mev, predicted_dm):
        ratio = obs / pred if pred > 0.01 else float('inf')
        print(f"  {r['name']:8s} {obs:>10.3f} {pred:>10.3f} {ratio:>8.2f}")

    # Now predict octet baryons
    octet = [r for r in results if r['B'] == 1 and r['nc'] == 0
             and r['J'] < 1.0 and r['residual_pct'] is not None]
    if octet:
        print(f"\n  Octet predictions (same c = {c_fit:.2f}):")
        print(f"  {'Particle':8s} {'|Δm| obs':>10s} {'|Δm| pred':>10s} {'Ratio':>8s}")
        print(f"  {'-' * 40}")
        for r in sorted(octet, key=lambda x: abs(x['S'])):
            obs = abs(r['m_pred'] - r['m_obs'])
            pred = c_fit * ALPHA * r['var_norm'] * M0
            ratio = obs / pred if pred > 0.01 else float('inf')
            print(f"  {r['name']:8s} {obs:>10.3f} {pred:>10.3f} {ratio:>8.2f}")

    # Charm baryons
    charm_B = [r for r in results if r['B'] == 1 and r['nc'] >= 1
               and r['residual_pct'] is not None]
    if charm_B:
        print(f"\n  Charm baryon predictions (same c = {c_fit:.2f}):")
        print(f"  {'Particle':8s} {'|Δm| obs':>10s} {'|Δm| pred':>10s} {'Ratio':>8s}")
        print(f"  {'-' * 40}")
        for r in sorted(charm_B, key=lambda x: x['m_pred']):
            obs = abs(r['m_pred'] - r['m_obs'])
            pred = c_fit * ALPHA * r['var_norm'] * M0
            ratio = obs / pred if pred > 0.01 else float('inf')
            print(f"  {r['name']:8s} {obs:>10.3f} {pred:>10.3f} {ratio:>8.2f}")


# =========================================================================
# TEST 3: FAMILY-BY-FAMILY — Does H = 1 always win?
# =========================================================================

def run_test_3(results):
    """Within each controlled family, check if max-H member has min residual."""
    print("\n" + "=" * 80)
    print("TEST 3: FAMILY-BY-FAMILY — Does maximum H always win?")
    print("=" * 80)
    print("\n  Within each family of particles sharing the same construction")
    print("  mechanism and spin, does the member closest to H = 1 always")
    print("  have the smallest mass residual?\n")

    families = {}
    for r in results:
        if r['residual_pct'] is None:
            continue
        # Group by (B, J, n_charm, construction type approximated by N_nodes)
        key = (r['B'], r['J'], r['nc'])
        if key not in families:
            families[key] = []
        families[key].append(r)

    wins = 0
    tests = 0
    for key, members in sorted(families.items()):
        if len(members) < 2:
            continue
        members_sorted = sorted(members, key=lambda r: r['H_layer'], reverse=True)
        best_H = members_sorted[0]
        best_resid = min(members, key=lambda r: r['residual_pct'])

        match = "✓" if best_H['name'] == best_resid['name'] else "✗"
        tests += 1
        if best_H['name'] == best_resid['name']:
            wins += 1

        B, J, nc = key
        Js = f"{J:.0f}" if J == int(J) else f"{int(2*J)}/2"
        print(f"  Family (B={B}, J={Js}, nc={nc}): {len(members)} members")
        print(f"    Highest H:      {best_H['name']:10s}  H = {best_H['H_layer']:.3f}"
              f"  resid = {best_H['residual_pct']:.4f}%")
        print(f"    Lowest residual: {best_resid['name']:10s}  H = {best_resid['H_layer']:.3f}"
              f"  resid = {best_resid['residual_pct']:.4f}%")
        print(f"    Match: {match}")
        print()

    print(f"  Score: {wins}/{tests} families have max-H = min-residual")
    pct = wins / tests * 100 if tests > 0 else 0
    print(f"  ({pct:.0f}%)")

    # Random baseline: if residual ranking is independent of H ranking,
    # expected fraction of max-H = min-residual is 1/n for a family of size n
    expected = sum(1.0 / len(m) for m in families.values() if len(m) >= 2)
    n_fam = sum(1 for m in families.values() if len(m) >= 2)
    print(f"\n  Random expectation: {expected:.1f}/{n_fam}"
          f" ({expected/n_fam*100:.0f}%)")


# =========================================================================
# COMBINED SUMMARY
# =========================================================================

def summarise(results):
    """Print combined assessment."""
    print("\n" + "=" * 80)
    print("COMBINED ASSESSMENT")
    print("=" * 80)

    # Universal confinement
    baryons = [r for r in results if r['B'] >= 1 and r['n_temporal'] > 0]
    all_confined = all(r['n_components'] > 1 for r in baryons)
    print(f"\n  Confinement (temporal binding):")
    print(f"    {len(baryons)} baryons tested, ALL fragment without temporal bonds: {all_confined}")

    # Decuplet correlation
    dec = [r for r in results if r['B'] == 1 and r['nc'] == 0
           and r['J'] >= 1.4 and r['residual_pct'] is not None]
    dec.sort(key=lambda r: r['H_layer'])
    resid_dec = [r['residual_pct'] for r in dec]
    H_dec = [r['H_layer'] for r in dec]
    # Check if residuals are monotonically decreasing with H
    monotone = all(resid_dec[i] >= resid_dec[i+1] for i in range(len(resid_dec)-1))
    rho_dec, p_dec = spearmanr([1-h for h in H_dec], resid_dec) if len(dec) >= 3 else (0, 1)
    print(f"\n  Decuplet correlation:")
    print(f"    Spearman ρ = {rho_dec:+.3f} (p = {p_dec:.4f})")
    print(f"    Monotonically decreasing residual with H: {monotone}")
    print(f"    Ω⁻ (H=1.000) residual: {dec[-1]['residual_pct']:.4f}%"
          f" — {'BEST' if dec[-1] == min(dec, key=lambda r: r['residual_pct']) else 'not best'}")

    # Single-layer strangeness
    kaon = [r for r in results if r['name'] in ('K±', 'K*') or 'K' in r['name']]
    zero_temp = [r for r in kaon if r['n_temporal'] == 0]
    print(f"\n  Kaon selection rule:")
    print(f"    Kaon-like particles with zero temporal bonds: {len(zero_temp)}")
    print(f"    → ΔS = 0 for graph-mediated interactions (strong/EM)")


# =========================================================================
# ENTRY POINT
# =========================================================================

if __name__ == '__main__':
    np.random.seed(42)
    print("Building full hadron catalogue...")
    results = build_catalogue()

    run_test_1(results)
    run_test_2(results)
    run_test_3(results)
    summarise(results)
