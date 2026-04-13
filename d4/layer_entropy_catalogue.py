#!/usr/bin/env python3
"""
layer_entropy_catalogue.py — D4 Layer Decomposition of the Full Hadron Catalogue
=================================================================================
Mitchell A. Cox, University of the Witwatersrand

Computes the D4 stacking-layer decomposition for every hadron whose defect
graph is available in the constructive Cosserat mass calculator, and tests the
prediction that the mass-formula residual correlates with the temporal
asymmetry of the layer distribution.

Physics:
    In the D4 interpretation of the FCC lattice, the ABC stacking direction is
    a compact fourth dimension with circumference 3ℓ.  Each node of a defect
    graph occupies one of three stacking layers (A, B, C), determined by:

        layer(r) = round((x + y + z) / 2) mod 3

    where (x, y, z) are Cartesian coordinates in the FCC primitive-vector basis
    a₁ = (1,1,0), a₂ = (1,0,1), a₃ = (0,1,1).

    The normalised Shannon entropy of the layer distribution,

        H = -1/ln(3) × Σᵢ (nᵢ/N) ln(nᵢ/N),

    ranges from 0 (all nodes on one layer) to 1 (equal populations).  The
    variance of the layer populations,

        σ² = (1/3) Σᵢ (nᵢ - N/3)²,

    measures the elastic dipole moment in the compact direction.  The mass
    formula m = N m₀ + Q mₑ is a 3D projection that misses an NLO correction
    proportional to σ²/N.  Particles with H → 1 (balanced) should have smaller
    residuals; particles with H ≪ 1 (imbalanced) should have larger residuals.

    Reference: Sec. layer_decomposition of the monograph (decay chapter).

Usage:
    python layer_entropy_catalogue.py          # full analysis + correlation
    python layer_entropy_catalogue.py --plot    # also produce scatter plot

Dependencies:
    numpy, scipy (for Spearman correlation), matplotlib (optional, for plot)
    cosserat_graph (from parent directory)
"""

import sys, os, math
import numpy as np

# Add parent directory to path for cosserat_graph import
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from cosserat_graph import predict_with_defect, QN, ALPHA, ME, M0

# =========================================================================
# LAYER ASSIGNMENT
# =========================================================================

def assign_layers(defect):
    """Assign each node of a Defect to a stacking layer (0, 1, 2).

    The layer index is computed from the FCC Cartesian coordinates:
        layer = round((x + y + z) / 2) mod 3

    Returns dict: {node_id: layer_index} and the counts [n_A, n_B, n_C].
    """
    assignments = {}
    for node_id in defect.positioned_nodes():
        pos = defect.pos[node_id]
        s = (pos[0] + pos[1] + pos[2]) / 2.0
        layer = int(round(s)) % 3
        assignments[node_id] = layer
    counts = [0, 0, 0]
    for layer in assignments.values():
        counts[layer] += 1
    return assignments, counts


def layer_entropy(counts):
    """Normalised Shannon entropy of the layer distribution.

    H = -1/ln(3) × Σ (nᵢ/N) ln(nᵢ/N)

    Returns 0 if all nodes on one layer, 1 if perfectly balanced.
    """
    N = sum(counts)
    if N == 0:
        return 0.0
    probs = [n / N for n in counts if n > 0]
    H = -sum(p * math.log(p) for p in probs)
    return H / math.log(3)


def layer_variance(counts):
    """Variance of the layer populations: (1/3) Σ (nᵢ - N/3)².

    Measures the elastic dipole moment in the compact direction.
    Returns (variance, variance/N) where the latter is the normalised form.
    """
    N = sum(counts)
    if N == 0:
        return 0.0, 0.0
    mean = N / 3.0
    var = sum((n - mean) ** 2 for n in counts) / 3.0
    return var, var / N


def bond_classification(defect, layer_assignments):
    """Classify bonds as spatial (same layer) or temporal (different layers).

    Returns (n_spatial, n_temporal, spatial_edges, temporal_edges).
    """
    spatial, temporal = [], []
    for a, b in defect.nn_edge_list():
        if a not in layer_assignments or b not in layer_assignments:
            continue
        if layer_assignments[a] == layer_assignments[b]:
            spatial.append((a, b))
        else:
            temporal.append((a, b))
    return len(spatial), len(temporal), spatial, temporal


def spatial_components(defect, layer_assignments):
    """Count connected components of the spatial-only sub-graph.

    If removing temporal bonds disconnects the graph, the particle
    requires temporal binding (confinement).
    """
    import networkx as nx
    pn = defect.positioned_nodes()
    G = nx.Graph()
    G.add_nodes_from(pn)
    for a, b in defect.nn_edge_list():
        if a in layer_assignments and b in layer_assignments:
            if layer_assignments[a] == layer_assignments[b]:
                G.add_edge(a, b)
    return nx.number_connected_components(G)


# =========================================================================
# PDG OBSERVED MASSES (for residual calculation)
# =========================================================================
# Sources: PDG 2024, supplemented by lattice QCD for unmeasured states.
# Key: (B, S, I, J, n_charm, level) → (name, m_obs_MeV)

PDG_MASSES = {
    # Light pseudoscalar mesons (B=0, J=0)
    (0,  0, 1, 0, 0, 1): ("π±",     139.570),
    (0,  0, 0, 0, 0, 1): ("π⁰/η",   547.862),   # η for I=0
    (0,  1, 0.5, 0, 0, 1): ("K±",   493.677),
    (0,  0, 0, 0, 0, 2): ("η'",     957.78),
    # Light vector mesons (B=0, J=1)
    (0,  0, 1, 1, 0, 1): ("ρ",       775.26),
    (0,  0, 0, 1, 0, 1): ("ω",       782.66),
    (0,  1, 0.5, 1, 0, 1): ("K*",   895.55),
    (0,  0, 0, 1, 0, 2): ("φ",      1019.461),
    # Tensor meson (B=0, J=2)
    (0,  0, 0, 2, 0, 1): ("f₂",    1275.5),
    # Octet baryons (B=1, J=1/2)
    (1,  0, 0.5, 0.5, 0, 1): ("p/n",    938.272),
    (1, -1, 0, 0.5, 0, 1): ("Λ",    1115.683),
    (1, -1, 1, 0.5, 0, 1): ("Σ",    1189.37),
    (1, -2, 0.5, 0.5, 0, 1): ("Ξ",    1314.86),
    # Decuplet baryons (B=1, J=3/2)
    (1,  0, 1.5, 1.5, 0, 1): ("Δ",    1232.0),
    (1, -1, 1, 1.5, 0, 1): ("Σ*",   1383.7),
    (1, -2, 0.5, 1.5, 0, 1): ("Ξ*",   1531.80),
    (1, -3, 0, 1.5, 0, 1): ("Ω⁻",   1672.45),
    # Charm mesons (B=0, J=0, n_charm)
    (0,  0, 0.5, 0, 1, 1): ("D",      1869.66),
    (0,  1, 0, 0, 1, 1): ("Ds",     1968.35),
    (0,  0, 0, 0, 2, 1): ("ηc",     2983.9),
    # Charm mesons (B=0, J=1, n_charm)
    (0,  0, 0.5, 1, 1, 1): ("D*",     2010.26),
    (0,  1, 0, 1, 1, 1): ("Ds*",    2112.2),
    (0,  0, 0, 1, 2, 1): ("J/ψ",    3096.9),
    # Charm baryons (B=1, n_charm=1)
    (1,  0, 0, 0.5, 1, 1): ("Λc",    2286.46),
    (1,  0, 1, 0.5, 1, 1): ("Σc",    2452.65),
    (1, -1, 0.5, 0.5, 1, 1): ("Ξc",    2467.71),
    (1, -2, 0, 0.5, 1, 1): ("Ωc",    2695.2),
    # Charm baryons (B=1, n_charm=2)
    (1,  0, 0.5, 0.5, 2, 1): ("Ξcc",   3621.2),
    # Charm baryons (B=1, n_charm=3)
    (1,  0, 0, 0.5, 3, 1): ("Ωccc",  4796.0),  # lattice QCD prediction
}


# =========================================================================
# MAIN CATALOGUE SCAN
# =========================================================================

def build_catalogue():
    """Build every particle, compute layer decomposition, return results.

    Returns list of dicts with all computed quantities.
    """
    results = []

    # All light quantum number combinations (mesons + baryons)
    light_states = []

    # Light mesons: (B, |S|, I, J, P, level)
    for absS in [0, 1]:
        for I in [0, 0.5, 1]:
            for J, P in [(0, -1), (1, -1), (2, +1)]:
                if absS == 1 and I != 0.5: continue
                if absS == 0 and I not in [0, 1]: continue
                if J == 2 and (I != 0 or absS != 0): continue
                for lv in [1, 2]:
                    if lv == 2 and (absS > 0 or I > 0 or J == 2): continue
                    light_states.append((0, absS, I, J, P, lv, 0))

    # Light baryons: (B, |S|, I, J, P, level, n_charm)
    BARYON_SU3 = [(0, 0.5), (0, 1.5), (1, 0), (1, 1), (2, 0.5), (3, 0)]
    DEC_SU3 = {(0, 1.5), (1, 1.0), (2, 0.5), (3, 0.0)}
    for absS, I in BARYON_SU3:
        for J in [0.5, 1.5]:
            if J >= 1.5 and (absS, I) not in DEC_SU3: continue
            light_states.append((1, absS, I, J, +1, 1, 0))

    # Charm mesons
    for nc in [1, 2]:
        for absS in [0, 1]:
            if nc == 2 and absS > 0: continue
            for I in [0, 0.5]:
                if absS == 0 and nc == 1 and I != 0.5: continue
                if absS == 0 and nc == 2 and I != 0: continue
                if absS == 1 and I != 0: continue
                for J in [0, 1]:
                    light_states.append((0, absS, I, J, -1, 1, nc))

    # Charm baryons
    for nc in [1, 2, 3]:
        for absS in range(4 - nc):
            n_light = 3 - nc
            if absS > n_light: continue
            n_ud = n_light - absS
            I_min = (n_ud % 2) * 0.5
            I_max = n_ud / 2.0
            Iv = I_min
            while Iv <= I_max + 0.01:
                light_states.append((1, absS, Iv, 0.5, +1, 1, nc))
                Iv += 1.0

    # Process each state
    seen_masses = set()
    for B, absS, I, J, P, lv, nc in light_states:
        S = -absS if B > 0 else absS
        I3 = I  # use maximal I3 (isospin partner with same graph)
        qn = QN(B=B, S=S, I=I, I3=I3, J=J, P=P, level=lv, n_charm=nc)

        try:
            result, defect, cluster_str = predict_with_defect(qn)
        except Exception:
            continue

        if result.mass <= 0 or result.cluster.startswith('forbidden'):
            continue
        if defect is None:
            continue

        # Skip duplicate masses (isospin partners have same graph)
        mass_key = round(result.mass, 1)
        if mass_key in seen_masses:
            continue
        seen_masses.add(mass_key)

        # Get positioned nodes
        pn = defect.positioned_nodes()
        N_nodes = len(pn)
        if N_nodes < 2:
            continue

        # Layer decomposition
        layer_assign, counts = assign_layers(defect)
        H = layer_entropy(counts)
        var, var_norm = layer_variance(counts)

        # Bond classification
        n_spat, n_temp, _, _ = bond_classification(defect, layer_assign)
        n_comp = spatial_components(defect, layer_assign)

        # Find observed mass for residual
        lookup_key = (B, S, I, J, nc, lv)
        pdg = PDG_MASSES.get(lookup_key, None)

        name = pdg[0] if pdg else result.cluster
        m_obs = pdg[1] if pdg else None
        m_pred = result.mass

        residual_pct = abs(m_pred - m_obs) / m_obs * 100 if m_obs else None

        entry = {
            'name': name,
            'B': B, 'S': S, 'I': I, 'J': J, 'nc': nc, 'level': lv,
            'N_nodes': N_nodes,
            'layers': counts,
            'H_layer': H,
            'var_layer': var,
            'var_norm': var_norm,
            'n_spatial': n_spat,
            'n_temporal': n_temp,
            'n_components': n_comp,
            'm_pred': m_pred,
            'm_obs': m_obs,
            'residual_pct': residual_pct,
            'cluster': result.cluster,
            'N_formula': result.N,
            'Q': result.Q,
        }
        results.append(entry)

    results.sort(key=lambda r: r['m_pred'])
    return results


def print_table(results):
    """Print the full layer decomposition table."""
    print("=" * 120)
    print("D4 LAYER DECOMPOSITION — FULL HADRON CATALOGUE")
    print("=" * 120)
    print(f"  {'Particle':10s} {'B':>1s} {'S':>2s} {'I':>4s} {'J':>4s}"
          f"  {'Nodes':>5s} {'Layers':>10s} {'H':>5s} {'σ²/N':>5s}"
          f"  {'Spat':>4s} {'Temp':>4s} {'Comp':>4s}"
          f"  {'m_pred':>8s} {'m_obs':>8s} {'Resid%':>7s}")
    print(f"  {'-' * 110}")

    for r in results:
        layers_str = f"{r['layers'][0]}:{r['layers'][1]}:{r['layers'][2]}"
        m_obs_str = f"{r['m_obs']:.1f}" if r['m_obs'] else "---"
        resid_str = f"{r['residual_pct']:.4f}" if r['residual_pct'] is not None else "---"
        Js = f"{r['J']:.0f}" if r['J'] == int(r['J']) else f"{int(2*r['J'])}/2"
        print(f"  {r['name']:10s} {r['B']:>1d} {r['S']:>+2d} {r['I']:>4.1f} {Js:>4s}"
              f"  {r['N_nodes']:>5d} {layers_str:>10s} {r['H_layer']:>5.3f} {r['var_norm']:>5.2f}"
              f"  {r['n_spatial']:>4d} {r['n_temporal']:>4d} {r['n_components']:>4d}"
              f"  {r['m_pred']:>8.1f} {m_obs_str:>8s} {resid_str:>7s}")


def correlation_analysis(results):
    """Compute and print correlation between layer entropy and residuals."""
    from scipy.stats import spearmanr, pearsonr

    print("\n" + "=" * 80)
    print("CORRELATION ANALYSIS: Layer Entropy vs Mass Residual")
    print("=" * 80)

    # Filter to particles with known observed masses
    with_obs = [r for r in results if r['residual_pct'] is not None]

    # Separate multi-layer (baryons, vector mesons) from single-layer
    multi = [r for r in with_obs if r['n_temporal'] > 0]
    single = [r for r in with_obs if r['n_temporal'] == 0]

    print(f"\n  Total particles with graphs: {len(results)}")
    print(f"  With observed masses:        {len(with_obs)}")
    print(f"  Multi-layer (temporal bonds): {len(multi)}")
    print(f"  Single-layer (no temporal):   {len(single)}")

    if len(multi) >= 4:
        asym = [1 - r['H_layer'] for r in multi]
        resid = [r['residual_pct'] for r in multi]
        rho_s, p_s = spearmanr(asym, resid)
        rho_p, p_p = pearsonr(asym, resid)
        print(f"\n  Multi-layer particles:")
        print(f"    Spearman  ρ = {rho_s:+.3f}  (p = {p_s:.4f})")
        print(f"    Pearson   r = {rho_p:+.3f}  (p = {p_p:.4f})")

    if len(with_obs) >= 4:
        # Also compute with variance/N as predictor
        var_n = [r['var_norm'] for r in with_obs]
        resid_all = [r['residual_pct'] for r in with_obs]
        rho_v, p_v = spearmanr(var_n, resid_all)
        print(f"\n  All particles (σ²/N vs residual):")
        print(f"    Spearman  ρ = {rho_v:+.3f}  (p = {p_v:.4f})")

    # Key observations
    print("\n  Key observations:")
    best = min(with_obs, key=lambda r: r['residual_pct'])
    worst = max(with_obs, key=lambda r: r['residual_pct'])
    print(f"    Best prediction:  {best['name']:10s}  H = {best['H_layer']:.3f}"
          f"  residual = {best['residual_pct']:.4f}%")
    print(f"    Worst prediction: {worst['name']:10s}  H = {worst['H_layer']:.3f}"
          f"  residual = {worst['residual_pct']:.4f}%")

    # Confinement check
    print("\n  Confinement check (spatial components without temporal bonds):")
    for r in results:
        if r['n_temporal'] > 0 and r['B'] >= 1:
            confined = "YES" if r['n_components'] > 1 else "no"
            print(f"    {r['name']:10s}  {r['n_components']} spatial components"
                  f"  → confined: {confined}")


def make_plot(results):
    """Scatter plot of layer entropy vs mass residual."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    with_obs = [r for r in results if r['residual_pct'] is not None]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Panel (a): H_layer vs residual
    for r in with_obs:
        color = '#2563eb' if r['B'] >= 1 else '#dc2626' if r['n_temporal'] > 0 else '#6b7280'
        marker = 's' if r['B'] >= 1 else 'o'
        ax1.scatter(r['H_layer'], r['residual_pct'], c=color, marker=marker,
                    s=80, edgecolors='k', linewidths=0.5, zorder=3)
        ax1.annotate(r['name'], (r['H_layer'], r['residual_pct']),
                     fontsize=7, ha='left', va='bottom', xytext=(3, 3),
                     textcoords='offset points')

    ax1.set_xlabel('Layer entropy $H_{\\mathrm{layer}}$', fontsize=12)
    ax1.set_ylabel('Mass residual |Δm/m| (%)', fontsize=12)
    ax1.set_title('(a) Layer entropy vs prediction accuracy', fontsize=13)
    ax1.set_xlim(-0.05, 1.1)
    ax1.grid(True, alpha=0.3)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#2563eb',
               markersize=8, label='Baryons'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#dc2626',
               markersize=8, label='Multi-layer mesons'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#6b7280',
               markersize=8, label='Single-layer mesons'),
    ]
    ax1.legend(handles=legend_elements, fontsize=9, loc='upper left')

    # Panel (b): σ²/N vs residual
    for r in with_obs:
        color = '#2563eb' if r['B'] >= 1 else '#dc2626' if r['n_temporal'] > 0 else '#6b7280'
        marker = 's' if r['B'] >= 1 else 'o'
        ax2.scatter(r['var_norm'], r['residual_pct'], c=color, marker=marker,
                    s=80, edgecolors='k', linewidths=0.5, zorder=3)
        ax2.annotate(r['name'], (r['var_norm'], r['residual_pct']),
                     fontsize=7, ha='left', va='bottom', xytext=(3, 3),
                     textcoords='offset points')

    ax2.set_xlabel('Normalised layer variance $\\sigma^2_{\\mathrm{layer}}/N$', fontsize=12)
    ax2.set_ylabel('Mass residual |Δm/m| (%)', fontsize=12)
    ax2.set_title('(b) Temporal asymmetry vs prediction accuracy', fontsize=13)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    outpath = os.path.join(os.path.dirname(__file__), 'layer_entropy_vs_residual.png')
    fig.savefig(outpath, dpi=200, bbox_inches='tight')
    print(f"\n  Plot saved to {outpath}")

    # Also save PDF for the monograph
    pdf_path = outpath.replace('.png', '.pdf')
    fig.savefig(pdf_path, bbox_inches='tight')
    print(f"  PDF saved to {pdf_path}")


# =========================================================================
# ENTRY POINT
# =========================================================================

if __name__ == '__main__':
    print("Building full hadron catalogue with D4 layer decomposition...")
    results = build_catalogue()
    print_table(results)
    correlation_analysis(results)

    if '--plot' in sys.argv:
        make_plot(results)
    else:
        print("\n  (Run with --plot to generate scatter plot)")
