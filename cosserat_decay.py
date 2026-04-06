#!/usr/bin/env python3
"""
cosserat_decay.py — Decay engine for the Cosserat supersolid lattice framework.

Computes decay widths and stability classifications from the defect graph
constructed by cosserat_graph.py.  No physics is duplicated: all lattice
constants, assembly logic, and Q computations live in cosserat_graph.

Tier 1 (pure graph theory):
    - Laplacian spectrum → Fiedler value λ₂
    - Stability classification: prime (irreducible) vs composite
    - Fiedler partition → decay products
    - Factorisation theorem: Γ = 2(|∂E|·mₑ + |∂V|·m₀/π)

Tier 2 (graph + physics rules from monograph):
    - Decuplet void deactivation: Γ ∝ n_free² / (N_c²−1)
    - Crossed-fault P-wave healing: Γ = p³ / (12π √N_c · m₀²)
    - Molecular docking → composite graph → factorisation width

Reference: Monograph Ch. "Decay rates and lifetimes", §§11.1–11.8.

Usage:
    python3 cosserat_decay.py
"""

import math
import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, List, Dict

from cosserat_graph import (
    QN, Result, Defect, FCCLattice,
    predict, predict_with_defect, predict_molecular,
    _lat, _select_docking_mode, _DOCK_MODES,
    ME, M0, ALPHA,
)

# ================================================================
# CONSTANTS (derived in monograph, not fitted)
# ================================================================
N_c   = _lat.N_c          # 3
Z1    = _lat.Z1            # 12
N_BL  = _lat.N_bilayer     # 8
N_CB  = _lat.N_CB          # 5
N_SUR = _lat.N_surplus     # 27 = N_c³

# Decuplet anchor: Δ(1232) width (Mode I, adjoint screening)
# Derived: Γ_Δ = n_free² / (N_c²-1) × (universal rate)
# The universal rate Γ₀ is fixed by the Δ observation.
GAMMA_DELTA_OBS = 117.0    # MeV, PDG 2024

# Hex-cap shadow factor (Eq. monograph_decay:shadow_factor)
# σ = (1/N_ring) × exp(-d/(3w))  where d=4ℓ/3, w=ℓ/π, N_ring=6
_d_over_w = 4.0 / (3.0 * (1.0 / math.pi))  # = 4π/3
SHADOW = (1.0 / 6.0) * math.exp(-_d_over_w)  # ≈ 0.030

# Factorisation theorem (Eq. monograph_decay:factorisation_width)
# Γ = 2 × (|∂E| × mₑ  +  |∂V| × m₀/π)
# Factor of 2 from optical theorem.


# ================================================================
# DATA: decay result
# ================================================================
@dataclass
class DecayResult:
    """Output of the decay analysis for one state."""
    name: str = ''
    mass: float = 0.0
    cluster: str = ''
    # Graph classification
    n_positioned: int = 0
    n_edges: int = 0
    fiedler: float = 0.0
    classification: str = ''     # 'prime', 'composite', 'void_unstable'
    # Decay channel (if applicable)
    channel: str = ''            # e.g. 'Δ→Nπ', 'ρ→ππ'
    mechanism: str = ''          # 'topologically_stable', 'void_deact', 'fault_healing', 'factorisation'
    width: float = 0.0           # MeV (0 = stable)
    width_obs: float = 0.0       # MeV observed (for comparison)
    notes: list = field(default_factory=list)


# ================================================================
# TIER 1: GRAPH SPECTRAL CLASSIFICATION
# ================================================================

def classify(defect: Defect) -> str:
    """Classify a defect graph as prime or composite.
    
    From the monograph (§parsing_results):
      λ₂ ≥ 2.0  →  'prime' (irreducible, topologically stable)
      λ₂ < 2.0  →  needs further analysis
    
    Void-containing states with free voids are 'void_unstable' — they decay
    by void deactivation (Mode I), not by graph fragmentation.
    """
    fv = defect.free_void_count()
    if fv > 0:
        return 'void_unstable'
    if defect.has_voids and defect.void_count() > 0:
        return 'void_unstable'  # has voids even if extensions block some
    
    lam2, _, _ = defect.fiedler()
    if lam2 >= 2.0 - 1e-6:
        return 'prime'
    else:
        return 'narrow'


def fiedler_partition(defect: Defect) -> Tuple[set, set, int, int]:
    """Partition the defect graph by the Fiedler vector.
    
    Returns:
        (set_A, set_B, n_cut_edges, n_shared_nodes)
    where set_A and set_B are sets of node IDs,
    n_cut_edges is the number of edges crossing the partition,
    n_shared_nodes is the number of nodes on the boundary (v₂ ≈ 0).
    """
    lam2, v2, nodes = defect.fiedler()
    if len(nodes) < 2:
        return set(nodes), set(), 0, 0
    
    set_A = {nodes[i] for i in range(len(nodes)) if v2[i] >= 0}
    set_B = {nodes[i] for i in range(len(nodes)) if v2[i] < 0}
    
    # Count edges crossing the partition
    n_cut = 0
    for a, b in defect.nn_edge_list():
        if (a in set_A and b in set_B) or (a in set_B and b in set_A):
            n_cut += 1
    
    # Shared nodes: those with |v₂| < threshold (sitting on the boundary)
    threshold = 0.05 * max(abs(v2))
    n_shared = sum(1 for v in v2 if abs(v) < threshold)
    
    return set_A, set_B, n_cut, n_shared


def factorisation_width(n_bond: int, n_node: int = 0) -> float:
    """Factorisation theorem width (Eq. monograph §11.3).
    
    Γ = 2 × (n_bond × mₑ  +  n_node × m₀/π)
    
    For molecular exotics:
      Mode F:  n_bond=6,  n_node=0  → Γ = 6.1 MeV
      Mode 2F: n_bond=12, n_node=0  → Γ = 12.3 MeV
      Mode V:  n_bond=27, n_node=1  → Γ = 72.2 MeV
    """
    return 2.0 * (n_bond * ME + n_node * M0 / math.pi)


# ================================================================
# TIER 2a: DECUPLET VOID DEACTIVATION (Mode I)
# ================================================================

def decuplet_width(qn: QN) -> DecayResult:
    """Compute decuplet baryon width from void deactivation.
    
    Physics (monograph §decuplet_widths):
      The rate has three factors:
        1. Void coherence: amplitude ∝ n_free, rate ∝ n_free²
        2. Adjoint screening: × 1/(N_c²−1)
        3. P-wave phase space: × p_CM³/m²
      All anchored to the Δ(1232) observation.
      NLO: hex-cap shadow σ per strange arm.
    """
    res, defect, cl = predict_with_defect(qn)
    if defect is None:
        return DecayResult(name='?', classification='no_defect')
    
    absS = abs(qn.S)
    n_free = max(0, _lat.n_voids - absS)
    
    if n_free == 0:
        return DecayResult(
            name=f'|S|={absS}', mass=res.mass, cluster=cl,
            classification='topologically_stable',
            mechanism='topologically_stable',
            channel='stable (no free voids)',
            notes=['n_free=0: all voids blocked']
        )
    
    m_parent = res.mass
    m_pi = predict(QN(B=0, S=0, I=1, I3=0, J=0, P=-1)).mass
    
    # Identify decay channels: decuplet(I,|S|) → octet(I',|S|) + π(I=1)
    # Isospin selection: I' ranges from |I-1| to I+1 in integer steps, 
    # but the daughter must be a valid octet baryon.
    # Octet baryons: (|S|,I) = (0,1/2), (1,0), (1,1), (2,1/2)
    octet_states = {0: [0.5], 1: [0, 1], 2: [0.5]}
    
    channels = []
    valid_I_daughters = octet_states.get(absS, [])
    for I_d in valid_I_daughters:
        # Check isospin coupling: |I_parent - 1| ≤ I_d ≤ I_parent + 1
        if abs(qn.I - 1) - 0.01 <= I_d <= qn.I + 1 + 0.01:
            I3_d = min(I_d, qn.I3)  # pick a valid I3
            m_oct = predict(QN(B=1, S=qn.S, I=I_d, I3=I3_d, J=0.5, P=1)).mass
            if m_oct > 0 and m_parent > m_oct + m_pi:
                p2 = ((m_parent**2 - (m_oct+m_pi)**2) *
                      (m_parent**2 - (m_oct-m_pi)**2)) / (4*m_parent**2)
                p_cm = math.sqrt(max(0, p2))
                # Isospin sum: CG²=1 per channel when summing over 
                # all charge states of daughter + pion (completeness).
                channels.append((f'→I={I_d}+π', p_cm, 1.0))
    
    if not channels:
        return DecayResult(name=f'|S|={absS}', mass=m_parent, cluster=cl,
                          classification='kinematically_forbidden')
    
    # Anchor: Δ→Nπ (the ONLY Δ channel, CG²=1)
    m_delta = predict(QN(B=1, S=0, I=1.5, I3=1.5, J=1.5, P=1)).mass
    m_N = predict(QN(B=1, S=0, I=0.5, I3=0.5, J=0.5, P=1)).mass
    p_delta = math.sqrt(((m_delta**2-(m_N+m_pi)**2)*(m_delta**2-(m_N-m_pi)**2))
                        / (4*m_delta**2))
    n_free_delta = _lat.n_voids  # = 4
    
    # Rate constant: Γ_Δ = R × n_free_Δ² × CG²_Δ × p_Δ³ / m_Δ²
    # CG²_Δ = 1 (only one channel), so:
    R_univ = GAMMA_DELTA_OBS * m_delta**2 / (n_free_delta**2 * p_delta**3)
    
    # Sum over channels with CG² weights
    gamma_lo = 0.0
    for ch_name, p_ch, cg2 in channels:
        gamma_lo += R_univ * n_free**2 * cg2 * p_ch**3 / m_parent**2
    
    # NLO: hex-cap shadow
    gamma_nlo = gamma_lo * (1.0 - absS * SHADOW)
    
    channels_str = {0: 'Δ→Nπ', 1: 'Σ*→Λπ/Σπ', 2: 'Ξ*→Ξπ'}
    channel = channels_str.get(absS, f'dec(|S|={absS})')
    obs = {0: 117.0, 1: 36.0, 2: 9.1}
    
    return DecayResult(
        name=channel.split('→')[0],
        mass=m_parent, cluster=cl,
        n_positioned=len(defect.positioned_nodes()),
        n_edges=len(defect.nn_edge_list()),
        fiedler=defect.fiedler()[0],
        classification='void_unstable',
        channel=channel,
        mechanism='void_deactivation',
        width=round(gamma_nlo, 1),
        width_obs=obs.get(absS, 0),
        notes=[
            f'n_free={n_free}',
            f'channels: {", ".join(f"{c[0]}(p={c[1]:.0f},CG²={c[2]:.2f})" for c in channels)}',
            f'LO={gamma_lo:.1f}, NLO={gamma_nlo:.1f} MeV'
        ]
    )


# ================================================================
# TIER 2b: CROSSED-FAULT P-WAVE HEALING (Mode I)
# ================================================================

def crossed_fault_width(qn: QN, daughters: list = None) -> DecayResult:
    """Compute crossed-fault meson width from the universal P-wave formula.
    
    Physics (monograph §rho_ab_initio, §N_cancellation):
      Γ = p³_CM / (12π √N_c · m₀²)
    
    The node count N cancels identically between coupling and mass.
    This formula applies to ANY crossed-fault defect decaying to two
    cell pairs (or hex caps for strange mesons).
    
    daughters: list of QN for the two daughter mesons.
               Defaults to [π⁺, π⁻] for ρ, [K⁺, K⁻] for φ, etc.
    """
    res, defect, cl = predict_with_defect(qn)
    if defect is None or cl not in ('crossed_fault', 'crossed_fault_I0',
                                      'bilayer_pair', 'strange_vector'):
        return DecayResult(
            name='?', mass=res.mass if res else 0, cluster=cl or '',
            classification='not_crossed_fault',
            notes=['This formula applies only to crossed-fault mesons.']
        )
    
    m_parent = res.mass
    
    # Default daughters based on quantum numbers
    absS = abs(qn.S)
    if daughters is None:
        if absS == 0 and qn.I >= 1:
            # ρ → π⁺π⁻
            daughters = [
                QN(B=0, S=0, I=1, I3=+1, J=0, P=-1),
                QN(B=0, S=0, I=1, I3=-1, J=0, P=-1)
            ]
        elif absS == 0 and qn.I == 0:
            # ω/φ → special channels
            return DecayResult(
                name='ω/φ', mass=m_parent, cluster=cl,
                classification='isoscalar',
                notes=['Isoscalar: needs OZI/mixing analysis, not simple P-wave.']
            )
        elif absS >= 1:
            # K* → Kπ
            daughters = [
                QN(B=0, S=1, I=0.5, I3=0.5, J=0, P=-1),
                QN(B=0, S=0, I=1, I3=-1, J=0, P=-1)
            ]
    
    m1 = predict(daughters[0]).mass
    m2 = predict(daughters[1]).mass
    
    if m_parent < m1 + m2:
        return DecayResult(
            name='subthreshold', mass=m_parent, cluster=cl,
            classification='kinematically_forbidden'
        )
    
    # CM momentum
    p_cm = math.sqrt(
        ((m_parent**2 - (m1 + m2)**2) *
         (m_parent**2 - (m1 - m2)**2)) / (4 * m_parent**2)
    )
    
    # Universal P-wave width (Eq. universal_Pwave)
    gamma = p_cm**3 / (12.0 * math.pi * math.sqrt(N_c) * M0**2)
    
    # For strange mesons: single active fault plane → ½ coupling²
    n_active_planes = 2 - min(absS, 1)  # 2 for ρ, 1 for K*
    plane_factor = n_active_planes / 2.0
    gamma *= plane_factor
    
    # Observed widths
    obs_map = {
        (0, 1): 147.4,    # ρ
        (1, 0.5): 49.1,   # K*
    }
    w_obs = obs_map.get((absS, qn.I), 0)
    
    name_map = {(0, 1): 'ρ', (1, 0.5): 'K*'}
    name = name_map.get((absS, qn.I), f'V(|S|={absS})')
    
    return DecayResult(
        name=name, mass=m_parent, cluster=cl,
        n_positioned=len(defect.positioned_nodes()),
        n_edges=len(defect.nn_edge_list()),
        fiedler=defect.fiedler()[0],
        classification='crossed_fault',
        channel=f'{name}→daughters',
        mechanism='fault_healing',
        width=round(gamma, 1),
        width_obs=w_obs,
        notes=[
            f'p_CM={p_cm:.1f} MeV',
            f'active_planes={n_active_planes}',
            f'Γ = p³/(12π√N_c·m₀²) × {plane_factor}'
        ]
    )


# ================================================================
# TIER 1+2: MOLECULAR COMPOSITE ANALYSIS
# ================================================================

def molecular_decay(qn_A: QN, qn_B: QN) -> DecayResult:
    """Analyse a molecular exotic state built from two docked hadrons.
    
    1. Determine docking mode from constituent structure
    2. Compute factorisation-theorem width from (n_bond, n_node)
    3. Report the graph-spectral diagnostics
    """
    mode = _select_docking_mode(qn_A, qn_B)
    n_node, n_bond = _DOCK_MODES[mode]
    
    r_A = predict(qn_A)
    r_B = predict(qn_B)
    
    if r_A.mass <= 0 or r_B.mass <= 0:
        return DecayResult(classification='constituent_fail')
    
    m_composite = r_A.mass + r_B.mass - n_node * M0 - n_bond * ME
    gamma = factorisation_width(n_bond, n_node)
    
    return DecayResult(
        name=f'{r_A.cluster}+{r_B.cluster}',
        mass=m_composite,
        cluster=f'molecular_{mode}',
        classification='composite',
        channel=f'mode {mode} → {r_A.cluster} + {r_B.cluster}',
        mechanism='factorisation',
        width=round(gamma, 1),
        notes=[
            f'mode={mode}',
            f'n_bond={n_bond}, n_node={n_node}',
            f'm_A={r_A.mass:.1f}, m_B={r_B.mass:.1f}',
            f'binding={n_node*M0 + n_bond*ME:.1f} MeV'
        ]
    )


# ================================================================
# FULL CATALOGUE SCAN
# ================================================================

# States to analyse, with their QN and known decay info
CATALOGUE = [
    # (name, QN, observed_width_MeV, dominant_channel, decay_type)
    # --- Mesons ---
    ('π⁺',    QN(B=0,S=0,I=1,I3=1,J=0,P=-1),      0,      'stable (weak)', 'stable'),
    ('π⁰',    QN(B=0,S=0,I=1,I3=0,J=0,P=-1),      7.8e-6, 'γγ',            'mode_II'),
    ('K⁺',    QN(B=0,S=1,I=0.5,I3=0.5,J=0,P=-1),  0,      'stable (weak)', 'stable'),
    ('η',     QN(B=0,S=0,I=0,J=0,P=-1),            1.3e-3, 'γγ/3π',         'mode_II'),
    ('ρ',     QN(B=0,S=0,I=1,I3=0,J=1,P=-1),       147.4,  'ππ',            'mode_I'),
    ('ω',     QN(B=0,S=0,I=0,J=1,P=-1),             8.68,  '3π',            'mode_I'),
    ('K*⁺',   QN(B=0,S=1,I=0.5,I3=0.5,J=1,P=-1),  49.1,   'Kπ',            'mode_I'),
    ('η\'',   QN(B=0,S=0,I=0,J=0,P=-1,level=2),    0.196,  'ηππ/ργ',        'mode_I'),
    # --- Baryon octet ---
    ('p',     QN(B=1,S=0,I=0.5,I3=0.5,J=0.5,P=1),  0,     'stable',        'stable'),
    ('n',     QN(B=1,S=0,I=0.5,I3=-0.5,J=0.5,P=1), 7.5e-28,'pe⁻ν̄',        'mode_III'),
    ('Λ',     QN(B=1,S=-1,I=0,I3=0,J=0.5,P=1),     0,      'stable (weak)', 'stable'),
    ('Σ⁺',    QN(B=1,S=-1,I=1,I3=1,J=0.5,P=1),     0,      'stable (weak)', 'stable'),
    ('Ξ⁰',    QN(B=1,S=-2,I=0.5,I3=0.5,J=0.5,P=1), 0,     'stable (weak)', 'stable'),
    # --- Baryon decuplet ---
    ('Δ',     QN(B=1,S=0,I=1.5,I3=1.5,J=1.5,P=1),  117.0, 'Nπ',            'mode_I'),
    ('Σ*⁺',   QN(B=1,S=-1,I=1,I3=1,J=1.5,P=1),     36.0,  'Λπ/Σπ',         'mode_I'),
    ('Ξ*⁰',   QN(B=1,S=-2,I=0.5,I3=0.5,J=1.5,P=1), 9.1,   'Ξπ',            'mode_I'),
    ('Ω⁻',    QN(B=1,S=-3,I=0,I3=0,J=1.5,P=1),     0,     'stable (weak)',  'stable'),
]

# Molecular exotics to analyse
MOLECULAR = [
    # (name, qn_A, qn_B, obs_width, obs_mass)
    ('d*(2380)',
     QN(B=1,S=0,I=1.5,I3=1.5,J=1.5,P=1),  # Δ
     QN(B=1,S=0,I=1.5,I3=-0.5,J=1.5,P=1), # Δ
     70.0, 2380.0),
    ('Pc(4457)',
     QN(B=1,S=0,I=0.5,I3=0.5,J=0.5,P=1,n_charm=1),  # Σ_c (J=1/2)
     QN(B=0,S=0,I=0.5,I3=-0.5,J=1,P=-1,n_charm=1),   # D̄*
     6.4, 4457.0),
    ('Pc(4312)',
     QN(B=1,S=0,I=0.5,I3=0.5,J=0.5,P=1,n_charm=1),  # Σ_c
     QN(B=0,S=0,I=0.5,I3=-0.5,J=0,P=-1,n_charm=1),   # D̄
     9.8, 4312.0),
    ('Pc(4440)',
     QN(B=1,S=0,I=1,I3=1,J=1.5,P=1,n_charm=1),       # Σ_c* (J=3/2)
     QN(B=0,S=0,I=0.5,I3=-0.5,J=1,P=-1,n_charm=1),   # D̄*
     20.6, 4440.0),
]


def analyse_single(name: str, qn: QN, obs_width: float, obs_channel: str,
                    obs_type: str) -> DecayResult:
    """Full analysis of a single-hadron state."""
    res, defect, cl = predict_with_defect(qn)
    
    if defect is None:
        return DecayResult(name=name, mass=res.mass if res else 0,
                          cluster=cl or '', classification='no_graph')
    
    pn = defect.positioned_nodes()
    ne = len(defect.nn_edge_list())
    lam2, _, _ = defect.fiedler()
    cls = classify(defect)
    
    dr = DecayResult(
        name=name, mass=res.mass, cluster=cl,
        n_positioned=len(pn), n_edges=ne,
        fiedler=round(lam2, 4),
        classification=cls,
        width_obs=obs_width,
    )
    
    # Dispatch to appropriate width calculation
    absS = abs(qn.S)
    is_decuplet = qn.B == 1 and qn.J >= 1.5 and qn.P == 1
    n_free_dec = max(0, _lat.n_voids - absS)
    
    if is_decuplet and n_free_dec > 0:
        # Decuplet with free voids/extensions: void deactivation width
        dec = decuplet_width(qn)
        dr.width = dec.width
        dr.channel = dec.channel
        dr.mechanism = 'void_deactivation'
        dr.notes = dec.notes
        dr.classification = 'void_unstable'
    
    elif cl in ('crossed_fault',) and qn.I >= 1:
        # Crossed fault meson: P-wave healing
        cf = crossed_fault_width(qn)
        dr.width = cf.width
        dr.channel = cf.channel
        dr.mechanism = 'fault_healing'
        dr.notes = cf.notes
    
    elif cl in ('strange_vector',):
        # K*: strange crossed-fault healing
        cf = crossed_fault_width(qn)
        dr.width = cf.width
        dr.channel = cf.channel
        dr.mechanism = 'fault_healing'
        dr.notes = cf.notes
    
    elif cls == 'prime':
        dr.mechanism = 'topologically_stable'
        dr.channel = 'stable (graph-prime)'
        dr.width = 0.0
    
    else:
        dr.mechanism = obs_type
        dr.channel = obs_channel
    
    return dr


def run_catalogue():
    """Run the full decay analysis catalogue and print results."""
    
    print("=" * 110)
    print("COSSERAT DECAY ENGINE — graph-spectral analysis + Tier 1/2 width calculations")
    print("=" * 110)
    
    # ── Single hadrons ──
    print(f"\n{'State':<8} {'mass':>7} {'|V|':>4} {'|E|':>4} {'λ₂':>7} "
          f"{'class':>12} {'mechanism':>18} {'Γ_pred':>8} {'Γ_obs':>8} {'Δ':>7} {'channel'}")
    print("-" * 110)
    
    for name, qn, obs_w, obs_ch, obs_type in CATALOGUE:
        dr = analyse_single(name, qn, obs_w, obs_ch, obs_type)
        
        # Residual
        if dr.width > 0 and dr.width_obs > 0:
            delta = f'{(dr.width - dr.width_obs)/dr.width_obs*100:+.1f}%'
        elif dr.width == 0 and dr.width_obs == 0:
            delta = '—'
        else:
            delta = '—'
        
        w_pred = f'{dr.width:.1f}' if dr.width > 0.01 else ('stable' if dr.width == 0 else f'{dr.width:.1e}')
        w_obs = f'{dr.width_obs:.1f}' if dr.width_obs > 0.01 else ('stable' if dr.width_obs == 0 else f'{dr.width_obs:.1e}')
        
        print(f"{name:<8} {dr.mass:>7.1f} {dr.n_positioned:>4} {dr.n_edges:>4} "
              f"{dr.fiedler:>7.3f} {dr.classification:>12} {dr.mechanism:>18} "
              f"{w_pred:>8} {w_obs:>8} {delta:>7} {dr.channel}")
    
    # ── Molecular exotics ──
    print(f"\n{'='*110}")
    print("MOLECULAR EXOTICS — composite graph analysis")
    print(f"{'='*110}")
    print(f"\n{'State':<12} {'m_pred':>7} {'m_obs':>7} {'mode':>6} "
          f"{'n_bond':>6} {'n_node':>6} {'Γ_pred':>8} {'Γ_obs':>8} {'Δ':>7}")
    print("-" * 85)
    
    for name, qnA, qnB, obs_w, obs_m in MOLECULAR:
        md = molecular_decay(qnA, qnB)
        
        mode = md.cluster.split('_')[-1] if '_' in md.cluster else '?'
        nb = 0; nn = 0
        for note in md.notes:
            if 'n_bond=' in note:
                parts = note.split(',')
                for p in parts:
                    if 'n_bond=' in p: nb = int(p.split('=')[1].strip())
                    if 'n_node=' in p: nn = int(p.split('=')[1].strip())
        
        if md.width > 0 and obs_w > 0:
            delta = f'{(md.width - obs_w)/obs_w*100:+.1f}%'
        else:
            delta = '—'
        
        print(f"{name:<12} {md.mass:>7.1f} {obs_m:>7.1f} {mode:>6} "
              f"{nb:>6} {nn:>6} {md.width:>8.1f} {obs_w:>8.1f} {delta:>7}")
    
    # ── Summary ──
    print(f"\n{'='*110}")
    print("SUMMARY")
    print(f"{'='*110}")
    
    # Collect all results with predictions
    results = []
    for name, qn, obs_w, obs_ch, obs_type in CATALOGUE:
        dr = analyse_single(name, qn, obs_w, obs_ch, obs_type)
        if dr.width > 0 and dr.width_obs > 0:
            results.append((name, dr.width, dr.width_obs))
    
    for name, qnA, qnB, obs_w, obs_m in MOLECULAR:
        md = molecular_decay(qnA, qnB)
        if md.width > 0 and obs_w > 0:
            results.append((name, md.width, obs_w))
    
    print(f"\n  Predictions with observed widths: {len(results)}")
    if results:
        residuals = [(p - o) / o * 100 for _, p, o in results]
        print(f"  Mean |residual|: {np.mean(np.abs(residuals)):.1f}%")
        print(f"  Max  |residual|: {np.max(np.abs(residuals)):.1f}%")
        for name, pred, obs in results:
            r = (pred - obs) / obs * 100
            print(f"    {name:<12}  pred={pred:>8.1f}  obs={obs:>8.1f}  Δ={r:+.1f}%")


# ================================================================
# GRAPH INVARIANTS TABLE (for the monograph)
# ================================================================

def print_graph_invariants():
    """Print the Laplacian spectrum and Fiedler values for all single hadrons.
    Extends monograph Table parsing_results from 6 to full catalogue."""
    
    print("\n" + "=" * 90)
    print("GRAPH INVARIANTS — Laplacian spectrum for every assembled defect")
    print("=" * 90)
    print(f"\n{'State':<10} {'|V|':>4} {'|E|':>4} {'λ₂':>7} {'λ_max':>7} "
          f"{'h(G)':>7} {'prime?':>6} {'cluster'}")
    print("-" * 70)
    
    all_states = [
        ('π⁺',    QN(B=0,S=0,I=1,I3=1,J=0,P=-1)),
        ('K⁺',    QN(B=0,S=1,I=0.5,I3=0.5,J=0,P=-1)),
        ('η',     QN(B=0,S=0,I=0,J=0,P=-1)),
        ('η\'',   QN(B=0,S=0,I=0,J=0,P=-1,level=2)),
        ('ρ',     QN(B=0,S=0,I=1,I3=0,J=1,P=-1)),
        ('ω',     QN(B=0,S=0,I=0,J=1,P=-1)),
        ('K*',    QN(B=0,S=1,I=0.5,I3=0.5,J=1,P=-1)),
        ('φ',     QN(B=0,S=0,I=0,J=1,P=-1,level=2)),
        ('p',     QN(B=1,S=0,I=0.5,I3=0.5,J=0.5,P=1)),
        ('n',     QN(B=1,S=0,I=0.5,I3=-0.5,J=0.5,P=1)),
        ('Λ',     QN(B=1,S=-1,I=0,I3=0,J=0.5,P=1)),
        ('Σ⁺',    QN(B=1,S=-1,I=1,I3=1,J=0.5,P=1)),
        ('Ξ⁰',    QN(B=1,S=-2,I=0.5,I3=0.5,J=0.5,P=1)),
        ('Δ⁺⁺',   QN(B=1,S=0,I=1.5,I3=1.5,J=1.5,P=1)),
        ('Σ*⁺',   QN(B=1,S=-1,I=1,I3=1,J=1.5,P=1)),
        ('Ξ*⁰',   QN(B=1,S=-2,I=0.5,I3=0.5,J=1.5,P=1)),
        ('Ω⁻',    QN(B=1,S=-3,I=0,I3=0,J=1.5,P=1)),
        ('Λc⁺',   QN(B=1,S=0,I=0,I3=0,J=0.5,P=1,n_charm=1)),
        ('Σc⁺⁺',  QN(B=1,S=0,I=1,I3=1,J=0.5,P=1,n_charm=1)),
        ('Ωc⁰',   QN(B=1,S=-2,I=0,I3=0,J=0.5,P=1,n_charm=1)),
    ]
    
    for name, qn in all_states:
        res, defect, cl = predict_with_defect(qn)
        if defect is None:
            print(f"{name:<10} — (no defect graph)")
            continue
        
        pn = defect.positioned_nodes()
        ne = len(defect.nn_edge_list())
        evals = defect.laplacian_spectrum()
        lam2, _, _ = defect.fiedler()
        lam_max = float(evals[-1]) if len(evals) > 0 else 0
        cls = classify(defect)
        is_prime = 'Y' if cls == 'prime' else 'N'
        
        # Cheeger lower bound: λ₂/2
        h_lower = lam2 / 2.0
        h_str = f'{h_lower:.3f}+'
        
        print(f"{name:<10} {len(pn):>4} {ne:>4} {lam2:>7.3f} {lam_max:>7.3f} "
              f"{h_str:>7} {is_prime:>6} {cl}")


# ================================================================
# MAIN
# ================================================================

if __name__ == '__main__':
    run_catalogue()
    print_graph_invariants()
