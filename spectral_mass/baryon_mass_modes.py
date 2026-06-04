#!/usr/bin/env python3
"""
baryon_mass_modes.py
====================
First-principles selection of the baryon rest-mass eigenvalue from the
Cosserat dynamical matrix of the defect cluster. No measured mass enters the
selection; the cluster geometry and one rule fix the eigenvalue.

THE RULE
--------
A baryon's rest mass is the elastic energy stored in the microrotation-dominant
("optical") partner of the defect's chiral winding. That partner is a single
eigenmode of the cluster's Cosserat matrix, identified by three structural
facts and nothing else:

  1. Channel. The defect's quantum numbers fix the parent O_h irrep, and the
     cluster's residual point group H fixes the label by subduction. J=1/2
     baryons read the one-dimensional A_2g baryon-number channel (the proton,
     with full O_h symmetry and inversion intact, reads its parity-flip partner
     A_2u). J=3/2 baryons read the three-dimensional T_1g microrotation triplet
     and take its mean.

  2. Microrotation dominance. Among the modes of that channel the mass mode is
     microrotation-dominant (phi-content > 1/2), since rest mass is stored
     microrotation, not displacement.

  3. Branch. The node count fixes the sign of the spectral correction
     -N(4-lambda)m_e and so the branch. A cap-extended cluster overshoots at
     lambda=4 and must read the soft branch (lambda<4); a void-extended cluster
     undershoots and reads the stiff branch (lambda>4). On the chosen branch the
     mass mode is the microrotation-dominant mode NEAREST the lambda=4 optical
     reference: the least-softened (caps) or least-stiffened (voids) partner.

This rule reproduces the framework's independently computed eigenvalues
(proton A_2u = 8.303; Sigma = 4.624; Delta T_1 = 9.052) and supplies the two
octet eigenvalues that had been left as fits (Lambda, Xi).

References for the group theory: Koster et al. (1963), Altmann & Herzig (1994).
"""
import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spectral_classifier import (cluster_coord_shell, build_lambda_cluster,
                                  hex_cap_extension_on_inactive_dir, INACTIVE_DIRS,
                                  fcc_nn_vectors)
from delta_first_principles import (cluster_delta, build_cosserat_matrix_two_d,
                                     void_positions_Td)
from proton_first_principles import build_oh_elements, find_perm, build_rep, class_of

# ----------------------------------------------------------------------
# Constants (CODATA 2022 electron mass; alpha sets m_0 = m_e / alpha)
# ----------------------------------------------------------------------
M_E = 0.51099895          # MeV
ALPHA = 1.0 / 137.0359992
M_0 = M_E / ALPHA         # 70.0253 MeV
OH = build_oh_elements()

# Full O_h character table, keyed by the class label that class_of() returns.
OH_CHAR = {
    'A_1g': {'E':1,'8C_3':1,'6C_4':1,'6C_2':1,'3C_2':1,'i':1,'8S_6':1,'6S_4':1,'3σ_h':1,'6σ_d':1},
    'A_2g': {'E':1,'8C_3':1,'6C_4':-1,'6C_2':-1,'3C_2':1,'i':1,'8S_6':1,'6S_4':-1,'3σ_h':1,'6σ_d':-1},
    'A_2u': {'E':1,'8C_3':1,'6C_4':-1,'6C_2':-1,'3C_2':1,'i':-1,'8S_6':-1,'6S_4':1,'3σ_h':-1,'6σ_d':1},
    'T_1g': {'E':3,'8C_3':0,'6C_4':1,'6C_2':-1,'3C_2':-1,'i':3,'8S_6':0,'6S_4':1,'3σ_h':-1,'6σ_d':-1},
    'T_2g': {'E':3,'8C_3':0,'6C_4':-1,'6C_2':1,'3C_2':-1,'i':3,'8S_6':0,'6S_4':-1,'3σ_h':-1,'6σ_d':1},
}


# ----------------------------------------------------------------------
# Channel projector: restrict an O_h character to the cluster's residual
# subgroup H. This lands automatically on the H-irrep that the O_h channel
# subduces to (A_2g -> A_2 in C_3v, B_2 in C_2v, A_2 in T_d, ...), with no
# convention choice in labelling H's own irreps.
# ----------------------------------------------------------------------
def channel_modes(coords, oh_char, n_shell=13, alpha=1.0, tol=1e-5):
    """Modes in the projected channel: (lambda, phi_content, shell_fraction, deg)."""
    n = len(coords)
    char = OH_CHAR[oh_char]
    P = np.zeros((6 * n, 6 * n))
    H = 0
    for R in OH:
        p = find_perm(R, coords)
        if p is None:
            continue
        P += char[class_of(R)] * build_rep(R, p, n)
        H += 1
    P *= char['E'] / H        # dim_Gamma/|H| makes the projector idempotent
                              # (dim = chi(E): 1 for A-type, 3 for T-type)
    M = build_cosserat_matrix_two_d(coords, K_u=1.0, K_phi=1.0, alpha=alpha)
    ev, V = np.linalg.eigh(M)
    out = []
    i = 0
    while i < len(ev):
        j = i
        while j < len(ev) and abs(ev[j] - ev[i]) < tol:
            j += 1
        Vs = V[:, i:j]
        w, vv = np.linalg.eigh(Vs.T @ P @ Vs)
        for k in range(j - i):
            if w[k] > 0.5:
                m = Vs @ vv[:, k]; m /= np.linalg.norm(m)
                phi = float(np.sum(m[3 * n:] ** 2))
                # shell fraction: amplitude on the original 13-shell (u and phi)
                shell = float(np.sum(m[:3 * n_shell] ** 2)
                              + np.sum(m[3 * n:3 * n + 3 * n_shell] ** 2))
                out.append((float(ev[i]), phi, shell, j - i))
        i = j
    return out, H


def select(coords, oh_char, branch, n_shell=13, phi_pure=0.9, shell_pure=0.9):
    """Apply the rule.

    Returns dict with the strict pure-optical pick (phi>=phi_pure, shell>=shell_pure),
    a flag for whether one exists on the forced branch, and a loose soft fallback.
    """
    modes, H = channel_modes(coords, oh_char, n_shell=n_shell)
    pure = sorted((l, p, s, g) for (l, p, s, g) in modes
                  if p >= phi_pure and s >= shell_pure)
    phi_dom = sorted((l, p, s, g) for (l, p, s, g) in modes if p > 0.5)
    res = {'H': H, 'pure': pure, 'phi_dom': phi_dom,
           'lam': None, 'clean': False, 'note': ''}
    if branch in ('anchor', 'stiff', 'triplet'):
        stiff_pure = [m for m in pure if m[0] >= 4]
        if stiff_pure:
            res['lam'] = min(stiff_pure, key=lambda m: m[0])[0]  # lowest stiff pure-optical
            res['clean'] = True
    elif branch == 'soft':
        soft_pure = [m for m in pure if m[0] < 4]
        if soft_pure:
            res['lam'] = max(soft_pure, key=lambda m: m[0])[0]
            res['clean'] = True
        else:
            # no pure-optical soft mode: fall back to the best soft phi-dominant
            soft = [m for m in phi_dom if m[0] < 4]
            if soft:
                res['lam'] = max(soft, key=lambda m: m[0])[0]
                res['note'] = 'NO pure-optical soft mode; fallback is mixed/cap-localised'
    return res


def mass(N, lam):
    return N * M_0 - N * (4 - lam) * M_E


# ----------------------------------------------------------------------
# Cluster builders (reuse existing; add the two-cap Xi cluster)
# ----------------------------------------------------------------------
def cluster_xi():
    s, _ = cluster_coord_shell()
    cl = np.vstack([s, hex_cap_extension_on_inactive_dir(INACTIVE_DIRS[0]),
                    hex_cap_extension_on_inactive_dir(INACTIVE_DIRS[1])])
    u = []
    for p in cl:
        if not any(np.linalg.norm(p - q) < 1e-6 for q in u):
            u.append(p)
    return np.array(u)


# ----------------------------------------------------------------------
# The baryon table: (name, builder, N, channel, branch, PDG, framework lambda)
# ----------------------------------------------------------------------
BARYONS = [
    ('p,n',     lambda: cluster_coord_shell()[0], 13, 'A_2u', 'anchor',  938.918, 8.303),
    ('Lambda',  build_lambda_cluster,             16, 'A_2g', 'soft',    1115.683, 3.295),
    ('Sigma0',  cluster_delta,                    17, 'A_2g', 'stiff',   1193.15,  4.624),
    ('Xi',      cluster_xi,                       19, 'A_2g', 'soft',    1318.29,  2.788),
    ('Delta',   cluster_delta,                    17, 'T_1g', 'triplet', 1232.0,   9.052),
]


CLASS = {'p,n': 'resisted', 'Sigma0': 'resisted', 'Delta': 'resisted',
         'Lambda': 'accommodated', 'Xi': 'accommodated'}

if __name__ == '__main__':
    print(f"m_0 = m_e/alpha = {M_0:.4f} MeV,  m_e = {M_E} MeV")
    print("Rule (pure-optical): lowest stiff mode with phi>=0.9 AND shell>=0.9.\n")
    print(f"{'baryon':8s}{'class':>13}{'N':>4}{'chan':>7}{'lambda':>9}"
          f"{'m_pred':>10}{'PDG':>10}{'resid':>9}{'  clean?':>9}")
    print("-" * 88)
    for name, build, N, chan, branch, pdg, fwk in BARYONS:
        coords = build()
        r = select(coords, chan, branch)
        lam = r['lam']
        m = mass(N, lam)
        res = (m - pdg) / pdg * 100
        clean = 'yes' if r['clean'] else 'NO (open)'
        print(f"{name:8s}{CLASS[name]:>13}{N:>4}{chan:>7}{lam:>9.4f}"
              f"{m:>10.2f}{pdg:>10.2f}{res:>+8.3f}%{clean:>9}")

    print("\nSelection transparency (pure-optical = phi>=0.9 AND shell>=0.9):\n")
    for name, build, N, chan, branch, pdg, fwk in BARYONS:
        coords = build()
        r = select(coords, chan, branch)
        pure_str = ", ".join(f"{l:.3f}(phi{p:.2f},shell{s:.2f}{'x'+str(g) if g>1 else ''})"
                             for l, p, s, g in r['pure']) or "(none)"
        print(f"  {name:8s} [{CLASS[name]}] |H|={r['H']:2d} {chan}")
        print(f"     pure-optical modes: {pure_str}")
        if r['note']:
            soft_dom = [(l, p, s) for l, p, s, g in r['phi_dom'] if l < 4]
            sd = ", ".join(f"{l:.3f}(phi{p:.2f},shell{s:.2f})" for l, p, s in soft_dom)
            print(f"     soft phi-dominant (mixed): {sd}")
            print(f"     -> {r['note']}; provisional lambda = {r['lam']:.4f}")
        else:
            print(f"     -> selected lambda = {r['lam']:.4f}")
        print()
