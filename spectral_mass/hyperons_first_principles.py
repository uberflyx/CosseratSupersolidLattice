#!/usr/bin/env python3
"""hyperons_first_principles.py -- rebuilt 2026-07 (original lost from repo).

Void-swap decuplet clusters and their T_1-daughter mass eigenvalues, plus
the octet daughter clusters (parent minus voids). Verified locks:
  scalar-graph invariants (|E|, lambda_2): Delta (58, 3.511),
  Sigma* (60, 1.292), Xi* (63, 1.087)  [tab:void_swap_ladder]
  daughter means: Sigma* 9.3108, Xi* 9.3307  [decuplet table]
  proton daughter mode 8.3028 on the 13-node shell.
Selection rule (monograph): the three lowest eigenvalues >= 85% microrotation
and >= 90% shell-concentrated; the Delta reads its exact T_d triplet 9.0515.
Note: for the Delta the generic rule also passes phi-dominant modes of other
irreps (4.62, 7.09); use the exact triplet there. Octet Lambda/Xi table
eigenvalues (3.295/2.788) live on the identified-hex-cap geometry
(sec:strange_hex_cap), not on these void-swap daughters: see open_items.
"""
import numpy as np, sys
sys.path.insert(0, '/home/claude/CosseratSupersolidLattice/spectral_mass')
import delta_first_principles as dfp

ELL = dfp.ELL
VOID_D = np.sqrt(6.0)/4.0*ELL

def cluster_voidswap(ns):
    """Centre + 12 shell + (4-ns) voids + ns FCC-continuation cap triangles.
    Cap k replaces void k: three nodes at A_i + A_j (i<j) over face k,
    with {A_1,A_2,A_3} the face vertices of the vacated {111} face."""
    base = dfp.cluster_delta()          # 17 coords: centre, shell, voids
    centre, shell, voids = base[0], base[1:13], base[13:]
    coords = [centre] + list(shell)
    roles = ['centre'] + ['shell']*12
    # face vertices of void k: the 3 shell nodes at VOID_D from it
    for vi in range(ns, 4):
        coords.append(voids[vi]); roles.append('void')
    for vi in range(ns):
        face = [s for s in shell if abs(np.linalg.norm(s - voids[vi]) - VOID_D) < 1e-6]
        assert len(face) == 3, len(face)
        for i in range(3):
            for j in range(i+1, 3):
                coords.append(face[i] + face[j] - centre); roles.append('cap')
    return np.array(coords), roles

def scalar_invariants(coords):
    n = len(coords)
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d = np.linalg.norm(coords[i]-coords[j])
            if abs(d-ELL) < 1e-6 or abs(d-VOID_D) < 1e-6:
                A[i, j] = A[j, i] = 1.0
    L = np.diag(A.sum(1)) - A
    w = np.linalg.eigvalsh(L)
    return int(A.sum()/2), w[1]

def daughters(ns):
    coords, roles = cluster_voidswap(ns)
    M = dfp.build_cosserat_matrix_two_d(coords)
    n = len(coords)
    w, v = np.linalg.eigh(M)
    shell_ids = [i for i, r in enumerate(roles) if r in ('centre', 'shell')]
    picks = []
    for k in range(len(w)):
        vec = v[:, k]
        pf = float(np.sum(vec[3*n:]**2))
        if pf < 0.85:
            continue
        u = vec[:3*n].reshape(n, 3); ph = vec[3*n:].reshape(n, 3)
        sf = float(sum(u[i]@u[i] + ph[i]@ph[i] for i in shell_ids))
        if sf < 0.90:
            continue
        if w[k] < 4.0:      # stiff branch only (mirrors the Delta selection)
            continue
        picks.append((w[k], k, pf, sf))
    picks.sort()
    return coords, roles, w, v, picks[:3]

for ns, nm, E_t, l2_t, lam_t in ((0, "Delta ", 58, 3.511, 9.0515),
                                  (1, "Sigma*", 60, 1.292, 9.3108),
                                  (2, "Xi*   ", 21, 1.087, 9.3307)):
    coords, roles = cluster_voidswap(ns)
    E, l2 = scalar_invariants(coords)
    coords, roles, w, v, tri = daughters(ns)
    lams = [t[0] for t in tri]
    print(f"{nm}: N={len(coords)}  |E|={E} (target {E_t if ns<2 else 63})  "
          f"lam2={l2:.3f} (target {l2_t})")
    print(f"   daughters: {[f'{l:.4f}' for l in lams]}  mean = {np.mean(lams):.4f} "
          f"(target {lam_t});  phi/shell: "
          f"{[(f'{t[2]:.2f}',f'{t[3]:.2f}') for t in tri]}")

# ================= daughter (octet) clusters and their mass modes ==========
def cluster_octet(ncaps):
    """Daughter clusters: parent minus all voids. 13 + 3*ncaps nodes."""
    coords, roles = cluster_voidswap(ncaps)
    keep = [i for i, r in enumerate(roles) if r != 'void']
    return coords[keep], [roles[i] for i in keep]

print("\n=== daughter-mode locks ===")
for ncaps, nm, lam_t, ref in ((0, "p (A2u)", 8.303, 'monograph'),
                               (1, "Lambda ", 3.295, 'HYP table'),
                               (2, "Xi     ", 2.788, 'HYP table')):
    coords, roles = cluster_octet(ncaps)
    M = dfp.build_cosserat_matrix_two_d(coords)
    n = len(coords)
    w, v = np.linalg.eigh(M)
    near = [k for k in range(len(w)) if abs(w[k]-lam_t) < 0.02]
    pf = [f"{float(np.sum(v[:,k][3*n:]**2)):.2f}" for k in near]
    print(f"{nm}: N={n}; modes near {lam_t}:", [f"{w[k]:.4f}" for k in near],
          " phi:", pf)

print("\n=== Lambda/Xi cluster spectra: phi-dominant soft modes (lam < 6) ===")
for ncaps, nm in ((1, "Lambda(16)"), (2, "Xi(19)")):
    coords, roles = cluster_octet(ncaps)
    M = dfp.build_cosserat_matrix_two_d(coords)
    n = len(coords)
    w, v = np.linalg.eigh(M)
    out = []
    for k in range(len(w)):
        if w[k] < 0.05 or w[k] > 6.0:
            continue
        pf = float(np.sum(v[:, k][3*n:]**2))
        if pf >= 0.70:
            out.append(f"{w[k]:.4f}({pf:.2f})")
    print(f"{nm}:", " ".join(out[:14]))

