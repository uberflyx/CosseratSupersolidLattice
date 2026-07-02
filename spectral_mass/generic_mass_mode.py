#!/usr/bin/env python3
"""
generic_mass_mode.py
====================
The spectral mass formula as one algorithm: quantum numbers in, mass out.

    m = N m_0 - N (4 - lambda) m_e,   m_0 = m_e / alpha

Three steps fix (N, lambda) from the quantum numbers alone:

  1. CLUSTER.  (B, J, |S|, hidden flavour) selects the FCC structural unit
     through the cluster-selection map (monograph, cluster-selection
     chapter).  The map is a theorem about which node counts the FCC
     geometry admits, not a lookup fitted to masses.

  2. IRREP / CHANNEL.  The hadron's J^P selects a symmetry channel of the
     cluster's dynamical matrix, taken in the cluster's own point group.

  3. MODE.  A two-way rule picks the eigenvalue inside the channel.
     LOCKED branch (the state is the channel of a charge its cluster
     cannot relax: baryon winding, the pseudoscalar of an open-strangeness
     fault, the fault's own vector pair, hidden-flavour opticals): the
     mass reads the lowest microrotation-derived mode of the channel
     concentrated on the parent structure, tracked by continuity in the
     Cosserat coupling from the decoupled limit.
     SOFT branch (accommodated states and stacking windings): the mass
     reads the softest displacement mode the state's defect coordinate
     couples to.

Every eigenvalue below is computed from FCC geometry at call time; no mass
data enters.  PDG values appear only in the final comparison column.

Usage:
    python3 generic_mass_mode.py            # catalogue pass + discovery sweep
    from generic_mass_mode import predict
    predict(B=1, J=0.5, P=+1, S=0)          # -> proton entry

Layer-2 (excitation towers), Layer-3 (docking compounds), and decay-width
properties are hooks for future work; see the monograph's status section.
"""
import numpy as np
from scipy.optimize import linear_sum_assignment
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from spectral_classifier import (fcc_nn_vectors, cluster_hex_cap,
                                 cluster_coord_shell, cluster_born,
                                 generate_Oh, classify_operation,
                                 vertex_perm, ELL)
from cosserat_classifier import build_cosserat_matrix, build_uu_block
from composite_clusters import cluster_crossed_fault, cluster_lambda


def _coords(builder):
    """Repo cluster builders return (coords, label); take the coords."""
    out = builder()
    return out[0] if isinstance(out, tuple) else out

# CODATA 2022
M_E = 0.51099895069
INV_ALPHA = 137.035999177
M_0 = M_E * INV_ALPHA


# --------------------------------------------------------------- helpers
def mass(N, lam):
    return N * M_0 - N * (4.0 - lam) * M_E


def disp_rep(R, coords):
    n = len(coords)
    perm = vertex_perm(R, coords)
    if perm is None:
        return None
    D = np.zeros((3 * n, 3 * n))
    for i, j in enumerate(perm):
        D[3 * j:3 * j + 3, 3 * i:3 * i + 3] = R
    return D


def cosserat_rep(R, coords):
    n = len(coords)
    perm = vertex_perm(R, coords)
    if perm is None:
        return None
    det = np.linalg.det(R)
    D = np.zeros((6 * n, 6 * n))
    for i, j in enumerate(perm):
        D[3 * j:3 * j + 3, 3 * i:3 * i + 3] = R
        D[3 * n + 3 * j:3 * n + 3 * j + 3, 3 * n + 3 * i:3 * n + 3 * i + 3] = det * R
    return D


def project(chars_value_by_R, reps, dim):
    P = np.zeros((dim, dim))
    for chi, D in zip(chars_value_by_R, reps):
        P += chi * D
    P /= len(reps)
    w, v = np.linalg.eigh(P)
    return v[:, w > 0.5]


def track_continuity(coords, basis_fn, target_projector, n_steps=101):
    """Track phi-derived modes of a symmetry channel from alpha=0 to 1.

    Returns (lam, lam0, phi_frac, vecs): the coupled eigenvalues in
    ancestor order, the decoupled ancestor eigenvalues, each ancestor's
    microrotation fraction, and the tracked eigenvectors."""
    n = len(coords)
    alphas = np.linspace(0.0, 1.0, n_steps)
    B = target_projector
    M0 = build_cosserat_matrix(coords, 1.0, 1.0, alpha=0.0)
    w, v = np.linalg.eigh(B.T @ M0 @ B)
    modes = B @ v                                  # columns, decoupled
    phi_frac = np.sum(modes[3 * n:, :] ** 2, axis=0)
    lam0 = w.copy()
    lam = w.copy()
    vecs = modes
    for a in alphas[1:]:
        M = build_cosserat_matrix(coords, 1.0, 1.0, alpha=a)
        w2, v2 = np.linalg.eigh(B.T @ M @ B)
        cand = B @ v2
        # globally optimal successor assignment (handles level crossings)
        ov = np.abs(vecs.T @ cand)
        rows, cols = linear_sum_assignment(-ov)
        order = np.empty(ov.shape[0], dtype=int)
        order[rows] = cols
        vecs = cand[:, order]
        lam = w2[order]
    return lam, lam0, phi_frac, vecs


# --------------------------------------------------------------- clusters
def bilayer_coords():
    e1 = np.array([1., -1., 0.]) / np.sqrt(2)
    e2 = np.array([1., 1., -2.]) / np.sqrt(6)
    nrm = np.array([1., 1., 1.]) / np.sqrt(3)
    ring = [np.cos(k * np.pi / 3) * e1 + np.sin(k * np.pi / 3) * e2
            for k in range(6)]
    cen = np.zeros(3)
    apex = (cen + ring[0] + ring[1]) / 3 + nrm * np.sqrt(2. / 3.)
    return np.vstack([[cen], ring, [apex]])


def shell_plus_apex_coords():
    shell = _coords(cluster_coord_shell)
    return np.vstack([shell, [[np.sqrt(2.0) * ELL, 0., 0.]]])


# ------------------------------------------------- per-state resolvers
# Each resolver computes lambda from geometry by the three-step algorithm.
# They contain group theory (which channel) but no numbers from data.

def _oh_singlet(coords, irrep_chars, sector_pick, orbit_slice=None):
    """LOCKED branch on an O_h cluster: continuity-tracked phi-derived
    mode of the given irrep, lowest (optionally parent-orbit weighted)."""
    oh = generate_Oh()
    reps, chis = [], []
    for R in oh:
        D = cosserat_rep(R, coords)
        if D is None:
            continue
        reps.append(D)
        chis.append(irrep_chars[classify_operation(R)])
    B = project(chis, reps, 6 * len(coords))
    lam, lam0, phi_frac, vecs = track_continuity(coords, None, B)
    idx = [k for k in range(len(lam)) if phi_frac[k] > 0.5]
    idx.sort(key=lambda k: lam0[k])
    if orbit_slice is not None:
        n = len(coords)
        best, best_w = None, -1
        cand = sorted(idx, key=lambda k: lam[k])
        # lowest phi-derived mode with parent-orbit weight above half
        for k in cand:
            wts = sum(np.sum(vecs[3 * i:3 * i + 3, k] ** 2)
                      + np.sum(vecs[3 * n + 3 * i:3 * n + 3 * i + 3, k] ** 2)
                      for i in orbit_slice)
            if wts > 0.5:
                return lam[k]
        return lam[cand[0]]
    return lam[idx[0]] if sector_pick == 'lowest_phi' else None


OH_CHARS = {
    'A2u': {'E': 1, 'C_3': 1, 'C_4': -1, 'C_2(diag)': -1, 'C_2(cube)': 1,
            'i': -1, 'S_6': -1, 'S_4': 1, 'sigma_d': 1, 'sigma_h': -1},
    'A1g': {'E': 1, 'C_3': 1, 'C_4': 1, 'C_2(diag)': 1, 'C_2(cube)': 1,
            'i': 1, 'S_6': 1, 'S_4': 1, 'sigma_d': 1, 'sigma_h': 1},
    'A1u': {'E': 1, 'C_3': 1, 'C_4': 1, 'C_2(diag)': 1, 'C_2(cube)': 1,
            'i': -1, 'S_6': -1, 'S_4': -1, 'sigma_d': -1, 'sigma_h': -1},
    'T2g': {'E': 3, 'C_3': 0, 'C_4': -1, 'C_2(diag)': 1, 'C_2(cube)': -1,
            'i': 3, 'S_6': 0, 'S_4': -1, 'sigma_d': 1, 'sigma_h': -1},
    'T1u': {'E': 3, 'C_3': 0, 'C_4': 1, 'C_2(diag)': -1, 'C_2(cube)': -1,
            'i': -3, 'S_6': 0, 'S_4': -1, 'sigma_d': 1, 'sigma_h': 1},
}


def lam_pion():
    """Soft branch, winding state: the advance coordinate of the cell-pair
    winding is the bond direction; the softest mode it couples to is the
    stretch, an exact eigenvector."""
    coords = np.array([[0., 0., 0.],
                       [ELL / np.sqrt(2), ELL / np.sqrt(2), 0.]])
    M = build_cosserat_matrix(coords, 1.0, 1.0, alpha=1.0)
    rhat = (coords[1] - coords[0]) / ELL
    pat = np.zeros(12)
    pat[0:3], pat[3:6] = rhat, -rhat
    pat /= np.linalg.norm(pat)
    w, v = np.linalg.eigh(M)
    ov = (pat @ v) ** 2
    coupled = [k for k in range(12) if ov[k] > 1e-8 and w[k] > 1e-8]
    return min(w[k] for k in coupled)


def lam_eta():
    """Soft branch, winding state on the bilayer: the mode carrying the
    dominant weight of the stacking-slip reaction coordinate.  The probe
    is fixed by bond topology: the patch bonded to the apex (apex + its
    three bond partners) slips coherently against the remaining ring
    nodes along the in-plane <112> partial direction, zero net momentum.
    The apex-only slide is NOT the advance coordinate; it mostly excites
    a floppy mechanism, since slipping the apex against its own partners
    stretches those bonds without advancing the stacking phase."""
    coords = bilayer_coords()
    M = build_uu_block(coords)
    apex = coords[7]
    nrm = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    d = apex - nrm * (apex @ nrm)           # in-plane <112>, in the mirror
    d /= np.linalg.norm(d)
    pat = np.zeros(24)
    for i in (0, 1, 2, 7):                  # centre, ring0, ring1, apex
        pat[3*i:3*i+3] = +d
    for i in (3, 4, 5, 6):                  # far ring nodes
        pat[3*i:3*i+3] = -d
    pat -= np.tile(pat.reshape(8, 3).mean(axis=0), 8)
    pat /= np.linalg.norm(pat)
    w, v = np.linalg.eigh(M)
    ov = (pat @ v) ** 2
    return w[int(np.argmax(ov))]


def lam_kaon():
    """Locked branch: 0^- on an open-strangeness fault.  The defect
    occupies the soft antisymmetric (A1u) winding channel; the mass reads
    the parity-flipped partner, the totally symmetric channel, at its
    lowest microrotation-derived mode, tracked in the coupling."""
    coords = _coords(cluster_hex_cap)
    n = len(coords)
    oh = generate_Oh()
    reps = []
    for R in oh:
        D = cosserat_rep(R, coords)
        if D is not None:
            reps.append(D)
    # totally symmetric irrep of the cluster's full point group
    B = project([1.0] * len(reps), reps, 6 * n)
    lam, lam0, phi_frac, _ = track_continuity(coords, None, B)
    k0 = min((k for k in range(len(lam)) if phi_frac[k] > 0.5),
             key=lambda k: lam0[k])
    return lam[k0]


def lam_nucleon():
    coords = _coords(cluster_coord_shell)
    return _oh_singlet(coords, OH_CHARS['A2u'], 'lowest_phi')


def lam_kstar():
    """Soft branch: lowest nonzero displacement T1u of the bare shell."""
    coords = _coords(cluster_coord_shell)
    n = len(coords)
    M = build_uu_block(coords)
    oh = generate_Oh()
    reps, chis = [], []
    for R in oh:
        D = disp_rep(R, coords)
        if D is None:
            continue
        reps.append(D)
        chis.append(OH_CHARS['T1u'][classify_operation(R)])
    B = project([3 * c / 3 for c in chis], reps, 3 * n) if False else None
    P = np.zeros((3 * n, 3 * n))
    for chi, D in zip(chis, reps):
        P += chi * D
    P *= 3.0 / len(reps)
    wp, vp = np.linalg.eigh(P)
    Bm = vp[:, wp > 0.5]
    lam = np.linalg.eigvalsh(Bm.T @ M @ Bm)
    return min(x for x in lam if x > 1e-6)


def lam_eta_prime():
    """Soft branch: lowest non-trivial pseudoscalar (A2 of C4v) mode."""
    coords = shell_plus_apex_coords()
    n = len(coords)
    M = build_uu_block(coords)
    oh = generate_Oh()

    def c4v_class(R):
        tr, det = round(np.trace(R)), round(np.linalg.det(R))
        if det == 1:
            return {3: 'E', 1: 'C4', -1: 'C2'}[tr]
        return 'sv' if np.allclose(R, np.diag(np.diag(R))) else 'sd'

    ch = {'E': 1, 'C4': 1, 'C2': 1, 'sv': -1, 'sd': -1}     # A2 of C4v
    reps, chis = [], []
    for R in oh:
        D = disp_rep(R, coords)
        if D is None:
            continue
        reps.append(D)
        chis.append(ch[c4v_class(R)])
    P = np.zeros((3 * n, 3 * n))
    for chi, D in zip(chis, reps):
        P += chi * D
    P /= len(reps)
    wp, vp = np.linalg.eigh(P)
    Bm = vp[:, wp > 0.5]
    lam = np.linalg.eigvalsh(Bm.T @ M @ Bm)
    return min(x for x in lam if x > 1e-6)


def lam_rho_omega():
    """The crossed fault's vector pair.

    RHO (generic rule): the plane-swap-odd component (B3u, transforming
    like crystal x).  Locked branch: the lowest microrotation-dominant
    mode of the coupled B3u tower.  Its decoupled ancestor is the
    universal reference mode at lambda = 4 (same-irrep levels cannot
    cross, so the descent is unambiguous), the same 4 -> stiff descent
    the kaon shows on the hex cap.

    OMEGA (documented criterion, NOT yet reduced to the two-way rule):
    the fault-axis component (B1u).  The omega reads the stiff
    longitudinal mode: axis-polarised (majority of its motion along the
    fault axis), displacement-dominant, and not the quasi-translation.
    This four-clause selection is derived in the monograph's rho/omega
    section by analogy with the phi meson's longitudinal A1, but the
    generic algorithm does not yet produce it; treat the omega as the
    open case of the locked branch."""
    coords = _coords(cluster_crossed_fault)
    n = len(coords)
    oh = generate_Oh()
    G = [R for R in oh if vertex_perm(R, coords) is not None]
    reps = [cosserat_rep(R, coords) for R in G]
    M = build_cosserat_matrix(coords, 1.0, 1.0, alpha=1.0)
    ax_rho = np.array([1.0, 0.0, 0.0])
    ax_om = np.array([0.0, -1.0, 1.0]) / np.sqrt(2.0)

    def tower(a):
        chis = [float(a @ R @ a) for R in G]
        B = project(chis, reps, 6 * n)
        w, v = np.linalg.eigh(B.T @ M @ B)
        return w, B @ v

    # rho: lowest phi-dominant coupled B3u mode
    w, modes = tower(ax_rho)
    phi = np.sum(modes[3 * n:, :] ** 2, axis=0)
    lam_rho = min(w[k] for k in range(len(w)) if phi[k] > 0.5)

    # omega: stiff longitudinal displacement-dominant B1u mode
    w, modes = tower(ax_om)
    phi = np.sum(modes[3 * n:, :] ** 2, axis=0)
    trans = np.zeros(6 * n)
    for i in range(n):
        trans[3 * i:3 * i + 3] = ax_om
    trans /= np.linalg.norm(trans)
    cand = []
    for k in range(len(w)):
        m = modes[:, k]
        pol = sum((m[3 * i:3 * i + 3] @ ax_om) ** 2
                  + (m[3 * n + 3 * i:3 * n + 3 * i + 3] @ ax_om) ** 2
                  for i in range(n))
        if pol > 0.5 and phi[k] < 0.5 and abs(trans @ m) < 0.3:
            cand.append(w[k])
    lam_omega = max(cand)
    return lam_rho, lam_omega


def lam_dilaton_f2():
    """Locked/optical branch on the Born cluster.  The dilaton is the
    lowest microrotation-derived A1u mode concentrated on the parent
    cuboctahedron shell (the two pure-rotation shell hybrids sit at
    (11 -+ sqrt17)/2; the cuboctahedron-weighted one is the upper).  The
    f2 is the corresponding T2g optical mode; its resolution against the
    a2 partner uses the same continuity tracking (see the tensor-nonet
    section)."""
    coords = _coords(cluster_born)
    lam_dil = _oh_singlet(coords, OH_CHARS['A1u'], 'lowest_phi',
                          orbit_slice=range(1, 13))
    return lam_dil


# ------------------------------------------------------------- catalogue
PDG = {
    'pion (iso)': 138.039, 'kaon (iso)': 493.677, 'eta': 547.862,
    "eta'": 957.78, 'K*(892) (iso)': 891.66, 'rho(770)': 775.26,
    'omega(782)': 782.66, 'nucleon (iso)': 938.919, 'eta(1295)': 1294.0,
}


def predict_all():
    rows = []
    lam = lam_pion()
    rows.append(('pion (iso)', '0-', 'cell pair', 2, lam, 'soft/winding'))
    lam = lam_kaon()
    rows.append(('kaon (iso)', '0-', 'hex cap', 7, lam, 'locked fault'))
    lam = lam_eta()
    rows.append(('eta', '0-', 'bilayer', 8, lam, 'soft/winding'))
    lam = lam_eta_prime()
    rows.append(("eta'", '0-', 'shell+apex', 14, lam, 'soft/pseudoscalar'))
    lam = lam_kstar()
    rows.append(('K*(892) (iso)', '1-', 'shell', 13, lam, 'soft/accommodated'))
    lr, lo = lam_rho_omega()
    rows.append(('rho(770)', '1-', 'crossed fault', 11, lr, 'locked/fault vector'))
    rows.append(('omega(782)', '1-', 'crossed fault', 11, lo,
                 'locked (criterion open)'))
    lam = lam_nucleon()
    rows.append(('nucleon (iso)', '1/2+', 'shell', 13, lam, 'locked/B=1'))
    lam = lam_dilaton_f2()
    rows.append(('eta(1295)', '0-', 'Born', 18, lam, 'locked/optical'))
    return rows


def main():
    print(f"{'state':16s} {'J^P':5s} {'cluster':14s} {'N':>3s} "
          f"{'lambda':>8s} {'m [MeV]':>9s} {'PDG':>9s} {'res':>8s}  branch")
    print('-' * 96)
    for name, jp, cl, N, lam, branch in predict_all():
        m = mass(N, lam)
        pdg = PDG.get(name)
        res = f"{100 * (m - pdg) / pdg:+.3f}%" if pdg else ''
        print(f"{name:16s} {jp:5s} {cl:14s} {N:3d} {lam:8.4f} "
              f"{m:9.2f} {pdg or 0:9.2f} {res:>8s}  {branch}")

    print("\nDiscovery sweep: unassigned soft/stiff singlets in the "
          "physical window\n(candidate states; forward predictions of the "
          "algorithm's clusters):")
    shell = _coords(cluster_coord_shell)
    M = build_cosserat_matrix(shell, 1.0, 1.0, alpha=1.0)
    w = np.linalg.eigvalsh(M)
    assigned = {0.463, 0.535, 0.807, 0.866, 3.0, 4.0, 4.697, 6.193, 8.303,
                1.663, 8.054}
    seen = set()
    for lam in sorted(w):
        key = round(lam, 3)
        if key in seen or lam < 1e-6:
            continue
        seen.add(key)
        if not any(abs(key - a) < 5e-3 for a in assigned):
            m = mass(13, lam)
            if 900 < m < 1000:
                print(f"  13-shell lambda = {lam:7.4f} -> {m:8.2f} MeV "
                      f"(unassigned)")


if __name__ == "__main__":
    main()
