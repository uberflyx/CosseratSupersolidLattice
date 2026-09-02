#!/usr/bin/env python3
"""
meson_baryon_docking.py
=======================
Meson-baryon and meson-decuplet compounds under the docking grammar.

Question: does a compound built from a light meson cluster docked on a baryon
cluster reach the mass of a state sitting tens of MeV below the two-body
threshold, as the kaon-nucleon reading of the Lambda(1405), the kaon-Delta
reading of the Sigma(1670) and the rho-Delta reading of the Delta(1920) would
require?

Method. The guest cluster (seven-node hex cap for the kaon, eleven-node crossed
fault for the rho) is placed at every FCC translate d with |d|^2 <= 12 (units
(a/2)^2) and every one of the 48 O_h orientations relative to the host (the
thirteen-node coordination shell for the nucleon, the seventeen-node void
cluster for the Delta). Sites the two clusters share are kept once in the
geometry; the node count follows the dimer rule, N = N_host + N_guest with
nothing subtracted because neither guest carries a void. For each docking the
two-length Cosserat matrix (unit stiffness on every bond, the discrete curl
weighted by 1/2d per bond) is diagonalised and the compound mode of best
geometric-mean local consistency against the two parents' mass eigenvectors is
reported, together with the additive-rule eigenvalue

    lambda_eff = 4 - [N_h (4 - lambda_h) + N_g (4 - lambda_g)] / (N_h + N_g),

which is the constituent sum, i.e. the threshold.

Result. No docking yields a compound mode below the host parent's own
eigenvalue: the kaon on the shell reads lambda = 8.26 to 10.5 (1444 to 1467 MeV
at N = 20, threshold 1432.7), the kaon on the Delta 9.07 to 11.9 (1743 to 1777
MeV at N = 24, threshold 1728.1), the rho on the Delta 9.0 to 12.8 (2033 to
2086 MeV at N = 28, threshold 2009.6). A compound with no void pair sits at or
above its threshold, so the N m_0 counts (1400.5, 1680.6, 1960.7 MeV) are
resets to the reference eigenvalue and not closures.

Depends on delta_first_principles, composite_clusters, spectral_classifier and
proton_first_principles in this directory.
"""
import numpy as np
from delta_first_principles import cluster_delta, build_cosserat_matrix_two_d
from composite_clusters import cluster_crossed_fault
from spectral_classifier import cluster_coord_shell, fcc_nn_vectors
from proton_first_principles import build_oh_elements

M_E = 0.51099895
ALPHA = 1.0 / 137.035999084
M_0 = M_E / ALPHA
A_LAT = np.sqrt(2.0)
VOID_D = np.sqrt(6.0) / 4.0
OH = build_oh_elements()


def spectrum(coords):
    return np.linalg.eigh(build_cosserat_matrix_two_d(coords))


def parent(coords, target, tol=2e-3):
    """Eigenvectors of the parent's mass mode (a degenerate multiplet is returned whole)."""
    w, V = spectrum(coords)
    idx = np.where(np.abs(w - target) < tol)[0]
    assert len(idx) > 0, f"no parent mode near {target}"
    return V[:, idx], w[idx[0]]


def hex_cap():
    n = np.array([1.0, 1.0, 1.0])
    ring = [v for v in fcc_nn_vectors() if abs(v @ n) < 1e-6]
    return np.vstack([[0.0, 0.0, 0.0]] + ring)


def site_key(v):
    return tuple(np.round(v * 4).astype(int))


def is_bond(d):
    return abs(d - 1.0) < 1e-6 or abs(d - VOID_D) < 1e-6


def fcc_translates(rmax2):
    out = []
    for h in range(-4, 5):
        for k in range(-4, 5):
            for l in range(-4, 5):
                if (h + k + l) % 2 == 0 and h * h + k * k + l * l <= rmax2:
                    out.append(np.array([h, k, l]) * A_LAT / 2)
    return out


def restrict(vec, n, idx):
    u = vec[:3 * n].reshape(n, 3)[idx].ravel()
    p = vec[3 * n:].reshape(n, 3)[idx].ravel()
    return np.concatenate([u, p])


def dock(host, P_host, lam_host, N_host, guest, P_guest, lam_guest, N_guest, label):
    N = N_host + N_guest
    lam_eff = 4 - (N_host * (4 - lam_host) + N_guest * (4 - lam_guest)) / N
    m_eff = N * M_0 - N * (4 - lam_eff) * M_E
    print(f"\n=== {label}: N = {N_host} + {N_guest} = {N}")
    print(f"additive rule: lambda_eff = {lam_eff:.3f}, m = {m_eff:.1f} MeV (the constituent sum)")
    hk = {site_key(v): i for i, v in enumerate(host)}
    nh = len(host)
    results, seen = {}, set()
    for d in fcc_translates(12):
        for R in OH:
            g = (R @ guest.T).T + d
            sk = frozenset(site_key(v) for v in g)
            if sk in seen:
                continue
            seen.add(sk)
            gkeys = [site_key(v) for v in g]
            shared = [k for k in gkeys if k in hk]
            new = [v for v, k in zip(g, gkeys) if k not in hk]
            if not new:
                continue
            ifb = sum(is_bond(np.linalg.norm(x - y)) for x in host for y in new)
            if not shared and ifb == 0:
                continue
            coords = np.vstack([host, np.array(new)])
            n = len(coords)
            newkeys = [site_key(v) for v in new]
            idx_g = [hk[k] if k in hk else nh + newkeys.index(k) for k in gkeys]
            w, V = spectrum(coords)
            best = None
            for k in range(len(w)):
                vh = restrict(V[:, k], n, list(range(nh)))
                vg = restrict(V[:, k], n, idx_g)
                ch = np.linalg.norm(P_host.T @ vh) ** 2 / max(vh @ vh, 1e-12)
                cg = np.linalg.norm(P_guest.T @ vg) ** 2 / max(vg @ vg, 1e-12)
                gm = np.sqrt(ch * cg)
                if best is None or gm > best[0]:
                    best = (gm, w[k], ch, cg)
            results.setdefault((int(round(2 * d @ d)), len(shared), ifb), []).append((best, n))
    for key in sorted(results):
        (gm, lam, ch, cg), n = max(results[key], key=lambda t: t[0][0])
        m = N * M_0 - N * (4 - lam) * M_E
        print(f"|d|^2={key[0]:2d} shared={key[1]} interface bonds={key[2]:2d} sites={n:2d}  "
              f"best GM={gm:.2f} (host {ch:.2f}, guest {cg:.2f})  lambda={lam:.3f}  m={m:.1f} MeV")


if __name__ == "__main__":
    shell, _ = cluster_coord_shell()
    delta = cluster_delta()
    cap = hex_cap()
    rho, _ = cluster_crossed_fault()
    P_N, l_N = parent(shell, 8.3028)
    P_D, l_D = parent(delta, 9.0515)
    P_K, l_K = parent(cap, 5.0, tol=1e-4)
    P_R, l_R = parent(rho, 4.8907)
    print(f"parents: nucleon A2u {l_N:.4f}, Delta T1 {l_D:.4f}, kaon A1g {l_K:.4f}, rho B3u {l_R:.4f}")
    dock(shell, P_N, l_N, 13, cap, P_K, l_K, 7, "kaon hex cap on nucleon shell (Lambda(1405) candidate)")
    dock(delta, P_D, l_D, 17, cap, P_K, l_K, 7, "kaon hex cap on Delta void cluster (Sigma(1670) candidate)")
    dock(delta, P_D, l_D, 17, rho, P_R, l_R, 11, "rho crossed fault on Delta void cluster (Delta(1920) candidate)")
