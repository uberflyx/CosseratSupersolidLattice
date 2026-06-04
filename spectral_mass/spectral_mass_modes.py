#!/usr/bin/env python3
"""
spectral_mass_modes.py
======================
One engine, run across the spectrum. For every hadron whose cluster builder
is wired in, this builds the cluster from FCC coordinates, assembles the same
6n x 6n Cosserat dynamical matrix, projects the quantum-number-selected
channel, weights every eigenmode by microrotation content and shell
concentration, and applies the channel's selection rule. No measured mass
enters the selection.

It then checks, for each hadron, whether the engine reproduces the eigenvalue
the monograph quotes. The point is the audit: where one engine suffices, and
where a sector still uses a bespoke construction.

Channels and selection types used here:
  anchor      proton: the lone microrotation-dominant A_2u root
  stiff       resisted J=1/2: lowest stiff pure-optical root in the channel
  soft        accommodated J=1/2: highest soft microrotation-dominant root
  triplet     J=3/2: lowest stiff pure-optical T_1g triplet (3-fold)
  tensor_lo   tensor meson f2: lowest stiff T_2g root (optical-mix)
  tensor_hi   tensor meson a2: purest stiff T_2g root
  shell_stiff excited baryon: lowest stiff pure-optical mode, full spectrum
              (no irrep projection; the proton-winding-derived modes)
"""
import os
import sys
import numpy as np
# resolve sibling modules (delta_first_principles, proton_first_principles, ...)
# from this script's own directory, so the driver runs from anywhere.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from delta_first_principles import build_cosserat_matrix_two_d, cluster_delta
from proton_first_principles import build_oh_elements, find_perm, build_rep, class_of
from spectral_classifier import (cluster_coord_shell, build_lambda_cluster,
                                  hex_cap_extension_on_inactive_dir, INACTIVE_DIRS,
                                  cluster_born, cluster_hex_cap)

M_E = 0.51099895
ALPHA = 1.0 / 137.0359992
M_0 = M_E / ALPHA
OH = build_oh_elements()

OH_CHAR = {
    'A_1g': {'E':1,'8C_3':1,'6C_4':1,'6C_2':1,'3C_2':1,'i':1,'8S_6':1,'6S_4':1,'3σ_h':1,'6σ_d':1},
    'A_2g': {'E':1,'8C_3':1,'6C_4':-1,'6C_2':-1,'3C_2':1,'i':1,'8S_6':1,'6S_4':-1,'3σ_h':1,'6σ_d':-1},
    'A_2u': {'E':1,'8C_3':1,'6C_4':-1,'6C_2':-1,'3C_2':1,'i':-1,'8S_6':-1,'6S_4':1,'3σ_h':-1,'6σ_d':1},
    'T_1g': {'E':3,'8C_3':0,'6C_4':1,'6C_2':-1,'3C_2':-1,'i':3,'8S_6':0,'6S_4':1,'3σ_h':-1,'6σ_d':-1},
    'T_2g': {'E':3,'8C_3':0,'6C_4':-1,'6C_2':1,'3C_2':-1,'i':3,'8S_6':0,'6S_4':-1,'3σ_h':-1,'6σ_d':1},
}


def mass(N, lam):
    return N * M_0 - N * (4 - lam) * M_E


def channel_modes(coords, oh_char, n_shell, tol=1e-5):
    """(lambda, phi, shell, deg) for modes in the projected O_h-irrep channel."""
    n = len(coords)
    char = OH_CHAR[oh_char]
    P = np.zeros((6 * n, 6 * n)); H = 0
    for R in OH:
        p = find_perm(R, coords)
        if p is None:
            continue
        P += char[class_of(R)] * build_rep(R, p, n); H += 1
    P *= char['E'] / H
    M = build_cosserat_matrix_two_d(coords, 1.0, 1.0, 1.0)
    ev, V = np.linalg.eigh(M)
    out = []; i = 0
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
                shell = float(np.sum(m[:3 * n_shell] ** 2)
                              + np.sum(m[3 * n:3 * n + 3 * n_shell] ** 2))
                out.append((float(ev[i]), phi, shell, j - i))
        i = j
    return out


def full_modes(coords, n_shell, tol=1e-5):
    """(lambda, phi, shell) for every mode, no projection."""
    n = len(coords)
    M = build_cosserat_matrix_two_d(coords, 1.0, 1.0, 1.0)
    ev, V = np.linalg.eigh(M)
    out = []
    for k in range(len(ev)):
        m = V[:, k] / np.linalg.norm(V[:, k])
        phi = float(np.sum(m[3 * n:] ** 2))
        shell = float(np.sum(m[:3 * n_shell] ** 2)
                      + np.sum(m[3 * n:3 * n + 3 * n_shell] ** 2))
        out.append((float(ev[k]), phi, shell))
    return out


def proton_psi():
    """The proton A_2u winding eigenvector (lambda=8.303) on the bare 13-shell."""
    shell = c_shell(); n = len(shell)
    M = build_cosserat_matrix_two_d(shell, 1.0, 1.0, 1.0)
    ev, V = np.linalg.eigh(M)
    k = int(np.argmin(np.abs(ev - 8.303)))
    return V[:, k], n


def embed(psi, n_src, n_dst):
    """Pad a source-cluster mode onto the first n_src nodes of a larger cluster."""
    e = np.zeros(6 * n_dst)
    e[:3 * n_src] = psi[:3 * n_src]
    e[3 * n_dst:3 * n_dst + 3 * n_src] = psi[3 * n_src:6 * n_src]
    return e / np.linalg.norm(e)


def select(coords, n_shell, kind, channel=None, phi_pure=0.9, shell_pure=0.9):
    if kind == 'shell_stiff':
        # excited baryons: the monograph uses the stricter >=95% purity cut
        cand = sorted((l, p, s) for (l, p, s) in full_modes(coords, n_shell)
                      if p >= 0.95 and s >= 0.95 and l > 4)
        return cand[0][0] if cand else None
    if kind == 'overlap':
        # radial excitation: the mode with maximum overlap with the embedded
        # proton winding (the Roper rule). No mass, just eigenvector geometry.
        psi, ns = proton_psi()
        pe = embed(psi, ns, len(coords))
        M = build_cosserat_matrix_two_d(coords, 1.0, 1.0, 1.0)
        ev, V = np.linalg.eigh(M)
        k = int(np.argmax([abs(V[:, j] @ pe) for j in range(len(ev))]))
        return float(ev[k])
    if kind == 'bilayer_triplet':
        # J=3/2 negative-parity: the bilayer-localised stiff phi-dominant
        # triplet, complementary to the shell-localised J=1/2 singlet.
        n = len(coords)
        M = build_cosserat_matrix_two_d(coords, 1.0, 1.0, 1.0)
        ev, V = np.linalg.eigh(M)
        trip = []
        for k in range(len(ev)):
            m = V[:, k] / np.linalg.norm(V[:, k])
            phi = float(np.sum(m[3 * n:] ** 2))
            ext = 1.0 - float(np.sum(m[:3 * n_shell] ** 2)
                              + np.sum(m[3 * n:3 * n + 3 * n_shell] ** 2))
            if ev[k] > 4 and phi >= 0.95 and ext >= 0.80:
                trip.append(float(ev[k]))
        trip = sorted(set(round(x, 3) for x in trip))
        return sum(trip) / len(trip) if trip else None
    modes = channel_modes(coords, channel, n_shell)
    pure = sorted((l, p, s, g) for (l, p, s, g) in modes
                  if p >= phi_pure and s >= shell_pure)
    phidom = sorted((l, p, s, g) for (l, p, s, g) in modes if p > 0.5)
    if kind == 'anchor':
        c = [m for m in pure if m[0] >= 4]; return max(c, key=lambda m: m[0])[0] if c else None
    if kind == 'stiff':
        c = [m for m in pure if m[0] >= 4]; return min(c, key=lambda m: m[0])[0] if c else None
    if kind == 'triplet':
        c = [m for m in pure if m[0] >= 4]; return min(c, key=lambda m: m[0])[0] if c else None
    if kind == 'soft':
        c = [m for m in phidom if m[0] < 4]; return max(c, key=lambda m: m[0])[0] if c else None
    if kind == 'tensor_lo':                      # f2: lowest stiff T_2g root
        c = sorted(m for m in phidom if m[0] >= 4); return c[0][0] if c else None
    if kind == 'tensor_hi':                      # a2: purest stiff T_2g root
        c = [m for m in phidom if m[0] >= 4]; return max(c, key=lambda m: m[1])[0] if c else None
    return None


# ---- cluster builders ----
def c_shell():
    s, _ = cluster_coord_shell(); return s

def c_xi():
    s, _ = cluster_coord_shell()
    cl = np.vstack([s, hex_cap_extension_on_inactive_dir(INACTIVE_DIRS[0]),
                    hex_cap_extension_on_inactive_dir(INACTIVE_DIRS[1])])
    u = []
    for p in cl:
        if not any(np.linalg.norm(p - q) < 1e-6 for q in u):
            u.append(p)
    return np.array(u)

def c_born():
    out = cluster_born(); return out[0] if isinstance(out, tuple) else out

def c_n1535():
    from n1535_first_principles import build_n1535_cluster; return build_n1535_cluster()

def c_secondshell():
    from backward_trace import add_second_shell
    s, _ = cluster_coord_shell(); return add_second_shell(s)

def c_dualorbit():
    from delta1600_dual_orbit import build_dual_orbit_cluster
    out = build_dual_orbit_cluster(); return out[0] if isinstance(out, tuple) else out

def c_hexcap():
    out = cluster_hex_cap(); return out[0] if isinstance(out, tuple) else out


# name, builder, N_mass, n_shell, kind, channel, monograph_lambda, obs_mass, sector
SPEC = [
    ('p,n',       c_shell,      13, 13, 'anchor',  'A_2u', 8.303,  938.92,  'octet'),
    ('Sigma0',    cluster_delta,17, 13, 'stiff',   'A_2g', 4.624,  1193.15, 'octet'),
    ('Lambda',    build_lambda_cluster, 16, 13, 'soft', 'A_2g', 3.204, 1115.68, 'octet'),
    ('Xi',        c_xi,         19, 13, 'soft',    'A_2g', 3.055,  1318.29, 'octet'),
    ('Delta',     cluster_delta,17, 13, 'triplet', 'T_1g', 9.052,  1232.0,  'decuplet'),
    ('f2(1270)',  c_born,       18, 19, 'tensor_lo','T_2g', 5.580, 1275.4,  'tensor'),
    ('a2(1320)',  c_born,       18, 19, 'tensor_hi','T_2g', 10.480,1318.2,  'tensor'),
    ('f1(1285)',  c_born,       18, 19, 'stiff',   'T_1g', 7.006,  1281.86, 'axial'),
    ('a1/b1',     c_born,       18, 19, 'soft',    'T_1g', 0.416,  1230.0,  'axial'),
    ('N(1535)',   c_n1535,      21, 13, 'shell_stiff',   None, 8.046, 1510.0, 'excited-N'),
    ('N(1520)',   c_n1535,      21, 13, 'bilayer_triplet',None, 8.325, 1510.0, 'excited-N'),
    ('N(1440)',   c_secondshell,19, 13, 'overlap',       None, 10.193, 1390.0,'excited-N'),
    ('Delta(1600)',c_dualorbit, 21, 13, 'shell_stiff',   None, 10.216, 1570.0,'excited-N'),
    ('K(494)',    c_hexcap,      7,  7, 'shell_stiff',   None, 5.000,  493.7, 'pseudoscalar'),
]

if __name__ == '__main__':
    print(f"m_0 = {M_0:.4f} MeV\n")
    print(f"{'hadron':12}{'N':>4}{'engine_lam':>11}{'m_pred':>10}{'obs':>9}"
          f"{'resid':>9}{'monograph':>11}{'match?':>8}")
    print("-" * 82)
    for name, build, N, ns, kind, chan, fwk, obs, sec in SPEC:
        try:
            coords = build()
            lam = select(coords, ns, kind, chan)
            if lam is None:
                print(f"{name:12}{N:>4}{'(no mode)':>11}"); continue
            m = mass(N, lam); res = (m - obs) / obs * 100
            match = 'YES' if abs(lam - fwk) < 0.05 else 'no'
            print(f"{name:12}{N:>4}{lam:>11.3f}{m:>10.1f}{obs:>9.1f}"
                  f"{res:>+8.2f}%{fwk:>11.3f}{match:>8}")
        except Exception as e:
            print(f"{name:12}{N:>4}  ERROR: {e}")
