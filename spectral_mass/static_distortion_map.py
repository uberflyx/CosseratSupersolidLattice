#!/usr/bin/env python3
"""
static_distortion_map.py
========================
Two routes to a baryon's static distortion u0, and what each gives when fed to
the mode sum m = N m_0 + N m_e (<u0|Phi|u0> - 4).

1. Eigenstrain relaxation. A Shockley partial writes the offset r_ij . b on
   every bond crossing its fault patch; the compatible relaxation is the
   least-squares solution u0 = B^+ d with B the compatibility (bond-stretch)
   matrix. For three partials at 120 degrees on one {111} facet the offsets sum
   to zero on every bond, so d = 0 and u0 vanishes: colour closure is zero
   eigenstrain. A single strange offset on the Lambda's cap patch relaxes to
   Rayleigh quotients far from the Lambda's 3.204 and is reported as a negative.

2. Mechanisms. With no compatible part, u0 lies in the kernel of B less the six
   rigid motions. The shell, the void cluster, and the one- and two-triangle
   hosts each carry exactly one mechanism. Re-carried into the axial sector the
   mechanism's Rayleigh quotient on the coupled Cosserat matrix is exactly 8 on
   the shell (0.92 of its weight on the proton's 8.303 root) and exactly 9 on
   the void cluster (the Delta's counted rung, 1233.9 MeV); on the cap hosts it
   sits at 7.60 and 7.31, three per cent above the Lambda and Xi, so the
   accommodated hyperons read something else.

Depends on delta_first_principles and spectral_classifier in this directory.
"""
import numpy as np
from delta_first_principles import build_cosserat_matrix_two_d, cluster_delta
from spectral_classifier import (cluster_coord_shell, build_lambda_cluster,
                                 hex_cap_extension_on_inactive_dir, INACTIVE_DIRS)

M_E = 0.51099895
ALPHA = 1.0 / 137.035999084
M_0 = M_E / ALPHA
VOID_D = np.sqrt(6.0) / 4.0


def bond_list(c):
    n = len(c)
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(c[i] - c[j])
            if abs(d - 1.0) < 1e-6 or abs(d - VOID_D) < 1e-6:
                out.append((i, j, d))
    return out


def compat(c, bl):
    B = np.zeros((len(bl), 3 * len(c)))
    for k, (i, j, d) in enumerate(bl):
        r = (c[j] - c[i]) / d
        B[k, 3 * j:3 * j + 3] = r
        B[k, 3 * i:3 * i + 3] = -r
    return B


def rigid_basis(c):
    n = len(c)
    R = []
    for a in range(3):
        v = np.zeros((n, 3)); v[:, a] = 1.0; R.append(v.ravel())
    for a in range(3):
        ax = np.zeros(3); ax[a] = 1.0; R.append(np.cross(ax, c - c.mean(0)).ravel())
    Q, _ = np.linalg.qr(np.array(R).T)
    return Q


def mechanisms(c):
    B = compat(c, bond_list(c))
    _, sv, vt = np.linalg.svd(B)
    null = vt[np.sum(sv > 1e-9):]
    Q = rigid_basis(c)
    null = null - (null @ Q) @ Q.T
    if null.shape[0] == 0:
        return np.zeros((0, 3 * len(c)))
    _, S, Vt = np.linalg.svd(null)
    return Vt[:np.sum(S > 1e-9)]


def mass(N, lam):
    return N * M_0 - N * (4.0 - lam) * M_E


def rayleigh(M, vec):
    v = vec / np.linalg.norm(vec)
    return v @ M @ v, v


def report_mechanisms(c, label, N):
    M = build_cosserat_matrix_two_d(c)
    w, V = np.linalg.eigh(M)
    n = len(c)
    mech = mechanisms(c)
    print(f"--- {label}: {len(mech)} mechanism(s)")
    for m in mech:
        lam, v = rayleigh(M, np.concatenate([np.zeros(3 * n), m]))
        wts = (V.T @ v) ** 2
        top = np.argsort(wts)[::-1][:3]
        print(f"    axial re-carrier: Rayleigh quotient {lam:.4f}, m = {mass(N, lam):.1f} MeV; "
              f"weights {[(round(w[k], 3), round(wts[k], 2)) for k in top]}")


def eigenstrain(c, far_face, bvec):
    """Least-squares relaxation of Shockley offsets on the bonds crossing the cut
    between far_face (slipped by bvec) and the rest of the cluster."""
    bl = bond_list(c)
    B = compat(c, bl)
    d = np.zeros(len(bl))
    for k, (i, j, dd) in enumerate(bl):
        r = (c[j] - c[i]) / dd
        if (j in far_face) != (i in far_face):
            d[k] = (r @ bvec) if j in far_face else -(r @ bvec)
    return np.linalg.pinv(B) @ d, int(np.sum(np.abs(d) > 1e-12))


def inplane_directions(normal):
    nrm = normal / np.linalg.norm(normal)
    out = []
    for ax in np.eye(3):
        s = ax - (ax @ nrm) * nrm
        out.append(s / np.linalg.norm(s))
    return out


if __name__ == "__main__":
    b = 1.0 / np.sqrt(3.0)   # Shockley step in bond-length units
    shell, _ = cluster_coord_shell()
    lam_cl = build_lambda_cluster()
    xi_cl = np.vstack([shell, hex_cap_extension_on_inactive_dir(INACTIVE_DIRS[0]),
                       hex_cap_extension_on_inactive_dir(INACTIVE_DIRS[1])])
    delta = cluster_delta()

    print("== Eigenstrain on the shell's top facet ==")
    top = [i for i in range(1, 13) if shell[i] @ np.array([1.0, 1.0, 1.0]) > 0.9]
    dirs = inplane_directions(np.array([1.0, 1.0, 1.0]))
    M = build_cosserat_matrix_two_d(shell)
    u1, ncross = eigenstrain(shell, top, b * dirs[0])
    lam_u, _ = rayleigh(M, np.concatenate([u1, np.zeros(39)]))
    print(f"one partial: {ncross} crossing bonds, |u0| = {np.linalg.norm(u1):.4f}, displacement quotient {lam_u:.4f}")
    utot = sum(eigenstrain(shell, top, b * s)[0] for s in dirs)
    print(f"three partials at 120 degrees: |u0| = {np.linalg.norm(utot):.2e}  (colour closure: zero eigenstrain)")

    print("\n== Eigenstrain of one strange partial on the Lambda's cap patch ==")
    ext = [i for i in range(len(lam_cl)) if np.linalg.norm(lam_cl[i]) > 1.2]
    M = build_cosserat_matrix_two_d(lam_cl)
    n = len(lam_cl)
    for s in inplane_directions(np.array(INACTIVE_DIRS[0], float)):
        u, ncross = eigenstrain(lam_cl, ext, b * s)
        lu, _ = rayleigh(M, np.concatenate([u, np.zeros(3 * n)]))
        lp, _ = rayleigh(M, np.concatenate([np.zeros(3 * n), u]))
        print(f"slip along {np.round(s, 3)}: displacement quotient {lu:.4f} ({mass(16, lu):.1f} MeV), "
              f"axial re-carrier {lp:.4f} ({mass(16, lp):.1f} MeV); target 3.204")

    print("\n== Host mechanisms and their axial re-carriers ==")
    report_mechanisms(shell, "coordination shell (proton, target 8.303)", 13)
    report_mechanisms(delta, "shell + 4 voids (Delta, target 9.052)", 17)
    report_mechanisms(lam_cl, "shell + 1 triangle (Lambda, target 3.204)", 16)
    report_mechanisms(xi_cl, "shell + 2 triangles (Xi, target 3.055)", 19)
