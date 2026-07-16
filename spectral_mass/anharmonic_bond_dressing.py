#!/usr/bin/env python3
"""Anharmonic bond dressing of cluster eigenvalues: K* and phi first test.

Monograph: spectral chapter, "The anharmonic residual band" (eq:bond_dressing,
eq:dressed_eigenvalue_shift). Status: mechanism sized (right band), K* sign
right at one-fifth size, phi sign wrong under the modelled single-partial
field; the derived static-distortion map is the missing input.

Mechanism: the strange partial's static PN field strains each cluster bond;
bond stiffness dresses as f_ij = 1 + g3*eps_ij + 0.5*g4*eps_ij^2 (both the
central spring and the rotational spring; the Cosserat coupling alpha is
left undressed as a first approximation). Eigenvalues shift; masses follow
m = N m0 - N(4-lam) me, so dm = +N me dlam.
"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spectral_classifier import fcc_nn_vectors, cluster_coord_shell, ELL
from cosserat_classifier import build_coupling_block

def coords_shell():
    out = cluster_coord_shell()
    return out[0] if isinstance(out, tuple) else out

def coords_shell_apex():
    return np.vstack([coords_shell(), [[np.sqrt(2.0)*ELL, 0., 0.]]])

def bonds(coords, tol=1e-6):
    bl = []
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):
            if abs(np.linalg.norm(coords[j]-coords[i]) - ELL) < tol:
                bl.append((i, j))
    return bl

def build_dressed(coords, f, alpha=1.0):
    """6n x 6n Cosserat matrix with per-bond stiffness factors f[(i,j)]."""
    n = len(coords)
    UU = np.zeros((3*n, 3*n)); PP = np.zeros((3*n, 3*n))
    for (i, j) in bonds(coords):
        r = coords[j]-coords[i]; rh = r/np.linalg.norm(r)
        D = np.outer(rh, rh)*f[(i, j)]
        UU[3*i:3*i+3, 3*i:3*i+3] += D; UU[3*j:3*j+3, 3*j:3*j+3] += D
        UU[3*i:3*i+3, 3*j:3*j+3] -= D; UU[3*j:3*j+3, 3*i:3*i+3] -= D
        I = np.eye(3)*f[(i, j)]
        PP[3*i:3*i+3, 3*i:3*i+3] += I; PP[3*j:3*j+3, 3*j:3*j+3] += I
        PP[3*i:3*i+3, 3*j:3*j+3] -= I; PP[3*j:3*j+3, 3*i:3*i+3] -= I
    M = np.zeros((6*n, 6*n))
    M[:3*n, :3*n] = UU; M[3*n:, 3*n:] = PP
    uu_c, pp_c, up_c, pu_c = build_coupling_block(coords, alpha=alpha)
    M[:3*n, :3*n] += uu_c; M[3*n:, 3*n:] += pp_c
    M[:3*n, 3*n:] += up_c; M[3*n:, :3*n] += pu_c
    return M

# --- static PN field of one Shockley partial through the cluster centre ---
B_PART = ELL/np.sqrt(3.0)          # Shockley Burgers magnitude
W_PN   = 0.783*(ELL/np.sqrt(3.0))  # PN half-width, w/d = ln(1/alpha)/2pi
NHAT = np.array([1., 1., 1.])/np.sqrt(3.)     # {111} glide-plane normal
BHAT = np.array([1., 1., -2.])/np.sqrt(6.)    # <112> Burgers direction

def u_static(p, scale=1.0):
    """PN slip field: relative jump -> b behind the core, 0 ahead."""
    xn = p @ NHAT                   # distance from glide plane
    xb = p @ BHAT                   # in-plane coord along Burgers
    mag = (scale*B_PART/(2*np.pi))*(np.pi/2 + np.arctan(xb/(W_PN + abs(xn))))
    return 0.5*np.sign(xn if abs(xn) > 1e-12 else 1.0)*mag*BHAT

def strains(coords, scale=1.0):
    eps = {}
    for (i, j) in bonds(coords):
        du = u_static(coords[j], scale) - u_static(coords[i], scale)
        rh = (coords[j]-coords[i])/ELL
        eps[(i, j)] = float(du @ rh)/ELL*ELL  # du.rhat / ell, ell=1 units safe
    return eps

def track(M0, M1, k0):
    """Follow eigenvector k0 of M0 into M1 by maximal overlap."""
    w0, v0 = np.linalg.eigh(M0); w1, v1 = np.linalg.eigh(M1)
    ov = np.abs(v1.T @ v0[:, k0])
    k1 = int(np.argmax(ov))
    return w0[k0], w1[k1], ov[k1]

def run(name, coords, lam_target, N, dm_needed, g3, g4, scale=1.0):
    f0 = {b: 1.0 for b in bonds(coords)}
    M0 = build_dressed(coords, f0)
    w0 = np.linalg.eigh(M0)[0]
    k0 = int(np.argmin(np.abs(w0 - lam_target)))
    eps = strains(coords, scale)
    f1 = {b: 1.0 + g3*e + 0.5*g4*e*e for b, e in eps.items()}
    M1 = build_dressed(coords, f1)
    lam0, lam1, ov = track(M0, M1, k0)
    dlam = lam1 - lam0
    ME = 0.51099895069
    print(f"{name}: lam0={lam0:.4f} (target {lam_target}), overlap={ov:.3f}")
    print(f"  bond strains: max {max(abs(e) for e in eps.values()):+.4f}, "
          f"rms {np.sqrt(np.mean([e*e for e in eps.values()])):.4f}")
    print(f"  dlam = {dlam:+.4f}  ->  dm = {N*ME*dlam:+.2f} MeV "
          f"(needed {dm_needed:+.2f} MeV, dlam {dm_needed/(N*ME):+.3f})\n")

if __name__ == "__main__":
    for g3, g4, tag in ((-6.4, 23.3, "dressed"), (-7.0, 38.1, "bare")):
        print(f"=== gamma3={g3}, gamma4={g4} ({tag}) ===")
        run("K*(892) shell T1u-low", coords_shell(), 1.663, 13, -3.14, g3, g4)
        run("phi(1020) shell+apex", coords_shell_apex(), 8.441, 14, +7.33, g3, g4)
        run("phi doubled-strain", coords_shell_apex(), 8.441, 14, +7.33, g3, g4, scale=2.0)
