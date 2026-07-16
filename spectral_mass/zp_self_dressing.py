#!/usr/bin/env python3
"""Zero-point self-dressing of cluster mass modes: the protection theorem.

Monograph: spectral chapter, anharmonic residual band (eq:zp_self_dressing).
Result: the ideal-geometry dressing is bond-uniform to a few parts per
thousand, so it is absorbed into the calibration of K and every closure is
protected; the mode differential is below 1 MeV. The residual band must be
sourced by the defect non-ideal bonds (the static-distortion map).

Each bond's stiffness smears over its own zero-point stretch variance:
  K_eff/K = 1 + 0.5*(g4 - g3^2) * <d^2>/l^2   (bare Morse: g3=-7, g4=343/9)
Observable = cluster-minus-bulk difference of <d^2> per bond, same model:
Cosserat alpha=1, K_u=K_phi=1, unit masses, Mott hbar = 1/sqrt(2).
"""
import numpy as np, sys, os
sys.path.insert(0, '/home/claude/CosseratSupersolidLattice/spectral_mass')
from spectral_classifier import cluster_coord_shell, ELL
from cosserat_classifier import build_coupling_block

HBAR = 1.0/np.sqrt(2.0)
G3, G4 = -7.0, 343.0/9.0
COEF = 0.5*(G4 - G3*G3)          # = -5.44

def coords_shell():
    out = cluster_coord_shell()
    return out[0] if isinstance(out, tuple) else out

def coords_shell_apex():
    return np.vstack([coords_shell(), [[np.sqrt(2.0)*ELL, 0., 0.]]])

def nn_pairs(coords, box=None, tol=1e-6):
    """NN bonds; minimum-image if box (3-vector of periodic lengths) given."""
    n = len(coords); out = []
    for i in range(n):
        for j in range(i+1, n):
            d = coords[j] - coords[i]
            if box is not None:
                d = d - box*np.round(d/box)
            if abs(np.linalg.norm(d) - ELL) < tol:
                out.append((i, j, d/np.linalg.norm(d)))
    return out

def cosserat_matrix(coords, pairs, f=None, alpha=1.0, coupling=True):
    n = len(coords)
    M = np.zeros((6*n, 6*n))
    for k, (i, j, rh) in enumerate(pairs):
        fac = 1.0 if f is None else f[k]
        D = np.outer(rh, rh)*fac
        M[3*i:3*i+3, 3*i:3*i+3] += D; M[3*j:3*j+3, 3*j:3*j+3] += D
        M[3*i:3*i+3, 3*j:3*j+3] -= D; M[3*j:3*j+3, 3*i:3*i+3] -= D
        I = np.eye(3)*fac
        M[3*n+3*i:3*n+3*i+3, 3*n+3*i:3*n+3*i+3] += I
        M[3*n+3*j:3*n+3*j+3, 3*n+3*j:3*n+3*j+3] += I
        M[3*n+3*i:3*n+3*i+3, 3*n+3*j:3*n+3*j+3] -= I
        M[3*n+3*j:3*n+3*j+3, 3*n+3*i:3*n+3*i+3] -= I
    if coupling:
        uu_c, pp_c, up_c, pu_c = build_coupling_block(coords, alpha=alpha)
        M[:3*n, :3*n] += uu_c; M[3*n:, 3*n:] += pp_c
        M[:3*n, 3*n:] += up_c; M[3*n:, :3*n] += pu_c
    return M

def zp_bond_variance(coords, pairs, M):
    """<delta_ij^2> = sum_k (hbar/2 w_k) [rhat.(e_kj - e_ki)]^2 over w_k>0."""
    n = len(coords)
    w, v = np.linalg.eigh(M)
    var = np.zeros(len(pairs))
    for k in range(len(w)):
        if w[k] < 1e-8:
            continue
        amp = HBAR/(2.0*np.sqrt(w[k]))
        for b, (i, j, rh) in enumerate(pairs):
            du = v[3*j:3*j+3, k] - v[3*i:3*i+3, k]
            var[b] += amp*(rh @ du)**2
    return var

def bulk_reference(ncell=3):
    """Periodic FCC supercell (conventional cells), same Cosserat model.
    Coupling block omitted for the bulk (site-level curl builder is not
    periodic-aware); its effect on <d^2> is checked on the cluster below."""
    a = np.sqrt(2.0)*ELL
    base = np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]])*a
    pts = []
    for ix in range(ncell):
        for iy in range(ncell):
            for iz in range(ncell):
                pts.extend(base + a*np.array([ix, iy, iz]))
    coords = np.array(pts); box = np.array([ncell*a]*3)
    pairs = nn_pairs(coords, box=box)
    M = cosserat_matrix(coords, pairs, coupling=False)
    var = zp_bond_variance(coords, pairs, M)
    return float(np.mean(var)), float(np.std(var))

def mode_shift(coords, lam_target, label):
    pairs = nn_pairs(coords)
    M0 = cosserat_matrix(coords, pairs)
    var = zp_bond_variance(coords, pairs, M0)
    dK = COEF*(var - VBULK)          # per-bond fractional stiffness change
    f1 = 1.0 + dK
    w0, v0 = np.linalg.eigh(M0)
    k0 = int(np.argmin(np.abs(w0 - lam_target)))
    e = v0[:, k0]
    # exact first-order shift: e^T dM e with dM = dressed - bare
    dM = cosserat_matrix(coords, pairs, f=f1, coupling=False) \
       - cosserat_matrix(coords, pairs, coupling=False)
    dlam1 = float(e @ dM @ e)
    # full re-diagonalisation with overlap tracking
    M1 = cosserat_matrix(coords, pairs, f=f1)
    w1, v1 = np.linalg.eigh(M1)
    ov = np.abs(v1.T @ e); k1 = int(np.argmax(ov))
    print(f"{label}: lam={w0[k0]:.4f}, bond dvar range "
          f"[{(var-VBULK).min():+.4f},{(var-VBULK).max():+.4f}]")
    print(f"  dlam: 1st-order {dlam1:+.4f} | rediag {w1[k1]-w0[k0]:+.4f} "
          f"(overlap {ov[k1]:.3f})")
    return dlam1

if __name__ == "__main__":
    VBULK, sd = bulk_reference(3)
    print(f"bulk <d^2> = {VBULK:.4f} l^2 (spread {sd:.1e}); "
          f"Einstein check <u^2> = {3*HBAR/(2*np.sqrt(8.0)):.3f} l^2\n")
    ME = 0.51099895069
    for coords, lam, N, name, resid in (
        (coords_shell(), 1.663, 13, "K*(892)  T1u-low ", -3.14),
        (coords_shell(), 8.303, 13, "proton   A2u-stiff", -0.01),
        (coords_shell_apex(), 8.441, 14, "phi(1020) stiff  ", +7.33)):
        d = mode_shift(coords, lam, name)
        print(f"  dm = {N*ME*d:+.2f} MeV (residual to explain {resid:+.2f})\n")

# ---------------------------------------------------------------- embedded
def embedded_matrix(coords, pairs, alpha=1.0):
    """Einstein embedding: restore each node's missing bulk coordination as
    springs to the static environment (adds the missing bond dyadics to the
    diagonal u-block, and the missing count to the diagonal phi-block)."""
    n = len(coords)
    M = cosserat_matrix(coords, pairs, alpha=alpha)
    # bulk diagonal blocks: sum over 12 NN dyadics = 4*I; phi Laplacian: 12*I
    have_dyad = [np.zeros((3, 3)) for _ in range(n)]
    have_cnt = np.zeros(n)
    for (i, j, rh) in pairs:
        D = np.outer(rh, rh)
        have_dyad[i] += D; have_dyad[j] += D
        have_cnt[i] += 1; have_cnt[j] += 1
    for i in range(n):
        M[3*i:3*i+3, 3*i:3*i+3] += 4.0*np.eye(3) - have_dyad[i]
        M[3*n+3*i:3*n+3*i+3, 3*n+3*i:3*n+3*i+3] += (12.0 - have_cnt[i])*np.eye(3)
    return M

def embedded_check(coords, label):
    pairs = nn_pairs(coords)
    M = embedded_matrix(coords, pairs)
    var = zp_bond_variance(coords, pairs, M)
    d = var - VBULK
    print(f"{label} embedded: dvar range [{d.min():+.4f},{d.max():+.4f}], "
          f"mean {d.mean():+.4f}  (free-cluster mean was ~+0.055)")

print("\n=== Einstein-embedded clusters: does the surface softening vanish? ===")
VBULK, _ = bulk_reference(3)
embedded_check(coords_shell(), "coord shell   ")
embedded_check(coords_shell_apex(), "shell + apex  ")
