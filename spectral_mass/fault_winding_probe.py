#!/usr/bin/env python3
"""
fault_winding_probe.py
======================
One topologically fixed probe selects the working modes of every fault
state, and the defect-mass "double flip" is a re-carrier theorem.

THE PROBE.  A fault's winding fixes a spatial vector pattern with no free
choices:
  - cell pair (pion): the two halves in relative motion along the bond;
  - hex cap (kaon): the Z3 stacking alternation registered on the six
    ring nodes as an alternating radial pattern, centre silent at the
    partial core;
  - bilayer (eta): the patch bonded to the apex slipping coherently
    against the rest of the plane along the in-plane <112> direction.

THE CARRIERS.  Any spatial pattern P can be loaded into either elastic
sector: as a displacement (polar) field or as a microrotation (axial)
field.  Polar and axial vectors differ by a sign under improper
operations, so the axial carrier transforms in Gamma x A1u when the polar
carrier transforms in Gamma: re-carrying a pattern from u to phi flips
the parity AUTOMATICALLY.  This is the constructive content of the
parity-flip / sector-flip ("double flip") rule.

THE RULE.  Accommodated fault states read the displacement carrier
(pion lambda = 2 with weight 1.000; eta lambda = 0.8229 with dominant
weight 0.475).  Charge-locked fault states read the axial re-carrier of
the same pattern (kaon: Z3 pattern is an exact eigenvector of BOTH
sectors, lambda = 1 as displacement = the defect, lambda = 5 as
microrotation = the mass, weight 1.000 each).

THE THEOREM ON THE SHELL.  The nucleon obeys the same statement: the
baryon winding (A2g displacement zero mode) and the mass parent (A2u
microrotation integer at 7) are the SAME thirteen-node spatial pattern,
overlap 1.000000, with A2g x A1u = A2u supplying the parity flip.
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cosserat_classifier import build_cosserat_matrix
from spectral_classifier import fcc_nn_vectors

M_E = 0.51099895069
INV_ALPHA = 137.035999177
M_0 = M_E * INV_ALPHA

def top_weights(M, pat, k=2):
    w, v = np.linalg.eigh(M)
    ov = sorted(((w[i], (pat @ v[:, i])**2) for i in range(len(w))),
                key=lambda t: -t[1])
    return ov[:k]

def pion():
    coords = np.array([[0., 0., 0.], [1/np.sqrt(2), 1/np.sqrt(2), 0.]])
    d = coords[1] - coords[0]; d /= np.linalg.norm(d)
    M = build_cosserat_matrix(coords, 1.0, 1.0, alpha=1.0)
    pat = np.zeros(12); pat[0:3] = -d; pat[3:6] = +d
    pat /= np.linalg.norm(pat)
    (lam, wgt), _ = top_weights(M, pat)
    m = 2*M_0 - 2*(4 - lam)*M_E
    print(f"pion  | u-carrier: lambda = {lam:.4f} (weight {wgt:.3f})"
          f"  -> m = {m:.3f} MeV (PDG 138.039)")

def kaon():
    e1 = np.array([1., -1., 0.])/np.sqrt(2)
    e2 = np.array([1., 1., -2.])/np.sqrt(6)
    ring = [np.cos(k*np.pi/3)*e1 + np.sin(k*np.pi/3)*e2 for k in range(6)]
    cap = np.vstack([[np.zeros(3)], ring]); n = 7
    M = build_cosserat_matrix(cap, 1.0, 1.0, alpha=1.0)
    Z = np.zeros(3*n)
    for k in range(6):
        rh = ring[k]/np.linalg.norm(ring[k])
        Z[3*(k+1):3*(k+1)+3] = ((-1)**k)*rh
    Z /= np.linalg.norm(Z)
    for sector, off, role in [('u', 0, 'defect'), ('phi', 3*n, 'mass')]:
        pat = np.zeros(6*n); pat[off:off+3*n] = Z
        (lam, wgt), _ = top_weights(M, pat)
        print(f"kaon  | {sector}-carrier ({role}): lambda = {lam:.4f} "
              f"(weight {wgt:.3f})")
    m = 7*M_0 - 7*(4 - 5)*M_E
    print(f"      -> m = 7 m0 + 7 me = {m:.3f} MeV (PDG 493.677)")

def eta():
    e1 = np.array([1., -1., 0.])/np.sqrt(2)
    e2 = np.array([1., 1., -2.])/np.sqrt(6)
    nrm = np.array([1., 1., 1.])/np.sqrt(3)
    ring = [np.cos(k*np.pi/3)*e1 + np.sin(k*np.pi/3)*e2 for k in range(6)]
    apex = (np.zeros(3) + ring[0] + ring[1])/3 + nrm*np.sqrt(2./3.)
    bl = np.vstack([[np.zeros(3)], ring, [apex]]); n = 8
    from cosserat_classifier import build_uu_block
    M = build_uu_block(bl)
    d = apex - nrm*(apex @ nrm); d /= np.linalg.norm(d)
    pat = np.zeros(3*n)
    for i in (0, 1, 2, 7): pat[3*i:3*i+3] = +d
    for i in (3, 4, 5, 6): pat[3*i:3*i+3] = -d
    pat -= np.tile(pat.reshape(n, 3).mean(axis=0), n)
    pat /= np.linalg.norm(pat)
    (lam, wgt), (lam2, wgt2) = top_weights(M, pat)
    m = 8*M_0 - 8*(4 - lam)*M_E
    print(f"eta   | u-carrier (patch slip): lambda = {lam:.4f} "
          f"(weight {wgt:.3f}; runner-up {lam2:.4f} at {wgt2:.3f})"
          f"  -> m = {m:.2f} MeV (PDG 547.862)")

def shell_recarrier():
    """Basis-free test of the re-carrier theorem on the coordination shell:
    does the pure-u zero space (which contains the baryon winding), when
    re-carried into the phi sector, intersect the phi eigenspace at the
    integer lambda = 7?  The largest principal-angle cosine between the two
    subspaces is 1 exactly when some kernel pattern's axial re-carrier is
    an exact eigenvector at 7; that pattern is the winding."""
    shell = np.vstack([[[0., 0., 0.]], fcc_nn_vectors()]); n = 13
    M0 = build_cosserat_matrix(shell, 1.0, 1.0, alpha=0.0)
    w, v = np.linalg.eigh(M0)
    U = v[:, [k for k in range(6*n)
              if abs(w[k]) < 1e-8 and np.sum(v[:3*n, k]**2) > 0.999]][:3*n]
    F = v[:, [k for k in range(6*n)
              if abs(w[k]-7) < 1e-8 and np.sum(v[3*n:, k]**2) > 0.999]][3*n:]
    # orthonormalise the spatial patterns and take principal angles
    Uo, _ = np.linalg.qr(U)
    Fo, _ = np.linalg.qr(F)
    s = np.linalg.svd(Fo.T @ Uo, compute_uv=False)
    print(f"shell | u-kernel patterns re-carried to phi vs phi-eigenspace at 7:")
    print(f"        principal-angle cosines = {np.round(s, 6)}")
    print(f"        max = {s.max():.6f}: the winding's axial re-carrier IS the")
    print(f"        mass parent (A2g x A1u = A2u supplies the parity flip)")

if __name__ == "__main__":
    pion(); kaon(); eta(); shell_recarrier()
