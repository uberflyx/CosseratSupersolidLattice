#!/usr/bin/env python3
"""
kk_cluster_matrix.py
====================
The Kaluza-Klein reduced four-dimensional Cosserat cluster matrix.

Each node carries ten degrees of freedom: the displacement u (3), the compact
displacement u_4, the spatial microrotation phi (3, the axial vector of phi_ij)
and the mixed rotations beta_i = phi_i4 (3). The compact direction enters as a
Bloch phase exp(i k_4 x_4) with k_4 d_111 = 2 pi n/3, n in {0, +1, -1}: the
compact derivative is the two-neighbour difference i sin(k_4 d_111) and the
compact Laplacian is 2(1 - cos k_4 d_111). Single-scale convention: unit spring
on every quadratic invariant (central springs on u and on u_4, bond-graph plus
compact Laplacian on every rotation, coupling alpha |phi - omega|^2 on all six
planes with the clamped-frame curl and gradient).

Checks and results:
  * at n = 0 the (u, phi) sector reproduces the 6n cluster spectrum exactly
    (the factorisation theorem in cluster form);
  * on the wound rung the uniform compact translation has Rayleigh quotient
    3 + (1/16n) sum_i |sum_j r_ij|^2, the compact Laplacian plus a boundary
    count, and the u_4-dominant band sits at 3.0 to 3.5 on every host;
  * the accommodated hyperons (Lambda 3.204, Xi 3.055, Omega 3.369) all read
    inside that band, so at leading order m = N (m_0 - m_e): Lambda -0.31 %,
    Xi +0.19 %, Omega -0.25 %; the nucleon and kaon do not (-3.75 %, -1.43 %);
  * the same rung carried to the heavy baryons, core plus dressing plus
    N (m_0 - m_e), gives seven states at mean -0.14 % and rms 0.50 % with no
    eigenvalue used anywhere (Lambda_c to Omega_b and the Xi_cc);
  * the mixed rotations phi_i4 carry a second band whose floor is
    3 + alpha_Cos = 4 plus a boundary dressing, unoccupied in the catalogue.

Depends on spectral_classifier, delta_first_principles and omega_triple_bilayer
in this directory.
"""
import numpy as np
from spectral_classifier import (cluster_coord_shell, build_lambda_cluster, fcc_nn_vectors,
                                 hex_cap_extension_on_inactive_dir, INACTIVE_DIRS)
from delta_first_principles import build_cosserat_matrix_two_d
from omega_triple_bilayer import build_barrel
M_E = 0.51099895
ALPHA = 1.0 / 137.035999084
M_0 = M_E / ALPHA
EPS = np.zeros((3, 3, 3)); EPS[0,1,2] = EPS[1,2,0] = EPS[2,0,1] = 1; EPS[0,2,1] = EPS[2,1,0] = EPS[1,0,2] = -1
ME, M0 = M_E, M_0

def build_kk(coords,n4,alpha=1.0):
    n=len(coords); k4d=2*np.pi*n4/3; s4=np.sin(k4d); L4=2*(1-np.cos(k4d))
    bl=[(i,j) for i in range(n) for j in range(i+1,n) if abs(np.linalg.norm(coords[i]-coords[j])-1)<1e-6]
    # blocks: U (3n), U4 (n), PHI (3n), BETA (3n)
    oU=0; o4=3*n; oP=4*n; oB=7*n; dim=10*n
    M=np.zeros((dim,dim),dtype=complex)
    Lap=np.zeros((n,n))
    for i,j in bl:
        r=(coords[j]-coords[i]); r/=np.linalg.norm(r); o=np.outer(r,r)
        for p,q in [(i,j),(j,i)]:
            M[oU+3*p:oU+3*p+3,oU+3*q:oU+3*q+3]-=o; M[oU+3*p:oU+3*p+3,oU+3*p:oU+3*p+3]+=o
        Lap[i,j]-=1; Lap[j,i]-=1; Lap[i,i]+=1; Lap[j,j]+=1
    # compact central spring on u4, compact Laplacian on all rotations
    M[o4:o4+n,o4:o4+n]+=L4*np.eye(n)
    M[oP:oP+3*n,oP:oP+3*n]+=np.kron(Lap,np.eye(3))+L4*np.eye(3*n)
    M[oB:oB+3*n,oB:oB+3*n]+=np.kron(Lap,np.eye(3))+L4*np.eye(3*n)
    # spatial curl C: omega_k = (1/2) sum_j eps r_ij x u_j  (clamped)
    C=np.zeros((3*n,3*n))
    G=np.zeros((3*n,n))      # gradient of a scalar: (D f)_i,a = (1/2) sum_j r_ij,a f_j
    for i,j in bl:
        r=(coords[j]-coords[i]); r/=np.linalg.norm(r)
        for p,q,sg in [(i,j,1),(j,i,-1)]:
            rh=sg*r
            C[3*p:3*p+3,3*q:3*q+3]+=np.einsum('abc,b->ac',EPS,rh)/2.0
            G[3*p:3*p+3,q]+=rh/2.0
    # omega_i4 = 1/2 (D_i u4 - d4 u_i), d4 -> i s4 : W = (1/2)[G | -i s4 I]
    W=np.zeros((3*n,dim),dtype=complex); W[:,o4:o4+n]=G/2; W[:,oU:oU+3*n]=-1j*s4*np.eye(3*n)/2
    # coupling energies alpha |phi - C u|^2 and alpha |beta - W v|^2
    E1=np.zeros((3*n,dim),dtype=complex); E1[:,oP:oP+3*n]=np.eye(3*n); E1[:,oU:oU+3*n]=-C
    E2=np.zeros((3*n,dim),dtype=complex); E2[:,oB:oB+3*n]=np.eye(3*n); E2-=W
    M+=alpha*(E1.conj().T@E1)+alpha*(E2.conj().T@E2)
    return M
def spectrum(coords,n4):
    M=build_kk(coords,n4); w,V=np.linalg.eigh(M); return w,V,M
def sector_weights(V,n):
    oU=0;o4=3*n;oP=4*n;oB=7*n
    return {'u':np.sum(abs(V[oU:oU+3*n])**2,0),'u4':np.sum(abs(V[o4:o4+n])**2,0),'phi':np.sum(abs(V[oP:oP+3*n])**2,0),'beta':np.sum(abs(V[oB:oB+3*n])**2,0)}
def report(coords, label, N, targets):
    n=len(coords)
    w0,V0,M0_=spectrum(coords,0)
    w6=np.linalg.eigvalsh(build_cosserat_matrix_two_d(coords))
    # check: at n4=0 the (u,phi) sector should reproduce the 6n spectrum
    sw=sector_weights(V0,n); mask=(sw['u']+sw['phi'])>0.999
    print(f"=== {label}: n4=0 factorises? (u,phi)-only eigenvalues match 6n spectrum: {np.allclose(np.sort(w0[mask]),np.sort(w6),atol=1e-8)}  [{mask.sum()} vs {len(w6)}]")
    w1,V1,_=spectrum(coords,1); sw1=sector_weights(V1,n)
    print(f"  wound rung n4=1: {len(w1)} modes, range {w1.min():.3f}..{w1.max():.3f}")
    for t in targets:
        k=np.argmin(abs(w1-t)); print(f"    nearest to target {t}: lambda={w1[k]:.4f} (u {sw1['u'][k]:.2f}, u4 {sw1['u4'][k]:.2f}, phi {sw1['phi'][k]:.2f}, beta {sw1['beta'][k]:.2f})  m={N*M0-N*(4-w1[k])*ME:.1f}")
    # beta-dominant modes (the mixed-rotation sector) at n4=1
    bd=[(round(w1[k],3),round(sw1['beta'][k],2)) for k in range(len(w1)) if sw1['beta'][k]>0.6]
    print(f"    beta-dominant (>0.6) modes at n4=1: {bd[:16]}")
    ud=[(round(w1[k],3),round(sw1['u4'][k],2)) for k in range(len(w1)) if sw1['u4'][k]>0.6]
    print(f"    u4-dominant (>0.6) modes at n4=1: {ud[:12]}")
    return w1,V1,sw1


def boundary_count(c):
    n = len(c)
    tot = 0.0
    for i in range(n):
        v = sum((c[j] - c[i]) / np.linalg.norm(c[j] - c[i]) for j in range(n)
                if abs(np.linalg.norm(c[i] - c[j]) - 1) < 1e-6)
        tot += float(np.dot(v, v))
    return tot / (16 * n)


if __name__ == "__main__":
    cap = np.vstack([[0, 0, 0]] + [v for v in fcc_nn_vectors() if abs(v @ np.array([1, 1, 1.])) < 1e-6])
    shell, _ = cluster_coord_shell()
    lam = build_lambda_cluster()
    xi = np.vstack([shell, hex_cap_extension_on_inactive_dir(INACTIVE_DIRS[0]),
                    hex_cap_extension_on_inactive_dir(INACTIVE_DIRS[1])])
    omega = build_barrel()
    report(cap, "hex cap (kaon)", 7, [5.0])
    report(shell, "coordination shell (nucleon)", 13, [8.303])
    report(lam, "Lambda", 16, [3.204])
    print("\n== the compact rung on the accommodated hosts ==")
    for label, c, N, mobs in [("Lambda", lam, 16, 1115.683), ("Xi", xi, 19, 1318.29),
                              ("Omega barrel", omega, 24, 1672.45),
                              ("shell (control)", shell, 13, 938.919), ("hex cap (control)", cap, 7, 493.677)]:
        n = len(c)
        w, V, M = spectrum(c, 1)
        sw = sector_weights(V, n)
        band = [w[k] for k in range(len(w)) if sw['u4'][k] > 0.5]
        v = np.zeros(10 * n, dtype=complex); v[3 * n:4 * n] = 1 / np.sqrt(n)
        rq = (v.conj() @ M @ v).real
        lreq = 4 + (mobs - N * M_0) / (N * M_E)
        print(f"{label:18s} N={N:2d}: u4 band {min(band):.3f}..{max(band):.3f}; uniform seed {rq:.4f} = 3 + {boundary_count(c):.4f}; "
              f"lambda_req {lreq:.3f}; N(m0-me) = {N*(M_0-M_E):.1f} vs {mobs} ({100*(N*(M_0-M_E)/mobs-1):+.2f}%)")
    print("\n== the rung carried to the heavy baryons ==")
    CHARM, BOT, DRESS = 18 * M_0, 6 * np.pi ** 2 * M_0, 5 / 3 * M_0
    heavy = [("Lambda_c", CHARM + DRESS, 13, 2286.46), ("Xi_c", CHARM + DRESS, 16, 2469.08),
             ("Omega_c", CHARM + DRESS, 19, 2695.2), ("Lambda_b", BOT + DRESS, 19, 5619.57),
             ("Xi_b", BOT + DRESS, 22, 5794.35), ("Omega_b", BOT + DRESS, 25, 6045.8),
             ("Xi_cc", 2 * CHARM + 8 / 3 * M_0, 13, 3621.6)]
    res = []
    for lab, core, N, obs in heavy:
        pred = core + N * (M_0 - M_E)
        res.append(100 * (pred / obs - 1))
        print(f"{lab:10s} N={N:2d}: {pred:8.1f} MeV  obs {obs:8.2f}  {res[-1]:+6.2f}%")
    res = np.array(res)
    print(f"mean {res.mean():+.2f}%, rms {np.sqrt((res ** 2).mean()):.2f}%")
