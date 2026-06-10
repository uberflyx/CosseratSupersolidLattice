#!/usr/bin/env python3
"""
hyperon_tensor_admixture.py
===========================
Does the quadrupole (couple-stress) admixture of the baryon mass mode grow
with strangeness? The octet ladder needs the answer: the per-swap linear
costs after bilinear and magnetic removal are 3.73 (N), 7.52 +- 0.07
(Sigma), 11.46 +- 0.41 (Xi) in units of m_e, against the targets 4, 8, 12
of an exact 4(1+n_s) reading. That reading survives only if the
spin-independent admixture per swap reaches about -0.25 MeV in BOTH
hyperons, three times the nucleon's, equal across two different hostings
(arm-riding strangeness in the Sigma, cap-riding in the Xi). This script
computes the admixture and closes the question.

Method
------
Second-order perturbation theory on each cluster's mass mode. The mode
(lambda_0, psi_0) is selected by the baryon_mass_modes rule. A charge
configuration perturbs the matrix as V = g * sum_i q_i C_i with C_i a unit
coupling operator at quark site i. Exactly in g:

  lambda = lambda_0 + g (a . q) + g^2 (q . S . q) + O(g^3),
  a_i  = <0|C_i|0>,
  S_ij = sum_{k!=0} <0|C_i|k><k|C_j|0> / (lambda_0 - lambda_k).

The splitting's second order is the bilinear form difference between the
two members; per swap it is the spin-independent admixture in lambda units
per g^2, and N m_e maps it to MeV (the spectral mass differential). The PT
is validated against exact diagonalisation over three decades of g.

Coupling models (pre-registered, physically motivated)
------------------------------------------------------
  u_iso   : isotropic displacement stiffness g q I on the site (the choice
            of nucleon_tensor_magnitude.py; charge as a displacement-sector
            stiffness defect).
  phi_iso : isotropic microrotation stiffness g q I on the site (charge as
            a rotational stiffness defect).
  chiral  : local modulation of the Cosserat coupling alpha (phi - omega)^2
            by g q at the site. A screw converts translation into rotation,
            and charge IS screw content (q = cos beta), so this is the most
            faithful model.
  (A central-force bond-stretch model annihilates the mass modes exactly:
   the phi-dominant modes carry no axial bond displacement, so charge does
   not couple as a stretch. Established and excluded.)

Charge placements follow the realised hostings of the isospin ledger
(hadrons/isospin_isotensor_count.py): light quarks on the mutually-NN
{111} face sites; the Sigma's strange on the third face site of the bare
17-node shell+voids cluster; the Xi's two strange charges on the two
cap sites adjacent to the light arm of the 19-node two-cap cluster.

Result
------
No coupling produces the uniform enhancement. Ratios of N*d2 per swap to
the nucleon: u_iso (+0.60, +2.35), phi_iso (-0.15, -0.02), chiral
(-0.13, +0.24). A single-scale fit to the required (-0.136, -0.245,
-0.276) MeV is rejected at chi^2 >= 19 on three points for every model.
The obstruction is structural: the nucleon mode sits at the TOP of its
channel, so level repulsion only pushes it down (sign-protected); the
Sigma (4.624) and Xi (3.055) modes sit mid-spectrum, where contributions
from above and below cancel, leaving admixtures that are small and of
unstable sign. Anchored to the nucleon's -0.08 MeV through the chiral
coupling (g = 0.24), the hyperon terms are below 0.02 MeV either way.

Verdict: the ladder step is NOT a quadrupole dressing. After bilinear and
magnetic removal the step stands at 3.79 +- 0.07 m_e, three standard
deviations below the territory base of four; its origin is the strange
structure's own linear (Peierls-Nabarro) sector, which remains the open
derivation.

Accompanies: M. A. Cox, "The Cosserat Supersolid" (2026).
https://doi.org/10.5281/zenodo.18636501

Author: Mitchell A. Cox, University of the Witwatersrand
"""
import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spectral_classifier import cluster_coord_shell
from delta_first_principles import cluster_delta, build_cosserat_matrix_two_d
from baryon_mass_modes import select, OH, OH_CHAR, cluster_xi, M_E
from proton_first_principles import find_perm, build_rep, class_of

QU, QD, QS = 2.0/3.0, -1.0/3.0, -1.0/3.0
VOID_D = np.sqrt(6.0) / 4.0

EPS = np.zeros((3, 3, 3))
for a, b, c in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    EPS[a, b, c] = 1.0
for a, b, c in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
    EPS[a, b, c] = -1.0


# ---------------------------------------------------------------- mass mode
def mode_vec(coords, chan, branch):
    """Mass mode (lambda_0, psi_0) plus the full spectrum, selected by the
    baryon_mass_modes rule with the eigenvector recovered through the
    channel projector."""
    n = len(coords)
    r = select(coords, chan, branch)
    lam0 = r['lam']
    char = OH_CHAR[chan]
    P = np.zeros((6*n, 6*n))
    for R in OH:
        p = find_perm(R, coords)
        if p is None:
            continue
        P += char[class_of(R)] * build_rep(R, p, n)
    P *= char['E'] / r['H']
    M0 = build_cosserat_matrix_two_d(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    ev, V = np.linalg.eigh(M0)
    idx = [k for k in range(len(ev)) if abs(ev[k] - lam0) < 1e-5]
    Vs = V[:, idx]
    w, vv = np.linalg.eigh(Vs.T @ P @ Vs)
    psi0 = Vs @ vv[:, -1]
    psi0 /= np.linalg.norm(psi0)
    return lam0, psi0, ev, V, M0


# ---------------------------------------------------------- coupling models
def curl_row(coords, i, tol=1e-6):
    """3 x 3n discrete-curl operator at site i: omega_i = sum over bonds of
    (1/2d) rhat x u_j, over both bond lengths of the cluster."""
    n = len(coords)
    R = np.zeros((3, 3*n))
    for j in range(n):
        if j == i:
            continue
        d = np.linalg.norm(coords[j] - coords[i])
        if abs(d - 1.0) < tol or abs(d - VOID_D) < tol:
            rh = (coords[j] - coords[i]) / d
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        if EPS[a, b, c]:
                            R[a, 3*j + c] += EPS[a, b, c] * rh[b] / (2*d)
    return R


def site_coupling(coords, i, model):
    """Unit coupling operator C_i for one charge at site i."""
    n = len(coords)
    C = np.zeros((6*n, 6*n))
    if model == 'u_iso':
        C[3*i:3*i+3, 3*i:3*i+3] = np.eye(3)
    elif model == 'phi_iso':
        C[3*n+3*i:3*n+3*i+3, 3*n+3*i:3*n+3*i+3] = np.eye(3)
    elif model == 'chiral':
        # d/d(alpha_i) of the Cosserat term alpha_i (phi_i - omega_i)^2
        R = curl_row(coords, i)
        C[:3*n, :3*n] += R.T @ R
        C[3*n+3*i:3*n+3*i+3, 3*n+3*i:3*n+3*i+3] += np.eye(3)
        C[3*n+3*i:3*n+3*i+3, :3*n] -= R
        C[:3*n, 3*n+3*i:3*n+3*i+3] -= R.T
    else:
        raise ValueError(model)
    return C


# ------------------------------------------------------------- perturbation
def S_form(coords, sites, chan, branch, model):
    """First-order vector a and second-order bilinear S over the charge
    sites, in lambda units per g and per g^2."""
    lam0, psi0, ev, V, _ = mode_vec(coords, chan, branch)
    Cs = [site_coupling(coords, i, model) for i in sites]
    a = np.array([psi0 @ C @ psi0 for C in Cs])
    T = np.array([V.T @ (C @ psi0) for C in Cs])
    mask = np.abs(ev - lam0) > 1e-8
    S = np.einsum('im,jm->ij', T[:, mask] / (lam0 - ev[mask]), T[:, mask])
    return lam0, a, S


def validate_pt(coords, sites, q_hi, q_lo, chan, branch, model, g=1e-2):
    """Exact-diagonalisation check of the second-order splitting via the
    even finite difference [lam(g) + lam(-g) - 2 lam0] / (2 g^2), with
    the perturbed mode tracked by maximal overlap with psi_0."""
    lam0, psi0, ev, V, M0 = mode_vec(coords, chan, branch)
    n = len(coords)
    out = {}
    for tag, q in (('hi', q_hi), ('lo', q_lo)):
        Bs = []
        for s in (+g, -g):
            Vq = sum(s * qi * site_coupling(coords, i, model)
                     for qi, i in zip(q, sites))
            e2, V2 = np.linalg.eigh(M0 + Vq)
            k = int(np.argmax((V2.T @ psi0) ** 2))
            Bs.append(e2[k])
        out[tag] = (Bs[0] + Bs[1] - 2*lam0) / (2*g*g)
    return out['hi'] - out['lo']


# --------------------------------------------------------------------- main
def main():
    shell, _ = cluster_coord_shell()
    sig = cluster_delta()
    xi = cluster_xi()
    cap_idx = [i for i in range(13, len(xi))
               if abs(np.linalg.norm(xi[i] - xi[1]) - 1.0) < 1e-6
               and abs(xi[i] @ np.ones(3) - 2*np.sqrt(2)) < 1e-6]

    runs = [
        ('N    ', shell, 13, [1, 2, 3], [QU, QD, QD], [QU, QU, QD],
         'A_2u', 'anchor', 1),
        ('Sigma', sig, 17, [1, 2, 3], [QD, QD, QS], [QU, QU, QS],
         'A_2g', 'stiff', 2),
        ('Xi   ', xi, 19, [1] + cap_idx, [QD, QS, QS], [QU, QS, QS],
         'A_2g', 'soft', 1),
    ]

    print("=" * 74)
    print("1. PT validation (nucleon, u_iso): analytic d2 vs exact FD")
    print("=" * 74)
    lam0, a, S = S_form(shell, [1, 2, 3], 'A_2u', 'anchor', 'u_iso')
    qh, ql = np.array([QU, QD, QD]), np.array([QU, QU, QD])
    d2_pt = qh @ S @ qh - ql @ S @ ql
    for g in (1e-3, 1e-2, 1e-1):
        d2_fd = validate_pt(shell, [1, 2, 3], qh, ql, 'A_2u', 'anchor',
                            'u_iso', g=g)
        print(f"   g = {g:.0e}:  FD {d2_fd:+.6f}   PT {d2_pt:+.6f}")

    print("\n" + "=" * 74)
    print("2. The admixture per swap, three couplings, three clusters")
    print("=" * 74)
    results = {}
    for model in ('u_iso', 'phi_iso', 'chiral'):
        print(f"\n  model: {model}")
        print(f"  {'baryon':7}{'lam0':>8}{'d1/sw [g]':>12}"
              f"{'d2/sw [g^2]':>13}{'N*d2/sw':>11}")
        for name, co, N, sites, q_hi, q_lo, ch, br, nsw in runs:
            lam0, a, S = S_form(co, sites, ch, br, model)
            q_hi, q_lo = np.array(q_hi), np.array(q_lo)
            d1 = (a @ q_hi - a @ q_lo) / nsw
            d2 = (q_hi @ S @ q_hi - q_lo @ S @ q_lo) / nsw
            results[(model, name)] = N * d2
            print(f"  {name:7}{lam0:>8.3f}{d1:>12.5f}{d2:>13.6f}"
                  f"{N*d2:>11.5f}")
        rN = results[(model, 'N    ')]
        print(f"  ratios to nucleon:  Sigma "
              f"{results[(model,'Sigma')]/rN:+.3f}   "
              f"Xi {results[(model,'Xi   ')]/rN:+.3f}")

    print("\n" + "=" * 74)
    print("3. Scale scan: can ANY g^2 deliver the 4(1+n_s) requirement?")
    print("=" * 74)
    need = np.array([-0.136, -0.245, -0.276])     # MeV per swap
    err = np.array([0.005, 0.038, 0.210])
    print("  required cs per swap: N -0.136, Sigma -0.245, Xi -0.276 MeV")
    for model in ('u_iso', 'phi_iso', 'chiral'):
        pm = np.array([results[(model, n)]
                       for n in ('N    ', 'Sigma', 'Xi   ')]) * M_E
        x = np.sum(need * pm / err**2) / np.sum(pm**2 / err**2)
        chi2 = np.sum(((need - x * pm) / err) ** 2)
        print(f"  {model:8}: best x = {x:8.3f}, chi2 = {chi2:6.1f} (3 pts) "
              f"-> N {x*pm[0]:+.3f}, Sig {x*pm[1]:+.3f}, "
              f"Xi {x*pm[2]:+.3f} MeV")

    print("\n" + "=" * 74)
    print("4. Anchor check and verdict")
    print("=" * 74)
    pm = results[('chiral', 'N    ')] * M_E
    g = np.sqrt(0.08 / abs(pm))
    s = results[('chiral', 'Sigma')] * M_E * g * g
    xq = results[('chiral', 'Xi   ')] * M_E * g * g
    print(f"  chiral coupling anchored to nucleon -0.08 MeV: g = {g:.3f}")
    print(f"  hyperon admixtures: Sigma {s:+.4f}, Xi {xq:+.4f} MeV per swap")
    print("""
  Verdict: no coupling, at any scale, supplies the uniform -0.25 MeV the
  exact-4(1+n_s) reading needs. The nucleon mode (top of channel) is
  sign-protected; the mid-spectrum hyperon modes are not. The ladder
  step stands at 3.79 +- 0.07 m_e: not a quadrupole dressing, and its
  derivation belongs to the strange structures' linear (PN) sector.""")


if __name__ == "__main__":
    main()
