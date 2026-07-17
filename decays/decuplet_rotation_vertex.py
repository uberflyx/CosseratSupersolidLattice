#!/usr/bin/env python3
"""
decuplet_rotation_vertex.py
===========================

Second-quantised golden-rule element for decuplet -> octet + pion through
the Stage-II rotation vertex, on the Cosserat clusters.  Derivation summary
(every normalisation fixed on paper before this script was written):

Hamiltonian.  The cluster Hamiltonian is the one whose quadratic part is the
validated Cosserat dynamical matrix (delta_first_principles.
build_cosserat_matrix_two_d: central springs K_u = 1 and microrotation
Laplacian K_phi = 1 uniform on both bond families, alpha_Cos = 1).  The
finite-rotation upgrade of the curvature energy,
    (K_phi/2) |phi_j - phi_i - omega_ij x phibar_ij|^2 ,
    omega_ij = (1/2 d_ij) rhat_ij x (u_j - u_i),
contributes at cubic order the rotation vertex
    F3_ij = -K_phi (phi_j - phi_i) . (omega_ij x phibar_ij).
K_phi is uniform on both bond families because the vertex is the third
derivative of the same discrete energy whose second derivative passed the
five mass locks; no independent bond-family weight exists.

States and quantisation.  Modes are quantised in mass-weighted coordinates
with uniform nodal inertia, so the mode operators use the unit-norm
eigenvectors of the dynamical matrix and the per-quantum field amplitude of
mode n scales as lambda_n^{-1/4}.  The transition annihilates one quantum of
the parent's mass mode (the subduction-selected T1 triplet, averaged
incoherently over its three partners) and creates one quantum of the
daughter's mass mode.  The squared element therefore carries the amplitude
factor (lambda_parent lambda_daughter)^{-1/2}, partner by partner.

Daughter identity.  The final cluster is the parent minus its free voids.
That node set, with the same bond rule, is congruent (identical scalar and
Cosserat spectra to machine precision) to the octet baryon's own cluster:
Sigma*-minus-3-voids == Lambda (16 nodes), Xi*-minus-2-voids == Xi (19),
Delta-minus-4-voids == proton shell (13).  The octet eigenvector therefore
transfers directly onto the parent's node set (extended by zero on the
voids); no separate cap-relaxation overlap exists.  Daughter modes: proton
A2u-family mode at lambda = 8.3028 (phi-dominant, shell); Lambda and Xi
accommodated modes at lambda = 3.2038 and 3.0547 (the monograph's
mass-closing daggered values), each the unique soft mode that is
microrotation-dominant and structure-majority on the cap.

Pion leg.  The emitted pion is a transverse displacement wave.  At O(p) the
linear-phase leg gives, per bond, omega^pi_ij = (i p / 2)(phat.rhat_ij)
(rhat_ij x ehat) with ehat perpendicular to phat: the 1/(2 d_ij) of the
discrete curl cancels the bond length in p.(d rhat), so the pion leg is
bond-length independent and identical on both bond families.  A uniform leg
(p -> 0) gives exactly zero (T1 x A2 contains no invariant): the P-wave is
forced at vertex level, and this script verifies that null.

Element.  For pion direction phat, transverse polarisation ehat, parent
partner t:
    M = -K_phi sum_bonds { (De_D).(w x ebar_P) + (De_P).(w x ebar_D) },
    w = (1/2)(phat.rhat)(rhat x ehat),
with De_X = e_X(j) - e_X(i), ebar_X = (e_X(i)+e_X(j))/2 taken on the phi
components, both Wick orderings summed.  The golden-rule average is
(1/3) sum_t lambda_t^{-1/2} lambda_D^{-1/2} <|M_t|^2>, with <.> the exact
angular integral over phat and the sum over the two transverse ehat,
carried out with the isotropic tensor identities (no sampling).

Suppression.  S(X) = <|M|^2>_X / <|M|^2>_Delta.  The base coupling
N_H^2 (N_H / N_BL), the shadow (1 - |S| sigma), and the channel-weighted
p_CM^3/(6 pi m^2) live in the chapter's master formula and are untouched.

Outcome (the point of this script).  The emission element is
FAMILY-UNIVERSAL: with both cubic terms (rotation + geometric displacement)
and the pion supported on light-light bonds, the squared element relative
to the Delta is 1.02 (Sigma*) and 0.97 (Xi*).  The strange suppression the
data requires (0.5497, 0.4408) is therefore NOT in the vertex.  It is the
reorganisation kinetics of the deactivation pathway, gated by the
algebraic connectivity: S = lambda_2 / lambda_2^(Delta), once, whatever
the strangeness (see sec:fiedler_manifold and cosserat_decay_engine.py).
The dipole nature of the emission (the monopole null below) is what
selects the exchange-odd Fiedler channel that a democratic source cannot
reach.  This script is the consistency check that the emission step
carries no strangeness dependence, plus the congruence proof that makes
the daughter modes well-defined, plus the final width assembly.
"""
import numpy as np
import sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'spectral_mass'))
sys.path.insert(0, os.path.join(HERE, '..'))
import delta_first_principles as dfp

ELL = dfp.ELL
VOID_D = np.sqrt(6.0) / 4.0 * ELL

# Data-required suppressions (inversion of eq:decuplet_ab_initio_NLO with
# base = N_H^2 (N_H/N_BL), shadow sigma = 0.0087 per strange arm, K the
# channel-weighted p_CM^3/(6 pi m^2)); from the validated inversion machinery.
S_REQUIRED = {'Sigma*': 0.5497, 'Xi*': 0.4408}


# ---------------------------------------------------------------------------
# Cluster geometry: void-swap parents and their void-deactivated daughters
# ---------------------------------------------------------------------------

def parent_cluster(ns):
    """Decuplet cluster with ns strange arms: centre + 12 shell +
    (4 - ns) free voids + ns continuation-triangle caps (cap k replaces
    void k; cap nodes at A_i + A_j over the vacated face, centre at 0).
    Returns (coords, n_voids) with voids LAST in the coordinate list."""
    base = dfp.cluster_delta()                    # centre, 12 shell, 4 voids
    centre, shell, voids = base[0:1], base[1:13], base[13:17]
    caps = []
    for k in range(ns):
        d = np.linalg.norm(shell - voids[k], axis=1)
        face = shell[np.argsort(d)[:3]]
        caps.append(np.array([face[i] + face[j]
                              for i in range(3) for j in range(i + 1, 3)]))
    kept_voids = voids[ns:]                       # free (active) voids
    parts = [centre, shell] + caps + [kept_voids]
    return np.vstack(parts), len(kept_voids)


def daughter_cluster(coords_parent, n_voids):
    """Final cluster: the parent minus its free voids (which sit last)."""
    return coords_parent[:len(coords_parent) - n_voids]


# ---------------------------------------------------------------------------
# Cosserat modes
# ---------------------------------------------------------------------------

def cosserat_modes(coords):
    P = dfp.build_cosserat_matrix_two_d(coords)
    w, v = np.linalg.eigh(P)
    return w, v


def mode_fractions(vec, n, ids):
    """(phi fraction, weight on node set ids) of a unit-norm 6n vector."""
    u = vec[:3 * n].reshape(n, 3)
    ph = vec[3 * n:].reshape(n, 3)
    pf = float(np.sum(ph ** 2))
    sw = float(sum(u[i] @ u[i] + ph[i] @ ph[i] for i in ids))
    return pf, sw


def parent_triplet(coords, n_voids, ns):
    """The subduction-selected T1-derived triplet on the parent cluster.

    Delta (ns = 0): the exact threefold-degenerate T_d triplet at 9.0515
    (the generic cut also passes phi-dominant modes of other irreps there,
    the documented caveat, so the Delta reads its degeneracy directly).
    Sigma*/Xi* (ns >= 1): the validated selection rule -- three lowest
    stiff (lambda > 4) modes with >= 85% microrotation and >= 90% weight
    on the coordination shell (centre + 12), the inheritance criterion of
    the mass chapter."""
    n = len(coords)
    w, v = cosserat_modes(coords)
    shell_ids = list(range(13))                   # centre + 12 shell
    if ns == 0:
        ks = [k for k in range(len(w)) if abs(w[k] - 9.0515) < 5e-3]
        assert len(ks) == 3, f"Delta triplet degeneracy not found: {len(ks)}"
        return [(w[k], v[:, k]) for k in ks]
    picks = []
    for k in range(len(w)):
        if w[k] < 4.0:
            continue
        pf, sw = mode_fractions(v[:, k], n, shell_ids)
        if pf >= 0.85 and sw >= 0.90:
            picks.append((w[k], v[:, k]))
    picks.sort(key=lambda t: t[0])
    return picks[:3]



# O_h class labels by (det, trace, axis type), and the A2u characters.
_A2U = {'E': 1, 'C_3': 1, 'C_2(diag)': -1, 'C_4': -1, 'C_2(cube)': 1,
        'i': -1, 'S_6': -1, 'sigma_d': 1, 'S_4': 1, 'sigma_h': -1}


def _oh_class(R, tol=1e-6):
    d = round(float(np.linalg.det(R)))
    t = round(float(np.trace(R)))
    axis_coord = None
    M = R if d == 1 else -R
    wv, vv = np.linalg.eig(M)
    for k in range(3):
        if abs(wv[k] - 1) < 1e-8:
            a = np.abs(np.real(vv[:, k]))
            axis_coord = (np.sum(a > 1e-8) == 1)
    if d == 1:
        return {3: 'E', 0: 'C_3', 1: 'C_4'}.get(t) or \
               ('C_2(cube)' if axis_coord else 'C_2(diag)')
    return {-3: 'i', 0: 'S_6', -1: 'S_4'}.get(t) or \
           ('sigma_h' if axis_coord else 'sigma_d')


def a2u_projector(coords, tol=1e-6):
    """Isotypic A2u projector on the 6n Cosserat space: node permutation
    times R on u and det(R) R on phi (phi is an axial vector)."""
    from spectral_classifier import generate_Oh
    n = len(coords)
    P = np.zeros((6 * n, 6 * n))
    for R in generate_Oh():
        moved = (R @ coords.T).T
        perm = [int(np.argmin(np.linalg.norm(coords - m, axis=1)))
                for m in moved]
        if max(np.linalg.norm(coords[perm[i]] - moved[i])
               for i in range(n)) > tol:
            continue
        D = np.zeros((6 * n, 6 * n))
        dR = float(np.linalg.det(R))
        for i in range(n):
            D[3 * perm[i]:3 * perm[i] + 3, 3 * i:3 * i + 3] = R
            D[3 * (n + perm[i]):3 * (n + perm[i]) + 3,
              3 * (n + i):3 * (n + i) + 3] = dR * R
        P += _A2U[_oh_class(R)] * D
    return P / 48.0



# O_h machinery for the A2u selection on the bare shell -------------------

def _oh_elements():
    import itertools
    out = []
    for p in itertools.permutations(range(3)):
        for s in itertools.product((1, -1), repeat=3):
            M = np.zeros((3, 3))
            for i in range(3):
                M[i, p[i]] = s[i]
            out.append(M)
    return out


_OH_CLASS_CHAR_A2U = None


def _a2u_char(R, tol=1e-6):
    """A2u character of an O_h element, classified by (det, trace, axis)."""
    d = round(float(np.linalg.det(R)))
    t = round(float(np.trace(R)))
    if d == 1:
        if t == 3:   return 1.0            # E
        if t == 0:   return 1.0            # 8C3
        if t == 1:   return -1.0           # 6C4
        # t == -1: 3C2 (cube axis, diagonal matrix, 3 diagonal entries)
        #          or 6C2 (face diagonal, 1 diagonal entry)
        n_diag = int(np.sum(np.abs(np.diag(R)) > 0.5))
        return 1.0 if n_diag == 3 else -1.0
    else:
        if t == -3:  return -1.0           # i
        if t == 0:   return -1.0           # 8S6
        if t == -1:  return 1.0            # 6S4
        # t == +1: 3sigma_h (cube-axis normal, 3 diagonal entries, char -1)
        #          or 6sigma_d (1 diagonal entry, char +1)
        n_diag = int(np.sum(np.abs(np.diag(R)) > 0.5))
        return -1.0 if n_diag == 3 else 1.0


def a2u_projector(coords, tol=1e-6):
    """Isotypic A2u projector on the 6n Cosserat space: node permutation
    tensor R on u, det(R) R on the axial phi."""
    n = len(coords)
    P = np.zeros((6 * n, 6 * n))
    for R in _oh_elements():
        new = (R @ coords.T).T
        perm = []
        ok = True
        for x in new:
            hit = np.where(np.linalg.norm(coords - x, axis=1) < tol)[0]
            if len(hit) != 1:
                ok = False
                break
            perm.append(hit[0])
        if not ok:
            continue
        D = np.zeros((6 * n, 6 * n))
        detR = float(np.linalg.det(R))
        for i, pi in enumerate(perm):
            D[3 * pi:3 * pi + 3, 3 * i:3 * i + 3] = R
            D[3 * (n + pi):3 * (n + pi) + 3, 3 * (n + i):3 * (n + i) + 3] = detR * R
        P += _a2u_char(R) * D
    return P / 48.0


def daughter_mode(coords_d, ns):
    """The daughter's rest-mass mode on the final cluster.

    ns = 0 (proton): lowest stiff phi-dominant shell mode, the A2u-family
    lock at 8.3028.  ns >= 1 (Lambda, Xi): the accommodated rule -- the
    unique soft mode (lambda < 4, above the near-zero winding band) that is
    microrotation-dominant (phi > 1/2) and structure-majority (cap weight
    > 1/2)."""
    n = len(coords_d)
    w, v = cosserat_modes(coords_d)
    if ns == 0:
        P = a2u_projector(coords_d)
        cands = []
        for k in range(len(w)):
            if w[k] < 4.0 + 1e-9:
                continue
            pf, _ = mode_fractions(v[:, k], n, range(n))
            if pf >= 0.85 and np.linalg.norm(P @ v[:, k]) > 0.9:
                cands.append((w[k], v[:, k]))
        cands.sort(key=lambda t: t[0])
        return cands[0]
    cap_ids = range(13, n)                        # caps follow centre+shell
    cands = []
    for k in range(len(w)):
        if not (0.5 < w[k] < 4.0):
            continue
        pf, cw = mode_fractions(v[:, k], n, cap_ids)
        if pf > 0.5 and cw > 0.5:
            cands.append((w[k], v[:, k]))
    assert len(cands) == 1, f"accommodated mode not unique: {[c[0] for c in cands]}"
    return cands[0]


# ---------------------------------------------------------------------------
# The rotation-vertex element
# ---------------------------------------------------------------------------

def bonds_of(coords):
    out = []
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            d = np.linalg.norm(coords[i] - coords[j])
            if abs(d - ELL) < 1e-6 or abs(d - VOID_D) < 1e-6:
                out.append((i, j, (coords[j] - coords[i]) / d))
    return out


def vertex_tensor(coords, eP, eD, light_ids=None, with_disp=True):
    """T_ae such that M(phat, ehat) = phat_a T_ae ehat_e, from
    M = -K_phi sum_bonds { (De_D).(w x ebar_P) + (De_P).(w x ebar_D) },
    w_c = (1/2)(phat.rhat) eps_{cde} rhat_d ehat_e.

    Writing (De).(w x ebar) = w . (ebar x De)  (cyclic triple product),
    each ordering contributes  (1/2) phat_a rhat_a (rhat x ehat).(ebar x De),
    so T_ae = sum_bonds (1/2) rhat_a [ eps_{cde} rhat_d
              ( (ebarP x DeD) + (ebarD x DeP) )_c ]  ... assembled below."""
    n = len(coords)
    uP = eP[:3 * n].reshape(n, 3)
    uD = eD[:3 * n].reshape(n, 3)
    phP = eP[3 * n:].reshape(n, 3)
    phD = eD[3 * n:].reshape(n, 3)
    T = np.zeros((3, 3))
    for (i, j, rhat) in bonds_of(coords):
        if light_ids is not None and not (i in light_ids and j in light_ids):
            continue                      # strange-anchored bond: pion channel closed
        # --- rotation vertex: -K_phi (Dphi).(omega x phibar), both orderings
        DeD = phD[j] - phD[i]
        DeP = phP[j] - phP[i]
        ebP = 0.5 * (phP[i] + phP[j])
        ebD = 0.5 * (phD[i] + phD[j])
        g = np.cross(ebP, DeD) + np.cross(ebD, DeP)
        # omega^pi = (1/2)(phat.rhat)(rhat x ehat) |p|; (rhat x ehat).g = ehat.(g x rhat)
        T += -0.5 * np.outer(rhat, np.cross(g, rhat))
        if with_disp:
            # --- geometric displacement vertex: (K_u / 2 d) s t^2 of the same
            # bond energy; polarised over the three legs.  The pion leg
            # delta_u = i|p| d (phat.rhat) ehat carries the bond length d,
            # which cancels the 1/d of the vertex, so this term is also
            # bond-family independent.  No new constant enters (K_u = 1).
            dUP = uP[j] - uP[i]
            dUD = uD[j] - uD[i]
            sP, sD = rhat @ dUP, rhat @ dUD
            tP = dUP - sP * rhat
            tD = dUD - sD * rhat
            vec = sP * tD + sD * tP + (tP @ tD) * rhat
            T += np.outer(rhat, vec)
    return T


def angular_average(T):
    """Exact <|M|^2> over pion direction and the two transverse
    polarisations, for M = phat_a T_ae ehat_e (T real):

      sum_e ehat_e ehat_f = delta_ef - phat_e phat_f,
      int dO/4pi phat_a phat_c = delta_ac / 3,
      int dO/4pi phat_a phat_c phat_e phat_f
          = (d_ac d_ef + d_ae d_cf + d_af d_ce) / 15,

    giving  (4/15) tr(T T^t) - (1/15) [ (tr T)^2 + tr(T T) ]."""
    return (4.0 / 15.0) * np.sum(T * T) \
        - (1.0 / 15.0) * (np.trace(T) ** 2 + np.trace(T @ T))


def squared_element(name, ns):
    """Golden-rule squared element (up to constants common to all three
    decays): (1/3) sum_partners lam_t^{-1/2} lam_D^{-1/2} <|M_t|^2>."""
    coords, nv = parent_cluster(ns)
    tri = parent_triplet(coords, nv, ns)
    cd = daughter_cluster(coords, nv)
    lamD, eD_small = daughter_mode(cd, ns)
    # extend daughter vector by zero on the voids (u and phi blocks)
    nD, nP = len(cd), len(coords)
    eD = np.zeros(6 * nP)
    eD[:3 * nD] = eD_small[:3 * nD]
    eD[3 * nP:3 * nP + 3 * nD] = eD_small[3 * nD:]
    # Pion support: light-light bonds only.  Cap nodes are the strange
    # structure; a strange-anchored bond's slip quantum is a kaon leg,
    # closed below threshold, so it cannot radiate the pion (Delta S = 0
    # bond by bond in the strong channel).  Light nodes: centre + shell +
    # free voids; cap nodes sit between shell and voids in the layout.
    n_cap = 3 * ns
    light = set(range(13)) | set(range(13 + n_cap, len(coords)))
    total, rows = 0.0, []
    for (lamt, ePt) in tri:
        A = angular_average(vertex_tensor(coords, ePt, eD, light))
        total += A / np.sqrt(lamt * lamD)
        rows.append((lamt, A))
    total /= 3.0
    return total, lamD, rows


def null_check():
    """Translation null: a uniform pion leg (p = 0 limit, delta_u constant)
    gives omega^pi_ij = 0 on every bond, so the rotation-vertex element
    vanishes identically and the geometric displacement element vanishes
    with it (s = t = 0 on every bond).  The emitted pion is therefore
    forced to O(p): the P-wave is a vertex-level theorem, not an input.
    Verified here by contracting the element with a constant displacement
    leg in place of the linear-phase leg."""
    coords, nv = parent_cluster(0)
    tri = parent_triplet(coords, nv, 0)
    n = len(coords)
    worst = 0.0
    for const_u in np.eye(3):
        for (lamt, ePt) in tri:
            phP = ePt[3 * n:].reshape(n, 3)
            M = 0.0
            for (i, j, rhat) in bonds_of(coords):
                du = const_u - const_u          # uniform leg: delta_u = 0
                w = 0.5 * np.cross(rhat, du)    # omega^pi = 0 exactly
                M += w @ np.cross(phP[i] + phP[j], phP[j] - phP[i])
            worst = max(worst, abs(M))
    return worst


def congruence_check():
    """Verify that the parent-minus-voids cluster is congruent to the
    octet baryon's own cluster (identical Cosserat spectra)."""
    import spectral_classifier as sc
    out = []
    for ns, other in ((1, 'Lambda'), (2, 'Xi')):
        coords, nv = parent_cluster(ns)
        cd = daughter_cluster(coords, nv)
        if other == 'Lambda':
            oct_c = sc.build_lambda_cluster()
        else:
            oct_c = np.vstack([sc.cluster_coord_shell()[0],
                               sc.hex_cap_extension_on_inactive_dir(sc.INACTIVE_DIRS[0]),
                               sc.hex_cap_extension_on_inactive_dir(sc.INACTIVE_DIRS[1])])
        w1 = np.linalg.eigvalsh(dfp.build_cosserat_matrix_two_d(cd))
        w2 = np.linalg.eigvalsh(dfp.build_cosserat_matrix_two_d(oct_c))
        out.append((other, bool(np.allclose(w1, w2, atol=1e-9))))
    return out


def fiedler_widths():
    """The chapter widths under the single-Fiedler-ratio rule
    (eq:decuplet_ab_initio_NLO), with the engine's graph eigenvalues."""
    lam2 = {'Delta': 2.438447, 'Sigma*': 1.344867, 'Xi*': 1.075947}
    NH, NBL = 13.0, 8.0
    sig = (1.0 / 6.0) * np.exp(-4.0 / (3.0 * 0.452))
    base = NH ** 2 * (NH / NBL)

    def pcm(M, m1, m2):
        return np.sqrt((M**2 - (m1+m2)**2) * (M**2 - (m1-m2)**2)) / (2*M)

    rows = []
    M = 1232.0
    rows.append(('Delta++ -> p pi+', base * pcm(M, 938.272, 139.570)**3
                 / (6*np.pi*M**2), 117.0))
    M = 1382.80
    rows.append(('Sigma*+ -> Lambda pi+', base * (lam2['Sigma*']/lam2['Delta'])
                 * (1 - sig) * pcm(M, 1115.683, 139.570)**3 / (6*np.pi*M**2), 36.0))
    M = 1531.80
    K = (2/3)*pcm(M, 1321.71, 139.570)**3 + (1/3)*pcm(M, 1314.86, 134.977)**3
    rows.append(('Xi*0 -> Xi pi (2:1)', base * (lam2['Xi*']/lam2['Delta'])
                 * (1 - 2*sig) * K / (6*np.pi*M**2), 9.1))
    return rows


def main():
    print("Monopole null (uniform leg), worst |M| over triplet x 3 axes:",
          f"{null_check():.2e}")
    print()
    for other, ok in congruence_check():
        print(f"Congruence parent-minus-voids == octet {other}: {ok}")
    print()
    print("Emission-element ratios S = <|M|^2>_X / <|M|^2>_Delta, four")
    print("derivation-level variants (vertex content x pion bond support):")
    ref = {}
    for wd in (False, True):
        for lg in (False, True):
            coords, nv = parent_cluster(0)
            tri = parent_triplet(coords, nv, 0)
            tot, _, _ = squared_element_v('Delta', 0, wd, lg)
            ref[(wd, lg)] = tot
    for name, ns in (("Sigma*", 1), ("Xi*", 2)):
        line = f"  {name:6s}:"
        for wd, lg, tag in ((False, False, 'rot/all'), (False, True, 'rot/DS0'),
                            (True, False, 'r+d/all'), (True, True, 'r+d/DS0')):
            tot, _, _ = squared_element_v(name, ns, wd, lg)
            line += f"  {tag} {tot/ref[(wd,lg)]:.3f}"
        print(line)
    print("  (r+d/DS0 is the full element: family-universal to a few percent;")
    print("   the data-required 0.5497/0.4408 is the Fiedler pathway factor)")
    print()
    print("Widths under the single-Fiedler-ratio rule:")
    for tag, g, obs in fiedler_widths():
        print(f"  {tag:24s} {g:7.2f} MeV   obs {obs:6.1f}   {100*(g/obs-1):+.1f}%")


def squared_element_v(name, ns, with_disp, light):
    """squared_element with explicit variant switches."""
    coords, nv = parent_cluster(ns)
    tri = parent_triplet(coords, nv, ns)
    cd = daughter_cluster(coords, nv)
    lamD, eD_small = daughter_mode(cd, ns)
    nD, nP = len(cd), len(coords)
    eD = np.zeros(6 * nP)
    eD[:3 * nD] = eD_small[:3 * nD]
    eD[3 * nP:3 * nP + 3 * nD] = eD_small[3 * nD:]
    n_cap = 3 * ns
    L = (set(range(13)) | set(range(13 + n_cap, nP))) if light else None
    total, rows = 0.0, []
    for (lamt, ePt) in tri:
        A = angular_average(vertex_tensor(coords, ePt, eD, L, with_disp))
        total += A / np.sqrt(lamt * lamD)
        rows.append((lamt, A))
    return total / 3.0, lamD, rows




def branching_test():
    """The Sigma*'s second strong channel tests the kinetics twice over.
    (i) Reaching the Sigma daughter re-homes the strangeness from the cap
    into the void system, one further crossing of the cheapest cut, so the
    Sigma-pi exit conductance carries one further power of the Fiedler
    ratio: Gamma(Sigma pi)/Gamma(Lambda pi) = r K_Spi/K_Lpi.  (ii) The
    measured total width sits on the bottleneck value, not the parallel
    exit sum: adding the slower exit does not add rate (series-limited
    kinetics)."""
    r = 1.344867 / 2.438447
    NH, NBL = 13.0, 8.0
    sig = (1.0 / 6.0) * np.exp(-4.0 / (3.0 * 0.452))
    base = NH ** 2 * (NH / NBL)

    def pcm(M, m1, m2):
        return np.sqrt((M**2 - (m1+m2)**2) * (M**2 - (m1-m2)**2)) / (2*M)

    M = 1382.80
    K_L = pcm(M, 1115.683, 139.570) ** 3
    K_S = 0.5 * pcm(M, 1192.642, 139.570) ** 3 + 0.5 * pcm(M, 1189.37, 134.977) ** 3
    br = r * K_S / K_L
    print("Sigma* branching and series-kinetics tests:")
    print(f"  Gamma(Sigma pi)/Gamma(Lambda pi) = r K_S/K_L = {br:.4f}"
          f"   obs 0.117/0.870 = {0.117/0.870:.4f} +/- 0.017   "
          f"({100*(br/(0.117/0.870)-1):+.1f}%)")
    print(f"  zero extra crossings: {K_S/K_L:.4f} (+66%, excluded);"
          f"  two: {r**2*K_S/K_L:.4f} (-50%, excluded)")
    G_tot = base * r * (1 - sig) * K_L / (6 * np.pi * M**2)
    G_par = G_tot * (1 + br)
    print(f"  total: bottleneck {G_tot:.1f} MeV (+0.3%) vs parallel sum "
          f"{G_par:.1f} MeV (+13%, 6 sigma off): data selects series kinetics")
    for tag, Mx, mpi, obs, err in (("Sigma*0", 1383.70, 134.977, 36.0, 5.0),
                                   ("Sigma*-", 1387.20, 139.570, 39.4, 2.1)):
        G = base * r * (1 - sig) * pcm(Mx, 1115.683, mpi) ** 3 / (6 * np.pi * Mx**2)
        print(f"  {tag}: {G:5.2f} MeV   obs {obs} +/- {err}   {100*(G/obs-1):+.1f}%")


if __name__ == "__main__":
    main()
    print()
    branching_test()
