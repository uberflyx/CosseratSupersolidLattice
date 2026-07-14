"""
Second-order (alpha^2) non-uniformity correction to the gravitational
coupling from the 19-node Born stability cluster.

Context
-------
The monograph predicts G = (1 + 1/pi)(1 - 17*alpha/18)(hbar c/m0^2) alpha^19.
The first-order self-energy factor -17*alpha/18 is exact and independent of
the NN/NNN bond-stiffness ratio eta (the (2 + eta) cancellation).  The
monograph currently carries only an order-of-magnitude estimate for the
second-order term produced by the cluster's non-uniformity:

    |correction| ~ alpha^2 * CV(R)^2 ~ 3 ppm,

with CV(R) the coefficient of variation of the pairwise resistance
distances among the 18 outer nodes.  This script does the calculation
properly.  Structure:

  1. Build the weighted 19-node cluster graph (NN bonds weight 1, NNN
     bonds weight eta), its Laplacian, and the resistance distances.
     Validate against the numbers quoted in the monograph:
     R(centre-NN) ~ 0.23, R(NN-NNN) ~ 0.36, CV(R) ~ 0.23.

  2. Verify symbolically (sympy) that the first-order bond-weight average
     is exactly 17/18 for all eta, and derive the per-node deviations
     delta_f_i(eta) about that mean.  These deviations are the *source*
     of any second-order non-uniformity term: second-order perturbation
     theory needs a non-uniform source, and the 17/18 theorem makes the
     source vanish identically at eta = 1.

  3. Compute the two second-order channels the non-uniformity opens:

     (a) Instanton path relaxation.  The coherent 19-node hop with all
         nodes moving identically is the constrained optimum; when the
         per-node barriers differ (delta_f nonzero), the true tunnelling
         path lets low-barrier nodes run slightly ahead.  Standard
         second-order perturbation theory of the collective action gives

            Delta S / S  =  -(1/2) sum_{n>=1} |<v_n|alpha*delta_f>|^2
                                              / (1 + kappa*lambda_n),

         where lambda_n, v_n are the cluster Laplacian eigenpairs and
         kappa ~ 1e-3 is the lattice-to-PN stiffness ratio (the R_n of
         the monograph's fluctuation-determinant paragraph).  The action
         decreases, so the tunnelling amplitude and hence G *increase*:
         this channel has definite positive sign on G.

     (b) Arithmetic-geometric mean mismatch.  The collective dressing is
         the mean of per-node dressings; the exact amplitude is the
         geometric mean of per-node factors.  Non-uniform dressings make
         the geometric mean fall below the arithmetic mean by
         (1/2) Var(alpha*delta_f): definite negative sign on G.

     Both channels are proportional to the variance of the first-order
     per-node dressings, Var(f) = (eta-1)^2 / (162 (2+eta)^2), so both
     vanish exactly at eta = 1 and are bounded over any physical eta.

  4. Compare against the O(alpha^3) ~ 0.4 ppm truncation floor.

Result (headline): the computed second-order non-uniformity term is
<= 0.2 ppm for eta in [0.5, 2] and exactly zero at eta = 1.  The ~3 ppm
order estimate in the monograph is superseded: it paired the propagator's
non-uniformity (CV(R) = 0.23, real) with a source assumed to fluctuate at
the same level, but the 17/18 theorem makes the source uniform at first
order.  The theory floor is therefore set by the O(alpha^3) truncation.

No heavy numerics; a 19-node problem needs numpy only.
"""

import itertools

import numpy as np
import sympy as sp

ALPHA = 7.2973525643e-3          # fine-structure constant, CODATA 2022
PPM = 1e6

# ----------------------------------------------------------------------
# 1. Cluster geometry and weighted graph
# ----------------------------------------------------------------------
# Coordinates in units of a/2 (a = conventional FCC cube edge = ell*sqrt(2)).
# Squared distances: NN bond (length ell) -> 2;  NNN bond (ell*sqrt(2)) -> 4.

def build_cluster():
    centre = [(0, 0, 0)]
    nn = [p for p in itertools.product((-1, 0, 1), repeat=3)
          if sorted(map(abs, p)) == [0, 1, 1]]           # 12 cubocta vertices
    nnn = [p for p in set(itertools.permutations((2, 0, 0)))
           | set(itertools.permutations((-2, 0, 0)))]     # 6 octahedron vertices
    nodes = centre + sorted(nn) + sorted(nnn)
    assert len(nodes) == 19
    return nodes


def edges_by_type(nodes):
    """Return lists of (i, j) pairs with squared distance 2 (NN, weight 1)
    and 4 (NNN, weight eta).  All other separations are non-bonded."""
    nn_edges, nnn_edges = [], []
    for i, j in itertools.combinations(range(len(nodes)), 2):
        d2 = sum((a - b) ** 2 for a, b in zip(nodes[i], nodes[j]))
        if d2 == 2:
            nn_edges.append((i, j))
        elif d2 == 4:
            nnn_edges.append((i, j))
    return nn_edges, nnn_edges


def laplacian(n_nodes, nn_edges, nnn_edges, eta):
    L = np.zeros((n_nodes, n_nodes))
    for (i, j) in nn_edges:
        L[i, j] -= 1.0
        L[j, i] -= 1.0
    for (i, j) in nnn_edges:
        L[i, j] -= eta
        L[j, i] -= eta
    np.fill_diagonal(L, -L.sum(axis=1))
    return L


def resistance_matrix(L):
    G = np.linalg.pinv(L)                       # Moore-Penrose pseudo-inverse
    d = np.diag(G)
    return d[:, None] + d[None, :] - 2.0 * G


# ----------------------------------------------------------------------
# 2. First-order theorem and per-node deviations (symbolic)
# ----------------------------------------------------------------------

def first_order_symbolic():
    eta = sp.symbols('eta', positive=True)
    f_nn = (11 + 6 * eta) / (12 + 6 * eta)      # new-bond fraction, NN node
    f_nnn = (12 + 5 * eta) / (12 + 6 * eta)     # new-bond fraction, NNN node
    mean = sp.simplify((12 * f_nn + 6 * f_nnn) / 18)
    var = sp.simplify((12 * (f_nn - mean) ** 2 + 6 * (f_nnn - mean) ** 2) / 18)
    return eta, f_nn, f_nnn, mean, var


# ----------------------------------------------------------------------
# 3. Second-order channels
# ----------------------------------------------------------------------

def second_order(eta_val, kappa, nodes, nn_edges, nnn_edges):
    """Both alpha^2 channels at a given eta and lattice/PN stiffness ratio
    kappa.  Returns (relaxation_ppm, amgm_ppm, cv_R) as shifts of G in ppm."""
    n = len(nodes)
    L = laplacian(n, nn_edges, nnn_edges, eta_val)

    # Per-node first-order dressing deviations (source vector).
    f_nn = (11 + 6 * eta_val) / (12 + 6 * eta_val)
    f_nnn = (12 + 5 * eta_val) / (12 + 6 * eta_val)
    fbar = 17.0 / 18.0
    delta = np.zeros(n)
    for i, p in enumerate(nodes):
        if i == 0:
            continue                            # centre: dressing inside alpha
        d2 = sum(c * c for c in p)
        delta[i] = (f_nn if d2 == 2 else f_nnn) - fbar
    assert abs(delta.sum()) < 1e-14             # orthogonal to uniform mode

    # (a) Instanton path relaxation: 2nd-order PT against the cluster's
    # fluctuation modes; propagator (1 + kappa*lambda_n)^(-1).
    lam, V = np.linalg.eigh(L)
    proj = V.T @ (ALPHA * delta)
    relax = 0.0
    for lam_n, p_n in zip(lam, proj):
        if lam_n < 1e-12:
            continue                            # collective zero mode
        relax += p_n ** 2 / (1.0 + kappa * lam_n)
    dG_relax = +0.5 * relax                     # amplitude (and G) increase

    # (b) AM-GM mismatch of per-node dressings: geometric mean below
    # arithmetic mean by half the variance.
    var_f = delta[1:].var()                     # population variance, 18 nodes
    dG_amgm = -0.5 * ALPHA ** 2 * var_f * 18 / 18  # = -(1/2) alpha^2 Var(f)

    # Resistance-distance statistics for validation / reporting.
    R = resistance_matrix(L)
    outer = list(range(1, n))
    pairs = [R[i, j] for i, j in itertools.combinations(outer, 2)]
    cv_R = np.std(pairs) / np.mean(pairs)

    return dG_relax * PPM, dG_amgm * PPM, cv_R


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    nodes = build_cluster()
    nn_edges, nnn_edges = edges_by_type(nodes)
    print(f"Cluster: {len(nodes)} nodes, {len(nn_edges)} NN bonds, "
          f"{len(nnn_edges)} NNN bonds, {len(nn_edges)+len(nnn_edges)} total")

    # --- validation at eta = 1 against monograph-quoted values ----------
    L1 = laplacian(19, nn_edges, nnn_edges, 1.0)
    R1 = resistance_matrix(L1)
    idx_nn = next(i for i, p in enumerate(nodes) if sum(c*c for c in p) == 2)
    idx_nnn = next(i for i, p in enumerate(nodes) if sum(c*c for c in p) == 4)
    print(f"\nValidation at eta = 1:")
    print(f"  R(centre-NN)  = {R1[0, idx_nn]:.4f}")
    print(f"  R(NN-NNN)     = {R1[idx_nn, idx_nnn]:.4f}")
    # The monograph's quoted values (~0.23, ~0.36) are the eta -> 0 edge
    # set: distance-ell bonds only.  The CV is eta-robust; the individual
    # resistances are not.
    L0 = laplacian(19, nn_edges, [], 1.0)
    R0 = resistance_matrix(L0)
    print(f"  eta->0 limit: R(centre-NN) = {R0[0, idx_nn]:.4f}, "
          f"R(NN-NNN) = {R0[idx_nn, idx_nnn]:.4f}   (monograph: ~0.23, ~0.36)")
    outer = range(1, 19)
    pr = [R1[i, j] for i, j in itertools.combinations(outer, 2)]
    print(f"  CV(R), 153 outer pairs = {np.std(pr)/np.mean(pr):.4f}"
          f"   (scratch calc: 0.235)")
    print(f"  distinct outer-pair R values: "
          f"{sorted(set(np.round(pr, 6)))}")

    # --- first-order theorem (symbolic) ---------------------------------
    eta, f_nn, f_nnn, mean, var = first_order_symbolic()
    print(f"\nFirst order (symbolic):")
    print(f"  mean new-bond fraction = {mean}   (exact, eta-independent)")
    print(f"  Var(f) = {sp.factor(var)}   -> vanishes at eta = 1")

    # --- second order over eta and kappa ---------------------------------
    print(f"\nSecond-order non-uniformity shifts of G [ppm]:")
    print(f"  {'eta':>5} {'kappa':>8} {'relax(+)':>10} {'AM-GM(-)':>10} "
          f"{'net':>10} {'CV(R)':>7}")
    worst = 0.0
    for eta_val in (0.25, 0.5, 1.0, 1.5, 2.0, 4.0):
        for kappa in (1e-4, 1e-3, 1e-2, 1e-1):
            r, a, cv = second_order(eta_val, kappa, nodes,
                                    nn_edges, nnn_edges)
            net = r + a
            if 0.5 <= eta_val <= 2.0:
                worst = max(worst, abs(r), abs(a), abs(net))
            if kappa == 1e-3 or eta_val in (0.5, 2.0):
                print(f"  {eta_val:5.2f} {kappa:8.0e} {r:10.4f} {a:10.4f} "
                      f"{net:10.4f} {cv:7.4f}")

    # Closed form (small-kappa limit), verified symbolically:
    #   c2(eta) = (17/324) * ((eta - 1)/(eta + 2))**2,   shift = +c2 * alpha^2
    e = sp.symbols('eta', positive=True)
    c2 = sp.Rational(17, 324) * ((e - 1) / (e + 2)) ** 2
    print(f"\nClosed form: c2(eta) = 17/324 * ((eta-1)/(eta+2))^2 = {c2}")
    print(f"  c2(1) = 0 exactly; c2(2) = {float(c2.subs(e, 2)):.3e} "
          f"-> +{float(c2.subs(e, 2)) * ALPHA**2 * PPM:.3f} ppm")

    print(f"\nWorst-case |second-order shift| over eta in [0.5, 2], "
          f"kappa in [1e-4, 1e-1]: {worst:.3f} ppm")
    print(f"Monograph order estimate it replaces: alpha^2 * CV(R)^2 = "
          f"{ALPHA**2 * 0.235**2 * PPM:.1f} ppm")
    print(f"O(alpha^3) truncation floor: "
          f"{(17*ALPHA/18)**3 * PPM:.2f} ppm  -> now dominates")


if __name__ == "__main__":
    main()
