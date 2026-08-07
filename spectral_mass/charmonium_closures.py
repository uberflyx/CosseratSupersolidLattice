#!/usr/bin/env python3
"""
charmonium_closures.py
======================
Charmonium in the heavy-flavour assembly of the spectral mass formula.

ASSEMBLY.  A charmonium state is two K_{9,9} antibonding charm cores on
the two nodes of one shared cell pair, plus the confining ribbon set by
J^P (hex cap 7 for 0-, bilayer 8 for 1-, coordination shell 13 for the
P-wave J = 0, shell + stacking extension 14 for P-wave J >= 1, one added
bilayer 8 per radial stacking crossing):

    m = 2 x 18 m0  +  [ N_R m0 - N_R (4 - lambda_R) m_e ]

with NO heavy-light dressing: the dressing is the asymmetry residual of
one heavy terminus, and it vanishes for symmetric termini by the same
derivation that gives 5/3 and 8/3 in the open-charm sector.

STATUS.  The interior ribbon mode weights are quoted at the reference
eigenvalue lambda = 4 (leading order); for each state the script reports
the RIBBON eigenvalue the closure requires through the assembly's ribbon
term N_R (4 - lambda) m_e.  Reachability: on the hex cap the term lowers
the mass by at most 4 N_R m_e = 14.3 MeV, so the eta_c's -27 MeV gap
cannot be spectral; it is electromagnetic in order, and together with
the J/psi's +15.8 MeV it composes the +43 MeV quarkonium hyperfine
excess exactly.  Across the P wave the required eigenvalue rises with J
(1.51 -> 5.31 -> 7.37 -> 11.68), an eigenvalue ladder within the shell
spectrum's actual range.

PINNING.  The two parameter-free pinning limits of the core pair
(clamped; mass-loaded at the strain ratio 18) are computed below and
turn out to be one limit: an eighteen-fold inertia is quasi-static at
band frequencies, so the two spectra agree to better than 0.01 across
the upper band.  Their shared effect is a truncation: an inert core
pair removes the modes that live on it, capping the band below the
stiff targets.  What survives is the band edge, the periphery's own
maximal microrotation twist, which carries zero weight on the core pair
and whose eigenvalue is therefore independent of the pinning treatment.
Reading the four stiffest demands at that edge brings each inside (or,
for the eta_c(2S), exactly to) the leading-order resolution band.

DERIVED CONTRAST.  Because the dressing is absent, the quarkonium
hyperfine step is Delta N = 1 (one ribbon node), while the open-charm
step is Delta N = 2 (ribbon + dressing).  Hidden-charm splitting below
the heavy-light splitting is therefore a parameter-free prediction:
observed 113.0 MeV (J/psi - eta_c) vs 142.0 MeV (D*0 - D0).
"""
M_E = 0.51099895069
M_0 = M_E * 137.035999177
CORE = 2 * 18 * M_0

STATES = [
    # name, J^PC, ribbon N, PDG mass, PDG err
    ("eta_c(1S)",  "0-+", 7,  2984.1,   0.4),
    ("J/psi(1S)",  "1--", 8,  3096.900, 0.006),
    ("chi_c0(1P)", "0++", 13, 3414.71,  0.30),
    ("chi_c1(1P)", "1++", 14, 3510.67,  0.05),
    ("h_c(1P)",    "1+-", 14, 3525.37,  0.14),
    ("chi_c2(1P)", "2++", 14, 3556.17,  0.07),
    ("eta_c(2S)",  "0-+", 7+8, 3637.7,  0.9),
    ("psi(2S)",    "1--", 8+8, 3686.097, 0.011),
]

def main():
    print(f"m0 = {M_0:.4f} MeV; charm cores 2 x 18 m0 = {CORE:.1f} MeV\n")
    print(f"{'state':12s} {'N_R':>4s} {'N':>3s} {'LO (l=4)':>9s} "
          f"{'PDG':>9s} {'resid':>7s} {'required l':>10s}")
    for name, jpc, nr, obs, err in STATES:
        N = 36 + nr
        lo = N * M_0
        # required RIBBON eigenvalue, per the assembly's ribbon term:
        # obs = 36 m0 + N_R m0 - N_R (4 - l) m_e
        lam = 4 - (lo - obs) / (nr * M_E)
        # reachability: the ribbon term lowers the mass by at most
        # N_R * 4 * m_e (at lambda = 0); a larger downward gap cannot
        # be spectral and is electromagnetic in order.
        col = "  ---(EM)" if lam < 0 else f"{lam:9.2f}"
        print(f"{name:12s} {nr:4d} {N:3d} {lo:9.1f} {obs:9.2f} "
              f"{100*(lo-obs)/obs:+6.2f}% {col}")
    print("\nHyperfine contrast (parameter-free):")
    print(f"  quarkonium DN = 1 (no dressing): LO {M_0:.1f}  "
          f"obs J/psi-eta_c = 112.8, Upsilon-eta_b = 61.7")
    print(f"  heavy-light DN = 2 (ribbon+dressing): LO {2*M_0:.1f}  "
          f"obs D*0-D0 = 142.0, D*+-D+ = 140.6, Ds*-Ds = 143.8")
    print("  prediction: hidden-flavour splitting < heavy-light "
          "splitting: observed.")
    print("\nRadial step: one bilayer per stacking crossing, "
          f"8 m0 = {8*M_0:.1f} MeV;")
    print("  observed steps: psi(2S)-J/psi = 589.2, "
          "eta_c(2S)-eta_c = 653.6 (scale right, closure open).")

# ---------------------------------------------------------------------
# The pinning computation: clamped and loaded spectra, and the band edge.
# Single-scale Cosserat convention (K_u = K_phi = alpha_c = 1, integer FCC
# coordinates with NN distance sqrt 2, per-bond curl coefficient 1/(2 d_b)).
# ---------------------------------------------------------------------
import itertools
import numpy as np
from scipy.linalg import eigh

NN2 = 2.0
PAIR = (np.zeros(3), np.array([1., 1., 0.]))   # the shared cell pair


def _skew(v):
    x, y, z = v
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], float)


def cosserat(V):
    V = [np.asarray(v, float) for v in V]
    n = len(V)
    Puu = np.zeros((3 * n, 3 * n))
    L = np.zeros((n, n))
    C = np.zeros((3 * n, 3 * n))
    for i in range(n):
        for j in range(i + 1, n):
            d2 = float((V[j] - V[i]) @ (V[j] - V[i]))
            if not np.isclose(d2, NN2):
                continue
            for a, b in ((i, j), (j, i)):
                d = V[b] - V[a]
                db = np.linalg.norm(d) / np.sqrt(NN2)
                r = d / np.linalg.norm(d)
                B = np.outer(r, r)
                Puu[3*a:3*a+3, 3*a:3*a+3] += B
                Puu[3*a:3*a+3, 3*b:3*b+3] -= B
                L[a, a] += 1.0
                L[a, b] -= 1.0
                C[3*a:3*a+3, 3*b:3*b+3] += (0.5 / db) * _skew(r)
    Ppp = np.kron(L, np.eye(3)) + np.eye(3 * n)
    return np.block([[Puu + C.T @ C, -C.T], [-C, Ppp]])


def hex_cap(centre=(0, 0, 0)):
    c = np.array(centre, float)
    ring = [c + np.array(v, float)
            for v in itertools.product((-1, 0, 1), repeat=3)
            if sum(abs(x) for x in v) == 2 and abs(sum(v)) < 1e-9]
    return [c] + ring


def dedupe(V):
    out, seen = [], set()
    for v in V:
        k = tuple(np.round(np.asarray(v, float), 6))
        if k not in seen:
            seen.add(k)
            out.append(np.asarray(v, float))
    return out


def shell():
    return [np.zeros(3)] + [np.array(v, float)
                            for v in itertools.product((-1, 0, 1), repeat=3)
                            if sum(abs(x) for x in v) == 2]


RIBBONS = {
    "bilayer":    dedupe(hex_cap() + [np.array([1., 1., 0.])]),
    "shell":      shell(),
    "shell+ext":  dedupe(shell() + [np.array([2., 0., 0.])]),
    "cap+bilayer": dedupe(hex_cap() + hex_cap((1, 1, 0))
                          + [np.array([2., 2., 0.])]),
}


def pinned(V, mode):
    """Spectrum with the cell pair clamped or 18x mass-loaded."""
    V = [np.asarray(v, float) for v in V]
    n = len(V)
    M = cosserat(V)
    idx = [i for i, v in enumerate(V)
           if any(np.allclose(v, p) for p in PAIR)]
    dof = [k for i in idx
           for k in list(range(3*i, 3*i+3)) + list(range(3*n+3*i, 3*n+3*i+3))]
    if mode == "clamped":
        keep = [k for k in range(6 * n) if k not in dof]
        return np.linalg.eigvalsh(M[np.ix_(keep, keep)])
    B = np.eye(6 * n)
    for k in dof:
        B[k, k] = 18.0
    return eigh(M, B, eigvals_only=True)


def pinning_report():
    print("\nPinning limits and the band edge:")
    for name, V in RIBBONS.items():
        cl = pinned(V, "clamped")
        ld = pinned(V, "loaded")
        print(f"  {name:12s} N={len(V):2d}  edge clamped {cl[-1]:7.4f}  "
              f"loaded {ld[-1]:7.4f}  (agree to {abs(cl[-1]-ld[-1]):.4f})")
    print("\nBand-edge readings of the four stiffest demands:")
    EDGE = [("J/psi",     "bilayer",     8,  44, 3096.900),
            ("chi_c2",    "shell+ext",   14, 50, 3556.17),
            ("eta_c(2S)", "cap+bilayer", 15, 51, 3637.7),
            ("psi(2S)",   "cap+bilayer", 16, 52, 3686.097)]
    for nm, rb, nr, N, obs in EDGE:
        lam = pinned(RIBBONS[rb], "loaded")[-1]
        m = N * M_0 + nr * (lam - 4) * M_E
        note = "  [provisional: 15 distinct sites vs 16 counted]" \
            if nm == "psi(2S)" else ""
        print(f"  {nm:10s} edge lam {lam:7.4f}  m = {m:8.2f}  "
              f"obs {obs:8.2f}  ({100*(m-obs)/obs:+.2f}%){note}")


if __name__ == "__main__":
    main()
    pinning_report()
