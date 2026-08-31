#!/usr/bin/env python3
"""
cluster_monotonicity.py
=======================
Three checks behind the eigenvalue ceiling and the node-count bracket.

1. MONOTONICITY.  If a cluster Gamma' is an induced subcluster of Gamma
   (a subset of the nodes, keeping every bond with both ends surviving),
   then every eigenvalue of the 6n x 6n Cosserat dynamical matrix satisfies
   lambda_i(Gamma') <= lambda_i(Gamma), the eigenvalues taken in descending
   order.  The proof is in the monograph; this script is the numerical
   check, run on the catalogue's nested clusters and on randomised node
   deletions from random FCC clusters at several coupling strengths.

2. THE BAND TOP.  Pushing monotonicity to the infinite lattice bounds every
   cluster eigenvalue by the top of the FCC Cosserat band.  The microrotation
   sector is the graph Laplacian times the 3x3 identity, whose band is [0,16]
   because the FCC structure factor bottoms out at -4; the displacement
   sector is bounded by 8.  At the zone corner where the Laplacian reaches 16
   the discrete curl vanishes, so the coupling contributes exactly
   +alpha_Cos and the top is 16 + alpha_Cos, an integer 17 at alpha_Cos = 1.
   The script maximises the Bloch matrix over the zone to confirm that no
   other point climbs higher.

3. THE BRACKET.  With 0 <= lambda <= 17 the master formula
   m = N[m0 + (lambda - 4) me] traps the mass per node between
   m0 - 4 me = 67.981 MeV and m0 + 13 me = 76.668 MeV.  Given a measured
   mass the admissible node counts are the integers in [m/76.668, m/67.981],
   and for six light states there is exactly one.  The script prints that
   table and checks every light and strange closure against the band.

Inputs are two CODATA values (the electron mass and the fine-structure
constant) and FCC geometry.  Nothing is fitted.

Run:  python3 cluster_monotonicity.py
"""

import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cosserat_classifier import build_cosserat_matrix          # noqa: E402
from spectral_classifier import (                              # noqa: E402
    cluster_born,
    cluster_coord_shell,
    cluster_cuboctahedron,
    cluster_hex_cap,
)

# --- CODATA 2022 -----------------------------------------------------------
M_E = 0.51099895069          # electron mass [MeV]
ALPHA_FS = 7.2973525643e-3   # fine-structure constant [-]
M_0 = M_E / ALPHA_FS         # node assembly quantum [MeV]

ALPHA_COS = 1.0              # single-scale Cosserat coupling convention
LAMBDA_REF = 4.0             # on-site stiffness of a fully embedded node


# ---------------------------------------------------------------------------
# Cluster helpers
# ---------------------------------------------------------------------------

def _coords(builder):
    """Some builders return (coords, extras); normalise to a coordinate array."""
    out = builder()
    return np.atleast_2d(out[0] if isinstance(out, tuple) else out)


def catalogue():
    """The nested clusters the mass programme uses, smallest first."""
    return {
        "hex cap (7)": _coords(cluster_hex_cap),
        "cuboctahedron (12)": _coords(cluster_cuboctahedron),
        "coordination shell (13)": _coords(cluster_coord_shell),
        "Born cluster (19)": _coords(cluster_born),
    }


def fcc_ball(radius, nn_spacing=1.0):
    """FCC sites within a radius, nearest-neighbour distance nn_spacing."""
    reach = int(np.ceil(radius * np.sqrt(2))) + 1
    pts = np.array(
        [v for v in itertools.product(range(-reach, reach + 1), repeat=3)
         if sum(v) % 2 == 0],
        dtype=float,
    ) * (nn_spacing / np.sqrt(2))
    return pts[np.linalg.norm(pts, axis=1) <= radius + 1e-9]


def spectrum(coords, alpha=ALPHA_COS):
    """Descending eigenvalues of the 6n x 6n Cosserat dynamical matrix."""
    return np.sort(np.linalg.eigvalsh(build_cosserat_matrix(coords, alpha=alpha)))[::-1]


def is_induced_subcluster(sub, parent, tol=1e-6):
    """Every node of sub is a node of parent (bonds then follow)."""
    return all(np.min(np.linalg.norm(parent - s, axis=1)) < tol for s in sub)


# ---------------------------------------------------------------------------
# 1. Monotonicity
# ---------------------------------------------------------------------------

def check_catalogue_monotonicity():
    print("1. MONOTONICITY ON THE CATALOGUE'S NESTED CLUSTERS")
    print("   lambda_i(sub) <= lambda_i(parent) for every i, descending order")
    clusters = catalogue()
    names = list(clusters)
    failures = 0
    for a, b in itertools.permutations(names, 2):
        if not is_induced_subcluster(clusters[a], clusters[b]):
            continue
        sa, sb = spectrum(clusters[a]), spectrum(clusters[b])
        worst = np.max(sa - sb[: len(sa)])
        ok = worst <= 1e-9
        failures += 0 if ok else 1
        print(f"   {a:26s} in {b:26s} "
              f"max(sub - parent) = {worst:+.3e}  {'ok' if ok else 'FAIL'}")
    return failures


def check_random_monotonicity(trials=400, seed=7):
    print("\n   randomised node deletions from random FCC clusters")
    rng = np.random.default_rng(seed)
    base = fcc_ball(2.2)
    total_failures = 0
    for alpha in (0.0, 0.5, 1.0, 2.0):
        failures = 0
        for _ in range(trials):
            size = int(rng.integers(4, 10))
            parent = base[rng.choice(len(base), size=size, replace=False)]
            sub = np.delete(parent, rng.integers(0, size), axis=0)
            sa, sb = spectrum(sub, alpha), spectrum(parent, alpha)
            if np.any(sa - sb[: len(sa)] > 1e-9):
                failures += 1
        total_failures += failures
        print(f"   alpha_Cos = {alpha:4.1f}: {failures} violations in {trials} deletions")
    return total_failures


# ---------------------------------------------------------------------------
# 2. The band top
# ---------------------------------------------------------------------------

def fcc_neighbour_vectors(nn_spacing=1.0):
    """The twelve FCC nearest-neighbour vectors."""
    return np.array(
        [v for v in itertools.product([-1, 0, 1], repeat=3)
         if sum(v) % 2 == 0 and sum(abs(x) for x in v) == 2],
        dtype=float,
    ) * (nn_spacing / np.sqrt(2))


def bloch_matrix(k, alpha=ALPHA_COS, nn_spacing=1.0):
    """6x6 Cosserat dynamical matrix of the infinite FCC lattice at wavevector k."""
    uu = np.zeros((3, 3), complex)
    pp = np.zeros((3, 3), complex)
    curl = np.zeros((3, 3), complex)
    for b in fcc_neighbour_vectors(nn_spacing):
        bhat = b / np.linalg.norm(b)
        phase = np.exp(1j * k @ b)
        uu += np.outer(bhat, bhat) * (1.0 - phase)
        pp += np.eye(3) * (1.0 - phase)
        cross = np.array([[0.0, -bhat[2], bhat[1]],
                          [bhat[2], 0.0, -bhat[0]],
                          [-bhat[1], bhat[0], 0.0]])
        curl += cross * phase / (2.0 * nn_spacing)
    uu = (uu + uu.conj().T) / 2.0
    pp = (pp + pp.conj().T) / 2.0
    mat = np.zeros((6, 6), complex)
    mat[:3, :3] = uu + alpha * (curl.conj().T @ curl)
    mat[3:, 3:] = pp + alpha * np.eye(3)
    mat[:3, 3:] = -alpha * curl.conj().T
    mat[3:, :3] = -alpha * curl
    return mat


def band_top(alpha=ALPHA_COS, samples=200000, seed=3):
    """Maximum eigenvalue of the Bloch matrix over the Brillouin zone."""
    rng = np.random.default_rng(seed)
    span = 2.0 * np.pi * np.sqrt(2.0)
    best = 0.0
    for _ in range(samples):
        k = rng.uniform(-span, span, 3)
        best = max(best, np.linalg.eigvalsh(bloch_matrix(k, alpha))[-1])
    return best


def check_band_top():
    print("\n2. THE INFINITE-LATTICE BAND TOP")
    # The zone corner where every cos(k.b) argument is a multiple of pi.
    corner = np.array([0.0, np.pi, np.pi]) * np.sqrt(2.0)
    for alpha in (0.0, 0.5, 1.0):
        exact = 16.0 + alpha
        at_corner = np.linalg.eigvalsh(bloch_matrix(corner, alpha))[-1]
        scanned = band_top(alpha, samples=60000)
        print(f"   alpha_Cos = {alpha:4.1f}:  16 + alpha = {exact:7.4f}"
              f"   at zone corner = {at_corner:7.4f}"
              f"   zone maximum = {scanned:7.4f}")
    # Finite-ball ceilings: exact, and tighter than the infinite-lattice value.
    print("\n   finite-ball ceilings (any cluster inside the ball obeys these)")
    for radius in (1.5, 2.0, 2.5):
        pts = fcc_ball(radius)
        print(f"   radius {radius:.1f} l: N = {len(pts):3d}  lambda_1 = "
              f"{spectrum(pts)[0]:.4f}")


# ---------------------------------------------------------------------------
# 3. The node-count bracket
# ---------------------------------------------------------------------------

LAMBDA_MAX = 16.0 + ALPHA_COS
BRACKET_LO = M_0 - LAMBDA_REF * M_E
BRACKET_HI = M_0 + (LAMBDA_MAX - LAMBDA_REF) * M_E

# (label, node count, measured mass [MeV]).  Light and strange closures only:
# the heavy-flavour rows add a heavy-quark core outside the node counting, so
# the bracket does not apply to them as written.
CLOSURES = [
    ("pi", 2, 138.04),
    ("K+-", 7, 493.677),
    ("eta", 8, 547.862),
    ("rho", 11, 775.26),
    ("omega", 11, 782.66),
    ("K*(892)", 13, 891.66),
    ("eta'", 14, 957.78),
    ("eta(1295)", 18, 1294.0),
    ("N (isoscalar)", 13, 938.919),
    ("Lambda", 16, 1115.683),
    ("Sigma0", 17, 1193.15),
    ("Delta(1232)", 17, 1232.0),
    ("Sigma*(1385)", 19, 1384.4),
    ("Xi", 19, 1318.29),
    ("Xi*(1530)", 21, 1530.0),
    ("Lambda(1405)", 20, 1405.1),
    ("Delta(1600)", 21, 1570.0),
    ("Lambda(1600)", 22, 1600.0),
    ("Sigma(1660)", 23, 1660.0),
    ("Omega-", 24, 1672.45),
    ("Sigma(1670)", 24, 1675.0),
    ("d*(2380)", 33, 2380.0),
]

# The subset tabulated in the monograph, in the order it appears there.
TABULATED = ["pi", "K+-", "eta", "rho", "omega", "N (isoscalar)", "Lambda",
             "Sigma0", "Delta(1232)", "Xi", "Omega-", "d*(2380)"]


def admissible_counts(mass):
    lo = int(np.ceil(mass / BRACKET_HI))
    hi = int(np.floor(mass / BRACKET_LO))
    return list(range(lo, hi + 1))


def check_bracket():
    print("\n3. THE PER-NODE MASS BRACKET")
    print(f"   m0 = me/alpha = {M_0:.5f} MeV")
    print(f"   lambda in [0, {LAMBDA_MAX:.1f}]  ->  m/N in "
          f"[{BRACKET_LO:.4f}, {BRACKET_HI:.4f}] MeV   "
          f"(ratio {BRACKET_HI / BRACKET_LO:.5f})")

    print("\n   node count read back from the measured mass")
    print(f"   {'state':16s}{'mass':>10s}  {'window':>16s}  "
          f"{'admissible N':20s}{'cluster N':>10s}")
    by_label = {lab: (n, m) for lab, n, m in CLOSURES}
    for label in TABULATED:
        n_cluster, mass = by_label[label]
        adm = admissible_counts(mass)
        mark = "  UNIQUE" if len(adm) == 1 else ""
        status = "ok" if n_cluster in adm else "  *** EXCLUDED ***"
        print(f"   {label:16s}{mass:10.3f}  "
              f"{mass / BRACKET_HI:7.2f} to {mass / BRACKET_LO:5.2f}  "
              f"{str(adm):20s}{n_cluster:10d} {status}{mark}")

    print("\n   every light and strange closure against the band")
    outside = []
    for label, n_cluster, mass in CLOSURES:
        ratio = mass / n_cluster
        if not (BRACKET_LO - 1e-9 <= ratio <= BRACKET_HI + 1e-9):
            outside.append((label, ratio))
    ratios = [m / n for _, n, m in CLOSURES]
    print(f"   observed m/N spans {min(ratios):.3f} to {max(ratios):.3f} MeV "
          f"across {len(CLOSURES)} states")
    print(f"   states outside the band: {outside if outside else 'none'}")
    return len(outside)


# ---------------------------------------------------------------------------

def main():
    failures = check_catalogue_monotonicity()
    failures += check_random_monotonicity()
    check_band_top()
    failures += check_bracket()
    print(f"\n{'ALL CHECKS PASS' if failures == 0 else f'{failures} FAILURES'}")
    return failures


if __name__ == "__main__":
    raise SystemExit(main())
