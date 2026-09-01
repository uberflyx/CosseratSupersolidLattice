#!/usr/bin/env python3
"""
continuation_certificate.py
===========================
Can the adiabatic continuations be certified rather than merely observed?

The chapter identifies which bare-cluster mode a hadron's mass mode descends
from by turning the extension on and following the eigenvector by maximum
overlap between adjacent steps.  The evidence offered is the overlap itself,
which is an observation about one discretisation and not a bound.

Two theorems are available and only one of them works here.

Davis-Kahan bounds how far an eigenvector rotates under a perturbation:
sin(theta) <= ||E||_2 / delta, with delta the separation from the rest of the
spectrum.  Chained along the path the bounds add, and the sum is roughly the
path's total variation divided by the smallest gap.  Refining the path does not
help, because halving each step doubles the number of steps.  On these clusters
the spectrum is dense with symmetry-degenerate levels, delta is tiny, and the
chained bound runs to thousands of degrees.  It certifies nothing, and the
failure is structural rather than numerical.

Kato's theorem on isolated spectral subsets does work, and it asks a
topological question instead of a metric one.  If a group of eigenvalues stays
separated from the rest of the spectrum along the whole path, its spectral
projector varies continuously and its continuation is unique.  No angle is
needed: what matters is only whether the gap ever closes.

Applied to the smallest well-separated block containing the mode of interest,
this separates the two cases the raw overlaps had only hinted at.

  Roper, the proton's A_2u across the second-shell turn-on.  The block has
      dimension 15 and its gap to the rest of the spectrum never falls below
      1.234.  The lineage is certified at block level.  What the certificate
      does not fix is which member of the block, and the block spans 0.444 in
      eigenvalue, worth 4.3 MeV at N = 19.

  N(1535), both candidate seeds across the bilayer turn-on.  The blocks'
      gaps close, to 4.0e-4 and 9.1e-6.  Kato does not apply, the lineage
      cannot be certified, and that is consistent with the path dependence
      already reported for this cluster.

So the tool confirms one lineage, refuses the other, and says which of the two
it is doing.  A refusal here is a statement about the path, not a denial of the
physical assignment.

Run:  python3 continuation_certificate.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from n1535_continuation import (                                    # noqa: E402
    bond_list, clusters, weighted_cosserat,
)
from roper_continuation import as_coords                            # noqa: E402
from spectral_classifier import cluster_born, cluster_coord_shell   # noqa: E402

ELL = 1.0
PROTON_LAMBDA = 8.3028      # the proton's A_2u stiff mass mode on the bare shell
NEG_PARITY_LAMBDA = 6.1926  # the A_2g singlet the parity-flip rule assigns to 1/2^-


def path_matrices(cluster, extension, steps):
    """Bond-scaling turn-on: every bond and coupling touching the extension."""
    bonds = bond_list(cluster)
    n = len(cluster)

    def build(tau):
        bw = {(i, j): (tau if (i in extension or j in extension) else 1.0)
              for (i, j) in bonds}
        sw = {i: (tau if i in extension else 1.0) for i in range(n)}
        return weighted_cosserat(cluster, bw, sw, alpha=1.0)

    return [build(t) for t in np.linspace(1e-4, 1.0, steps)]


def certify_block(matrices, seed_lambda, window=0.2):
    """Track the smallest well-separated block containing the seed.

    Returns the block dimension, the smallest gap between the block and the
    rest of the spectrum anywhere on the path, and the block's eigenvalues at
    both ends.  A gap that never closes is a Kato certificate: the block's
    spectral projector continues uniquely along the path.
    """
    vals, vecs = np.linalg.eigh(matrices[0])
    k = int(np.argmin(np.abs(vals - seed_lambda)))
    lo = hi = k
    while lo > 0 and vals[lo] - vals[lo - 1] < window:
        lo -= 1
    while hi < len(vals) - 1 and vals[hi + 1] - vals[hi] < window:
        hi += 1
    basis = vecs[:, lo:hi + 1]
    dim = hi - lo + 1
    start = np.sort(vals[lo:hi + 1])
    min_gap = np.inf

    for mat in matrices[1:]:
        vals, vecs = np.linalg.eigh(mat)
        overlap = np.array([float(np.sum((vecs[:, i] @ basis) ** 2))
                            for i in range(len(vals))])
        pick = np.sort(np.argsort(overlap)[::-1][:dim])
        others = np.delete(np.arange(len(vals)), pick)
        gap = float(np.min(np.abs(vals[others][None, :] - vals[pick][:, None])))
        min_gap = min(min_gap, gap)
        basis = vecs[:, pick]

    return dim, min_gap, start, np.sort(vals[pick])


M_E = 0.51099895069        # electron mass [MeV]
GAP_FLOOR = 1e-2           # below this the block is not separated in practice


def report(label, dim, min_gap, start, end, n_nodes):
    """Print one block's certificate."""
    certified = min_gap > GAP_FLOOR
    spread = float(end.max() - end.min())
    print(f"  {label}")
    print(f"    block dimension        {dim}")
    print(f"    smallest gap on path   {min_gap:.3e}"
          f"   -> {'never closes: CERTIFIED' if certified else 'closes: not certified'}")
    print(f"    block at the start     {np.round(start, 4)[:6]} ...")
    print(f"    block at the end       {np.round(end, 4)[:6]} ...")
    if certified:
        print(f"    residual freedom       {spread:.4f} in eigenvalue"
              f" = {n_nodes * spread * M_E:.2f} MeV at N = {n_nodes}")


def main():
    print("Kato certificates for the extension-cluster continuations\n")

    steps = 300
    print(f"path resolution: {steps} steps\n")

    shell = as_coords(cluster_coord_shell())
    born = as_coords(cluster_born())
    idx = [int(np.argmin(np.linalg.norm(born - p, axis=1))) for p in shell]
    ext = [j for j in range(len(born)) if j not in idx]
    mats = path_matrices(born, ext, steps)
    report("Roper, proton's A_2u across the second-shell turn-on",
           *certify_block(mats, PROTON_LAMBDA), 19)

    _shell, child, _shell_idx, bilayer_idx = clusters()
    mats = path_matrices(child, bilayer_idx, steps)
    print()
    report("N(1535), proton's A_2u across the bilayer turn-on",
           *certify_block(mats, PROTON_LAMBDA), 21)
    print()
    report("N(1535), the A_2g singlet the parity rule selects",
           *certify_block(mats, NEG_PARITY_LAMBDA), 21)

    print("\nA gap that never closes is a Kato certificate: the block's spectral")
    print("projector continues uniquely, so the lineage holds at block level.")
    print("A gap that closes means the path meets an avoided crossing, and the")
    print("refusal is a statement about the path rather than about the physics.")


if __name__ == "__main__":
    raise SystemExit(main())
