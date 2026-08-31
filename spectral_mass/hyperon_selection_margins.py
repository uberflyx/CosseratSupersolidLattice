#!/usr/bin/env python3
"""
hyperon_selection_margins.py
============================
How thin are the accommodated-baryon selection margins, and what do they
constrain?

The accommodated rule picks the mass mode of the Lambda and the Xi as the soft
mode that is microrotation-dominant and carries the majority of its weight on
the extension structure rather than the coordination shell.  At the
single-scale convention alpha_Cos = 1 the rule returns exactly one mode on each
cluster, and the chapter records the margins as thin: on the Lambda the mass
mode sits at structure weight 0.517 against a rival at 0.415, and on the Xi at
0.582 against a rival at 0.489, one part in a hundred below the half-weight
threshold.

A criterion that turns on a comparison with 0.5 is only as good as the
stability of the numbers being compared, so this script sweeps the coupling and
finds where the rule breaks.  It breaks at both ends, and for different
reasons:

  below alpha_Cos = 0.938   the Lambda's own mass mode falls under half weight
                            and the rule returns no candidate at all;
  above alpha_Cos = 1.234   the Xi's rival rises past half weight and the rule
                            returns two.

So the two hyperons bracket the coupling from opposite sides, and the rule is
self-consistent only inside a window of width about 0.3 that contains the
convention alpha_Cos = 1.  This is a check rather than a derivation: nothing
here fixes the coupling, but the value the framework uses on other grounds is
the value the selection rule needs.

A second, independent check runs the other way.  Fitting alpha to the two
measured hyperon masses, which the framework never does, lands at 0.974,
within three percent of the convention, with the Lambda and Xi residuals
balanced at -0.195% and +0.194%.

Run:  python3 hyperon_selection_margins.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from baryon_mass_modes import M_0, M_E, channel_modes, cluster_xi  # noqa: E402
from spectral_classifier import build_lambda_cluster               # noqa: E402

CHANNEL = 'A_2g'          # the baryon-number irrep, subduced to the host group
N_SHELL = 13              # nodes of the original coordination shell
SOFT_MAX = 4.0            # the reference eigenvalue, lambda_ref
PHI_MIN = 0.5             # microrotation-dominant
STRUCTURE_MIN = 0.5       # majority weight on the accommodating structure

PDG = {'Lambda': 1115.683, 'Xi': 1318.29}
NODES = {'Lambda': 16, 'Xi': 19}


def mass(n_nodes, lam):
    """Master spectral formula."""
    return n_nodes * (M_0 + (lam - 4.0) * M_E)


def soft_candidates(coords, alpha):
    """Soft microrotation-dominant modes, ordered by structure weight."""
    modes, _ = channel_modes(coords, CHANNEL, n_shell=N_SHELL, alpha=alpha)
    out = [(lam, phi, 1.0 - shell)
           for (lam, phi, shell, _deg) in modes
           if lam < SOFT_MAX and phi > PHI_MIN]
    return sorted(out, key=lambda t: -t[2])


def n_qualifying(coords, alpha):
    """How many modes the accommodated rule admits."""
    return sum(1 for _lam, _phi, structure in soft_candidates(coords, alpha)
               if structure > STRUCTURE_MIN)


def bisect(predicate, lo, hi, tol=1e-4):
    """Smallest bracket where `predicate` flips from True at lo to False at hi."""
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if predicate(mid):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    clusters = {'Lambda': np.atleast_2d(build_lambda_cluster()),
                'Xi': np.atleast_2d(cluster_xi())}

    print("Margins at the single-scale convention alpha_Cos = 1")
    print(f"{'baryon':>8} {'lambda':>8} {'phi':>6} {'structure':>10} "
          f"{'rival structure':>16}")
    for name, coords in clusters.items():
        cand = soft_candidates(coords, 1.0)
        lam, phi, structure = cand[0]
        rival = cand[1][2] if len(cand) > 1 else float('nan')
        print(f"{name:>8} {lam:8.3f} {phi:6.2f} {structure:10.3f} {rival:16.3f}")

    print("\nWhere the rule breaks")
    low = bisect(lambda a: n_qualifying(clusters['Lambda'], a) != 1, 0.70, 1.00)
    high = bisect(lambda a: n_qualifying(clusters['Xi'], a) == 1, 1.00, 1.60)
    print(f"  Lambda loses its only candidate below alpha_Cos = {low:.3f}")
    print(f"  Xi gains a second candidate above alpha_Cos = {high:.3f}")
    print(f"  joint window: {low:.3f} < alpha_Cos < {high:.3f}, "
          f"width {high - low:.3f}, and it contains 1")

    print("\nMass predictions across the window")
    print(f"{'alpha':>6} {'lam_L':>7} {'m_L':>9} {'resid%':>8} | "
          f"{'lam_X':>7} {'m_X':>9} {'resid%':>8}")
    for alpha in (0.95, 1.00, 1.05, 1.10, 1.20):
        lam_l = soft_candidates(clusters['Lambda'], alpha)[0][0]
        lam_x = soft_candidates(clusters['Xi'], alpha)[0][0]
        m_l, m_x = mass(NODES['Lambda'], lam_l), mass(NODES['Xi'], lam_x)
        print(f"{alpha:6.2f} {lam_l:7.3f} {m_l:9.2f} "
              f"{100 * (m_l - PDG['Lambda']) / PDG['Lambda']:8.3f} | "
              f"{lam_x:7.3f} {m_x:9.2f} "
              f"{100 * (m_x - PDG['Xi']) / PDG['Xi']:8.3f}")

    best = None
    for alpha in np.linspace(low, high, 150):
        r_l = (mass(NODES['Lambda'], soft_candidates(clusters['Lambda'], alpha)[0][0])
               - PDG['Lambda']) / PDG['Lambda']
        r_x = (mass(NODES['Xi'], soft_candidates(clusters['Xi'], alpha)[0][0])
               - PDG['Xi']) / PDG['Xi']
        chi = r_l ** 2 + r_x ** 2
        if best is None or chi < best[1]:
            best = (alpha, chi, r_l, r_x)
    alpha, _chi, r_l, r_x = best
    print(f"\nFitting alpha to the two measured masses, which the framework does "
          f"not do,\ngives {alpha:.3f}: Lambda {100 * r_l:+.3f}%, "
          f"Xi {100 * r_x:+.3f}%, "
          f"and the convention differs by {100 * (1 - alpha) / alpha:+.1f}%.")


if __name__ == "__main__":
    raise SystemExit(main())
