#!/usr/bin/env python3
"""
strange_activation.py
=====================
Does the wound (strange) sector pay a crossing cost at a character defect
that the static sector does not? The D4 counting argument says yes: the
swapped bond's cross-layer sheath holds four lobes, dormant for light
flavours by the matched-ends cancellation, activated by a strange winding
(hadrons/isospin_isotensor_count.py, section 10). This script gives that
claim its first dynamical test on the existing 20x20 transfer matrix
(d4_transfer_matrix_zeros.py): a three-period compact loop with one
period carrying a character defect, and the channel eigenvalue tracked
between the clean loop T^3 and the defective loop T^2 D.

Defect models (the charge/character couples chirally; the chiral model
won the delta_cs survey in spectral_mass/hyperon_tensor_admixture.py):
  chiral : kappa_c -> kappa_c (1 + g) inside the defect period
  barrier: V0      -> V0      (1 + g) inside the defect period

Winding values. Two conventions coexist in the d4 scripts: the zeros
script calls k4*ell = 2*pi/3 the first KK level (layer spacing ell),
while the chirality resummation uses k4 = 2*pi/(sqrt(6)*ell), i.e.
k4*ell = 2.566 (three layers on the sqrt(6)*ell loop, phase 2*pi/3 per
interface). Both are scanned; the conclusions below hold for either.

Findings
--------
1. Static slopes are clean and linear: +1.01 g per defect period for
   the chiral model, -2.28 g for the barrier model.
2. At the nominal k4*ell = 2*pi/3 there is NO differential: the wound
   slope equals the static slope within a percent in both models. At
   the physical winding of the sqrt(6) ell loop, k4 = 2*pi/(sqrt(6)
   ell), the differential is order unity and even in k4 in both models
   (chiral -1.24, barrier +1.03, in static units). The matched-ends
   activation switches on exactly at the strange's true winding, which
   settles the convention split between the two d4 scripts in favour
   of the resummation value.
3. The k4-even parity means anti-hyperons pay the same crossing cost,
   as CPT requires.
4. Caveat and next steps: at the physical winding the loop spectrum is
   crowded and the small-g extraction jitters; degenerate-aware
   tracking is the next refinement, followed by the calibration chain
   from d ln|lambda| per defect period to electron masses per link,
   which would turn the O(1) statement into the counted four.

Accompanies: M. A. Cox, "The Cosserat Supersolid" (2026).
https://doi.org/10.5281/zenodo.18636501

Author: Mitchell A. Cox, University of the Witwatersrand
"""
import contextlib
import io
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load_machinery():
    """Import d4_transfer_matrix_zeros quietly (it prints on import)."""
    spec = importlib.util.spec_from_file_location(
        "tmz", os.path.join(HERE, "d4_transfer_matrix_zeros.py"))
    tmz = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(tmz)
    return tmz


TMZ = load_machinery()
ALPHA = 1.0 / 137.0359992


def em_track(T_loop, v_ref):
    """EM-channel eigenvalue of a (non-symmetric) loop matrix, tracked by
    maximal eigenvector overlap with the single-period EM eigenvector
    (the loop and the clean period share eigenvectors, since T^3 v =
    lambda^3 v)."""
    ev, V = np.linalg.eig(T_loop)
    ov = np.abs(np.conj(v_ref) @ V) / np.linalg.norm(V, axis=0)
    k = int(np.argmax(ov))
    return ev[k], float(ov[k])


def em_anchor(T_period):
    """EM eigenpair of ONE period: |lambda| nearest alpha, the script's
    own calibration point."""
    ev, V = np.linalg.eig(T_period)
    k = int(np.argmin(np.abs(np.abs(ev) - ALPHA)))
    return ev[k], V[:, k] / np.linalg.norm(V[:, k])


def period(k4, model=None, g=0.0, n_steps=4000):
    """One-period transfer matrix, optionally with a defect model."""
    if model is None or g == 0.0:
        return TMZ.transfer_matrix(TMZ.M20, k4, n_steps=n_steps)
    if model == "chiral":
        kc0, mt0 = TMZ.kappa_c, TMZ.mu_tot
        TMZ.kappa_c = kc0 * (1.0 + g)
        TMZ.mu_tot = TMZ.mu + TMZ.kappa_c
        D = TMZ.transfer_matrix(TMZ.M20, k4, n_steps=n_steps)
        TMZ.kappa_c, TMZ.mu_tot = kc0, mt0
        return D
    if model == "barrier":
        v0 = TMZ.V0
        TMZ.V0 = v0 * (1.0 + g)
        D = TMZ.transfer_matrix(TMZ.M20, k4, n_steps=n_steps)
        TMZ.V0 = v0
        return D
    raise ValueError(model)


def crossing(k4, model, g, n_steps=4000):
    """d ln|lambda_EM| of the 3-period compact loop when one period
    carries the defect."""
    T = period(k4, n_steps=n_steps)
    _, v0 = em_anchor(T)
    lam0, ov0 = em_track(T @ T @ T, v0)
    D = period(k4, model, g, n_steps=n_steps)
    lam1, ov1 = em_track(T @ T @ D, v0)
    if min(ov0, ov1) < 0.5:
        return float('nan')
    return float(np.log(abs(lam1)) - np.log(abs(lam0)))


def main():
    k4_grid = [0.0, 2*np.pi/3, -2*np.pi/3,
               2*np.pi/np.sqrt(6), -2*np.pi/np.sqrt(6)]
    gs = (0.025, 0.05, 0.10)

    print("=" * 72)
    print("1. Crossing cost per defect period: d ln|lambda_EM|, 3-period loop")
    print("=" * 72)
    slopes = {}
    for model in ("chiral", "barrier"):
        print(f"\n  model: {model}")
        print(f"  {'k4*ell':>9}" + "".join(f"{f'g={g}':>13}" for g in gs)
              + f"{'slope':>11}")
        for k4 in k4_grid:
            vals = [crossing(k4, model, g) for g in gs]
            slope = np.polyfit(gs, vals, 1)[0]
            slopes[(model, round(k4, 4))] = slope
            print(f"  {k4:>9.4f}" + "".join(f"{v:>13.5e}" for v in vals)
                  + f"{slope:>11.4f}")

    print("\n" + "=" * 72)
    print("2. Parity, differential, and the activation ratio")
    print("=" * 72)
    for model in ("chiral", "barrier"):
        a0 = slopes[(model, 0.0)]
        for k4 in (2*np.pi/3, 2*np.pi/np.sqrt(6)):
            ap = slopes[(model, round(k4, 4))]
            am = slopes[(model, round(-k4, 4))]
            even = 0.5 * (ap + am)
            odd = 0.5 * (ap - am)
            print(f"  {model:8s} k4*ell={k4:.3f}: slope {ap:+.4f}, parity "
                  f"even {even:+.4f} / odd {odd:+.1e}; "
                  f"differential (wound - static) = {ap - a0:+.4f}; "
                  f"|diff/static| = {abs((ap - a0)/a0):.2f}")

    print("\n" + "=" * 72)
    print("3. Convergence (chiral, k4*ell = 2*pi/3, g = 0.05)")
    print("=" * 72)
    for n in (2000, 4000, 8000):
        print(f"  n_steps = {n}: {crossing(2*np.pi/3, 'chiral', 0.05, n):.6e}")

    print("""
Verdict: at the physical winding 2*pi/(sqrt(6) ell) the wound sector
engages a character defect with an order-unity, k4-even differential in
both defect models, while at the nominal 2*pi/3 it does not. The
activation behind the ladder step is dynamically real and sits at the
strange's true winding. Open: degenerate-aware tracking at the crowded
winding, then the d ln|lambda| -> m_e per link calibration chain.""")


if __name__ == "__main__":
    main()
