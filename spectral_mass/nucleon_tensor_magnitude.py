#!/usr/bin/env python3
"""
nucleon_tensor_magnitude.py
===========================
Quantitative test: does the couple-stress (tensor) channel supply the
~0.8 m_e that the proton-neutron splitting is missing?

Mechanism
---------
The rest mass of the nucleon sits in the A_2u microrotation mode at
lambda = 8.303 on the 13-node coordination shell. The charge configuration
(uud vs udd) perturbs the cluster. The perturbation shifts the mass-mode
eigenvalue, and via the master spectral formula a shift d(lambda) maps to
a mass shift N * m_e * d(lambda), N = 13.

The splitting is computed by direct diagonalisation: build the perturbed
dynamical matrix for the proton (uud) and the neutron (udd), find the
mass mode in each, and take the eigenvalue difference. This captures all
orders, and its second-order part is the tensor (T_2g) admixture, since
A_2u (x) T_1u = T_2g routes the charge break into that channel and nowhere
else.

Parameter policy
----------------
The perturbation has ONE overall strength. It is fixed, not fitted, by
requiring the FIRST-order (diagonal) part of the splitting to reproduce
the framework's strong-isospin character cost m_d - m_u = 5 m_e on the
cluster. With that single calibration spent, the SECOND-order tensor part
is a prediction. The script reports the first-order, second-order, and
total, and tests robustness against the shape of the charge coupling.

What is parameter-free regardless of model
------------------------------------------
  - the selection rule A_2u (x) T_1u = T_2g (pure group theory);
  - the tensor spectrum on the shell (pure geometry);
  - the parity of the result (odd in I_3, from the I_3*Y cross term).
The MAGNITUDE depends on the coupling shape; that dependence is measured
below so the reader sees how robust the number is.

Place in the budget (see nucleon_magnetic_selfenergy.py)
-------------------------------------------------------
The number this script returns (~-0.075 MeV) is the spin-INDEPENDENT
couple-stress admixture. Its spin-DEPENDENT partner is the nucleon magnetic
self-energy (~-0.11 MeV, nucleon_magnetic_selfenergy.py). The two are the
microrotation-sector share of the splitting; the rest of the gap is the
three-body charge overlap plus the inelastic Compton part.

Author: Mitchell A. Cox, with Claude (Anthropic).  License: MIT.
"""

import numpy as np
from cosserat_classifier import build_cosserat_matrix
from spectral_classifier import fcc_nn_vectors

M_E = 0.51099895069   # MeV
N_NODES = 13
ELL = 1.0             # NN distance in builder units

# quark electric charges (screw content q = cos beta)
QU, QD = 2.0/3.0, -1.0/3.0


def cluster():
    co = fcc_nn_vectors()
    return np.vstack([np.array([[0., 0., 0.]]), co])


def mass_mode(M, coords):
    """Return (eigenvalue, eigenvector) of the phi-derived A_2u mass mode,
    identified as the eigenvalue nearest 8.303 with >80% microrotation."""
    n = len(coords)
    ev, vec = np.linalg.eigh(M)
    best, best_d = None, 1e9
    for k in range(len(ev)):
        u_c = np.sum(vec[:3*n, k]**2)
        p_c = np.sum(vec[3*n:, k]**2)
        if p_c > 0.8 and abs(ev[k] - 8.303) < best_d:
            best_d, best = abs(ev[k] - 8.303), k
    return ev[best], vec[:, best]


def charge_perturbation(coords, charges_on_quarks, quark_idx, g, shape="radial"):
    """A displacement-sector stiffness perturbation localised on the quark
    sites, proportional to the site charge. This is the photon (T_1u)
    coupling: charge sources the displacement-shear field.

    shape='radial' : stiffness g*q along the site's radial direction r-hat.
    shape='iso'    : isotropic stiffness g*q on the site (g*q * I_3).
    """
    n = len(coords)
    V = np.zeros((6*n, 6*n))
    for q, i in zip(charges_on_quarks, quark_idx):
        r = coords[i]
        rhat = r / np.linalg.norm(r)
        if shape == "radial":
            block = np.outer(rhat, rhat)
        else:
            block = np.eye(3)
        V[3*i:3*i+3, 3*i:3*i+3] += g * q * block
    return V


def splitting(coords, quark_idx, g, shape):
    """Mass-mode eigenvalue for proton (uud) and neutron (udd), and the
    first- and second-order decomposition of their difference."""
    M0 = build_cosserat_matrix(coords, 1.0, 1.0, alpha=1.0)
    lam0, psi0 = mass_mode(M0, coords)

    # charge assignments: site order (1,2,3) -> (A,B,C)
    q_p = [QU, QU, QD]   # uud
    q_n = [QU, QD, QD]   # udd

    Vp = charge_perturbation(coords, q_p, quark_idx, g, shape)
    Vn = charge_perturbation(coords, q_n, quark_idx, g, shape)

    lam_p, _ = mass_mode(M0 + Vp, coords)
    lam_n, _ = mass_mode(M0 + Vn, coords)

    # First-order (diagonal) shift for each: <psi0|V|psi0>
    f1_p = psi0 @ Vp @ psi0
    f1_n = psi0 @ Vn @ psi0
    # Total shift from diagonalisation
    tot_p = lam_p - lam0
    tot_n = lam_n - lam0
    # Second-order (the rest): tensor admixture
    s2_p = tot_p - f1_p
    s2_n = tot_n - f1_n

    return dict(lam0=lam0,
                d1=(f1_n - f1_p),       # first-order splitting in lambda  (n - p)
                d2=(s2_n - s2_p),       # second-order splitting in lambda (n - p)
                dtot=(tot_n - tot_p))   # total splitting in lambda        (n - p)


def main():
    coords = cluster()
    quark_idx = [1, 2, 3]

    print("=" * 74)
    print("Calibration: fix g so the FIRST-ORDER splitting = strong 5 m_e")
    print("=" * 74)
    # 5 m_e in eigenvalue units: d(lambda) = 5 m_e / (N m_e) = 5/N
    target_d1_lambda = 5.0 / N_NODES
    print(f"  strong target: m_d - m_u = 5 m_e  ->  d(lambda)_1 = 5/{N_NODES} = {target_d1_lambda:.5f}")

    for shape in ("radial", "iso"):
        # first-order is linear in g: calibrate with g=1, then rescale
        r1 = splitting(coords, quark_idx, g=1.0, shape=shape)
        if abs(r1["d1"]) < 1e-12:
            print(f"\n  shape='{shape}': first-order splitting vanishes; "
                  f"cannot calibrate this shape.")
            continue
        g_cal = target_d1_lambda / r1["d1"]
        r = splitting(coords, quark_idx, g=g_cal, shape=shape)
        d1_mev = N_NODES * M_E * r["d1"]
        d2_mev = N_NODES * M_E * r["d2"]
        dtot_mev = N_NODES * M_E * r["dtot"]
        print(f"\n  shape = '{shape}'   (calibrated g = {g_cal:+.4f})")
        print(f"    first-order  (strong, calibrated): {d1_mev:+.4f} MeV  "
              f"= {d1_mev/M_E:+.3f} m_e")
        print(f"    second-order (tensor, predicted) : {d2_mev:+.4f} MeV  "
              f"= {d2_mev/M_E:+.3f} m_e")
        print(f"    total mass-mode splitting        : {dtot_mev:+.4f} MeV  "
              f"= {dtot_mev/M_E:+.3f} m_e")

    print("\n" + "=" * 74)
    print("Reference targets")
    print("=" * 74)
    print(f"  observed m_n - m_p                     = +1.2933 MeV = +2.531 m_e")
    print(f"  framework strong (first order, 5 m_e)  = +2.555 MeV = +5.000 m_e")
    print(f"  framework Coulomb (separate, EM)       = -0.852 MeV = -1.667 m_e")
    print(f"  -> strong + Coulomb = +1.703 MeV; gap to observed = -0.410 MeV "
          f"(-0.80 m_e)")
    print(f"\n  The tensor second-order is the candidate for that gap. Its SIGN")
    print(f"  must reduce the splitting (neutron-heavy excess removed), i.e. be")
    print(f"  negative in (n - p).")


if __name__ == "__main__":
    main()
