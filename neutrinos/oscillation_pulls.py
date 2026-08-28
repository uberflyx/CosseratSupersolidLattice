#!/usr/bin/env python3
"""
oscillation_pulls.py
====================
Verification of the framework's neutrino-oscillation confrontations
(Confrontation chapter and neutrino-oscillations chapter).

Recomputes every quoted pull against the verified mid-2026 references:

- NuFIT 6.1 (v61 parameter table), BOTH variants: IC23 without
  SK-atmospheric, and IC24 with SK-atmospheric. 6.1 already folds in the
  JUNO first measurement (Abusleme et al., arXiv 2511.14593; 59.1 days),
  so the solar pair is taken from the global fit rather than from JUNO
  as well; quoting both would enter the same measurement twice.
  The two variants agree on the theta_23 octant, returning 0.470
  (+0.017 -0.014) in each, where 6.0 put the no-SK best fit at 0.561 in
  the upper octant. One dataset swap did that, IC23 replacing IC19, so
  the agreement is a movement in the fit and not a resolved octant. The
  theta_23 likelihood is still bimodal and the 3-sigma ranges still span
  both octants (0.432-0.587 no SK, 0.435-0.584 with SK), which remains
  the sturdier statement. The variants do still differ on Delta m^2_3l
  and on delta_CP, and both are reported.
- Post-JUNO global fit (Esteban et al., JHEP 04 (2026) 089): Normal
  Ordering preferred at Delta chi^2 = 4.6 to 9.4, the two values being
  without and with the SK-atm and IceCube-24 samples. Note that the
  "without" case still carries IC23 atmospheric data; the contrast is
  SK-atm plus IC24 against IC23 alone, not atmospherics against none.
  The framework REQUIRES Normal Ordering.

Pulls are (prediction - measurement)/sigma with the sigma on the side the
prediction falls; asymmetric errors handled accordingly.
"""
import numpy as np

# ---- framework predictions (zero free parameters) ----
pred = {
    "Dm21 [1e-5 eV^2]":  7.47,
    "Dm31 [1e-3 eV^2]":  2.523,
    "sin2_th12":         0.3088,
    "sin2_th23":         0.4682,
    "sin2_th13":         0.0220,
    "deltaCP [deg]":     182.0,
}

# ---- verified measurements ----
# NuFIT 6.1 IC23 WITHOUT SK-atm, Normal Ordering (v61 parameter table).
# 6.1 already includes the JUNO first result, so the solar pair is taken from
# the same fit rather than quoted separately; quoting both would double-count it.
nufit_noSK = {
    "Dm21 [1e-5 eV^2]": (7.537, 0.094, 0.10),
    "Dm31 [1e-3 eV^2]": (2.521, 0.026, 0.018),
    "sin2_th12":        (0.3088, 0.0067, 0.0066),
    "sin2_th23":        (0.470, 0.017, 0.014),
    "sin2_th13":        (0.02249, 0.00057, 0.00057),
    "deltaCP [deg]":    (207.0, 23.0, 20.0),
}
# NuFIT 6.1 IC24 WITH SK-atm, Normal Ordering (v61 parameter table);
# theta_23 sits in the FIRST octant in this variant, on the prediction.
nufit_SK = {
    "Dm21 [1e-5 eV^2]": (7.537, 0.094, 0.10),
    "Dm31 [1e-3 eV^2]": (2.511, 0.021, 0.020),
    "sin2_th12":        (0.3088, 0.0067, 0.0066),
    "sin2_th23":        (0.470, 0.017, 0.014),
    "sin2_th13":        (0.02248, 0.00055, 0.00059),
    "deltaCP [deg]":    (212.0, 26.0, 36.0),
}


def pull(p, meas):
    c, up, dn = meas
    s = up if p > c else dn
    return (p - c) / s


def report(name, ref):
    chi2, n = 0.0, 0
    print(f"\n--- pulls vs {name} ---")
    for k, m in ref.items():
        z = pull(pred[k], m)
        chi2 += z**2; n += 1
        print(f"  {k:18s}: pred {pred[k]:8.4f}  meas {m[0]:8.4f} "
              f"(+{m[1]:.4f}/-{m[2]:.4f})  pull {z:+.2f} sigma")
    return chi2, n


if __name__ == "__main__":
    c3, n3 = report("NuFIT 6.1 WITH SK-atm (primary anchor)", nufit_SK)
    print(f"\nsix-observable chi^2 (NuFIT 6.1 with SK-atm, single anchor) "
          f"= {c3:.2f} for {n3} observables")

    c2, n2 = report("NuFIT 6.1 without SK-atm (contrast: same octant, different dm31 and dcp)", nufit_noSK)
    print(f"\nWith the SK-atm variant the theta_23 Gaussianised pull is "
          f"{pull(pred['sin2_th23'], nufit_SK['sin2_th23']):+.1f} sigma, but the "
          f"likelihood is bimodal:")
    print("  the 3-sigma ranges of BOTH variants span both octants, so the")
    print("  prediction remains inside 3 sigma of either fit. The octant is the")
    print("  one genuinely unresolved parameter (NuFit 6.0/6.1 abstract), the")
    print("  datasets split (T2K+NOvA joint: upper, Bayes factor 3.5; SK-atm:")
    print("  lower), and the framework stakes a sharp lower-octant claim:")
    print("  sin^2 th23 = (1 - m_mu/m_tau)/2 = 0.470, half the mu-to-tau mass")
    print("  ratio below maximal; the inter-family coupling is circulant and")
    print("  adds nothing.")
    print("\nMass ordering: framework REQUIRES Normal Ordering. Post-JUNO global")
    print("fit (JHEP 04 (2026) 089, Eq. 5.1): NO preferred at Delta chi^2 = 4.62")
    print("without and 9.41 with the SK-atm and IceCube-24 samples (about 2.1 to")
    print("3.1 sigma). The 'without' case still carries IC23 atmospheric data.")
