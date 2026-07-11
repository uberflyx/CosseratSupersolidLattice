#!/usr/bin/env python3
"""
oscillation_pulls.py
====================
Verification of the framework's neutrino-oscillation confrontations
(Confrontation chapter and neutrino-oscillations chapter).

Recomputes every quoted pull against the verified mid-2026 references:

- JUNO first measurement (Abusleme et al., arXiv 2511.14593; 59.1 days):
  world-leading Delta m^2_21 and sin^2 theta_12.
- NuFit 6.0 (Esteban et al., JHEP 12 (2024) 216), BOTH variants:
  IC19 without SK-atmospheric, and IC24 with SK-atmospheric. The two
  variants disagree on the theta_23 octant (0.561 upper vs 0.470 first),
  so the theta_23 pull is variant-dependent and is reported honestly for
  both. The theta_23 likelihood is bimodal; Gaussianised pulls against
  the disfavoured octant overstate the tension, so the 3-sigma ranges
  (which span both octants in both variants) are the sturdier statement.
- Post-JUNO global fit (Esteban et al., JHEP 04 (2026) 089): Normal
  Ordering preferred at Delta chi^2 = 4.6 (without atmospherics) to 9.4
  (with SK + IceCube-24). The framework REQUIRES Normal Ordering.

Pulls are (prediction - measurement)/sigma with the sigma on the side the
prediction falls; asymmetric errors handled accordingly.
"""
import numpy as np

# ---- framework predictions (zero free parameters) ----
pred = {
    "Dm21 [1e-5 eV^2]":  7.47,
    "Dm31 [1e-3 eV^2]":  2.523,
    "sin2_th12":         0.3148,
    "sin2_th23":         0.5549,
    "sin2_th13":         0.0223,
    "deltaCP [deg]":     182.0,
}

# ---- verified measurements ----
# JUNO first measurement (arXiv 2511.14593), Normal Ordering
juno = {
    "Dm21 [1e-5 eV^2]": (7.50, 0.12, 0.12),
    "sin2_th12":        (0.3092, 0.0087, 0.0087),
}
# NuFit 6.0 IC19 WITHOUT SK-atm, Normal Ordering (v60 parameter table)
nufit_noSK = {
    "Dm31 [1e-3 eV^2]": (2.534, 0.025, 0.023),
    "sin2_th23":        (0.561, 0.012, 0.015),
    "sin2_th13":        (0.02195, 0.00054, 0.00058),
    "deltaCP [deg]":    (177.0, 19.0, 20.0),
}
# NuFit 6.0 IC24 WITH SK-atm, Normal Ordering (v60 parameter table);
# theta_23 flips to the FIRST octant in this variant.
nufit_SK = {
    "Dm31 [1e-3 eV^2]": (2.513, 0.021, 0.019),
    "sin2_th23":        (0.470, 0.017, 0.013),
    "deltaCP [deg]":    (212.0, 26.0, 41.0),
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
    c1, n1 = report("JUNO (2511.14593)", juno)
    c2, n2 = report("NuFit 6.0 without SK-atm", nufit_noSK)
    print(f"\nsix-observable chi^2 (JUNO solar sector + NuFit-6.0-without-SK "
          f"atmospheric sector) = {c1 + c2:.2f} for {n1 + n2} observables")

    c3, n3 = report("NuFit 6.0 WITH SK-atm (theta_23 octant flips)", nufit_SK)
    print(f"\nWith the SK-atm variant the theta_23 Gaussianised pull is "
          f"{pull(pred['sin2_th23'], nufit_SK['sin2_th23']):+.1f} sigma, but the "
          f"likelihood is bimodal:")
    print("  the 3-sigma ranges of BOTH variants span both octants, so the")
    print("  prediction remains inside 3 sigma of either fit. The octant is the")
    print("  one genuinely unresolved parameter (NuFit 6.0 abstract), the")
    print("  datasets split (T2K+NOvA joint: upper, Bayes factor 3.5; SK-atm:")
    print("  lower), and the framework stakes a sharp upper-octant claim.")
    print("\nMass ordering: framework REQUIRES Normal Ordering. Post-JUNO global")
    print("fit (JHEP 04 (2026) 089): NO preferred at Delta chi^2 = 4.6 without")
    print("atmospherics, 9.4 with SK + IceCube-24 (about 2.1 to 3.1 sigma).")
