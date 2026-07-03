"""Confrontation of the neutrino-sector predictions with the mid-2026 data landscape.

Predictions (fully exponentiated generator; see pmns_construction.py and
neutrino_predictions.py) against the sharpest published anchor for each
observable:

  theta12, Dm2_21 : JUNO first reactor measurement (arXiv:2511.14593),
                    59 days, 1.6-1.8x better than all previous combined.
  theta13         : NuFit 6.0 (IC24, without SK atmospherics), nu-fit.org.
  theta23, Dm2_32 : T2K+NOvA joint analysis, Nature (2025),
                    doi:10.1038/s41586-025-09599-3 (arXiv:2510.19888).
  delta_CP        : NuFit 6.0 normal-ordering best fit (without SK).
  Sum m_nu        : DESI DR2 BAO + CMB (arXiv:2503.14738, 2503.14744).

Pulls use the side-of-prediction convention for asymmetric errors: the
uncertainty on the side of the best fit where the prediction falls.
"""

import numpy as np
from scipy.stats import chi2

# Predictions (zero adjustable parameters; inputs m_e, m_mu, m_tau, alpha)
PRED = {
    "s12": 0.31480,        # sin^2 theta_12
    "s23": 0.55493,        # sin^2 theta_23
    "s13": 0.02230,        # sin^2 theta_13
    "dcp": 182.4,          # delta_CP [deg]
    "dm21": 7.4711e-5,     # eV^2
    "dm31": 2.5234e-3,     # eV^2
}
PRED["dm32"] = PRED["dm31"] - PRED["dm21"]

# Anchors: (best fit, sigma_plus, sigma_minus)
DATA = {
    "s12":  ("JUNO 2025",        0.3092,  0.0087, 0.0087),
    "dm21": ("JUNO 2025",        7.50e-5, 0.12e-5, 0.12e-5),
    "s13":  ("NuFit 6.0",        0.02195, 0.00054, 0.00058),
    "s23":  ("T2K+NOvA joint",   0.56,    0.03,   0.05),
    "dm32": ("T2K+NOvA joint",   2.43e-3, 0.04e-3, 0.03e-3),
    "dcp":  ("NuFit 6.0 (NO)",   177.0,   19.0,   20.0),
}
# Cross-checks with the alternative atmospheric anchor
ALT = {
    "s23":  ("NuFit 6.0 (no SK)", 0.561,   0.012,  0.015),
    "dm31": ("NuFit 6.0 (no SK)", 2.534e-3, 0.025e-3, 0.023e-3),
}

def pull(pred, best, sp, sm):
    """Side-of-prediction pull: sigma above the best fit if pred > best."""
    return (pred - best) / (sp if pred > best else sm)

print(f"{'observable':<8}{'anchor':<20}{'prediction':>12}{'best fit':>12}{'pull':>8}")
chisq = 0.0
for key, (name, best, sp, sm) in DATA.items():
    p = pull(PRED[key], best, sp, sm)
    chisq += p * p
    print(f"{key:<8}{name:<20}{PRED[key]:>12.5g}{best:>12.5g}{p:>+8.2f}")

print(f"\njoint chi^2 = {chisq:.2f} for 6 observables, 0 parameters")
print(f"P(chi^2 <= {chisq:.2f} | 6 dof) = {chi2.cdf(chisq, 6):.3f}"
      f"  (the fit is 'too good' at the {100*chi2.cdf(chisq,6):.0f}% level)")

print("\ncross-checks against the alternative atmospheric anchor:")
for key, (name, best, sp, sm) in ALT.items():
    p = pull(PRED[key], best, sp, sm)
    print(f"{key:<8}{name:<20}{PRED[key]:>12.5g}{best:>12.5g}{p:>+8.2f}")

# Sum of masses vs cosmology
print("\nSum m_nu = 65.5 meV vs DESI DR2 + CMB:")
print("  LambdaCDM 95% limit: 64.2 meV  -> prediction 2% above the bound")
print("  (the same analysis's Feldman-Cousins limit, 53 meV, breaches even")
print("   the minimal normal-ordering sum of 58.8 meV)")
print("  w0wa dark energy 95% limit: 160 meV -> no tension")
