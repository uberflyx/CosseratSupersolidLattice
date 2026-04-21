"""
PDG comparison: compare the generated catalogue against observed hadron masses.

PDG values taken from the project's pdg.pdf reference (to be updated from file if needed).
This module encodes a subset of well-established ground-state hadrons with their
quantum numbers and measured masses, then matches them against the constructor's
predictions.

The catalogue is not meant to predict individual-state masses with better than
a few percent accuracy at leading order (since the Q-correction from the full
electromagnetic algorithm is not applied here). The goal is to verify that the
structural-coding framework produces the right state for each observed hadron.
"""
from __future__ import annotations

from fractions import Fraction
from dataclasses import dataclass

from coding import SMLabels, sm_to_structural
from topology import stage_a_topology_class
from mass import mass_prediction


@dataclass
class PDGEntry:
    name: str
    sm: SMLabels
    mass_obs_MeV: float


PDG_OBSERVATIONS: list[PDGEntry] = [
    PDGEntry("pi+/-", SMLabels(B=0, S=0, I=Fraction(1), I_3=Fraction(1), J=Fraction(0), P=-1), 139.570),
    PDGEntry("K+/-", SMLabels(B=0, S=1, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(0), P=-1), 493.677),
    PDGEntry("eta", SMLabels(B=0, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(0), P=-1), 547.862),
    PDGEntry("eta'", SMLabels(B=0, S=2, I=Fraction(0), I_3=Fraction(0), J=Fraction(0), P=-1), 957.78),
    PDGEntry("rho", SMLabels(B=0, S=0, I=Fraction(1), I_3=Fraction(0), J=Fraction(1), P=-1), 775.26),
    PDGEntry("omega", SMLabels(B=0, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(1), P=-1), 782.66),
    PDGEntry("K*", SMLabels(B=0, S=1, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(1), P=-1), 893.6),
    PDGEntry("phi", SMLabels(B=0, S=2, I=Fraction(0), I_3=Fraction(0), J=Fraction(1), P=-1), 1019.46),
    PDGEntry("f_2(1270)", SMLabels(B=0, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(2), P=+1), 1275.5),
    PDGEntry("p", SMLabels(B=1, S=0, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(1, 2), P=+1), 938.272),
    PDGEntry("n", SMLabels(B=1, S=0, I=Fraction(1, 2), I_3=Fraction(-1, 2), J=Fraction(1, 2), P=+1), 939.565),
    PDGEntry("Lambda", SMLabels(B=1, S=1, I=Fraction(0), I_3=Fraction(0), J=Fraction(1, 2), P=+1), 1115.683),
    PDGEntry("Sigma+", SMLabels(B=1, S=1, I=Fraction(1), I_3=Fraction(1), J=Fraction(1, 2), P=+1), 1189.37),
    PDGEntry("Sigma0", SMLabels(B=1, S=1, I=Fraction(1), I_3=Fraction(0), J=Fraction(1, 2), P=+1), 1192.64),
    PDGEntry("Sigma-", SMLabels(B=1, S=1, I=Fraction(1), I_3=Fraction(-1), J=Fraction(1, 2), P=+1), 1197.45),
    PDGEntry("Xi0", SMLabels(B=1, S=2, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(1, 2), P=+1), 1314.86),
    PDGEntry("Xi-", SMLabels(B=1, S=2, I=Fraction(1, 2), I_3=Fraction(-1, 2), J=Fraction(1, 2), P=+1), 1321.71),
    PDGEntry("Delta++", SMLabels(B=1, S=0, I=Fraction(3, 2), I_3=Fraction(3, 2), J=Fraction(3, 2), P=+1), 1232.0),
    PDGEntry("Sigma*+", SMLabels(B=1, S=1, I=Fraction(1), I_3=Fraction(1), J=Fraction(3, 2), P=+1), 1382.8),
    PDGEntry("Xi*0", SMLabels(B=1, S=2, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(3, 2), P=+1), 1531.80),
    PDGEntry("Omega-", SMLabels(B=1, S=3, I=Fraction(0), I_3=Fraction(0), J=Fraction(3, 2), P=+1), 1672.45),
    PDGEntry("D0", SMLabels(B=0, S=0, I=Fraction(1, 2), I_3=Fraction(-1, 2), J=Fraction(0), P=-1, n_c=1), 1864.84),
    PDGEntry("D+", SMLabels(B=0, S=0, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(0), P=-1, n_c=1), 1869.66),
    PDGEntry("Ds+", SMLabels(B=0, S=1, I=Fraction(0), I_3=Fraction(0), J=Fraction(0), P=-1, n_c=1), 1968.35),
    PDGEntry("D*0", SMLabels(B=0, S=0, I=Fraction(1, 2), I_3=Fraction(-1, 2), J=Fraction(1), P=-1, n_c=1), 2006.85),
    PDGEntry("eta_c", SMLabels(B=0, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(0), P=-1, n_c=2), 2983.9),
    PDGEntry("J/psi", SMLabels(B=0, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(1), P=-1, n_c=2), 3096.9),
    PDGEntry("chi_c0", SMLabels(B=0, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(0), P=+1, n_c=2), 3414.71),
    PDGEntry("chi_c1", SMLabels(B=0, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(1), P=+1, n_c=2), 3510.67),
    PDGEntry("Lambda_c+", SMLabels(B=1, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(1, 2), P=+1, n_c=1), 2286.46),
    PDGEntry("Sigma_c++", SMLabels(B=1, S=0, I=Fraction(1), I_3=Fraction(1), J=Fraction(1, 2), P=+1, n_c=1), 2453.97),
    PDGEntry("Xi_c+", SMLabels(B=1, S=1, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(1, 2), P=+1, n_c=1), 2467.71),
    PDGEntry("Omega_c0", SMLabels(B=1, S=2, I=Fraction(0), I_3=Fraction(0), J=Fraction(1, 2), P=+1, n_c=1), 2695.2),
    PDGEntry("Xi_cc++", SMLabels(B=1, S=0, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(1, 2), P=+1, n_c=2), 3621.2),
    PDGEntry("B0", SMLabels(B=0, S=0, I=Fraction(1, 2), I_3=Fraction(-1, 2), J=Fraction(0), P=-1, n_b=1), 5279.65),
    PDGEntry("B+", SMLabels(B=0, S=0, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(0), P=-1, n_b=1), 5279.34),
    PDGEntry("Bs0", SMLabels(B=0, S=1, I=Fraction(0), I_3=Fraction(0), J=Fraction(0), P=-1, n_b=1), 5366.92),
    PDGEntry("Bc+", SMLabels(B=0, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(0), P=-1, n_c=1, n_b=1), 6274.47),
    PDGEntry("Upsilon(1S)", SMLabels(B=0, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(1), P=-1, n_b=2), 9460.40),
    PDGEntry("Lambda_b0", SMLabels(B=1, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(1, 2), P=+1, n_b=1), 5619.60),
    PDGEntry("Sigma_b+", SMLabels(B=1, S=0, I=Fraction(1), I_3=Fraction(1), J=Fraction(1, 2), P=+1, n_b=1), 5810.56),
    PDGEntry("Xi_b0", SMLabels(B=1, S=1, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(1, 2), P=+1, n_b=1), 5791.9),
    PDGEntry("Omega_b-", SMLabels(B=1, S=2, I=Fraction(0), I_3=Fraction(0), J=Fraction(1, 2), P=+1, n_b=1), 6046.1),
]


def run_pdg_comparison() -> None:
    """For each PDG entry, compute the predicted mass and compare."""
    print(f"{'Name':<15} {'Predicted MeV':>15} {'Observed MeV':>15} {'Residual %':>12}  Topology")
    print("-" * 110)
    n_entries = 0
    n_within_1pct = 0
    n_within_5pct = 0
    max_residual_pct = 0.0
    total_abs_residual_pct = 0.0

    for entry in PDG_OBSERVATIONS:
        try:
            sc = sm_to_structural(entry.sm)
        except Exception as e:
            print(f"{entry.name:<15} {'SM error:':>15} {str(e)[:60]}")
            continue

        cls = stage_a_topology_class(sc)
        if cls is None:
            print(f"{entry.name:<15} {'no topology':>15} {entry.mass_obs_MeV:>15.2f}  (forbidden by constructor)")
            continue

        mass_pred, N_eff = mass_prediction(cls, J=entry.sm.J)
        residual_pct = (mass_pred - entry.mass_obs_MeV) / entry.mass_obs_MeV * 100
        n_entries += 1
        total_abs_residual_pct += abs(residual_pct)
        max_residual_pct = max(max_residual_pct, abs(residual_pct))
        if abs(residual_pct) < 1.0:
            n_within_1pct += 1
        if abs(residual_pct) < 5.0:
            n_within_5pct += 1

        N_str = str(cls.total_N.numerator) if cls.total_N.denominator == 1 else f"{cls.total_N.numerator}/{cls.total_N.denominator}"
        print(f"{entry.name:<15} {mass_pred:>15.2f} {entry.mass_obs_MeV:>15.3f} {residual_pct:>+11.2f}%  {cls.name} (N={N_str}, H={cls.residual.value})")

    print("-" * 110)
    print(f"Entries tested: {n_entries}")
    print(f"Within 1% (leading order only, no Q-correction): {n_within_1pct} ({n_within_1pct/n_entries*100:.0f}%)")
    print(f"Within 5%: {n_within_5pct} ({n_within_5pct/n_entries*100:.0f}%)")
    print(f"Mean |residual|: {total_abs_residual_pct/n_entries:.2f}%")
    print(f"Max |residual|: {max_residual_pct:.2f}%")
    print()
    print("Note: leading-order mass is m = N * m_0. The Q * m_e electromagnetic")
    print("correction (O(a) of total mass) is not applied here. With the full Q-algorithm,")
    print("residuals should reduce to O(0.1%) for ground states per the monograph.")


if __name__ == "__main__":
    run_pdg_comparison()
