"""
Leading-order mass formula for the constructor catalogue.

This module applies m = N * m_0 + Q * m_e to the topology-class output.

Only the leading-order (integer-N, zero-Q) mass is computed here. The full
Q-algorithm (bond-counting, colour-channel, surface, isospin) from monograph
Sec. 10.10 is a separate module and is not re-implemented here. The goal is
to verify that the leading-order N * m_0 reproduces the correct mass scale
for each catalogue entry, with Q providing the few-percent fine structure.

Spin corrections (fractional N values like 144/13 for rho, 27/2 for nucleons)
are applied as additive corrections to the leading-order N where applicable.
"""
from __future__ import annotations

from fractions import Fraction

from topology import TopologyClass


M_E = 0.5109989  # MeV, electron mass
ALPHA = 1.0 / 137.035999177  # fine structure constant
M_0 = M_E / ALPHA  # MeV, node mass


def leading_order_mass(cls: TopologyClass, spin_correction: Fraction = Fraction(0)) -> float:
    """Compute m = (N + spin_correction) * m_0 in MeV."""
    total_N = cls.total_N + spin_correction
    return float(total_N * M_0)


def spin_correction_from_flags(cls: TopologyClass, J: Fraction) -> Fraction:
    """Apply spin-correction rules to the topology class.

    These are the rules that produce the fractional-N results:
        - B=1 defects with the shell have winding correction +1/2.
        - M_vector_light (rho/omega on crossed fault) has delocalised spin giving
          total N = Z^2/(Z+1) = 144/13, i.e., correction = 144/13 - 11 = 1/13.
        - M_vector_strange_singlet (phi-like) has N = 29/2, correction = 1/2 on N=14.
    """
    if "winding_correction" in cls.flags:
        return Fraction(1, 2)
    if "delocalised_spin" in cls.flags:
        return Fraction(1, 13)
    if "pinned_bilayer_axial" in cls.flags:
        return Fraction(1, 2)
    return Fraction(0)


def mass_prediction(cls: TopologyClass, J: Fraction) -> tuple[float, Fraction]:
    """Returns (mass in MeV, effective N as a Fraction)."""
    correction = spin_correction_from_flags(cls, J)
    total_N_eff = cls.total_N + correction
    mass = float(total_N_eff * M_0)
    return mass, total_N_eff


def _self_test() -> None:
    from topology import BuildingBlock, PointGroup, TopologyClass
    from fractions import Fraction

    pi_cls = TopologyClass(
        name="M_pseudoscalar_light",
        blocks=(BuildingBlock.CELL_PAIR,),
        total_N=Fraction(2),
        residual=PointGroup.D_2h,
        flags=frozenset({"cell_pair", "translational_sector"}),
    )
    mass, N_eff = mass_prediction(pi_cls, J=Fraction(0))
    assert abs(mass - 140.05) < 1.0, f"pi mass {mass} not ~140 MeV"

    p_cls = TopologyClass(
        name="B_nucleon",
        blocks=(BuildingBlock.COORD_SHELL,),
        total_N=Fraction(13),
        residual=PointGroup.O_h,
        flags=frozenset({"shell", "winding_correction"}),
    )
    mass, N_eff = mass_prediction(p_cls, J=Fraction(1, 2))
    assert N_eff == Fraction(27, 2), f"Expected N=27/2, got {N_eff}"
    assert abs(mass - 945.3) < 1.0, f"p leading-order mass {mass} not ~945 MeV"

    print(f"Mass module self-test passed.")
    print(f"  pi mass (leading order): {mass_prediction(pi_cls, J=Fraction(0))[0]:.2f} MeV")
    print(f"  p mass (leading order):  {mass_prediction(p_cls, J=Fraction(1, 2))[0]:.2f} MeV")


if __name__ == "__main__":
    _self_test()
