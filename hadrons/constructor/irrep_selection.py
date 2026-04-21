"""
Stage B: irrep selection from (J, P, H).

Given angular momentum J, parity P, and residual point group H, this module
returns the target irrep(s) of H that the defect's lowest-lying eigenmode
must transform in.

Two-step process:
    1. Subduce D^{J,P} from SO(3) to O_h.
    2. Subduce from O_h to H (the residual group of the defect).

The lattice-specific fine point: the pseudoscalar J^P = 0^- is carried by the
A_{2u} irrep of O_h (basis function xyz), not the continuum-limit A_{1u}. This
is because the A_{1u} basis functions on the cubic lattice appear only at
degree 9 (too high for hadronic scale), while A_{2u}'s xyz is degree 3.

Reference: monograph Ch. 4 (spectral structure), line 2396 onwards.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from topology import PointGroup


# --------------------------------------------------------------------------
# O_h irreps
# --------------------------------------------------------------------------

class OhIrrep:
    A1g = "A_1g"
    A2g = "A_2g"
    Eg = "E_g"
    T1g = "T_1g"
    T2g = "T_2g"
    A1u = "A_1u"
    A2u = "A_2u"
    Eu = "E_u"
    T1u = "T_1u"
    T2u = "T_2u"
    E_half_g = "E_{1/2,g}"
    E_half_u = "E_{1/2,u}"
    F_32_g = "F_{3/2,g}"
    F_32_u = "F_{3/2,u}"
    E_52_g = "E_{5/2,g}"
    E_52_u = "E_{5/2,u}"


# --------------------------------------------------------------------------
# Step 1: SO(3) -> O_h subduction for (J, P)
# --------------------------------------------------------------------------

def subduce_so3_to_oh(J: Fraction, P: int, lattice_correction: bool = True) -> list[str]:
    """Decompose D^{(J, P)} of SO(3) into O_h irreps.

    If lattice_correction is True, pseudoscalars J=0, P=-1 return A_2u (the
    physically realised lattice pseudoscalar via the xyz harmonic), not A_1u.
    """
    parity_g = (P == +1)

    if J == 0:
        if parity_g:
            return [OhIrrep.A1g]
        else:
            if lattice_correction:
                return [OhIrrep.A2u]
            else:
                return [OhIrrep.A1u]

    if J == Fraction(1, 2):
        return [OhIrrep.E_half_g if parity_g else OhIrrep.E_half_u]

    if J == 1:
        return [OhIrrep.T1g if parity_g else OhIrrep.T1u]

    if J == Fraction(3, 2):
        if parity_g:
            return [OhIrrep.E_half_g, OhIrrep.F_32_g]
        else:
            return [OhIrrep.E_half_u, OhIrrep.F_32_u]

    if J == 2:
        if parity_g:
            return [OhIrrep.Eg, OhIrrep.T2g]
        else:
            return [OhIrrep.Eu, OhIrrep.T2u]

    if J == Fraction(5, 2):
        if parity_g:
            return [OhIrrep.E_52_g, OhIrrep.F_32_g]
        else:
            return [OhIrrep.E_52_u, OhIrrep.F_32_u]

    if J == 3:
        if parity_g:
            return [OhIrrep.A2g, OhIrrep.T1g, OhIrrep.T2g]
        else:
            return [OhIrrep.A2u, OhIrrep.T1u, OhIrrep.T2u]

    raise NotImplementedError(f"J={J} not yet supported")


# --------------------------------------------------------------------------
# Step 2: O_h -> H subduction
# --------------------------------------------------------------------------

OH_TO_TD: dict[str, list[str]] = {
    OhIrrep.A1g: ["A_1"],
    OhIrrep.A2g: ["A_2"],
    OhIrrep.Eg:  ["E"],
    OhIrrep.T1g: ["T_1"],
    OhIrrep.T2g: ["T_2"],
    OhIrrep.A1u: ["A_2"],
    OhIrrep.A2u: ["A_1"],
    OhIrrep.Eu:  ["E"],
    OhIrrep.T1u: ["T_2"],
    OhIrrep.T2u: ["T_1"],
    OhIrrep.E_half_g: ["E_{1/2}"],
    OhIrrep.E_half_u: ["E_{1/2}"],
    OhIrrep.F_32_g: ["F_{3/2}"],
    OhIrrep.F_32_u: ["F_{3/2}"],
}

OH_TO_D3D: dict[str, list[str]] = {
    OhIrrep.A1g: ["A_1g"],
    OhIrrep.A2g: ["A_2g"],
    OhIrrep.Eg:  ["E_g"],
    OhIrrep.T1g: ["A_2g", "E_g"],
    OhIrrep.T2g: ["A_1g", "E_g"],
    OhIrrep.A1u: ["A_1u"],
    OhIrrep.A2u: ["A_2u"],
    OhIrrep.Eu:  ["E_u"],
    OhIrrep.T1u: ["A_2u", "E_u"],
    OhIrrep.T2u: ["A_1u", "E_u"],
    OhIrrep.E_half_g: ["E_{1/2,g}"],
    OhIrrep.E_half_u: ["E_{1/2,u}"],
    OhIrrep.F_32_g: ["E_{1/2,g}", "E_{3/2,g}"],
    OhIrrep.F_32_u: ["E_{1/2,u}", "E_{3/2,u}"],
}

OH_TO_C3V: dict[str, list[str]] = {
    OhIrrep.A1g: ["A_1"],
    OhIrrep.A2g: ["A_2"],
    OhIrrep.Eg:  ["E"],
    OhIrrep.T1g: ["A_2", "E"],
    OhIrrep.T2g: ["A_1", "E"],
    OhIrrep.A1u: ["A_2"],
    OhIrrep.A2u: ["A_1"],
    OhIrrep.Eu:  ["E"],
    OhIrrep.T1u: ["A_1", "E"],
    OhIrrep.T2u: ["A_2", "E"],
    OhIrrep.E_half_g: ["E_{1/2}"],
    OhIrrep.E_half_u: ["E_{1/2}"],
    OhIrrep.F_32_g: ["E_{1/2}", "E_{3/2}"],
    OhIrrep.F_32_u: ["E_{1/2}", "E_{3/2}"],
}

OH_TO_D2H: dict[str, list[str]] = {
    OhIrrep.A1g: ["A_g"],
    OhIrrep.A2g: ["B_{1g}"],
    OhIrrep.Eg:  ["A_g", "B_{1g}"],
    OhIrrep.T1g: ["B_{1g}", "B_{2g}", "B_{3g}"],
    OhIrrep.T2g: ["A_g", "B_{2g}", "B_{3g}"],
    OhIrrep.A1u: ["A_u"],
    OhIrrep.A2u: ["B_{1u}"],
    OhIrrep.Eu:  ["A_u", "B_{1u}"],
    OhIrrep.T1u: ["B_{1u}", "B_{2u}", "B_{3u}"],
    OhIrrep.T2u: ["A_u", "B_{2u}", "B_{3u}"],
    OhIrrep.E_half_g: ["E_{1/2,g}"],
    OhIrrep.E_half_u: ["E_{1/2,u}"],
    OhIrrep.F_32_g: ["E_{1/2,g}", "E_{1/2,g}"],
    OhIrrep.F_32_u: ["E_{1/2,u}", "E_{1/2,u}"],
}

OH_TO_D3H: dict[str, list[str]] = {
    OhIrrep.A1g: ["A_1'"],
    OhIrrep.A2g: ["A_2'"],
    OhIrrep.Eg:  ["E'"],
    OhIrrep.T1g: ["A_2'", "E''"],
    OhIrrep.T2g: ["A_1'", "E''"],
    OhIrrep.A1u: ["A_1''"],
    OhIrrep.A2u: ["A_2''"],
    OhIrrep.Eu:  ["E''"],
    OhIrrep.T1u: ["A_2''", "E'"],
    OhIrrep.T2u: ["A_1''", "E'"],
    OhIrrep.E_half_g: ["E_{1/2}"],
    OhIrrep.E_half_u: ["E_{1/2}"],
    OhIrrep.F_32_g: ["E_{1/2}", "E_{3/2}"],
    OhIrrep.F_32_u: ["E_{1/2}", "E_{3/2}"],
}

OH_TO_C2V: dict[str, list[str]] = {
    OhIrrep.A1g: ["A_1"],
    OhIrrep.A2g: ["A_2"],
    OhIrrep.Eg:  ["A_1", "A_2"],
    OhIrrep.T1g: ["A_2", "B_1", "B_2"],
    OhIrrep.T2g: ["A_1", "B_1", "B_2"],
    OhIrrep.A1u: ["A_2"],
    OhIrrep.A2u: ["A_1"],
    OhIrrep.Eu:  ["A_1", "A_2"],
    OhIrrep.T1u: ["A_1", "B_1", "B_2"],
    OhIrrep.T2u: ["A_2", "B_1", "B_2"],
    OhIrrep.E_half_g: ["E_{1/2}"],
    OhIrrep.E_half_u: ["E_{1/2}"],
    OhIrrep.F_32_g: ["E_{1/2}", "E_{1/2}"],
    OhIrrep.F_32_u: ["E_{1/2}", "E_{1/2}"],
}


OH_TO_H: dict[PointGroup, dict[str, list[str]]] = {
    PointGroup.O_h: {k: [k] for k in [OhIrrep.A1g, OhIrrep.A2g, OhIrrep.Eg, OhIrrep.T1g, OhIrrep.T2g,
                                        OhIrrep.A1u, OhIrrep.A2u, OhIrrep.Eu, OhIrrep.T1u, OhIrrep.T2u,
                                        OhIrrep.E_half_g, OhIrrep.E_half_u, OhIrrep.F_32_g, OhIrrep.F_32_u]},
    PointGroup.T_d: OH_TO_TD,
    PointGroup.D_3d: OH_TO_D3D,
    PointGroup.C_3v: OH_TO_C3V,
    PointGroup.D_2h: OH_TO_D2H,
    PointGroup.D_3h: OH_TO_D3H,
    PointGroup.C_2v: OH_TO_C2V,
}


def subduce_oh_to_h(oh_irreps: list[str], H: PointGroup) -> list[str]:
    """Subduce O_h irreps to the residual subgroup H."""
    if H not in OH_TO_H:
        raise NotImplementedError(f"Subduction to H={H.value} not tabulated")
    out = []
    table = OH_TO_H[H]
    for gamma in oh_irreps:
        if gamma not in table:
            raise KeyError(f"O_h irrep {gamma} not in subduction table for H={H.value}")
        out.extend(table[gamma])
    return out


# --------------------------------------------------------------------------
# Full Stage B
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class IrrepSelection:
    """Result of Stage B: the irreps that the defect's lowest mode must carry."""
    oh_irreps: tuple[str, ...]
    h_irreps: tuple[str, ...]
    residual: PointGroup


def stage_b_irrep(J: Fraction, P: int, H: PointGroup) -> IrrepSelection:
    """Run Stage B: map (J, P, H) to the target irrep(s) of H.

    Returns an IrrepSelection listing the O_h irreps from SO(3) subduction
    and the H irreps after O_h -> H subduction.
    """
    oh_irreps = subduce_so3_to_oh(J, P)
    h_irreps = subduce_oh_to_h(oh_irreps, H)
    return IrrepSelection(
        oh_irreps=tuple(oh_irreps),
        h_irreps=tuple(h_irreps),
        residual=H,
    )


def _self_test() -> None:
    sel = stage_b_irrep(J=Fraction(0), P=-1, H=PointGroup.D_2h)
    assert OhIrrep.A2u in sel.oh_irreps
    assert "B_{1u}" in sel.h_irreps

    sel = stage_b_irrep(J=Fraction(1), P=-1, H=PointGroup.D_2h)
    assert OhIrrep.T1u in sel.oh_irreps
    assert all(x in sel.h_irreps for x in ["B_{1u}", "B_{2u}", "B_{3u}"])

    sel = stage_b_irrep(J=Fraction(1, 2), P=+1, H=PointGroup.O_h)
    assert sel.h_irreps == (OhIrrep.E_half_g,)

    sel = stage_b_irrep(J=Fraction(3, 2), P=+1, H=PointGroup.T_d)
    assert "F_{3/2}" in sel.h_irreps

    sel = stage_b_irrep(J=Fraction(0), P=-1, H=PointGroup.C_3v)
    assert "A_1" in sel.h_irreps

    print("Irrep selection module self-test passed.")


if __name__ == "__main__":
    _self_test()
