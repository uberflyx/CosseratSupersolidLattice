"""
Stage A: topology-class enumeration from Burgers-vector closure rules.

Given a structural coding tuple, this module determines:
    - The topology class (which building blocks compose the defect)
    - The residual point group H ⊆ O_h
    - The vertex count N and edge count |E| of the composed graph

The algorithm runs deterministically from the structural coding alone. No
particle names or mass lookups are consulted.

Building blocks:
    CELL_PAIR           N=2,  H=D_2h (dipole along <110>)
    HEX_CAP             N=7,  H=C_3v (one {111} stacking-fault cluster)
    BILAYER             N=8,  H=C_3v (hex cap + one adjacent layer node)
    COORD_SHELL         N=13, H=O_h (central node + 12 NN)
    VOID_ACTIVATED      +4,   H=T_d (four tetrahedral voids of one inversion class)
    HEX_CAP_EXTENSION   +3,   H=C_3v (three new nodes beyond the shell on a {111})
    CROSSED_FAULT       N=11, H=D_2h (two hex caps sharing 2 vertices and a centre)
    TRIPLE_BILAYER      N=24, H=D_3h (three bilayers on three {111} planes)

Residual symmetry hierarchy:
    O_h > T_d > D_3d > D_3h > C_3v
    O_h > D_4h > D_2h > C_2v > C_2

Heavy-quark insertions preserve the light-sector residual symmetry; they add
nodes inside the shell without breaking the point group.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from enum import Enum
from typing import Optional

from coding import StructuralCoding


# --------------------------------------------------------------------------
# Building blocks
# --------------------------------------------------------------------------

class BuildingBlock(Enum):
    CELL_PAIR = "cell_pair"
    HEX_CAP = "hex_cap"
    BILAYER = "bilayer"
    COORD_SHELL = "coord_shell"
    VOID_ACTIVATION = "void_activation"
    HEX_CAP_EXTENSION = "hex_cap_extension"
    CROSSED_FAULT = "crossed_fault"
    TRIPLE_BILAYER = "triple_bilayer"
    CHARM_BLOCK = "charm_block"
    BOTTOM_BLOCK = "bottom_block"
    RIBBON_PLUS_ONE = "ribbon_plus_one"
    STRANGE_CASIMIR = "strange_casimir"


NODES_PER_BLOCK: dict[BuildingBlock, Fraction] = {
    BuildingBlock.CELL_PAIR: Fraction(2),
    BuildingBlock.HEX_CAP: Fraction(7),
    BuildingBlock.BILAYER: Fraction(8),
    BuildingBlock.COORD_SHELL: Fraction(13),
    BuildingBlock.VOID_ACTIVATION: Fraction(4),
    BuildingBlock.HEX_CAP_EXTENSION: Fraction(3),
    BuildingBlock.CROSSED_FAULT: Fraction(11),
    BuildingBlock.TRIPLE_BILAYER: Fraction(24),
    BuildingBlock.CHARM_BLOCK: Fraction(18),
    BuildingBlock.BOTTOM_BLOCK: Fraction(66),
    BuildingBlock.RIBBON_PLUS_ONE: Fraction(1),
    BuildingBlock.STRANGE_CASIMIR: Fraction(4, 3),
}


# --------------------------------------------------------------------------
# Residual point groups
# --------------------------------------------------------------------------

class PointGroup(Enum):
    O_h = "O_h"
    T_d = "T_d"
    D_3d = "D_3d"
    D_3h = "D_3h"
    C_3v = "C_3v"
    D_4h = "D_4h"
    D_2h = "D_2h"
    C_2v = "C_2v"
    C_s = "C_s"
    C_1 = "C_1"


# --------------------------------------------------------------------------
# Topology class
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class TopologyClass:
    """A topology class with its building-block composition.

    blocks is an ordered list of building blocks used to construct the defect.
    total_N is the vertex count of the composed graph (may be a Fraction for
        states with transcendental or fractional quark contributions such as
        the strange Casimir 4/3 or bottom 6pi^2).
    residual is the residual point group.
    flags records derived properties for downstream use.
    """
    name: str
    blocks: tuple[BuildingBlock, ...]
    total_N: Fraction
    residual: PointGroup
    flags: frozenset[str]


def compute_total_N(blocks: tuple[BuildingBlock, ...]) -> Fraction:
    """Apply the intersection lemma to compute the net vertex count.

    Rules:
        - A single block contributes its NODES_PER_BLOCK value.
        - Two hex caps share exactly 2 ring sites plus 1 centre, so the second
          hex cap contributes 7 - 3 = 4 new nodes (6 ring - 2 shared - already
          counted centre).
        - Three hex caps fill the shell: total = 13 (not 7 + 4 + 4 = 15).
        - A hex cap extension adds 3 NNN nodes beyond the shell.
        - Void activation adds 4 tetrahedral-void nodes.
        - Charm and bottom blocks add their full node count (internal to shell).
        - Triple bilayer: explicit 3 * 8 = 24 (three bilayers on three planes,
          no sharing at the Y-junction level).
        - Strange Casimir adds 4/3 (for heavy-meson charmed-strange states).
    """
    if not blocks:
        return Fraction(0)

    n = Fraction(0)
    hex_caps_used = 0
    has_shell = False

    for block in blocks:
        if block == BuildingBlock.CELL_PAIR:
            n += Fraction(2)
        elif block == BuildingBlock.HEX_CAP:
            if hex_caps_used == 0:
                n += Fraction(7)
            elif hex_caps_used == 1:
                n += Fraction(4)
            elif hex_caps_used == 2:
                n += Fraction(2)
            else:
                raise ValueError("More than 3 hex caps not supported")
            hex_caps_used += 1
        elif block == BuildingBlock.BILAYER:
            n += Fraction(8)
        elif block == BuildingBlock.COORD_SHELL:
            if has_shell:
                raise ValueError("Cannot add a second coord shell")
            n += Fraction(13)
            has_shell = True
        elif block == BuildingBlock.VOID_ACTIVATION:
            n += Fraction(4)
        elif block == BuildingBlock.HEX_CAP_EXTENSION:
            n += Fraction(3)
        elif block == BuildingBlock.CROSSED_FAULT:
            n += Fraction(11)
        elif block == BuildingBlock.TRIPLE_BILAYER:
            n += Fraction(24)
        elif block == BuildingBlock.CHARM_BLOCK:
            n += Fraction(18)
        elif block == BuildingBlock.BOTTOM_BLOCK:
            n += Fraction(66)
        elif block == BuildingBlock.RIBBON_PLUS_ONE:
            n += Fraction(1)
        elif block == BuildingBlock.STRANGE_CASIMIR:
            n += Fraction(4, 3)
        else:
            raise ValueError(f"Unknown block {block}")

    return n


# --------------------------------------------------------------------------
# Stage A: topology class from structural coding
# --------------------------------------------------------------------------

def stage_a_topology_class(sc: StructuralCoding) -> Optional[TopologyClass]:
    """Run Stage A: derive the topology class from structural coding.

    Returns None if no valid class exists (forbidden tuple).
    """
    B = sc.junction_count
    S = sc.S
    J = sc.J
    P = sc.P
    n_c = sc.n_c
    n_b = sc.n_b

    if B not in (0, 1):
        return None

    if S + n_c + n_b > 3 if B == 1 else S + n_c + n_b > 2:
        return None

    if B == 1 and J.denominator == 1:
        return None
    if B == 0 and J.denominator != 1:
        return None

    if B == 0:
        cls = _meson_topology(sc)
    elif B == 1:
        cls = _baryon_topology(sc)
    else:
        return None

    if cls is None:
        return None

    if n_c > 0 or n_b > 0:
        cls = _insert_heavy(cls, n_c, n_b)

    return cls


def _meson_topology(sc: StructuralCoding) -> Optional[TopologyClass]:
    """Topology class for a meson (B=0, junction_count=0).

    Selection:
        J=0, P=-1 (pseudoscalar) -> cell pair / hex cap / bilayer / shell+1
        J=1, P=-1 (vector)       -> crossed fault / pinned shell / bilayer+axial
        J=0, P=+1 (scalar)       -> tensor channel, Born cluster
        J=1, P=+1 (axial)        -> tensor channel, Born cluster
        J=2, P=+1 (tensor f_2)   -> Born cluster

    Heavy mesons use N = N_quarks + N_ribbon(J^P).
    Isosinglet vs isovector distinguishes flavour-mixing states.
    """
    S = sc.S
    J = sc.J
    P = sc.P
    n_heavy = sc.n_c + sc.n_b
    isosinglet = sc.isosinglet_flag

    is_pseudoscalar = (J == 0 and P == -1)
    is_vector = (J == 1 and P == -1)
    is_scalar = (J == 0 and P == +1)
    is_axial = (J == 1 and P == +1)
    is_tensor = (J == 2 and P == +1)

    if not (is_pseudoscalar or is_vector or is_scalar or is_axial or is_tensor):
        return None

    if n_heavy == 0:
        if is_pseudoscalar and S == 0 and not isosinglet:
            return TopologyClass(
                name="M_pseudoscalar_light",
                blocks=(BuildingBlock.CELL_PAIR,),
                total_N=Fraction(2),
                residual=PointGroup.D_2h,
                flags=frozenset({"cell_pair", "translational_sector"}),
            )
        if is_pseudoscalar and S == 1:
            return TopologyClass(
                name="M_pseudoscalar_strange",
                blocks=(BuildingBlock.HEX_CAP,),
                total_N=Fraction(7),
                residual=PointGroup.C_3v,
                flags=frozenset({"hex_cap", "translational_sector"}),
            )
        if is_pseudoscalar and S == 0 and isosinglet:
            return TopologyClass(
                name="M_pseudoscalar_flavour_singlet",
                blocks=(BuildingBlock.BILAYER,),
                total_N=Fraction(8),
                residual=PointGroup.D_3d,
                flags=frozenset({"bilayer", "translational_sector", "flavour_singlet"}),
            )
        if is_pseudoscalar and S == 2 and isosinglet:
            return TopologyClass(
                name="M_pseudoscalar_strange_singlet",
                blocks=(BuildingBlock.COORD_SHELL, BuildingBlock.RIBBON_PLUS_ONE),
                total_N=Fraction(14),
                residual=PointGroup.C_3v,
                flags=frozenset({"shell_plus_one", "translational_sector", "flavour_singlet"}),
            )

        if is_vector and S == 0 and not isosinglet:
            return TopologyClass(
                name="M_vector_light",
                blocks=(BuildingBlock.CROSSED_FAULT,),
                total_N=Fraction(11),
                residual=PointGroup.D_2h,
                flags=frozenset({"crossed_fault", "translational_sector", "delocalised_spin"}),
            )
        if is_vector and S == 0 and isosinglet:
            return TopologyClass(
                name="M_vector_light_singlet",
                blocks=(BuildingBlock.CROSSED_FAULT,),
                total_N=Fraction(11),
                residual=PointGroup.D_2h,
                flags=frozenset({"crossed_fault", "translational_sector", "flavour_singlet"}),
            )
        if is_vector and S == 1:
            return TopologyClass(
                name="M_vector_strange",
                blocks=(BuildingBlock.COORD_SHELL,),
                total_N=Fraction(13),
                residual=PointGroup.C_3v,
                flags=frozenset({"pinned_shell", "translational_sector"}),
            )
        if is_vector and S == 2 and isosinglet:
            return TopologyClass(
                name="M_vector_strange_singlet",
                blocks=(BuildingBlock.COORD_SHELL, BuildingBlock.RIBBON_PLUS_ONE),
                total_N=Fraction(14),
                residual=PointGroup.D_3d,
                flags=frozenset({"pinned_bilayer_axial", "translational_sector", "flavour_singlet"}),
            )

        if is_scalar or is_axial or is_tensor:
            if isosinglet and is_tensor:
                return TopologyClass(
                    name="M_tensor_flavour_singlet",
                    blocks=(BuildingBlock.COORD_SHELL, BuildingBlock.HEX_CAP_EXTENSION, BuildingBlock.VOID_ACTIVATION),
                    total_N=Fraction(18),
                    residual=PointGroup.O_h,
                    flags=frozenset({"born_cluster", "rotational_sector", "flavour_singlet"}),
                )
            if isosinglet and is_axial:
                return TopologyClass(
                    name="M_axial_flavour_singlet",
                    blocks=(BuildingBlock.COORD_SHELL, BuildingBlock.HEX_CAP_EXTENSION, BuildingBlock.VOID_ACTIVATION),
                    total_N=Fraction(18),
                    residual=PointGroup.O_h,
                    flags=frozenset({"born_cluster", "rotational_sector", "flavour_singlet"}),
                )
            ribbon_blocks = (BuildingBlock.COORD_SHELL,) if is_scalar else (BuildingBlock.COORD_SHELL, BuildingBlock.RIBBON_PLUS_ONE)
            ribbon_N = 13 if is_scalar else 14
            return TopologyClass(
                name=f"M_{['scalar','axial','tensor'][is_axial + 2*is_tensor]}_light",
                blocks=ribbon_blocks,
                total_N=ribbon_N,
                residual=PointGroup.O_h if is_scalar else PointGroup.C_3v,
                flags=frozenset({"born_cluster", "rotational_sector"}),
            )
    else:
        if is_pseudoscalar:
            ribbon = (BuildingBlock.HEX_CAP,)
            ribbon_N = Fraction(7)
            residual = PointGroup.C_3v
        elif is_vector:
            ribbon = (BuildingBlock.BILAYER,)
            ribbon_N = Fraction(8)
            residual = PointGroup.C_3v
        elif is_scalar:
            ribbon = (BuildingBlock.COORD_SHELL,)
            ribbon_N = Fraction(13)
            residual = PointGroup.O_h
        elif is_axial or is_tensor:
            ribbon = (BuildingBlock.COORD_SHELL, BuildingBlock.RIBBON_PLUS_ONE)
            ribbon_N = Fraction(14)
            residual = PointGroup.C_3v
        else:
            return None

        base_flags = {"heavy_ribbon", "translational_sector" if P == -1 else "rotational_sector"}
        blocks_out = ribbon
        total_N = ribbon_N
        if S == 1:
            blocks_out = ribbon + (BuildingBlock.STRANGE_CASIMIR,)
            total_N = ribbon_N + Fraction(4, 3)
            base_flags.add("strange_casimir")

        return TopologyClass(
            name=f"M_heavy_{['pseudoscalar','vector','scalar','axial','tensor'][is_vector + 2*is_scalar + 3*is_axial + 4*is_tensor]}",
            blocks=blocks_out,
            total_N=total_N,
            residual=residual,
            flags=frozenset(base_flags),
        )

    return None


def _baryon_topology(sc: StructuralCoding) -> Optional[TopologyClass]:
    """Topology class for a baryon (B=1, junction_count=1).

    Selection is driven by:
        - S (number of strange arms)
        - J (1/2 or 3/2 for ground-state decuplet)
        - pauli_flag (whether voids or extensions are needed)
    """
    S = sc.S
    J = sc.J
    pauli = sc.pauli_flag

    if S == 0:
        if J == Fraction(1, 2) and not pauli:
            return TopologyClass(
                name="B_nucleon",
                blocks=(BuildingBlock.COORD_SHELL,),
                total_N=Fraction(13),
                residual=PointGroup.O_h,
                flags=frozenset({"shell", "winding_correction"}),
            )
        if J == Fraction(1, 2) and pauli:
            return TopologyClass(
                name="B_sigma_like",
                blocks=(BuildingBlock.COORD_SHELL, BuildingBlock.VOID_ACTIVATION),
                total_N=Fraction(17),
                residual=PointGroup.T_d,
                flags=frozenset({"shell", "void_activated", "winding_correction"}),
            )
        if J == Fraction(3, 2) and pauli:
            return TopologyClass(
                name="B_delta",
                blocks=(BuildingBlock.COORD_SHELL, BuildingBlock.VOID_ACTIVATION),
                total_N=Fraction(17),
                residual=PointGroup.T_d,
                flags=frozenset({"shell", "void_activated", "winding_correction"}),
            )
        return None

    if S == 1:
        if J == Fraction(1, 2) and not pauli:
            return TopologyClass(
                name="B_lambda",
                blocks=(BuildingBlock.COORD_SHELL, BuildingBlock.HEX_CAP_EXTENSION),
                total_N=Fraction(16),
                residual=PointGroup.C_3v,
                flags=frozenset({"shell", "hex_cap_extension"}),
            )
        if J == Fraction(1, 2) and pauli:
            return TopologyClass(
                name="B_sigma",
                blocks=(BuildingBlock.COORD_SHELL, BuildingBlock.VOID_ACTIVATION),
                total_N=Fraction(17),
                residual=PointGroup.C_3v,
                flags=frozenset({"shell", "void_activated"}),
            )
        if J == Fraction(3, 2) and pauli:
            return TopologyClass(
                name="B_sigma_star",
                blocks=(BuildingBlock.COORD_SHELL, BuildingBlock.VOID_ACTIVATION, BuildingBlock.HEX_CAP_EXTENSION),
                total_N=Fraction(20),
                residual=PointGroup.C_3v,
                flags=frozenset({"shell", "void_activated", "hex_cap_extension"}),
            )
        return None

    if S == 2:
        if J == Fraction(1, 2) and not pauli:
            return TopologyClass(
                name="B_xi",
                blocks=(BuildingBlock.COORD_SHELL, BuildingBlock.HEX_CAP_EXTENSION, BuildingBlock.HEX_CAP_EXTENSION),
                total_N=Fraction(19),
                residual=PointGroup.C_2v,
                flags=frozenset({"shell", "double_hex_cap"}),
            )
        if J == Fraction(3, 2) and pauli:
            return TopologyClass(
                name="B_xi_star",
                blocks=(BuildingBlock.COORD_SHELL, BuildingBlock.HEX_CAP_EXTENSION, BuildingBlock.HEX_CAP_EXTENSION, BuildingBlock.HEX_CAP_EXTENSION),
                total_N=Fraction(22),
                residual=PointGroup.C_3v,
                flags=frozenset({"shell", "triple_hex_cap"}),
            )
        return None

    if S == 3:
        if J == Fraction(3, 2):
            return TopologyClass(
                name="B_omega",
                blocks=(BuildingBlock.TRIPLE_BILAYER,),
                total_N=Fraction(24),
                residual=PointGroup.D_3h,
                flags=frozenset({"triple_bilayer"}),
            )
        return None

    return None


def _insert_heavy(cls: TopologyClass, n_c: int, n_b: int) -> TopologyClass:
    """Insert heavy-quark blocks into an existing topology class.

    Each charm quark adds an 18-node K_{9,9} antibonding block inside the shell.
    Each bottom quark adds a 66-node disclination monopole block on a {111} plane
    (for baryons; for mesons this should use the transcendental 6pi^2 value, but
    we use the integer approximation 66 here consistently).

    The residual symmetry of the light sector is preserved.
    """
    extra_blocks = (BuildingBlock.CHARM_BLOCK,) * n_c + (BuildingBlock.BOTTOM_BLOCK,) * n_b
    new_blocks = cls.blocks + extra_blocks
    new_N = cls.total_N + n_c * Fraction(18) + n_b * Fraction(66)
    new_name = cls.name + f"_heavy_c{n_c}b{n_b}"
    return TopologyClass(
        name=new_name,
        blocks=new_blocks,
        total_N=new_N,
        residual=cls.residual,
        flags=cls.flags | frozenset({f"heavy_c{n_c}_b{n_b}"}),
    )


# --------------------------------------------------------------------------
# Self-test
# --------------------------------------------------------------------------

def _self_test() -> None:
    from coding import sm_to_structural, SMLabels

    pi = SMLabels(B=0, S=0, I=Fraction(1), I_3=Fraction(0), J=Fraction(0), P=-1)
    cls = stage_a_topology_class(sm_to_structural(pi))
    assert cls.total_N == Fraction(2), f"pi should have N=2, got {cls.total_N}"
    assert cls.residual == PointGroup.D_2h

    k = SMLabels(B=0, S=1, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(0), P=-1)
    cls = stage_a_topology_class(sm_to_structural(k))
    assert cls.total_N == Fraction(7), f"K should have N=7, got {cls.total_N}"

    rho = SMLabels(B=0, S=0, I=Fraction(1), I_3=Fraction(0), J=Fraction(1), P=-1)
    cls = stage_a_topology_class(sm_to_structural(rho))
    assert cls.total_N == Fraction(11), f"rho should have N=11, got {cls.total_N}"

    kstar = SMLabels(B=0, S=1, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(1), P=-1)
    cls = stage_a_topology_class(sm_to_structural(kstar))
    assert cls.total_N == Fraction(13), f"K* should have N=13, got {cls.total_N}"

    p = SMLabels(B=1, S=0, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(1, 2), P=+1)
    cls = stage_a_topology_class(sm_to_structural(p))
    assert cls.total_N == Fraction(13), f"p should have N=13, got {cls.total_N}"
    assert cls.residual == PointGroup.O_h

    delta = SMLabels(B=1, S=0, I=Fraction(3, 2), I_3=Fraction(3, 2), J=Fraction(3, 2), P=+1)
    cls = stage_a_topology_class(sm_to_structural(delta))
    assert cls.total_N == Fraction(17), f"Delta should have N=17, got {cls.total_N}"
    assert cls.residual == PointGroup.T_d

    lam = SMLabels(B=1, S=1, I=Fraction(0), I_3=Fraction(0), J=Fraction(1, 2), P=+1)
    cls = stage_a_topology_class(sm_to_structural(lam))
    assert cls.total_N == Fraction(16), f"Lambda should have N=16, got {cls.total_N}"

    sig = SMLabels(B=1, S=1, I=Fraction(1), I_3=Fraction(1), J=Fraction(1, 2), P=+1)
    cls = stage_a_topology_class(sm_to_structural(sig))
    assert cls.total_N == Fraction(17), f"Sigma should have N=17, got {cls.total_N}"

    xi = SMLabels(B=1, S=2, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(1, 2), P=+1)
    cls = stage_a_topology_class(sm_to_structural(xi))
    assert cls.total_N == Fraction(19), f"Xi should have N=19, got {cls.total_N}"

    omega = SMLabels(B=1, S=3, I=Fraction(0), I_3=Fraction(0), J=Fraction(3, 2), P=+1)
    cls = stage_a_topology_class(sm_to_structural(omega))
    assert cls.total_N == Fraction(24), f"Omega should have N=24, got {cls.total_N}"

    lambda_c = SMLabels(B=1, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(1, 2), P=+1, n_c=1)
    cls = stage_a_topology_class(sm_to_structural(lambda_c))
    assert cls.total_N == Fraction(31), f"Lambda_c should have N=31, got {cls.total_N}"

    lambda_b = SMLabels(B=1, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(1, 2), P=+1, n_b=1)
    cls = stage_a_topology_class(sm_to_structural(lambda_b))
    assert cls.total_N == Fraction(79), f"Lambda_b should have N=79, got {cls.total_N}"

    d0 = SMLabels(B=0, S=0, I=Fraction(1, 2), I_3=Fraction(-1, 2), J=Fraction(0), P=-1, n_c=1)
    cls = stage_a_topology_class(sm_to_structural(d0))
    assert cls.total_N == Fraction(25), f"D0 should have N=25, got {cls.total_N}"

    ds = SMLabels(B=0, S=1, I=Fraction(0), I_3=Fraction(0), J=Fraction(0), P=-1, n_c=1)
    cls = stage_a_topology_class(sm_to_structural(ds))
    assert cls.total_N == Fraction(79, 3), f"D_s should have N=25+4/3=79/3, got {cls.total_N}"

    eta_c = SMLabels(B=0, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(0), P=-1, n_c=2)
    cls = stage_a_topology_class(sm_to_structural(eta_c))
    assert cls.total_N == Fraction(43), f"eta_c should have N=43, got {cls.total_N}"

    jpsi = SMLabels(B=0, S=0, I=Fraction(0), I_3=Fraction(0), J=Fraction(1), P=-1, n_c=2)
    cls = stage_a_topology_class(sm_to_structural(jpsi))
    assert cls.total_N == Fraction(44), f"J/psi should have N=44, got {cls.total_N}"

    print("Topology module self-test passed.")


if __name__ == "__main__":
    _self_test()
