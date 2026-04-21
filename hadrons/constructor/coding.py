"""
Structural coding for hadrons in the Cosserat supersolid lattice framework.

The Standard Model labels hadrons by (B, S, I, I_3, J, P, n_c, n_b) plus derived
labels (Q, C, G-parity). The lattice framework classifies hadrons structurally
by the topology of the FCC defect graph.

This module provides translations between the two codings.

Structural coding tuple:
    (junction_count, S, n_c, n_b, J, P, pauli_flag, isosinglet_flag)

where:
    junction_count  = B       (number of Y-junctions)
    S               = |S|     (number of {111}-plane extensions)
    n_c             = n_c     (antibonding K_{9,9} excitations, one per charm quark)
    n_b             = n_b     (disclination monopoles, one per bottom quark)
    J               = J       (total angular momentum)
    P               = P       (parity eigenvalue)
    pauli_flag      = derived (whether tetrahedral voids are activated)
    isosinglet_flag = derived (whether the state is an SU(3) flavour singlet)

The Pauli and isosinglet flags are derived from I and the quark content, which in
turn is derived from (B, S, I, n_c, n_b). The flags are stored explicitly so
downstream code does not need to re-derive them.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Optional


# --------------------------------------------------------------------------
# Standard Model coding
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class SMLabels:
    """Standard Model quantum numbers for a hadron.

    B is baryon number (0 for mesons, 1 for baryons, -1 for antibaryons).
    S is strangeness (negative for quark strangeness by PDG convention, but
    this module uses |S| internally; the sign is carried only by `S_sign`).
    I is isospin.
    I_3 is the third component of isospin.
    J is total angular momentum.
    P is parity (+1 or -1).
    n_c is number of charm quarks minus antiquarks.
    n_b is number of bottom quarks minus antiquarks.
    S_sign is +1 for baryons containing strange antiquarks, -1 for strange quarks,
        0 if no strange content. Only relevant for sign conventions; does not
        affect the structural coding.
    """
    B: int
    S: int
    I: Fraction
    I_3: Fraction
    J: Fraction
    P: int
    n_c: int = 0
    n_b: int = 0
    S_sign: int = -1

    def __post_init__(self) -> None:
        if self.P not in (+1, -1):
            raise ValueError(f"P must be +1 or -1, got {self.P}")
        if self.B not in (-1, 0, 1, 2):
            raise ValueError(f"B must be in {{-1, 0, 1, 2}}, got {self.B}")
        if self.S < 0:
            raise ValueError(f"S stored as |S| >= 0, got {self.S}")

    @property
    def Q(self) -> Fraction:
        """Electric charge via Gell-Mann-Nishijima."""
        return self.I_3 + Fraction(self.B + self.S_sign * self.S + self.n_c - self.n_b, 2)


# --------------------------------------------------------------------------
# Structural coding
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class StructuralCoding:
    """Structural coding of a hadron in the FCC lattice framework.

    junction_count counts genuine Y-junctions where three stacking-fault
        partials meet. B=0 mesons have 0 junctions; B=1 baryons have 1.
    S is the number of {111}-plane extensions (strange quarks).
    n_c is the number of antibonding K_{9,9} excitations (charm quarks).
    n_b is the number of disclination monopoles (bottom quarks).
    J is total angular momentum.
    P is parity.
    pauli_flag is True when identical light quarks require spatial antisymmetry
        via tetrahedral void activation or additional hex-cap extension.
    isosinglet_flag is True when the state is an SU(3) flavour singlet (e.g.,
        eta, eta', omega, phi with specific quark-mixing structure).
    """
    junction_count: int
    S: int
    n_c: int
    n_b: int
    J: Fraction
    P: int
    pauli_flag: bool
    isosinglet_flag: bool

    def __post_init__(self) -> None:
        if self.P not in (+1, -1):
            raise ValueError(f"P must be +1 or -1, got {self.P}")
        if self.junction_count < 0:
            raise ValueError(f"junction_count must be >= 0, got {self.junction_count}")


# --------------------------------------------------------------------------
# Quark-content derivation
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class QuarkContent:
    """Quark content of a hadron: counts of each flavour.

    For baryons, the counts are quark counts (antibaryons would have antiquark
    counts). For mesons, the counts are of a quark-antiquark pair, so only one
    flavour of each sign is present at a time in the simplest case.
    """
    n_u: int
    n_d: int
    n_s: int
    n_c: int
    n_b: int

    def n_light(self) -> int:
        """Number of light quarks (u + d)."""
        return self.n_u + self.n_d

    def n_identical_light(self) -> int:
        """Number of identical light quarks in the defect.

        Returns 2 if there are two identical light quarks (uu or dd), 3 if
        three (uuu or ddd), 0 otherwise. Used for Pauli filter.
        """
        if self.n_u >= 2 and self.n_d == 0:
            return self.n_u
        if self.n_d >= 2 and self.n_u == 0:
            return self.n_d
        if self.n_u >= 1 and self.n_d >= 1:
            # ud, uud, udd, uuud, etc: at most one pair of identical quarks
            if self.n_u == 2 or self.n_d == 2:
                return 2
            return 0
        return 0


def quark_content_from_sm(sm: SMLabels) -> QuarkContent:
    """Derive the quark content from SM labels.

    This works for the canonical ground-state hadrons. For baryons (B=1)
    with given (S, I, I_3, n_c, n_b), the quark content is determined up to
    convention; for mesons (B=0) it is determined up to the flavour-singlet
    mixing ambiguity.

    Returns a QuarkContent with quark-antiquark counts summed:
        n_u = (number of u quarks) + (number of u antiquarks)
    This is appropriate for structural coding, which cares about how many
    partial dislocations of each flavour are present in the defect graph,
    not their sign.
    """
    if sm.B == 1:
        n_heavy = sm.n_c + sm.n_b
        n_strange = sm.S
        n_light = 3 - n_strange - n_heavy
        if n_light < 0:
            raise ValueError(f"Cannot construct baryon: {n_heavy} heavy + {n_strange} strange > 3 quarks")
        n_u = Fraction(n_light, 2) + sm.I_3
        n_d = Fraction(n_light, 2) - sm.I_3
        if n_u.denominator != 1 or n_d.denominator != 1:
            raise ValueError(f"Non-integer quark content for SM labels {sm}")
        return QuarkContent(
            n_u=int(n_u),
            n_d=int(n_d),
            n_s=n_strange,
            n_c=sm.n_c,
            n_b=sm.n_b,
        )
    elif sm.B == 0:
        n_heavy = sm.n_c + sm.n_b
        if n_heavy == 2:
            return QuarkContent(n_u=0, n_d=0, n_s=0, n_c=sm.n_c, n_b=sm.n_b)
        elif n_heavy == 1 and sm.S == 1:
            n_u = Fraction(1, 2) + sm.I_3
            n_d = Fraction(1, 2) - sm.I_3
            return QuarkContent(
                n_u=int(n_u),
                n_d=int(n_d),
                n_s=1,
                n_c=sm.n_c,
                n_b=sm.n_b,
            )
        elif n_heavy == 1 and sm.S == 0:
            n_u = Fraction(1, 2) + sm.I_3
            n_d = Fraction(1, 2) - sm.I_3
            return QuarkContent(
                n_u=int(n_u),
                n_d=int(n_d),
                n_s=0,
                n_c=sm.n_c,
                n_b=sm.n_b,
            )
        elif sm.S == 2:
            return QuarkContent(n_u=0, n_d=0, n_s=2, n_c=0, n_b=0)
        elif sm.S == 1:
            n_u_light = Fraction(1, 2) + sm.I_3
            n_d_light = Fraction(1, 2) - sm.I_3
            return QuarkContent(
                n_u=int(n_u_light),
                n_d=int(n_d_light),
                n_s=1,
                n_c=0,
                n_b=0,
            )
        else:
            n_u = Fraction(2, 2) if sm.I_3 > 0 else (Fraction(0) if sm.I_3 < 0 else Fraction(1))
            n_d = Fraction(2, 2) - n_u + Fraction(1) if sm.I_3 == 0 else Fraction(2) - n_u
            return QuarkContent(
                n_u=int(n_u),
                n_d=int(n_d),
                n_s=0,
                n_c=0,
                n_b=0,
            )
    else:
        raise NotImplementedError(f"B={sm.B} not supported")


def is_isosinglet(sm: SMLabels) -> bool:
    """Check whether a state is an SU(3) flavour singlet.

    For mesons, this means I=0 with specific quark-antiquark mixing
    (e.g., eta, eta', omega, phi). For baryons, this means the flavour
    wavefunction is in the SU(3) singlet representation (rare for ground
    states; the Lambda is I=0 but flavour-octet, not singlet).
    """
    if sm.B == 0 and sm.I == 0 and sm.S == 0 and sm.n_c == 0 and sm.n_b == 0:
        return True
    if sm.B == 0 and sm.I == 0 and sm.S == 2 and sm.n_c == 0 and sm.n_b == 0:
        return True
    return False


# --------------------------------------------------------------------------
# Pauli filter
# --------------------------------------------------------------------------

def pauli_activates(sm: SMLabels, qc: Optional[QuarkContent] = None) -> bool:
    """Determine whether Pauli exclusion requires void activation or spatial extension.

    Rule (monograph Sec. 10.8):
        pauli activates <=> n_identical_light >= 2 AND (I >= 1 OR J >= 3/2)

    Strange quarks do not count towards n_identical_light because their
    {111}-plane extensions spatially separate them automatically.
    Two strange quarks (as in the Xi*) still trigger activation via the
    J >= 3/2 branch of the rule.
    """
    if qc is None:
        qc = quark_content_from_sm(sm)

    n_id_light = qc.n_identical_light()
    if n_id_light >= 2 and (sm.I >= 1 or sm.J >= Fraction(3, 2)):
        return True

    if sm.B == 1 and qc.n_s >= 2 and sm.J >= Fraction(3, 2):
        return True

    return False


# --------------------------------------------------------------------------
# Translation: SM -> structural
# --------------------------------------------------------------------------

def sm_to_structural(sm: SMLabels) -> StructuralCoding:
    """Translate SM quantum numbers to structural coding."""
    qc = quark_content_from_sm(sm)
    return StructuralCoding(
        junction_count=sm.B,
        S=sm.S,
        n_c=sm.n_c,
        n_b=sm.n_b,
        J=sm.J,
        P=sm.P,
        pauli_flag=pauli_activates(sm, qc),
        isosinglet_flag=is_isosinglet(sm),
    )


# --------------------------------------------------------------------------
# Translation: structural -> SM (partial: I_3 not recoverable)
# --------------------------------------------------------------------------

def structural_to_sm_family(sc: StructuralCoding) -> list[SMLabels]:
    """Translate structural coding to the isospin multiplet of SM labels.

    The structural coding does not distinguish isospin partners (e.g., pi+, pi0, pi-
    all map to the same structural code). This function returns the full isospin
    multiplet: the list of (I, I_3) pairs consistent with the structural inputs.
    """
    if sc.junction_count == 0:
        n_light_quark_pairs = 2 - sc.S - sc.n_c - sc.n_b
        if sc.isosinglet_flag:
            I_candidates = [Fraction(0)]
        elif n_light_quark_pairs >= 2:
            if sc.S == 0 and sc.n_c == 0 and sc.n_b == 0:
                I_candidates = [Fraction(1)]
            else:
                I_candidates = [Fraction(1, 2)]
        elif n_light_quark_pairs == 1:
            I_candidates = [Fraction(1, 2)]
        else:
            I_candidates = [Fraction(0)]
    elif sc.junction_count == 1:
        n_light = 3 - sc.S - sc.n_c - sc.n_b
        if n_light < 0:
            return []
        if sc.pauli_flag and n_light >= 2:
            I_candidates = [Fraction(n_light, 2)]
        elif n_light == 3:
            I_candidates = [Fraction(1, 2), Fraction(3, 2)]
        elif n_light == 2:
            I_candidates = [Fraction(0), Fraction(1)]
        elif n_light == 1:
            I_candidates = [Fraction(1, 2)]
        else:
            I_candidates = [Fraction(0)]
    else:
        raise NotImplementedError(f"junction_count={sc.junction_count} not supported")

    labels = []
    for I in I_candidates:
        num_steps = int(2 * I) + 1
        for step in range(num_steps):
            I_3 = I - step
            labels.append(SMLabels(
                B=sc.junction_count,
                S=sc.S,
                I=I,
                I_3=I_3,
                J=sc.J,
                P=sc.P,
                n_c=sc.n_c,
                n_b=sc.n_b,
            ))
    return labels


# --------------------------------------------------------------------------
# Quick roundtrip validation
# --------------------------------------------------------------------------

def _self_test() -> None:
    """Basic validation of the translation layer with a few known hadrons."""
    pi_plus = SMLabels(B=0, S=0, I=Fraction(1), I_3=Fraction(1), J=Fraction(0), P=-1)
    sc = sm_to_structural(pi_plus)
    assert sc.junction_count == 0
    assert sc.S == 0
    assert sc.J == 0
    assert sc.P == -1
    assert not sc.pauli_flag

    proton = SMLabels(B=1, S=0, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(1, 2), P=+1)
    sc = sm_to_structural(proton)
    assert sc.junction_count == 1
    assert sc.S == 0
    assert sc.J == Fraction(1, 2)
    assert not sc.pauli_flag

    delta_plusplus = SMLabels(B=1, S=0, I=Fraction(3, 2), I_3=Fraction(3, 2), J=Fraction(3, 2), P=+1)
    sc = sm_to_structural(delta_plusplus)
    assert sc.pauli_flag, "Delta should have Pauli activation"

    sigma_plus = SMLabels(B=1, S=1, I=Fraction(1), I_3=Fraction(1), J=Fraction(1, 2), P=+1)
    sc = sm_to_structural(sigma_plus)
    assert sc.pauli_flag, "Sigma should have Pauli activation"

    lambda0 = SMLabels(B=1, S=1, I=Fraction(0), I_3=Fraction(0), J=Fraction(1, 2), P=+1)
    sc = sm_to_structural(lambda0)
    assert not sc.pauli_flag, "Lambda should not have Pauli activation"

    xi_star = SMLabels(B=1, S=2, I=Fraction(1, 2), I_3=Fraction(1, 2), J=Fraction(3, 2), P=+1)
    sc = sm_to_structural(xi_star)
    assert sc.pauli_flag, "Xi* should trigger Pauli via J>=3/2 with two strange quarks"

    print("Coding module self-test passed.")


if __name__ == "__main__":
    _self_test()
