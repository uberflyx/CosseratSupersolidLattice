"""
Hadron catalogue generator.

Enumerates structural-coding tuples (B, S, I, J, P, n_c, n_b) over a specified
range and, for each tuple, runs Stages A and B to produce the topology class,
residual symmetry, target irrep, and predicted mass.

The output is a table of predicted hadrons that can be compared against the
PDG to test the framework.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Optional
import sys

from coding import SMLabels, sm_to_structural, StructuralCoding
from topology import stage_a_topology_class, TopologyClass, PointGroup
from irrep_selection import stage_b_irrep, IrrepSelection
from mass import mass_prediction


@dataclass
class CatalogueEntry:
    """One row of the predicted hadron catalogue."""
    B: int
    S: int
    I: Fraction
    I_3: Fraction
    J: Fraction
    P: int
    n_c: int
    n_b: int
    pauli_flag: bool
    isosinglet_flag: bool
    topology_name: str
    N_integer: Fraction
    N_effective: Fraction
    residual: str
    oh_irreps: tuple[str, ...]
    h_irreps: tuple[str, ...]
    mass_leading_MeV: float
    building_blocks: tuple[str, ...]
    flags: tuple[str, ...]

    def to_row(self) -> list:
        def fmt_frac(f: Fraction) -> str:
            if f.denominator == 1:
                return str(f.numerator)
            return f"{f.numerator}/{f.denominator}"

        return [
            self.B,
            self.S,
            fmt_frac(self.I),
            fmt_frac(self.I_3),
            fmt_frac(self.J),
            "+" if self.P == +1 else "-",
            self.n_c,
            self.n_b,
            "Y" if self.pauli_flag else "N",
            "Y" if self.isosinglet_flag else "N",
            self.topology_name,
            fmt_frac(self.N_integer),
            fmt_frac(self.N_effective),
            self.residual,
            ",".join(self.oh_irreps),
            ",".join(self.h_irreps),
            f"{self.mass_leading_MeV:.2f}",
            ",".join(self.building_blocks),
            ",".join(sorted(self.flags)),
        ]


HEADERS = [
    "B", "S", "I", "I_3", "J", "P", "n_c", "n_b",
    "pauli", "isosinglet",
    "topology", "N_int", "N_eff", "residual",
    "O_h_irreps", "H_irreps", "mass_MeV (leading)",
    "blocks", "flags",
]


def generate_catalogue(
    B_range: tuple[int, ...] = (0, 1),
    S_range: tuple[int, ...] = (0, 1, 2, 3),
    J_mesons: tuple[Fraction, ...] = (Fraction(0), Fraction(1), Fraction(2)),
    J_baryons: tuple[Fraction, ...] = (Fraction(1, 2), Fraction(3, 2)),
    P_range: tuple[int, ...] = (-1, +1),
    n_c_range: tuple[int, ...] = (0, 1, 2, 3),
    n_b_range: tuple[int, ...] = (0, 1, 2),
    include_isosinglet_candidates: bool = True,
) -> list[CatalogueEntry]:
    """Generate the catalogue of predicted hadrons.

    Iterates over the Cartesian product of the input ranges. For each tuple,
    attempts Stage A; if it succeeds, runs Stage B, applies the mass formula,
    and records the result. Skips forbidden tuples.
    """
    catalogue: list[CatalogueEntry] = []
    seen_topologies: set = set()

    for B in B_range:
        J_range = J_baryons if B == 1 else J_mesons
        for S in S_range:
            for n_c in n_c_range:
                for n_b in n_b_range:
                    if B == 1 and S + n_c + n_b > 3:
                        continue
                    if B == 0 and S + n_c + n_b > 2 and not (n_c == 2 or n_b == 2):
                        continue

                    if B == 1:
                        n_light = 3 - S - n_c - n_b
                        if n_light < 0:
                            continue
                        if n_light == 0:
                            I_candidates = [Fraction(0)]
                        elif n_light == 1:
                            I_candidates = [Fraction(1, 2)]
                        elif n_light == 2:
                            I_candidates = [Fraction(0), Fraction(1)]
                        else:
                            I_candidates = [Fraction(1, 2), Fraction(3, 2)]
                    elif B == 0:
                        if n_c + n_b == 2:
                            I_candidates = [Fraction(0)]
                        elif n_c + n_b == 1 and S == 1:
                            I_candidates = [Fraction(0)]
                        elif n_c + n_b == 1 and S == 0:
                            I_candidates = [Fraction(1, 2)]
                        elif S == 2:
                            I_candidates = [Fraction(0)]
                            if include_isosinglet_candidates:
                                pass
                        elif S == 1:
                            I_candidates = [Fraction(1, 2)]
                        else:
                            I_candidates = [Fraction(0), Fraction(1)]
                    else:
                        continue

                    for I in I_candidates:
                        I_3_values = [I - k for k in range(int(2 * I) + 1)]
                        I_3 = I_3_values[0]

                        for J in J_range:
                            for P in P_range:
                                if B == 0 and not ((J == 0 and P == -1) or (J == 1 and P == -1)
                                                   or (J == 0 and P == +1) or (J == 1 and P == +1)
                                                   or (J == 2 and P == +1)):
                                    continue
                                if B == 1 and P != +1:
                                    continue

                                try:
                                    sm = SMLabels(B=B, S=S, I=I, I_3=I_3, J=J, P=P, n_c=n_c, n_b=n_b)
                                    sc = sm_to_structural(sm)
                                except (ValueError, NotImplementedError):
                                    continue

                                cls = stage_a_topology_class(sc)
                                if cls is None:
                                    continue

                                try:
                                    sel = stage_b_irrep(J=J, P=P, H=cls.residual)
                                except (KeyError, NotImplementedError):
                                    sel = IrrepSelection(oh_irreps=(), h_irreps=(), residual=cls.residual)

                                mass, N_eff = mass_prediction(cls, J=J)

                                entry = CatalogueEntry(
                                    B=B, S=S, I=I, I_3=I_3, J=J, P=P, n_c=n_c, n_b=n_b,
                                    pauli_flag=sc.pauli_flag,
                                    isosinglet_flag=sc.isosinglet_flag,
                                    topology_name=cls.name,
                                    N_integer=cls.total_N,
                                    N_effective=N_eff,
                                    residual=cls.residual.value,
                                    oh_irreps=sel.oh_irreps,
                                    h_irreps=sel.h_irreps,
                                    mass_leading_MeV=mass,
                                    building_blocks=tuple(b.value for b in cls.blocks),
                                    flags=tuple(cls.flags),
                                )

                                key = (B, S, I, J, P, n_c, n_b, sc.isosinglet_flag)
                                if key not in seen_topologies:
                                    catalogue.append(entry)
                                    seen_topologies.add(key)

                        if B == 0 and I == Fraction(0) and include_isosinglet_candidates:
                            for J in J_range:
                                for P in P_range:
                                    if not ((J == 0 and P == -1) or (J == 1 and P == -1)):
                                        continue
                                    try:
                                        sc_alt = StructuralCoding(
                                            junction_count=0,
                                            S=S,
                                            n_c=n_c,
                                            n_b=n_b,
                                            J=J,
                                            P=P,
                                            pauli_flag=False,
                                            isosinglet_flag=True,
                                        )
                                    except ValueError:
                                        continue
                                    cls = stage_a_topology_class(sc_alt)
                                    if cls is None:
                                        continue
                                    key = (B, S, I, J, P, n_c, n_b, True)
                                    if key in seen_topologies:
                                        continue
                                    try:
                                        sel = stage_b_irrep(J=J, P=P, H=cls.residual)
                                    except (KeyError, NotImplementedError):
                                        sel = IrrepSelection(oh_irreps=(), h_irreps=(), residual=cls.residual)
                                    mass, N_eff = mass_prediction(cls, J=J)
                                    entry = CatalogueEntry(
                                        B=B, S=S, I=I, I_3=Fraction(0), J=J, P=P, n_c=n_c, n_b=n_b,
                                        pauli_flag=False, isosinglet_flag=True,
                                        topology_name=cls.name,
                                        N_integer=cls.total_N,
                                        N_effective=N_eff,
                                        residual=cls.residual.value,
                                        oh_irreps=sel.oh_irreps,
                                        h_irreps=sel.h_irreps,
                                        mass_leading_MeV=mass,
                                        building_blocks=tuple(b.value for b in cls.blocks),
                                        flags=tuple(cls.flags),
                                    )
                                    catalogue.append(entry)
                                    seen_topologies.add(key)

    catalogue.sort(key=lambda e: (e.B, e.S, e.n_c, e.n_b, e.J, e.P, e.I))
    return catalogue


def write_catalogue_csv(catalogue: list[CatalogueEntry], path: str) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(HEADERS)
        for entry in catalogue:
            writer.writerow(entry.to_row())


def write_catalogue_markdown(catalogue: list[CatalogueEntry], path: str) -> None:
    with open(path, "w") as f:
        f.write("# Predicted hadron catalogue\n\n")
        f.write(f"Total entries: {len(catalogue)}\n\n")

        light_mesons = [e for e in catalogue if e.B == 0 and e.n_c == 0 and e.n_b == 0]
        light_baryons = [e for e in catalogue if e.B == 1 and e.n_c == 0 and e.n_b == 0]
        charm_mesons = [e for e in catalogue if e.B == 0 and e.n_c >= 1 and e.n_b == 0]
        charm_baryons = [e for e in catalogue if e.B == 1 and e.n_c >= 1 and e.n_b == 0]
        bottom_mesons = [e for e in catalogue if e.B == 0 and e.n_b >= 1]
        bottom_baryons = [e for e in catalogue if e.B == 1 and e.n_b >= 1]

        sections = [
            ("Light mesons", light_mesons),
            ("Light baryons", light_baryons),
            ("Charm mesons", charm_mesons),
            ("Charm baryons", charm_baryons),
            ("Bottom mesons", bottom_mesons),
            ("Bottom baryons", bottom_baryons),
        ]

        for title, entries in sections:
            if not entries:
                continue
            f.write(f"\n## {title} ({len(entries)} entries)\n\n")
            f.write("| " + " | ".join(HEADERS) + " |\n")
            f.write("| " + " | ".join("---" for _ in HEADERS) + " |\n")
            for entry in entries:
                row = entry.to_row()
                f.write("| " + " | ".join(str(x) for x in row) + " |\n")


def compare_with_pdg_light_sector(catalogue: list[CatalogueEntry]) -> str:
    """Produce a short readable summary of the catalogue."""
    lines = []
    lines.append(f"Catalogue size: {len(catalogue)} distinct structural tuples.")
    lines.append("")

    light_mesons = [e for e in catalogue if e.B == 0 and e.n_c == 0 and e.n_b == 0]
    light_baryons = [e for e in catalogue if e.B == 1 and e.n_c == 0 and e.n_b == 0]
    charm_mesons = [e for e in catalogue if e.B == 0 and e.n_c >= 1 and e.n_b == 0]
    charm_baryons = [e for e in catalogue if e.B == 1 and e.n_c >= 1 and e.n_b == 0]
    bottom_mesons = [e for e in catalogue if e.B == 0 and e.n_b >= 1]
    bottom_baryons = [e for e in catalogue if e.B == 1 and e.n_b >= 1]

    lines.append(f"Light mesons: {len(light_mesons)}")
    lines.append(f"Light baryons: {len(light_baryons)}")
    lines.append(f"Charm mesons: {len(charm_mesons)}")
    lines.append(f"Charm baryons: {len(charm_baryons)}")
    lines.append(f"Bottom mesons: {len(bottom_mesons)}")
    lines.append(f"Bottom baryons: {len(bottom_baryons)}")
    lines.append("")

    lines.append("Light mesons (leading-order mass):")
    for e in sorted(light_mesons, key=lambda e: e.mass_leading_MeV):
        j_str = f"{e.J.numerator}" if e.J.denominator == 1 else f"{e.J.numerator}/{e.J.denominator}"
        N_int_str = str(e.N_integer.numerator) if e.N_integer.denominator == 1 else f"{e.N_integer.numerator}/{e.N_integer.denominator}"
        lines.append(f"  ({e.B},{e.S},{e.I},{j_str}^{'+' if e.P == +1 else '-'}) {'[sing]' if e.isosinglet_flag else '     '}  "
                    f"N={N_int_str:>5} ({e.N_effective})  H={e.residual:<6}  {e.mass_leading_MeV:7.1f} MeV  "
                    f"{e.topology_name}")
    lines.append("")

    lines.append("Light baryons (leading-order mass):")
    for e in sorted(light_baryons, key=lambda e: e.mass_leading_MeV):
        j_str = f"{e.J.numerator}/{e.J.denominator}"
        N_int_str = str(e.N_integer.numerator) if e.N_integer.denominator == 1 else f"{e.N_integer.numerator}/{e.N_integer.denominator}"
        lines.append(f"  ({e.B},{e.S},{e.I},{j_str}^+)  "
                    f"N={N_int_str:>5} ({e.N_effective})  H={e.residual:<6}  {e.mass_leading_MeV:7.1f} MeV  "
                    f"{e.topology_name}")
    lines.append("")

    return "\n".join(lines)


if __name__ == "__main__":
    catalogue = generate_catalogue()
    print(compare_with_pdg_light_sector(catalogue))
    write_catalogue_csv(catalogue, "./catalogue.csv")
    write_catalogue_markdown(catalogue, "./catalogue.md")
    print()
    print(f"Full catalogue written to catalogue.csv and catalogue.md")
