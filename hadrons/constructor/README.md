# Hadron constructor

A deterministic pipeline that takes Standard Model quantum numbers `(B, S, I, I_3, J, P, n_c, n_b)` and produces the predicted hadron mass through a staged algorithm with no particle-name lookups.

**Status: working.** Stages A–C (topology, irrep selection, leading-order mass) are complete and validated against 43 ground-state hadrons from the PDG. Stage D (electromagnetic `Q` correction) is a prototype.

## Pipeline

The constructor runs four stages in sequence:

1. **Stage A (`topology.py`)** selects a building-block composition from the structural coding. Given `(B, S, J, P, n_c, n_b, pauli_flag, isosinglet_flag)`, it returns a `TopologyClass` specifying which blocks compose the defect (cell pair, hex cap, coordination shell, voids, extensions, charm block, bottom block, strange Casimir) and its residual point group `H ⊆ O_h`.

2. **Stage B (`irrep_selection.py`)** subduces the spin-parity `D^{(J, P)}` of SO(3) down through O_h to `H`, returning the target irrep of the defect's lowest-lying mode. The subduction uses the lattice-specific `A_{2u}` assignment for pseudoscalars (via the `xyz` basis function) rather than the continuum-limit `A_{1u}`.

3. **Stage C (`mass.py`)** applies `m = N * m_0 + Q * m_e`. At leading order it computes `N * m_0` with the appropriate fractional-`N` correction: winding-mode `+1/2` for baryons, delocalised-spin `+1/13` for the `rho`-family vectors, pinned-bilayer-axial `+1/2` for the `phi`-like strange singlets.

4. **Stage D (`fcc_graph.py`, prototype)** is the electromagnetic correction `Q`. Currently implemented as a simplicial-complex face-enumeration prototype that produces correct `Q` for the proton and the pion directly from vertex-type pairs and activated faces, but requires physics-derived dispatches for other baryons. A first-principles Cayley-graph circulation-integral formulation (eliminating the dispatches) is in development.

## Files

| File | Purpose |
|---|---|
| `coding.py` | SM ↔ structural coding translation. `SMLabels` and `StructuralCoding` dataclasses, the `sm_to_structural()` and `structural_to_sm_family()` maps, Pauli-void and isosinglet flag derivation. |
| `topology.py` | Stage A. Building-block enum, intersection lemma for shared vertices (two hex caps share 2 ring sites plus a centre), residual-symmetry assignment, heavy-quark block insertion. |
| `irrep_selection.py` | Stage B. SO(3) → O_h subduction with lattice pseudoscalar correction, plus O_h → `H` subduction tables for `T_d`, `D_3d`, `C_3v`, `D_4h`, `D_2h`, `D_3h`, `C_2v`. |
| `mass.py` | Stage C. `leading_order_mass()` and `mass_prediction()` functions with flag-driven spin corrections. |
| `generate_catalogue.py` | Iterates over structural-coding tuples and emits a full catalogue of predicted hadrons as CSV and Markdown. 93 distinct entries across light, charm, and bottom sectors. |
| `pdg_comparison.py` | Compares the leading-order catalogue against 43 observed PDG hadrons. |
| `fcc_graph.py` | Stage D prototype. Simplicial complex builder with vertex-type and face-activation tagging. |

## Usage

Run any module directly to execute its self-tests:

```bash
python coding.py                # Translation layer self-test
python topology.py              # Stage A self-test on 14 hadrons
python irrep_selection.py       # Stage B subduction self-test
python mass.py                  # Stage C mass-formula self-test
python generate_catalogue.py    # Produce the full catalogue (CSV + Markdown)
python pdg_comparison.py        # Run the 43-hadron PDG benchmark
python fcc_graph.py             # Stage D prototype on 6 ground states
```

To use as a library:

```python
from fractions import Fraction
from coding import SMLabels, sm_to_structural
from topology import stage_a_topology_class
from irrep_selection import stage_b_irrep
from mass import mass_prediction

sm = SMLabels(B=1, S=0, I=Fraction(1, 2), I_3=Fraction(1, 2),
              J=Fraction(1, 2), P=+1)
sc = sm_to_structural(sm)
cls = stage_a_topology_class(sc)
sel = stage_b_irrep(J=sm.J, P=sm.P, H=cls.residual)
mass, N_eff = mass_prediction(cls, J=sm.J)
print(f"{cls.name}: N={cls.total_N}, H={cls.residual.value}, m={mass:.2f} MeV")
```

## PDG benchmark

At leading order (no electromagnetic correction):

* 20 of 43 hadrons within 1% of observed mass (47%).
* 38 of 43 within 5% (88%).
* Mean absolute residual 2.08%, maximum 9.28%.
* All `N` values match the monograph's derivations exactly.

Residuals are dominated by the missing `Q * m_e` correction (order `α ≈ 0.7%` of total mass for light hadrons, up to a few percent for heavy mesons). With the full `Q`-algorithm the residuals are expected to drop to the monograph's sub-percent level for ground states.

## Stage D status

The `Q`-algorithm is under active development. Three formulations have been explored:

1. The monograph's four-term decomposition `Q = Q_bond + Q_col + Q_surf + Q_iso`. Diagnostic and physically motivated but requires particle-by-particle dispatching.

2. A simplicial-complex face-enumeration formulation (implemented here as a prototype). Faces of the defect's 2-complex are tagged by activation state; each activated face contributes `N_c^(dim - 1)` to `Q` via its colour multiplicity. Handles proton and pion directly from geometry; other baryons require isospin-dependent void rules and plane-coupling rules that currently must be supplied.

3. A Cayley-graph circulation-integral formulation (in development). `Q` as the circulation of the Cosserat microrotation field `phi(v)` around the defect's boundary. Expected to generalise automatically to any hadron without dispatching, by having Burgers-vector allocations at Y-junction vertices encode all the structural information currently supplied separately. The equivalence with the four-term decomposition is the subject of the monograph's reformulation theorem.

## Related files in the repository

* `../verify_fcc_geometry.py` — FCC coordinate checks used by the constructor.
* `../excited_baryons.py` — excited-state mass calculations, complementary to the ground-state catalogue.
* `../exotic_catalogue.py` — exotic-hadron classification; the constructor is the ground-state counterpart.
* `../pdg_comparison.md` — narrative comparison report; this subdirectory's `pdg_comparison.py` is the automated harness.
* `../../cosserat_graph.py` — the full constructive calculator referenced throughout the monograph.
