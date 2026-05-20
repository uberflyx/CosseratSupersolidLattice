# Spectral mass construction

Hadron masses from the elastic spectrum of an FCC defect cluster.  A
cluster is built from FCC geometry, its Cosserat dynamical matrix is
diagonalised, the mass mode is identified by its symmetry, and the mass
follows from one law:

```
m = N * m_0  -  N * (4 - lambda) * m_e
```

`N` is the integer node count of the cluster, `lambda` is the
dimensionless eigenvalue of the mass mode (in units of the bond
stiffness), `m_0 = m_e / alpha` is the node mass, and `m_e` is the
electron mass.  At the structural value `lambda = 4` the mass is exactly
`N m_0`; the elastic correction moves it off that line.

This replaces the earlier node-counting-plus-`Q` phenomenology.  Nothing
here is fitted: every cluster, every eigenvalue, and every irrep label
comes from the lattice geometry before any mass is consulted.

## How a mass is obtained

1. Build the cluster from FCC sites (a coordination shell, a bilayer
   extension, a second-shell octahedron, a dual-orbit void set, and so
   on).
2. Assemble the Cosserat dynamical matrix on those sites
   (`cosserat_classifier.build_cosserat_matrix`).
3. Diagonalise and decompose each eigenspace under the cluster's point
   group, a subgroup of `O_h`, using the character formula.
4. Pick the mass mode by the parity-flip rule: its spatial irrep has
   parity opposite to the hadron's spinor parity.
5. Read the mass from the eigenvalue.

## Core modules

- `spectral_classifier.py` -- FCC cluster builders, the displacement-only
  dynamical matrix, the `O_h` character table, and the irrep
  decomposition.  Also holds the Lambda ground-state builder.
- `cosserat_classifier.py` -- the full Cosserat (displacement plus
  microrotation) dynamical matrix `build_cosserat_matrix`.
- `proton_first_principles.py` -- the proton's `A_2u` mass mode, the
  `lambda = 8.303` reference for `J^P = 1/2+`.
- `delta_first_principles.py` -- Delta cluster builders and the
  two-displacement-field matrix used by several scripts.
- `n1535_first_principles.py` -- the `N(1535)` bilayer cluster builder.
- `hadron_spectral_mass.py` -- the mass law `m = N m_0 - N(4-lambda) m_e`
  and shared helpers.

## Reference spectrum and rules

- `bare_shell_irreps.py` -- the bare 13-node shell spectrum, every
  eigenvalue decomposed by `O_h` irrep.  This is the lookup table the
  other scripts draw on.  Stiff phi-dominated references: `A_2u` at
  8.303 (`1/2+`), `T_1u` at 8.054 (`3/2+`), `A_2g` at 6.193 (`1/2-`),
  `T_1g` at 8.071 (`3/2-`).
- `parity_flip_rule.py` -- verifies the parity-flip rule mapping each
  `J^P` to its parent irrep on the bare shell.

## Excited and resonance baryons

- `retrospective_adiabatic.py` -- the adiabatic audit.  Each excited
  state's mass mode is the adiabatic continuation of its parent's mass
  mode as the cluster is built.  Covers `N(1535)`, the Roper `N(1440)`,
  and `Lambda(1600)`.
- `delta1600_dual_orbit.py` -- the `Delta(1600)` dual-orbit cluster
  (N = 21), its `O_h` decomposition, and the `T_1u` mass mode at
  lambda = 10.216, m = 1537 MeV.
- `delta_adiabatic.py` -- the adiabatic interpolation from `Delta(1232)`
  to `Delta(1600)`, tracking the `T_1u` mode by eigenvector overlap.
- `backward_trace.py` -- traces bilayer cluster modes back to the
  bare-shell parent irrep they descend from.
- `delta1600_irrep_check.py` -- confirms the `Delta(1600)` mass mode is
  `T_1u`, not `T_1g`.

## Usage

Run any script from inside this directory; the modules import their
siblings by name.

```
cd spectral_mass
python bare_shell_irreps.py        # the bare-shell reference spectrum
python parity_flip_rule.py         # the J^P -> parent irrep map
python retrospective_adiabatic.py  # the excited-baryon audit
python delta1600_dual_orbit.py     # the Delta(1600) closure
```

## Three radial-excitation mechanisms

The excited baryons fall into three structural families, each a
different way of adding nodes to the 13-shell:

1. Bilayer extension (`N(1535)`, `N(1520)`): 13 + 8 = 21, `C_s`, a
   parity-flipping mode from a gerade bare-shell irrep.
2. Second-shell octahedron (Roper `N(1440)`, `Lambda(1600)`): 13 + 6,
   `O_h` preserved.
3. Dual-orbit voids (`Delta(1600)`): 13 + 8 = 21, `O_h` restored from
   the `T_d` of the ground-state Delta.
