# D4 Interpretation

Scripts supporting the four-dimensional (D4) interpretation of the Cosserat
supersolid lattice.

The D4 lattice (densest 4D packing, 24-cell nearest-neighbour shell) has one
compact direction of circumference L₄ = 3ℓ ≈ 8.5 fm identified as Euclidean
time. Every quantitative prediction is identical in FCC (3+1) and D4 (4+0)
below the Kaluza–Klein scale.

## Scripts

| Script | Purpose |
|--------|---------|
| `d4_lattice_sums.py` | Elastic constants: verifies A = 1 (exact isotropy), S₂₂₂ = 0, and cross-index invariance for both 12-vector (FCC) and 24-vector (24-cell) shells. |
| `d4_mu_prime.py` | Acoustoelastic check: verifies μ' = 2 identically in FCC and D4 after KK reduction. |
| `d4_scales.py` | Temperature hierarchy and the exact identity G × K_sf = c⁴/ℓ². |
| `d4_chirality.py` | **20-branch spectrum**, chirality splitting, and the imaginary phase of the compact-direction self-energy. |

## Key results

### Elastic isotropy (d4_lattice_sums.py)
The Zener ratio A = 1 exactly in D4, versus A ≠ 1 in FCC. The inter-layer
bonds stiffen the ⟨100⟩ directions without affecting the shear modulus.

### Chirality splitting (d4_chirality.py)
The 4D Cosserat theory has 20 phonon branches (not 12), because SO(4) has
6 rotation generators. At k₄ = 0, the 20 branches factorise exactly into
System A (the familiar 12) and System B (8 massive KK modes).

At k₄ ≠ 0, the two EM polarisations **split**: the self-dual component
tunnels 1.2% more easily than the anti-self-dual at the first KK level.
The origin is the sign correlation between the curl coupling (opposite for
the two polarisations) and the compact-direction mixing (same sign for both).

The compact-rotation component of the EM eigenvector has phase exactly π/2
— the coupling Σ₁₃ is **purely imaginary**, making it time-reversal-odd.
This is geometric CP violation from the stacking chirality A→B→C.

The effect exists only in the thermal window 23 MeV < T < 156 MeV (between
the geometric temperature and the deconfinement transition). It is active
during hadronisation in heavy-ion collisions and in neutron star mergers.
