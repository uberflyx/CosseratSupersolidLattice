# D4 Interpretation

Scripts supporting the four-dimensional (D4) interpretation of the Cosserat
supersolid lattice (Sec. 2.8 of the monograph).

The D4 interpretation treats the vacuum as a four-dimensional FCC crystal (the
D4 root lattice / 24-cell) with one compact direction of circumference
L₄ = 3ℓ ≈ 8.5 fm identified as Euclidean time. Every quantitative prediction
is identical in D3 (3+1) and D4 (4+0); the scripts verify this and compute
the structural advantages of D4.

## Scripts

| Script | Purpose |
|--------|---------|
| `d4_lattice_sums.py` | Computes all SOEC/TOEC lattice sums for the 12-vector (FCC) and 24-vector (24-cell) shells. Verifies A = 1, S₂₂₂ = 0, and the cross-index invariance. |
| `d4_mu_prime.py` | Verifies μ' = (3−ξ)/5 = 2.000 identically in D3 and D4 after KK reduction. Resolves the apparent discrepancy from the naive 4D formula. |
| `d4_scales.py` | Computes the temperature hierarchy (T_geom, m_KK, T_c, Λ_QCD, T_D) and verifies the exact identity G × K_sf = c⁴/ℓ². |

## Key results

- **Exact isotropy**: The Zener ratio A = 2C₄₄/(C₁₁−C₁₂) = 1 in D4, versus
  A = 2 in D3. The inter-layer bonds stiffen the ⟨100⟩ directions without
  affecting the shear modulus.

- **μ' = 2 in both D3 and D4**: The inter-layer bonds contribute to C₁₁ but
  not C₄₄ (each has only one nonzero 3D component, so R₁²R₂² = 0). The ratio
  μ_V/(3K) = 1/5 is invariant under KK reduction.

- **S₂₂₂ = 0 in both D3 and D4**: The C₄₅₆ = 0 selection rule that protects
  the α derivation holds in D4.

- **G × K_sf = c⁴/ℓ²**: Newton's constant times the superfluid bulk modulus
  is an exact algebraic identity of the bootstrap equations.

## References

- Conway & Sloane, *Sphere Packings, Lattices and Groups* (Springer, 1999):
  D4 root lattice and the 24-cell.
- Kaluza, Sitzungsber. Preuss. Akad. Wiss. (1921) 966; Klein, Z. Phys. 37,
  895 (1926): Kaluza–Klein compactification.
- Gross, Pisarski & Yaffe, Rev. Mod. Phys. 53, 43 (1981): deconfinement as
  a KK phase transition.
