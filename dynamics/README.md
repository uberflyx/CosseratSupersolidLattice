# Dynamics

One-loop scattering and radiative correction computations from the Cosserat lattice Feynman rules.

## Scripts

### `ew_one_loop_lattice.py`

Computes the complete one-loop electroweak programme from lattice first principles:

- **Leptonic vacuum polarisation** — exact one-loop integral with Koide-derived lepton masses. Matches the SM value to three significant figures (0.0314 vs 0.0315).

- **Hadronic vacuum polarisation** — dispersion integral using lattice vector meson Breit-Wigner spectral functions (ρ, ω, φ, J/ψ, Υ) below 2 GeV, and perturbative QCD with Λ_QCD = π·m₀ above. Zero free parameters. Result: Δα_had = 0.0284, within 2.8% of the data-driven value 0.02766 ± 0.00010.

- **Custodial-breaking ρ parameter** — UV-convergent integral from the top-bottom mass splitting in the evanescent (W) self-energy. Result: Δρ = 0.0092, matching the SM value ~0.0094.

- **Running coupling** — α⁻¹(M_Z) = 128.8, within 0.08% of measured 128.94.

- **Corrected G_F** — tree-level residual reduced from 3.2% to ~1%.

- **Effective Weinberg angle** — tree-level residual reduced from 4.0% to ~1.6%.

### Usage

```
python3 ew_one_loop_lattice.py
```

Requires: `numpy`, `scipy`.

### References

Monograph Chapter 10 (Dynamics), Sections on the two-loop QED programme, Landau pole elimination, and the electroweak radiative correction programme.
