# Dynamics

One-loop scattering and radiative correction computations from the Cosserat lattice Feynman rules.

## Scripts

### `ew_one_loop_lattice.py`

Computes the complete one-loop electroweak programme from lattice first principles:

- **Leptonic vacuum polarisation** — exact one-loop integral with Koide-derived lepton masses. Matches the SM value to three significant figures (0.0314 vs 0.0315).

- **Hadronic vacuum polarisation** — once-subtracted dispersion integral over the lattice vector meson spectral function (ρ, ω, φ, J/ψ, Υ(1S)), with perturbative QCD at Λ_QCD = π·m₀ above the quark-hadron duality matching point. Zero free parameters. Result: Δα_had = 0.0273 ± 0.0003 against the data-driven 0.02766 ± 0.00010. The radial vector tower widths are derived, not fitted: the Cornell potential built from the framework's own string tension and α_s, solved by radial diagonalisation, reproduces both measured bottomonium radial ratios to 1–3% where the pure-linear limit was off by factors of 1.9 and 2.4. The uncertainty is now dominated by the duality matching point and the string tension.

- **Custodial-breaking ρ parameter** — UV-convergent integral from the top-bottom mass splitting in the evanescent (W) self-energy. Result: Δρ = 0.0092, matching the SM value ~0.0094.

- **Running coupling** — α⁻¹(M_Z) = 128.98 ± 0.04 against the measured 128.943 ± 0.014, an agreement at 0.9σ.

- **Corrected G_F** — Δr = Δα − (c²/s²)Δρ = 0.0265 reduces the tree-level residual from 3.2% to 0.6%.

- **Effective Weinberg angle** — tree-level residual reduced from 4.0% to 1.7%.

The script carries five analytic checks that run before any physics is reported: the
dispersion integral against its closed form for constant R, the Breit-Wigner spectral
area against the arctangent, the narrow-resonance pole formula against area times weight,
the KSRF self-consistency of the pion form factor (the Breit-Wigner assembled from
the KSRF coupling returns the tree-level KSRF leptonic width, 4.91 keV, to five digits),
and the relation |ψ(0)|² = (μ/2π)⟨V'⟩ against hydrogen, where both sides are known.

### Usage

```
python3 ew_one_loop_lattice.py
```

Requires: `numpy`, `scipy`.

### References

Monograph Chapter 10 (Dynamics), Sections on the two-loop QED programme, Landau pole elimination, and the electroweak radiative correction programme.
