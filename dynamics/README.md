# Dynamics

One-loop scattering and radiative correction computations from the Cosserat lattice Feynman rules.

## Scripts

### `ew_one_loop_lattice.py`

Computes the complete one-loop electroweak programme from lattice first principles:

- **Leptonic vacuum polarisation** — exact one-loop integral with Koide-derived lepton masses. Matches the SM value to three significant figures (0.0314 vs 0.0315).

- **Hadronic vacuum polarisation** — once-subtracted dispersion integral over the lattice vector meson spectral function (ρ, ω, φ, J/ψ, Υ(1S)), with perturbative QCD at Λ_QCD = π·m₀ above the quark-hadron duality matching point. Zero free parameters. Result: Δα_had = 0.0273 ± 0.0005 against the data-driven 0.02766 ± 0.00010. The radial vector tower widths are derived, not fitted: the Cornell potential built from the framework's own string tension and α_s, solved by radial diagonalisation, reproduces both measured bottomonium radial ratios to 1–3% where the pure-linear limit was off by factors of 1.9 and 2.4. All three vector families are built explicitly on their own Regge trajectories, the continuum carries the standard QCD series through c₃ (verified against the PDG QCD review, eq. 9.9), and the matching point is quantised by the tower's own spacing: s_match(N) = m_ρ² + (N+½)·2πσ with N = 1..4 admissible. The quoted uncertainty is the spread over those four points, the spectral-completion band, plus the input budget.

- **Custodial-breaking ρ parameter** — UV-convergent integral from the top-bottom mass splitting in the evanescent (W) self-energy. Result: Δρ = 0.0092, matching the SM value ~0.0094.

- **Running coupling** — α⁻¹(M_Z) = 128.98 ± 0.07 against the measured 128.943 ± 0.014, an agreement at 0.5σ.

- **Corrected G_F** — Δr = Δα − (c²/s²)Δρ = 0.0265 reduces the tree-level residual from 3.2% to 0.6%.

- **Effective Weinberg angle** — tree-level residual reduced from 4.0% to 1.7%.

- **Muon anomaly** — a_mu^HVP = (733 ± 10) × 10⁻¹⁰ against WP25's 713 ± 6, 1.7σ. The ρ leptonic width carries a derived finite-pion-mass correction: both Weinberg sum rules are chiral-limit statements, and giving the pion its mass moves its pole off the origin, adding f_π²m_π² to the first moment and multiplying the coupling by (1 − m_π²/m_a1²) = −1.08%. The factor is small because the pole enters the *first* moment, so it is weighed against m_a1 rather than m_ρ.

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

### `rho_sector_audit.py`

Audits the rho sector, which supplies the single most consequential number in the
hadronic programme. Establishes what is in good order (mass 0.05%, KSRF width 0.23%,
KSRF coupling tested against the measured width 0.15%) and localises what is not:

- **Gamma_ee(rho) is +4.6%**, and is the outlier among the three vector leptonic widths
  the framework computes (omega +0.5%, phi -1.1%); no single rescaling fits rho and phi
  together, so the defect is specific to the rho's coupling to the photon.
- **The two routes to that width are not independent.** They agree only through the
  numerical near-identity 8 pi^3 ~ 2 (Z^2/(Z+1))^2, true to 1.1%, which makes
  m_a1^2 = 3 m_rho^2 and the Weinberg factor 3/2. Their agreement checks neither.
- **The a_1 mass cannot rescue it.** Lowering m_a1 raises the sum-rule factor, so the
  observed 1230 MeV gives 8.15 keV; the measured width would need m_a1 = 1409, further
  from observation than the framework's 1348.
- **The Regge trajectory does not separate orbital from radial excitation.** The observed
  slopes differ by 1.69 (0.91e6 MeV^2 to the a_1 against 1.55e6 to the rho(1450)), while
  the framework uses one step for both, which is why the same formula puts the a_1 10%
  high and the rho(1450) 8% low.
