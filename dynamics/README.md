# Dynamics

One-loop scattering and radiative correction computations from the Cosserat lattice Feynman rules.

## Scripts

### `ew_one_loop_lattice.py`

Computes the complete one-loop electroweak programme from lattice first principles:

- **Leptonic vacuum polarisation** — exact one-loop integral with Koide-derived lepton masses. Matches the SM value to three significant figures (0.0314 vs 0.0315).

- **Hadronic vacuum polarisation** — once-subtracted dispersion integral over the lattice vector meson spectral function (ρ, ω, φ, J/ψ, Υ(1S)), with perturbative QCD at Λ_QCD = π·m₀ above the quark-hadron duality matching point. Zero free parameters. Result: Δα_had = 0.027526 ± 0.000066 against the data-driven 0.027609 ± 0.000112, a residual of −0.30% or 0.6σ. Against the PDG electroweak review's α_s-driven 0.02783 ± 0.00006 the residual is −1.09%, which the band turns into 3.4σ, so the construction excludes one of the two established determinations rather than sitting between them. The radial vector tower widths are derived, not fitted: the Cornell potential built from the framework's own string tension and α_s, solved by radial diagonalisation, reproduces both measured bottomonium radial ratios to 1–3% where the pure-linear limit was off by factors of 1.9 and 2.4. All three vector families are built explicitly on their own Regge trajectories, the continuum carries the standard QCD series through c₃ (verified against the PDG QCD review, eq. 9.9), and the matching point is quantised by the tower's own spacing: s_match(N) = m_ρ² + (N+½)·2πσ. The resolution condition Γ_n < Δm_n admits N = 1 and N = 2 only, the second recurrence being the last resolved state and marginal at Γ_n/Δm_n = 0.95. The quoted uncertainty is the spread over those two points (±0.000042), the stated Cornell coupling (±0.000047) and the constituent mass (±0.000020), in quadrature.

- **Custodial-breaking ρ parameter** — UV-convergent integral from the top-bottom mass splitting in the evanescent (W) self-energy. Result: Δρ = 0.0092, matching the SM value ~0.0094.

- **Running coupling** — α⁻¹(M_Z) = 128.957 ± 0.009 against the measured 128.946 ± 0.015, an agreement at 0.7σ. Read backwards, the same agreement bounds any charged state the framework's spectrum does not contain: below the Z, N_c Q² < 0.038 at 131 MeV rising to 0.733 at 45 GeV, or Q < 0.2 e at the light end for a colour singlet.

- **Corrected G_F** — Δr = Δα − (c²/s²)Δρ = 0.0265 reduces the tree-level residual from 3.2% to 0.6%.

- **Effective Weinberg angle** — tree-level residual reduced from 4.0% to 1.7%.

- **Muon anomaly** — a_mu^HVP,LO = 700 (+4/−3) × 10⁻¹⁰, one per cent above the data-driven 692.78 ± 2.42 and two per cent below WP25's 713 ± 6. The π π channel uses a two-resonance Gounaris-Sakurai form factor with F_π(0) = 1 imposed exactly; a rescaled Breit-Wigner violates that constraint by 48%, inflating the low-s region the muon kernel weights most. c_0 is fixed by the **pole residue**, not by the channel's spectral area: Γ_ee sets the peak of the form factor, and reading the narrow-resonance relation at s = m_ρ² gives c_0 |A_0(m_ρ²)| = sqrt(36 Γ_ee / (α² β³ Γ_ρ)), hence c_0 = 1.132. The area rule returns c_0 = 1.199 on the same inputs, high by eight per cent in c_0 and sixteen in the channel; applied to *measured* inputs the residue rule returns 1.1101 against the 1.1124 the dispersive two-pion determination independently demands, an agreement to a fifth of a per cent, which is what settles the choice.

  The dominant systematic is the F_π truncation, and its sign follows from that choice. c_2 is the product of a photon coupling the Cornell solve gives and a pion coupling the framework does not, so it is bounded rather than predicted: |c_2| ≤ 0.112 on the loosest of three residue conventions (the other two give 0.099 and 0.087). Because the peak is pinned, c_0 is re-solved whenever c_2 moves, so a third resonance *displaces the second* instead of rescaling the channel, and **both signs raise a_mu**, by at most 2.3 units. The two-resonance value is therefore a floor rather than a maximum, and the error is one-sided upward, toward WP25 rather than away from it. `formfactor_truncation_independent.py` checks this on a different line shape and returns the same sign and order; under the area rule the same sweep moves a_mu *downward* by an order of magnitude more, which is one more reason that rule is the wrong one.

  Γ_ee(ρ) carries a derived finite-pion-mass correction: both Weinberg sum rules are chiral-limit statements, and giving the pion its mass moves its pole off the origin, adding f_π²m_π² to the first moment and multiplying the coupling by (1 − m_π²/m_a1²). The factor is small, −1.05%, because the pole enters the *first* moment and is therefore weighed against m_a1 rather than m_ρ. The result, 7.288 keV, stands +3.5% above the pre-2023 average of 7.04 ± 0.06 keV. Whether that is a defect depends on which two-pion measurements survive: the construction lands within 0.1σ of CMD-3 and well above BaBar and KLOE, so it takes the high side of a disagreement that itself decides whether the muon anomaly survives. What is not a bet is the shape, and it is the sharper result: the predicted ρ is too narrow by about five per cent, read off a ratio of window integrals that cancels the normalisation entirely.

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

### `hvp_chapter_audit.py`

Recomputes every framework-derived decimal the monograph's hadronic-vacuum-polarisation
and running-coupling chapter prints, from alpha and m_e through the closed forms the
chapter itself states, and compares each with the printed value at the printed
precision.  It shares no code with `ew_one_loop_lattice.py`, so agreement is independent
evidence rather than a copy; that separation is the point, since a script and a text that
were copied from one another agree even when both are wrong.  Around ninety literals:
the node mass and the scales built on it, the light hadron sector, both routes to the rho
leptonic width and the near-identity that makes them agree, the chirally split doublet
and where the Regge trajectory stops working, the tower and the duality matching, the
leptonic running and the zone-boundary waypoint, the spectral-sensitivity table, the spin
read backwards from the measurement, the node form factor, and the closure census.
Exits non-zero if anything fails to reproduce.

Two conventions are set deliberately in the source and are easy to get wrong.  Every
evaluation of the one-loop shift is at **spacelike** momentum transfer, since the running
is a function of momentum transfer; evaluated timelike, the real part near a pair
threshold turns negative and the charge-census bound at 45 GeV comes out meaningless.
And the Weinberg sum rules take the **isospin-average** pion, being isospin statements
about the pion pole, where the two-pion threshold takes the charged mass.

```
python3 hvp_chapter_audit.py
```

### `formfactor_truncation_independent.py`

A second, deliberately unshared calculation of the same truncation systematic that
`formfactor_truncation.py` reports, built on a plain relativistic Breit-Wigner rather
than the Gounaris-Sakurai line shape the production calculation uses.  The sensitivity
being tested is a property of how the normalisation constraint propagates, so it has to
survive a change of line shape to mean anything.  Computes the systematic under both
candidate rules, area and peak, and asserts F_pi(0) = 1, the pole condition against the
leptonic width, and both analytic limits of the muon kernel before quoting a
sensitivity.  Slow: the area-rule solve nests a quadrature inside a root find.

```
python3 formfactor_truncation_independent.py
```

### References

Monograph Chapter 10 (Dynamics), Sections on the two-loop QED programme, Landau pole elimination, and the electroweak radiative correction programme.

### `regge_massive_endpoints.py`

Solves the classical rotating Nambu-Goto string with massive endpoints exactly, on the
framework's fixed inputs (sigma = (2 pi m_0)^2, m_q = 315 MeV), against the leading
light trajectory rho(770), a2(1320), rho3(1690), a4(1970) and the chiral doublet
centres.  Written to test the corpus claim that endpoint masses close the -7% Regge
slope gap; it disproves it.  Massive ends raise every M^2 step above 2 pi sigma and
push the fitted slope to 0.71-0.77 GeV^-2, away from the observed 0.92, while the
measured steps (1.137, 1.114, 1.017 GeV^2) sit below 2 pi sigma and shrink, a pattern
that tracks the growth of the widths along the trajectory (unitarity self-energy, sized
in the script at Gamma M / pi = 0.04-0.20 GeV^2 against step deficits of 0.08-0.20)
rather than any endpoint mass.  A labelled one-parameter diagnostic inversion is
included; it returns m = 0 as the best fit at fixed sigma, which is the same exclusion
read backwards.

```
python3 regge_massive_endpoints.py
```

### `rho_sector_audit.py`

Audits the rho sector, which supplies the single most consequential number in the
hadronic programme. Establishes what is in good order (mass 0.05%, KSRF width 0.23%,
KSRF coupling tested against the measured width 0.15%) and localises what is not:

- **Gamma_ee(rho) is the outlier** among the three vector leptonic widths the framework
  computes (omega +0.5%, phi -1.1%); no single rescaling fits rho and phi together, so
  the defect is specific to the rho's coupling to the photon.  The script reports +4.6%,
  which is the chiral-limit width against the pre-2023 average.  The derived value now
  carries the finite-pion-mass correction above and stands at +3.5%; the audit's
  conclusion is unchanged, since the correction is common to every route.
- **The two routes to that width are not independent.** They agree only through the
  numerical near-identity 8 pi^3 ~ 2 (m_rho/m_0)^2, 248.05 against 245.16 and so true
  to 1.2%, which makes
  m_a1^2 = 3 m_rho^2 and the Weinberg factor 3/2. Their agreement checks neither.
- **The a_1 mass cannot rescue it.** Lowering m_a1 raises the sum-rule factor, so the
  observed 1230 MeV gives 8.15 keV; the measured width would need m_a1 = 1409, further
  from observation than the framework's 1348.
- **The Regge trajectory does not separate orbital from radial excitation.** The observed
  slopes differ by 1.69 (0.91e6 MeV^2 to the a_1 against 1.55e6 to the rho(1450)), while
  the framework uses one step for both, which is why the same formula puts the a_1 10%
  high and the rho(1450) 8% low.
