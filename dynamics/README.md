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

### `chiral_ladder.py`

Chases the observation regge_massive_endpoints.py recorded: the Adler-zero intercept
alpha_rho(0) = 1/2 (Lovelace-Shapiro), combined with the framework's own slope
1/(2 pi sigma), anchors the chiral vector trajectory at M_L^2 = (2L+1) pi sigma with
nothing fitted, i.e. m_rho = 2 pi^(3/2) m_0 in the chiral limit.  Against the corpus's
own chiral-limit constructions: the L = 1 doublet centre lands at -0.14% (0.1 sigma),
the L = 2 centre at +3.3% (the unitarity bend, separately established).  Two corpus
puzzles close structurally: the 8 pi^3 = 2 (m_rho/m_0)^2 near-identity (1.2%) becomes
exact on the ladder, the 1.2% now measuring the spectral rho's offset from
sqrt(pi sigma), and the same offset is the 0.4% gap between the two Gamma_ee routes,
which merge identically once m_a1^2 = 3 m_rho^2 is exact.  Honest caveats are printed:
the physical-pion anchor overshoots the rho by 2.2%, empirical intercept fits sit at
0.42-0.48 rather than exactly 1/2, and the pion-trajectory daughter check is weak.

The sharpest test is the level structure rather than any single mass.  At level N the
amplitude's residue carries spins 0 to N, all degenerate, so the leading J = N state
and the J = 1 daughter must sit together; the J = 1 daughter is an axial/vector pair
split by chiral symmetry breaking, whose centre is the doublet construction the corpus
already uses.  At N = 3 the two land on top of each other, rho3(1690) at 1688.8 against
the a1(1640)/rho(1700) centre at 1687.8, and at N = 2 the a2 sits 34 MeV low.  Neither
figure survives being read alone, and the script's full-rung test says why: the
amplitude puts EVERY spin from 0 to N on rung N, so the test is whether all of them sit
together.  Rung 2 spreads over 236 MeV, 17.6% of its mean; rung 3 over 65 MeV, 3.8%.
The degeneracy is badly broken low and nearly exact high, improving by a factor of 4.6
in one step.  That cuts both ways: the a2's 34 MeV is one slice of a badly split rung
rather than a puzzle about the a2, and rung 3's 1 MeV is closer than its own rung
warrants.  What the ladder gets right is the trend, being a spin-blind construction
whose departures are spin-dependent forces that fall with excitation at a rate the
framework observes rather than predicts.  The ladder's own centre also tightens the
corpus's rho' prediction from -0.42% to -0.09% and its half-splitting from -1.8% to
+0.44%.

```
python3 chiral_ladder.py
```

### `pull_scorecard.py`

Recomputes every paper-relevant framework prediction and scores it, reporting the
residual, the pull in experimental sigma alone, and the pull against experiment and
the framework's own stated resolution combined.  The three columns are kept separate
because ranking on experimental sigma alone gets the answers backwards: m_pi at 0.024%
scores 3.2 sigma against a razor-sharp PDG average while Gamma_rho scores 0.4 sigma at
ten times the residual.  Currently 18 of 25 entries sit below one sigma.  The script
also records that the two a_mu benchmarks disagree with each other by 3.1 sigma, so no
prediction can sit within one sigma of both and the pair bounds what is achievable.

```
python3 pull_scorecard.py
```

### `wsr_tower_saturation.py`

Asks whether saturating the Weinberg sum rules with the framework's own derived tower
closes the +3.5% excess on Gamma_ee(rho), which is the largest single pull in the
scorecard.  It does not, and the script is the record of why.

Three findings, all negative or diagnostic:

- The excess is not the framework's inputs.  Feeding the corpus formula MEASURED f_pi
  and m_rho still returns 7.240 keV, +2.8%, so single-state saturation itself carries
  most of it and the derived inputs add 0.7%.
- The tower is the wrong tool.  Solving both sum rules with two vectors and two axials
  on the framework's own tower moves f_rho^2 by -16% to +14% depending on the axial
  coupling ratio, which the framework does not derive.  A 15% tool cannot fix a 3%
  gap.  Truncated sum rules are known to behave this way, since the two sums converge
  only in the difference and that difference is controlled by chiral restoration at
  high rung rather than by the first recurrence.
- The comparison value is not soft.  All six published Gounaris-Sakurai determinations
  are listed; they agree with each other at chi2/dof = 0.76, so the PDG error needs no
  inflating and the pull is real.

What the script does establish in the framework's favour: Weinberg's classic
m_a1^2 = 2 m_rho^2 misses the measured width by 39%, and the framework's own a1 at
m_a1^2 = 3 m_rho^2 brings that to 3.5%, a factor of eleven.  Closing the remainder
needs m_a1^2/m_rho^2 = 3.26, or m_a1 = 1399 MeV, which no route in the corpus produces.

```
python3 wsr_tower_saturation.py
```

### `vector_normalisation.py`

Asks whether the 3.5% excess on Gamma_ee(rho) belongs to the rho or to the whole light
vector sector.  The omega-phi mixing angle is fixed by scaling the isoscalars off the
MEASURED rho, and doing it instead off the framework's own f_pi makes the two isoscalars
disagree (4.5 degrees against 0.7); the whole of that disagreement is the rho's own
excess, so the isoscalars carry information about it.

Fitting one common normalisation k and one mixing angle to all three measured widths at
once returns k = 0.957 and delta = 3.08 degrees, chi2 = 3.3 on one degree of freedom,
bringing all three residuals below 1.1 sigma where ideal mixing leaves them at 4.1, 8.9
and 1.5.  What matters is not the fit quality but that the rho alone asks for k = 0.966
while the two isoscalars alone, fixing delta from their ratio with no normalisation
involved, then ask for k = 0.941: two determinations from different states and different
charge weights, agreeing to 2.5% with no chance to negotiate.

So the excess is a sector-wide normalisation rather than a rho-channel defect, about 4%
in the width and 2% in the coupling.  Three candidate sources are then tested and all
three are excluded, which is recorded so they are not tried again.  The chiral-limit
decay constant has the right direction and badly the wrong size, k wanting
F_0/f_pi = 0.978 where the lattice sits near 0.94, which would drive the rho to -8.5%.
The N_c exponent of the f_pi derivation would need N_c^0.230 against a topological
susceptibility that confirms N_c^0.249.  Giving each channel its own chiral partner
instead of a common Weinberg factor fails informatively: the rho lands at +4.2% while
the omega and phi both land at +34%, agreeing to half a per cent, so the isoscalar pair
shares a factor the rho does not, which is the opposite of a sector-wide k.

What survives is the measurement and not its explanation.  k is reported as a
diagnostic and the corpus widths are not rescaled by it.

```
python3 vector_normalisation.py
```

### `regge_bend_mechanisms.py`

Tests, against the spectrum rather than by argument, the three mechanisms proposed for
the bend in the chiral ladder, and then asks whether the reported Regge slope defect is
a defect at all.

Massive endpoints are excluded elsewhere (regge_massive_endpoints.py).  The corpus's
unitarity claim, that wider states are pulled further below their bare mass, is
falsified on the rungs: the broadest states on the trajectory are the a1(1260) and
rho(1450) at 400 MeV apiece and they sit on the rung the ladder places best, while the
narrower a1(1640) and rho(1700) at 250 MeV sit on the rung it misses.  On the leading
trajectory the widths do grow and a regularity survives for the a2, rho3 and a4 at
0.62, 0.69 and 0.61 of Gamma M before the rho breaks it by a factor of ten.  A 1/N_c
reading does not scale either, the rung with the largest Gamma/M showing the smallest
departure.  Nor does the square-root trajectory alpha(s) = 1/2 + a1(sqrt(s0) -
sqrt(s0-s)), the form analyticity motivates: with the framework's slope imposed, each
rung solves for s0 separately and returns 67, 29 and 13 GeV^2, three scales falling by a
factor of five rather than one scale.  Four mechanisms excluded by calculation.

The last section is the useful one.  A linear Regge fit has two parameters, and fitting
both to a bending trajectory spreads the bend across them, so part of the apparent
slope error is the intercept absorbing curvature.  Holding the intercept at the Adler
value of exactly 1/2 lets each rung state a tension of its own: the first returns
2 pi sigma = 1.2136 GeV^2 against the derived 1.2163, agreeing to 0.23%, while the
second and third fall to 1.1972 and 1.1460 as the trajectory bends.  So the derived
string tension is confirmed to a quarter of a per cent where the trajectory is straight,
and the "-7% Regge slope" is the bend counted twice rather than a separate defect.  No
single tension can absorb the bend in any case, since the M^2 deficits divided by
(N - 1/2) run 0.0027, 0.0191, 0.0703 and are not constant.

```
python3 regge_bend_mechanisms.py
```

### `omega_phi_mixing.py`

Fits the single omega-phi mixing angle that reconciles the two isoscalar leptonic
widths.  Ideal mixing puts the omega 21% high and the phi 6% low, opposite signs that no
rescaling can fix and which are the signature of a rotation.  One angle acting with
opposite effect on the two, delta = 3.09 degrees, brings both from 7.0 and 5.4 sigma to
below one, at chi2 = 1.77 on one degree of freedom, and the two isoscalars agree on the
angle separately (3.60 +- 0.52 and 2.64 +- 0.50) having had no chance to.  The check the
fit could not arrange: the resulting nonet angle is 35.26 + 3.09 = 38.36 degrees against
the standard phenomenological value near 39, and nothing in the two widths knows about
that number.

The isoscalars are scaled off the MEASURED rho rather than the framework's own f_pi, and
that is not a free choice: vector_normalisation.py shows the difference between the two
routes is the rho's own leptonic excess, and uses it to establish that the excess belongs
to the whole sector.  The angle is fitted, not derived; deriving it needs the off-diagonal
element of the isoscalar mass matrix.

```
python3 omega_phi_mixing.py
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
