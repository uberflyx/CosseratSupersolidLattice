# Cosmology

Scripts for the cosmological predictions of the Cosserat supersolid lattice framework.

## Files

### gw_spectrum_crystallisation.py
Computes the stochastic gravitational wave background from the vacuum lattice
crystallisation phase transition. Uses the standard fitting formulas (Caprini et al.
2016, 2020; Hindmarsh et al. 2017; Ellis & Lewicki 2020) with all four
thermodynamic parameters determined from FCC lattice mechanics:

- **T\*** ≈ 150–220 MeV (nucleation temperature, from the Debye temperature Θ_D = πm₀)
- **α** ≈ 0.03–1.0 (transition strength, principal uncertainty)
- **β/H\*** ≈ 5–8 (inverse duration, from CNT barrier ΔG\* ≈ 1.17 GeV)
- **v_w** = c/√3 (Jouguet detonation velocity)

**Key result**: the GW spectrum peaks at ~350 nHz (in the PTA band), with amplitude
h²Ω_peak ~ 10⁻¹⁵ to 10⁻⁹ depending on α. The prediction is testable by SKA.

Produces two figures:
- `gw_spectrum_crystallisation.png/pdf` — spectrum for four α scenarios vs NANOGrav data
- `gw_spectrum_decomposition.png/pdf` — sound wave / bubble / turbulence decomposition

**Dependencies**: numpy, matplotlib
