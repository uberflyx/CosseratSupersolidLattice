# Cosmology

Scripts for the cosmological predictions of the Cosserat supersolid lattice.

## Scripts

| Script | Purpose |
|--------|---------|
| `gw_spectrum_crystallisation.py` | Gravitational wave spectrum from the vacuum phase transition (fluid → FCC crystal). Sound-wave mechanism with all four thermodynamic parameters derived from lattice mechanics. Peak at ~354 nHz (testable by SKA). |
| `potts_fcc_mc.py` | Monte Carlo simulation of the 3-state Potts model on FCC, verifying the deconfinement temperature T_c = 156.1 MeV from the Polyakov loop mechanism. |

## Key results

### GW spectrum
- **Peak frequency**: 354 nHz (robust, set by nucleation geometry)
- **Peak amplitude**: h²Ω ≤ 2.6 × 10⁻¹⁰ (upper bound — see below)
- **Transition strength**: α = 0.94 (derived from bag model, not fitted)
- **Self-consistency**: B^(1/4) = 228 MeV matches Λ_QCD = 220 MeV to 3.8%

The peak amplitude is an **upper bound** because the standard sound-wave
efficiency κ_v was derived for a classical fluid pushed by expanding bubbles.
The vacuum lattice crystallises by material addition at the crystal front,
which may generate less anisotropic stress than the standard scenario assumes.
Bubble collisions are negligible (κ_wall ≈ 2 × 10⁻¹⁹) because the crystal
is incompressible and absorbs wall kinetic energy.

### Deconfinement temperature
T_c = 156.1 MeV from the D4 Polyakov loop mechanism, matching lattice QCD
(156 ± 3 MeV) to 0.06%. This is the first observable that distinguishes
FCC from D4 — FCC gives only the Lindemann estimate ~230 MeV (48% high).
