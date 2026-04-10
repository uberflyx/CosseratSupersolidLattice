# Cosmology

Scripts for the cosmological predictions of the Cosserat supersolid lattice framework.

## Files

### gw_spectrum_crystallisation.py
Gravitational wave spectrum from the vacuum lattice crystallisation phase transition.
Uses standard fitting formulas (Caprini et al. 2016; Hindmarsh et al. 2017) with
lattice-determined parameters. Peak at ~350 nHz (PTA/SKA band).

### potts_fcc_mc.py
Monte Carlo simulation of the 3D ferromagnetic 3-state Potts model on the FCC lattice.
Determines the critical coupling K_c(FCC) = 0.258 ± 0.003, which converts to the
QCD deconfinement temperature T_c = 156.1 MeV via the D4 Polyakov loop mechanism.

**Key result**: T_c = ε_SF / ((3+√3) × K_c) = 156.1 MeV, matching lattice QCD
(156 ± 3 MeV) to 0.06%. This is the first observable that distinguishes D3 from D4:
D3 gives T_melt ≈ 230 MeV (48% off), D4 gives 156.1 MeV (0.06% off).

Inputs:
- ε_SF = (π√3/2)m₀ = 190.5 MeV (stacking-fault energy per temporal plaquette)
- K_c = 0.258 (from MC on FCC lattice, L=8-12, numba-accelerated)
- Z_eff = 3 + √3 (geometric factor from {111} fault Burgers vector projection)

**Dependencies**: numpy, numba, scipy
