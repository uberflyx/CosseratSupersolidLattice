# Cosmology

Scripts for the cosmological predictions of the Cosserat supersolid lattice.

## Scripts

| Script | Purpose |
|--------|---------|
| `gw_spectrum_crystallisation.py` | Gravitational wave spectrum from the vacuum phase transition (fluid → FCC crystal). Sound-wave mechanism with all four thermodynamic parameters derived from lattice mechanics. Peak at ~354 nHz (testable by SKA). |
| `potts_fcc_mc.py` | Monte Carlo simulation of the 3-state Potts model on FCC, verifying the deconfinement temperature T_c = 156.1 MeV from the Polyakov loop mechanism. |
| `crystallisation_baryogenesis.py` | Order-of-magnitude estimate of the cosmic baryon asymmetry from the vacuum-lattice crystallisation transition. Geometric baryon-number violation (frozen stacking winding), first-order out-of-equilibrium, CP from compact-direction propagation chirality. Reports the CP-bias bracket, the required transition efficiency, and the radiation-era epoch (~24 us). |
| `vacuum_line_web.py` | Monte Carlo percolation of the partial-dislocation line network left by the random Z3 stacking choice at crystallisation (Kibble / Vachaspati-Vilenkin construction run natively on the framework's 3-state order parameter). Traces both the line network and the domain network of a random Z3 field, measures the infinite-string fraction and the loop-size law, and checks the triple-junction coincidence: every line lies on an edge where all three stacking domains meet (the line web is the filament skeleton of the foam). Generates the line-sector chapter's web figure. |
| `time_dislocation_calcs.py` | Tests the time-dislocation reading of the compact direction against the framework's own numbers (CODATA 2022 constants). Verifies the baryogenesis chain (eta_B = 8.5 theta_ch^2 reproduces Planck to 0.3%; the D4 route with its Bose factor agrees iff the survival efficiency ~ theta_ch); tests time-dislocation baryogenesis as a Boltzmann estimate (matter/antimatter = forward/backward windings; reproduces the structure but not the second power of theta_ch -- a consistency check, not a new result); and quantifies chronology protection (a CTC needs the light cone tilted past vertical, which saturates the lattice at ~5e32 J/m^3, ~10^42 x dark energy). |

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
T_c = 156 MeV from the D4 Polyakov loop mechanism. The rigorous, parameter-free
part is the pure-gauge scale: the bare three-bond count gives 246 MeV against the
pure-gauge lattice value 277 MeV (11%). Dynamical-quark (flavour) screening lowers
this to the full-QCD 156 MeV, an agreement at the ~10% level of the screening
estimate rather than an exact central-value match. This is the first observable that
distinguishes FCC from D4. FCC gives only the Lindemann estimate ~230 MeV (48% high).

### Baryon asymmetry
The deconfinement transition freezes net stacking winding into net baryon number,
with CP supplied by the compact-direction propagation chirality. The CP supply
exceeds the Standard Model by ~10^16, so the mechanism overshoots the observed
eta_B = 6.1 x 10^-10 at unit efficiency and requires a net transition efficiency
of ~10^-6 to ~10^-1, depending on the (open) resummation power of alpha. The sign
is a definite prediction: matter over antimatter, fixed by the stacking handedness
and correlated with the sign of the heavy-ion chiral-magnetic charge correlator.
The freeze-in is a QCD-epoch event at t ~ 24 us.

### Vacuum line web
The partial-dislocation lines written by the random Z3 stacking choice at the
freeze percolate. A single connected, system-spanning string carries about
three-quarters of all line length (infinite-string fraction settling toward the
0.75-0.80 Vachaspati-Vilenkin band as the box grows), with the remainder in
finite loops whose sizes follow the random-walk law n(l) ~ l^{-5/2}. Percolation
probability is 1 at every size tested, and the result is robust to lattice and
discretisation.

The lines also lie exactly on the edges of the stacking-domain foam: a Z3
winding requires all three registries, so 100% of line plaquettes touch all
three labels (triple-junction edges), and about two-thirds of triple-junction
edges carry a net line. In the foam-to-cosmic-web dictionary (cells=voids,
faces=walls, edges=filaments, vertices=nodes) the line web is the filament
skeleton, provided the stacking coherence length is grain-scale rather than
microscopic; Lorentz invariance disfavours the microscopic case. Whether the
filament lines persist across cosmic time is the open survival question.
