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
| `wall_template_flow.py` | Toy model of the peculiar-velocity field sourced by the pinned fossil-wall dark energy. A periodic Poisson-Voronoi foam (21.6 Mpc cells) in a 1280 Mpc box on a 384^3 grid carries the full dark-energy density on its faces with the active weight -2 that the fossil stress fixes (the walls repel); the perturbation Poisson equation is solved by FFT and velocities are built kinematically over a Hubble time (an upper estimate). Reports point velocities, bulk flows in spheres, and an uncorrected mock distance ladder, pooled over eight foam realisations. Exports `wall_template()` for the two scripts below. |
| `wall_rsd_bound.py` | Pairwise velocity of wall-resident tracers in the same template, against the window that redshift-space growth data leave, with a finer-grid resolution check; and the displacement the template builds at the BAO scale. |
| `sn_wall_residual.py` | The template as a low-redshift supernova systematic: the coherent sky-mean magnitude offset between a wall-resident observer and wall-resident hosts, per redshift bin. |
| `wall_foam_lensing.py` | Lensing and microwave-sky signals of the walls against Planck LCDM (CAMB): the foam's contrast spectrum, the walls' potential power relative to the matter's, the Limber convergence power for galaxy and CMB sources in the proper and comoving readings, the integrated Sachs-Wolfe temperature power, and a void-stack excess surface density. |

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

### Pinned-wall gravity
The fossil walls carry a stress trace of -3 times their energy density at every point, so in the static weak field their Newtonian source is -2 epsilon and their lensing source -epsilon: they push matter away and defocus light, with the active density of a cosmological constant gathered into sheets. The scripts above size the consequences on the kinematic upper estimate, with the full observed dark-energy density on the walls.

- Velocities (`wall_template_flow.py`, eight realisations): about 455 km/s rms; bulk flow medians 198, 113, 83 and 68 km/s in spheres of 30, 100, 150 and 200 Mpc, falling roughly as R^(-1/2) because the foam's density is white noise above the cell scale. That is far below the CosmicFlows-4 bulk flow. The uncorrected ladder bias has an observer-to-observer spread of 0.57 km/s/Mpc around a mean of +0.009, so the local-side wall channel on the Hubble tension is closed at about half a km/s/Mpc.
- Redshift-space consistency (`wall_rsd_bound.py`): the pairwise outflow of wall-resident tracers is 48 km/s at 10-15 Mpc (56 km/s on a 1.67 Mpc grid), 1.6 to 1.9 times the ~30 km/s window that growth data allow, so the full template is in tension with redshift-space distortions on this estimate. At 60-100 Mpc it is 2.4 km/s, a displacement of 0.018 Mpc at the acoustic scale, about 0.01% of the sound horizon: the BAO-plus-CMB neutrino bound is untouched.
- Supernova residual (`sn_wall_residual.py`): rms magnitude offset 0.039, 0.019 and 0.010 mag at z = 0.015, 0.025 and 0.04 for wall-resident hosts and observer, capped at about 0.02-0.025 mag at z = 0.015 by the redshift-space bound. Wall residence is an assumption, since the walls repel.
- Lensing and CMB (`wall_foam_lensing.py`): the walls carry 8.7% and 7.6% of the matter's convergence power at multipole 300 for sources at z = 0.6 and 1 (proper reading; up to a factor of two lower in the comoving reading), 1-2% at multipoles 30 and 3000, and below 1% for CMB lensing. Their integrated Sachs-Wolfe power is 0.2% of the Planck temperature power at multipole 30 in the comoving reading, a third of it at the quadrupole, below cosmic variance throughout. A void stack at the cell scale (mean cell radius 13 Mpc) gets +0.46 Msun/pc^2 from the walls near 0.85 R_v, against -0.17 for a matter void of the DES profile, so cell-sized voids would show a lensing profile of the opposite sign unless most of the matter sits on the walls.
