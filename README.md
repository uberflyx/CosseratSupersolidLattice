# Cosserat Supersolid Lattice

**Mitchell A. Cox, University of the Witwatersrand**

Companion code for [*The Cosserat Supersolid*](https://doi.org/10.5281/zenodo.18636501)

Full monograph: `https://doi.org/10.5281/zenodo.18636501`

---

## Background

Since I was a kid watching Star Trek, I've wanted to build a warp drive. Figuring out the physics was obviously the first step. That's probably why I eventually went into academia — I started with an MSc in high-energy physics and a member of ATLAS/CERN, then I did a PhD in photonics. My day job is essentially electrical engineering: telecoms, optical communications, structured light photonics. But the dream never went away. I've kept my mind on the exotic stuff: on what space actually *is*, for years. I have several physicist friends who have patiently listened to my hare-brained theories over the years. I appreciate you guys!

The core observation is simple. Materials scientists have spent a century classifying what can go wrong in a crystal: dislocations, vacancies, stacking faults, grain boundaries. Particle physicists have spent a century classifying what exists in the vacuum: electrons, quarks, photons, maybe gravitons, maybe dark matter and something we call dark energy. At the same time, physicists often do things like subtract infinities (renormalise), and the direct inspiration for this work was: what if we walked a parallel history to when we renormalised Dirac's particle sea, but used the math and physics we know nowadays to avoid doing that? Now, this isn't literally that, but it's what got me thinking.

The catalyst was LLMs. Whenever a random idea popped into my head — like the Dirac sea question or "what if the fine structure constant is a tunnelling amplitude?" or "what if gravity is just the compression channel of the same lattice?" — I could actually *try it out*, immediately, at 9pm, without needing to be an expert in every subfield the question touched. Modern LLMs have, in some practical sense, the sum total of human knowledge inside them. I felt that we collectively know enough to make this work — it just needed someone to sit down and try to bring it all together. This monograph is the result.

I don't want to trivialise it. This is extremely multidisciplinary — it draws on Cosserat continuum mechanics, crystallography, quantum field theory, general relativity, nuclear physics, and cosmology. Academia unfortunately encourages single-minded focus, and a cross-disciplinary approach is usually frowned upon. I'm aware that an electrical engineer claiming results in particle physics will raise eyebrows. But the derivations are all written out, the code reproduces every number, and the predictions are falsifiable. I'd rather put it out there and be wrong than sit on it and wonder.

I used AI (Claude Opus and a little bit of Gemini) as iterative discussion partners throughout — to stress-test arguments, check and derive mathematics, spot errors, and move much faster than I could alone. All responsibility for the physics reasoning, derivations, and conclusions is mine. Fortunately my general knowledge is good and my MSc helped me a lot, but I have spent quite some time brushing up on materials science so that I can smoke test this work!

The result appears to be on par with a theory of everything. I'm aware of how that sounds. But the numbers are what they are: ~60 predictions from three inputs (*c*, *ℏ*, *mₑ*), spanning from the fine structure constant to the Higgs mass, with no fitted parameters and nearly all retrodictions at pull < 1. The framework derives both α (to 0.003 ppb) and *G* (to 2 ppm) from the same lattice geometry. It predicts hadron masses, lepton masses, quark masses, mixing angles, neutrino mass splittings, dark matter, dark energy and the cosmological baryon asymmetry. It is not numerology — every integer that appears is derived from the FCC coordination geometry before any mass is consulted, and the monograph shows the working step by step.

The monograph is a work in progress but fairly complete. It has many falsifiable predictions (many are explicitly catalogued, several testable with existing data, particularly at ALICE). It also has open questions and calculations that I haven't finished. I'm working toward splitting it into peer-reviewable papers, but the claims are unconventional and that process will be very challenging. In the meantime, I'm putting everything out here — the full monograph, all the derivations, and the code to reproduce every number — because I think it deserves scrutiny, and I'd rather people can check the work now than wait.

Whether the leading-order agreement is evidence for the hypothesis or an elaborate coincidence is for the reader to judge. The purpose of the monograph is to present the derivations in sufficient detail that the question can be answered.

Please go read it: `https://doi.org/10.5281/zenodo.18636501`

---

## Summary

The vacuum is modelled as an FCC Cosserat supersolid — a face-centred cubic crystal that is simultaneously superfluid (rigid to shear, frictionless to translation). Particles are topological defects: screw dislocations (electrons), partial dislocations (quarks, with the FCC threefold stacking degeneracy providing colour), edge dislocations (neutrinos), and vacancies (dark matter). Forces are elastic waves: transverse shear (electromagnetism), stacking-fault elasticity (strong force), evanescent modes (weak force), and longitudinal compression (gravity).

The fine structure constant is the Boltzmann factor for tunnelling through the Peierls–Nabarro barrier. Newton's constant requires coherent tunnelling of all 19 nodes in the FCC coordination cluster, giving *G* ∝ α¹⁹ — the 10³⁸ hierarchy between electromagnetism and gravity is a counting problem (1 node vs 19). A universal mass formula *m* = *m*ₑ(*N*/α*ᵏ* + *Q*) yields the full particle spectrum from five geometric building blocks and a four-term electromagnetic correction, each derived from the defect geometry with no reference to experiment.

The obvious objection — that a crystal breaks Lorentz invariance — is addressed quantitatively in the monograph and in a [dedicated paper](https://doi.org/10.5281/zenodo.18739953) (submitted, under review). Two suppression mechanisms intrinsic to the lattice (the Peierls–Nabarro form factor and Debye–Waller smearing from quantum zero-point delocalisation) combine to push Bragg scattering five orders of magnitude below the Fermi-LAT gamma-ray bound. The lattice is there, but it hides well.

---

## Files

| Script | What it does |
|---|---|
| `cosserat_calculator.py` | Full prediction engine — every mass and coupling in the monograph |
| `pn_variational.py` | Derives α and *G* from the Peierls–Nabarro tunnelling calculation |
| `verify_fcc_geometry.py` | Verifies every geometric claim in Ch. 9 against explicit FCC coordinates |
| `exotic_catalogue.py` | Blind mass-formula test on all 32 known exotic hadrons |
| `exotic_filling_fraction.py` | Honesty check — proves the mass formula is trivially flexible without independent *N* derivation |

Requirements: Python 3, NumPy, SciPy.

---

## Output: `cosserat_calculator.py`

```
╔══════════════════════════════════════════════════════════════════════╗
║  THE COSSERAT SUPERSOLID: Complete Mass & Coupling Calculator      ║
║  Mitchell A. Cox, University of the Witwatersrand                  ║
║                                                                    ║
║  Inputs:  c, ħ, mₑ = 0.51100 MeV                                 ║
║  Everything else derived from FCC Cosserat geometry.               ║
╚══════════════════════════════════════════════════════════════════════╝

  FUNDAMENTAL CONSTANTS
  α⁻¹ (fine structure const.)              137.036      137.036          +0.000%   0.00
  G (Newton's constant)                6.67643e-11   6.6743e-11 m³/kg/s²   +0.032%   0.04

  QUARK MASSES
  m_u (up quark)                           2.18696         2.16 MeV      +1.248%   1.71
  m_d (down quark)                         4.74196         4.67 MeV      +1.541%   2.11
  m_s (strange)                             93.367         93.4 MeV      -0.035%   0.05
  m_c (charm)                               1269.65         1270 MeV      -0.027%   0.04
  m_b (bottom)                              4176.99         4180 MeV      -0.072%   0.10
  m_t (top)                                172.815       172.57 GeV      +0.142%   0.19

  HADRON MASSES
  π± (pion)                                 139.54       139.57 MeV      -0.022%   0.03
  K± (kaon)                                493.754      493.677 MeV      +0.016%   0.02
  η (eta)                                  547.938      547.862 MeV      +0.014%   0.02
  ρ (rho)                                  775.664       775.26 MeV      +0.052%   0.07
  φ (phi)                                  1019.45      1019.46 MeV      -0.001%   0.00
  p (proton)                               939.209      938.272 MeV      +0.100%   0.14
  Λ (Lambda)                               1115.81      1115.68 MeV      +0.011%   0.01
  Ω⁻ (Omega)                               1672.43      1672.45 MeV      -0.001%   0.00
  Δ(1232)                                  1233.11         1232 MeV      +0.090%   0.12
  d*(2380) dibaryon                        2382.39         2380 MeV      +0.100%   0.14

  ELECTROWEAK SECTOR
  m_H (Higgs boson)                        125.203       125.25 GeV      -0.038%   0.05
  M_W (W boson)                            80.3911        80.37 GeV      +0.026%   0.04
  M_Z (Z boson)                             91.155       91.188 GeV      -0.036%   0.05
  v (Higgs VEV)                            246.194       246.22 GeV      -0.011%   0.01

  CHARGED LEPTONS
  Σm_ℓ = 27m₀ − 15mₑ                       1883.02      1883.03 MeV      -0.001%   0.00
```

Full output includes neutrino parameters, mixing angles, cosmological predictions, Regge trajectories, and derivation notes for each particle.

---

## Output: `pn_variational.py`

```
  alpha_PN^-1  = 137.035999177348
  CODATA       = 137.035999177 +/- 2.1e-08
  Residual     = +0.003 ppb  (+0.017 sigma)

  G_PN         = 6.6743e-11 m^3 kg^-1 s^-2
  G_CODATA     = 6.6743e-11 +/- 1.50e-15
  Offset       = -2.1 ppm
```

---

## Output: `verify_fcc_geometry.py`

```
CHECK 1: Coordination shell (N = 13)         ✓
CHECK 2: {111} intersection lemma             ✓ (all 6 pairings: 2 shared)
CHECK 3: Hexagonal cap (N = 7)               ✓
CHECK 4: Tetrahedral voids (+4)              ✓
CHECK 5: Hexagonal bilayer (N = 8)           ✓ (2 ring + 1 centre = 3 bonds)
CHECK 6: Cell pair common-NN → factor 5      ✓
CHECK 7: Crossed fault N = 144/13            ✓
CHECK 8: Proton mass = 939.21 MeV            ✓ (+0.100%)
CHECK 9: Spin self-energy ratio = 6.5×       ✓
CHECK 10: Bond energy = αm₀ = mₑ            ✓

RESULTS: 21/21 passed, 0/21 failed
```

---

## Output: `exotic_filling_fraction.py`

This script exists because honesty matters more than a good scorecard.

```
  1. The mass formula m = Nm₀ + Qmₑ can fit ANY mass (100% filling fraction)
     → 32/32 PASS is mathematically trivial, NOT a physics result

  2. The REAL test is whether N decomposes into known constituent N values
     → Molecular decomposition works for threshold states
     → Far-from-threshold states remain unconstrained

  3. HONEST FRAMING:
     → The formula ACCOMMODATES all 32 exotics (consistency check)
     → But it does not PREDICT them (no a priori N derivation)
     → This is a necessary condition, not a sufficient one
```

---

## AI acknowledgement

Large language models (Claude Opus and Gemini) were used as iterative discussion partners throughout the development of both the monograph and this code: to stress-test physical arguments, identify logical vulnerabilities, check mathematical consistency, and improve clarity. All responsibility for the scientific content, physical reasoning, derivations, and conclusions is the author's.

## Licence

MIT

---

Warp drive doesn't jump out of this framework (unfortunately). But if the vacuum really is a crystal, then the universe is (now) an engineering problem — and that's a start.
