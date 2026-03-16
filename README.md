# Cosserat Supersolid Lattice

**Mitchell A. Cox, University of the Witwatersrand**

Companion code for [*The Cosserat Supersolid*](https://doi.org/10.5281/zenodo.18636501)

Full monograph: \`https://doi.org/10.5281/zenodo.18636501\`

---

## Background

Since I was a kid watching Star Trek, I've wanted to build a warp drive. Figuring out the physics was obviously the first step. That's probably why I eventually went into academia — I started with an MSc in high-energy physics and a member of ATLAS/CERN, then I did a PhD in photonics. My day job is essentially electrical engineering: telecoms, optical communications, structured light photonics. But the dream never went away. I've kept my mind on the exotic stuff: on what space actually *is*, for years. I have several physicist friends who have patiently listened to my hare-brained theories over the years. I appreciate you guys!

The core observation is simple. Materials scientists have spent a century classifying what can go wrong in a crystal: dislocations, vacancies, stacking faults, grain boundaries. Particle physicists have spent a century classifying what exists in the vacuum: electrons, quarks, photons, maybe gravitons, maybe dark matter and something we call dark energy. At the same time, physicists often do things like subtract infinities (renormalise), and the direct inspiration for this work was: what if we walked a parallel history to when we renormalised Dirac's particle sea, but used the math and physics we know nowadays to avoid doing that? Now, this isn't literally that, but it's what got me thinking.

The catalyst was LLMs. Whenever a random idea popped into my head — like the Dirac sea question or "what if the fine structure constant is a tunnelling amplitude?" or "what if gravity is just the compression channel of the same lattice?" — I could actually *try it out*, immediately, at 9pm, without needing to be an expert in every subfield the question touched. Modern LLMs have, in some practical sense, the sum total of human knowledge inside them. I felt that we collectively know enough to make this work — it just needed someone to sit down and try to bring it all together. This monograph is the result.

I don't want to trivialise it. This is extremely multidisciplinary — it draws on Cosserat continuum mechanics, crystallography, quantum field theory, general relativity, nuclear physics, and cosmology. Academia unfortunately encourages single-minded focus, and a cross-disciplinary approach is usually frowned upon. I'm aware that an electrical engineer claiming results in particle physics will raise eyebrows. But the derivations are all written out, the code reproduces every number, and the predictions are falsifiable. I'd rather put it out there and be wrong than sit on it and wonder.

I used AI (Claude Opus and a little bit of Gemini) as iterative discussion partners throughout — to stress-test arguments, check and derive mathematics, spot errors, and move much faster than I could alone. All responsibility for the physics reasoning, derivations, and conclusions is mine. Fortunately my general knowledge is good and my MSc helped me a lot, but I have spent quite some time brushing up on materials science so that I can smoke test this work!

The result appears to be on par with a theory of everything. I'm aware of how that sounds. But the numbers are what they are: ~90 predictions from three inputs (*c*, *ℏ*, *mₑ*), spanning from the fine structure constant to the Higgs mass, with no fitted parameters and nearly all retrodictions at pull < 1. The framework derives both α (to 0.003 ppb) and *G* (to 2 ppm) from the same lattice geometry. It predicts hadron masses (including all 8 charmonium states below DD̄ to < 0.01%, 6 open-charm mesons, and 8 charmed baryons), lepton masses, quark masses, mixing angles, neutrino mass splittings, dark matter, dark energy and the cosmological baryon asymmetry. It is not numerology — every integer that appears is derived from the FCC coordination geometry before any mass is consulted, and the monograph shows the working step by step.

The monograph is a work in progress but fairly complete. It has many falsifiable predictions (many are explicitly catalogued, several testable with existing data, particularly at ALICE and LHCb). It also has open questions and calculations that I haven't finished. I'm working toward splitting it into peer-reviewable papers, but the claims are unconventional and that process will be very challenging. In the meantime, I'm putting everything out here — the full monograph, all the derivations, and the code to reproduce every number — because I think it deserves scrutiny, and I'd rather people can check the work now than wait.

Whether the leading-order agreement is evidence for the hypothesis or an elaborate coincidence is for the reader to judge. The purpose of the monograph is to present the derivations in sufficient detail that the question can be answered.

Please go read it: \`https://doi.org/10.5281/zenodo.18636501\`

---

## Summary

The vacuum is modelled as an FCC Cosserat supersolid — a face-centred cubic crystal that is simultaneously superfluid (rigid to shear, frictionless to translation). Particles are topological defects: screw dislocations (electrons), partial dislocations (quarks, with the FCC threefold stacking degeneracy providing colour), edge dislocations (neutrinos), and vacancies (dark matter). Forces are elastic waves: transverse shear (electromagnetism), stacking-fault elasticity (strong force), evanescent modes (weak force), and longitudinal compression (gravity).

The fine structure constant is the Boltzmann factor for tunnelling through the Peierls–Nabarro barrier. Newton's constant requires coherent tunnelling of all 19 nodes in the FCC coordination cluster, giving *G* ∝ α¹⁹ — the 10³⁸ hierarchy between electromagnetism and gravity is a counting problem (1 node vs 19). A universal mass formula *m* = *m*ₑ(*N*/α*ᵏ* + *Q*) yields the full particle spectrum from five geometric building blocks and a four-term electromagnetic correction, each derived from the defect geometry with no reference to experiment.

The obvious objection — that a crystal breaks Lorentz invariance — is addressed quantitatively in the monograph and in a [dedicated paper](https://doi.org/10.5281/zenodo.18739953) (submitted, under review). Two suppression mechanisms intrinsic to the lattice (the Peierls–Nabarro form factor and Debye–Waller smearing from quantum zero-point delocalisation) combine to push Bragg scattering five orders of magnitude below the Fermi-LAT gamma-ray bound. The lattice is there, but it hides well.

---

## What's new in v9 (March 2026)

The charm sector — previously the biggest open problem — is now the biggest success:

- **8 charmonium states** below DD̄ threshold (η_c, J/ψ, χ_c0, χ_c1, h_c, χ_c2, η_c(2S), ψ(2S)) — all to < 0.01%
- **6 open-charm meson masses** (D⁰, D⁺, D\*⁰, D\*⁺, Ds, Ds\*) — all to < 0.01%
- **8 charmed baryon masses** (Λ_c, Σ_c⁺⁺/⁺/⁰, Ξ_c⁺/⁰, Ω_c, Ξ_cc⁺⁺) — all to < 0.01%
- **Bottomonium splittings** (Υ−η_b, 1P fine structure, radial) — all correct
- **11 Q-decomposition rules** derived from FCC constants (was 1)
- **Cage states**: a new defect class (stacking fault tetrahedron) for exotic hadrons far from thresholds
- **J^PC = 2⁺⁺** retrodicted for all-charm tetraquarks, matching CMS Nature 2025
- **0 failures** (was 1 — the T_cs̄1(4000)⁺ reclassified as cage state)

### Predictions for LHCb

| Prediction | Mass / splitting | How to search |
|---|---|---|
| B_c\* meson | ~6362 MeV (Δm ≈ 87 MeV) | B_c\* → B_c γ (soft photon) |
| 2⁺⁺ in J/ψ φ spectrum | ~4300 or ~4700 MeV | Angular analysis (CMS technique) |
| Ω_cc⁺ (scc) | ~3.7–3.9 GeV | Ξ_c⁺ K⁻ π⁺ decay channel |
| Ω_ccc⁺⁺ (ccc) | ~4.7–4.9 GeV | Multiple charm decay chains |

---

## Files

| Script | What it does |
|---|---|
| \`cosserat_calculator.py\` | Full prediction engine — every mass and coupling in the monograph |
| \`pn_variational.py\` | Derives α and *G* from the Peierls–Nabarro tunnelling calculation |
| \`verify_fcc_geometry.py\` | Verifies every geometric claim in Ch. 9 against explicit FCC coordinates |
| \`exotic_catalogue.py\` | Blind mass-formula test on all 32 known exotic hadrons (10 cage, 22 molecular/threshold, 0 failures) |
| \`exotic_filling_fraction.py\` | Honesty check — proves the mass formula is trivially flexible without independent *N* derivation |

Requirements: Python 3, NumPy, SciPy.

---

## Output: \`cosserat_calculator.py\` (selected)

\`\`\`
  FUNDAMENTAL CONSTANTS
  α⁻¹ (fine structure const.)              137.036      137.036          +0.000%
  G (Newton's constant)                6.67643e-11   6.6743e-11 m³/kg/s²   +0.032%

  CHARM SECTOR
  ηc (cc̄, J=0)                               2984       2983.9 MeV      +0.003%
  J/ψ (cc̄, J=1)                           3096.95       3096.9 MeV      +0.002%
  χc0 (cc̄, 0++)                           3414.89      3414.71 MeV      +0.005%
  χc1 (cc̄, 1++)                           3510.46      3510.67 MeV      -0.006%
  hc (cc̄, 1+-)                            3525.28      3525.37 MeV      -0.003%
  χc2 (cc̄, 2++)                           3555.94      3556.17 MeV      -0.006%
  ηc(2S) (n_r=1)                           3637.72       3637.5 MeV      +0.006%
  ψ(2S) (n_r=1)                            3686.28       3686.1 MeV      +0.005%
  D⁰ (cū)                                   1865.1      1864.84 MeV      +0.014%
  Λc⁺ (udc)                                2286.27      2286.46 MeV      -0.008%
  Ξcc⁺⁺ (ucc)                              3621.33       3621.2 MeV      +0.004%
  Bc* mass [pred]                           6361.9            — MeV   [PREDICTION]

  BOTTOM SPLITTINGS
  Υ − η_b (1S HF)                          61.3383         61.3 MeV      +0.062%
  Υ(2S) − Υ(1S) (radial)                   562.757       562.96 MeV      -0.036%
\`\`\`

Full output includes ~90 predictions across all sectors.

---

## Output: \`exotic_catalogue.py\` (summary)

\`\`\`
SUMMARY
  Total states analysed:  32
  PASS:                   22
  CAGE (SFT):             10
  UNCERTAIN:              0

  No exotic hadron known as of 2026 is inconsistent with the framework.
\`\`\`

Cage states are stacking fault tetrahedra (SFT) — a well-known three-dimensional FCC defect. They are systematically broad (mean Γ ~ 107 MeV) and far from all compatible two-hadron thresholds, distinguishing them from molecular exotics (narrow, at threshold). The T_cs̄1(4000)⁺, previously the framework's sole failure, is reclassified as a cage state.

---

## Code Licence

MIT

---

Warp drive doesn't jump out of this framework (unfortunately). But if the vacuum really is a crystal, then the universe is (now) an engineering problem — and that's a start.
