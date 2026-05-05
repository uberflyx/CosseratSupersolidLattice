# Cosserat Supersolid Lattice

## Repository structure

```
cosserat_calculator.py          # Legacy mass calculator (lookup-table version)
cosserat_graph.py               # ★ Constructive graph calculator — all masses from FCC edge counting
├── foundations/                 # Ch. 3–6: PN barrier, α derivation, symmetry channels
│   ├── pn_variational.py       #   PN tunnelling amplitude (Ch. 5)
│   ├── theta_ch_ab_initio_v2.py#   Chirality parameter θ_ch (Ch. 6)
│   ├── oh_irrep_overlaps.py    #   O_h irrep overlaps (Ch. 4)
│   ├── cosserat_transfer_matrix.py # Transfer matrix eigenvalues (Ch. 4)
│   ├── cosserat_dw.py          #   Debye-Waller factor (Ch. 3)
│   └── gamma4_scp.py           #   Quartic anharmonic correction (Ch. 11)
├── hadrons/                    # Ch. 9–10: mass formula, exotics, excited baryons
│   ├── verify_fcc_geometry.py  #   FCC coordinate checks (Ch. 9)
│   ├── excited_baryons.py      #   N(1440), N(1535), Λ(1670) masses (Ch. 12)
│   ├── exotic_catalogue.py     #   Exotic hadron classification
│   ├── exotic_filling_fraction.py
│   ├── three_mechanisms.py     #   Vacuum energy mechanisms (Ch. 13)
│   ├── probe_cosserat_graph.py #   PDG validation harness for cosserat_graph
│   ├── pdg_comparison.md       #   Full prediction-vs-PDG comparison report
│   └── fcc_defect_catalogue.py #   Composition catalogue of FCC defect cores
├── decays/                     # Ch. 12: decay rates and lifetimes
│   ├── cosserat_decay_engine.py#   Full decay rate calculator
│   └── test_decay_engine.py    #   Regression tests
├── neutrinos/                  # Ch. 8: neutrino sector
│   ├── neutrino_predictions.py #   ν mass and mixing predictions
│   └── neutrino_error_budget.py#   Error propagation
├── gravity/                    # Ch. 15–16: black holes, gravity
│   ├── bh_entropy_derivation.py#   Bekenstein-Hawking entropy
│   └── bh_statistical_mechanics.py # BH statistical mechanics
└── d4/                         # Sec. 2.8: four-dimensional interpretation
    ├── d4_lattice_sums.py      #   D3 vs D4 elastic lattice sums
    ├── d4_mu_prime.py          #   μ' = 2 self-consistency: KK reduction
    └── d4_scales.py            #   Temperature hierarchy, G × K_sf = c⁴/ℓ²
```


**Mitchell A. Cox, University of the Witwatersrand**

Companion code for [*The Cosserat Supersolid*](https://doi.org/10.5281/zenodo.18636501)

Full monograph: [https://doi.org/10.5281/zenodo.18636501](https://doi.org/10.5281/zenodo.18636501)

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

Please go read it: [https://doi.org/10.5281/zenodo.18636501](https://doi.org/10.5281/zenodo.18636501)

---

## Summary

The vacuum is modelled as an FCC Cosserat supersolid — a face-centred cubic crystal that is simultaneously superfluid (rigid to shear, frictionless to translation). Particles are topological defects: screw dislocations (electrons), partial dislocations (quarks, with the FCC threefold stacking degeneracy providing colour), edge dislocations (neutrinos), and vacancies (dark matter). Forces are elastic waves: transverse shear (electromagnetism), stacking-fault elasticity (strong force), evanescent modes (weak force), and longitudinal compression (gravity).

The fine structure constant is the Boltzmann factor for tunnelling through the Peierls–Nabarro barrier. Newton's constant requires coherent tunnelling of all 19 nodes in the FCC coordination cluster, giving *G* ∝ α¹⁹ — the 10³⁸ hierarchy between electromagnetism and gravity is a counting problem (1 node vs 19). A universal mass formula *m* = *m*ₑ(*N*/α*ᵏ* + *Q*) yields the full particle spectrum from geometric building blocks, each derived from the defect geometry with no reference to experiment.

The obvious objection — that a crystal breaks Lorentz invariance — is addressed quantitatively in the monograph and in a [dedicated paper](https://doi.org/10.5281/zenodo.18739953) (submitted, under review). Two suppression mechanisms intrinsic to the lattice (the Peierls–Nabarro form factor and Debye–Waller smearing from quantum zero-point delocalisation) combine to push Bragg scattering five orders of magnitude below the Fermi-LAT gamma-ray bound. The lattice is there, but it hides well.

---

## Code Licence

MIT

---

Warp drive doesn't jump out of this framework (unfortunately). But if the vacuum really is a crystal, then the universe is (now) an engineering problem — and that's a start. 

For the brave: I was wondering if the Lattice supports some sort of long-range propagation for information — like electromagnetism and gravity, but not those. Kinks on lattice dislocation lines might be something, with heavy nuclei as "transducers" for the interface.
