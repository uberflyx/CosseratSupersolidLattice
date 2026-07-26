# Cosserat Supersolid Lattice

**Mitchell A. Cox, University of the Witwatersrand**

Companion code for the monograph *The Cosserat Supersolid*. The scripts verify;
they do not derive. Every result is worked through in full in the monograph, so
a reader can rebuild any of it from the prose alone, and these files exist to
make that check cheap and to guard the numbers against regression.

**Monograph: [doi.org/10.5281/zenodo.18636501](https://doi.org/10.5281/zenodo.18636501)**

---

## Papers

Results are drawn out of the monograph into stand-alone papers as they reach a
state fit for external review.

**Published**

- M. A. Cox, "The Peierls–Nabarro amplitude as an emergent coupling constant: an
  explicit evaluation for the FCC Cosserat lattice", *Annals of Physics* **493**,
  170613 (2026). Open access:
  [doi.org/10.1016/j.aop.2026.170613](https://doi.org/10.1016/j.aop.2026.170613)
  → reproduced by [`pn_variational.py`](pn_variational.py), which computes the fine
  structure constant from the shear channel and Newton's constant from the
  compression channel.

Further papers are in preparation or under peer review at any given time, and
they are added here once a decision is in. The monograph remains the complete
and current statement of the framework.

---

## Repository structure

Scripts are sorted by topic, with one exception kept at the root: the script
behind the published paper, so a reader arriving from it finds it immediately.

```
pn_variational.py               # ★ Published paper: α and G from the PN amplitude
cosserat_calculator.py          # Legacy mass calculator (lookup-table version)
├── foundations/                # PN barrier, α derivation, symmetry channels
│   ├── theta_ch_ab_initio_v2.py#   Chirality parameter θ_ch
│   ├── oh_irrep_overlaps.py    #   O_h irrep overlaps
│   ├── cosserat_transfer_matrix.py # Transfer matrix eigenvalues
│   ├── cosserat_dw.py          #   Debye–Waller factor
│   └── koide_sigma_mott.py     #   Koide amplitude and quartic from the Mott point
├── hadrons/                    # Mass formula, exotics, excited baryons
│   ├── verify_fcc_geometry.py  #   FCC coordinate checks
│   ├── exotic_catalogue.py     #   Exotic hadron classification
│   ├── exotic_filling_fraction.py
│   ├── three_mechanisms.py     #   Vacuum energy mechanisms
│   ├── pdg_comparison.md       #   Full prediction-vs-PDG comparison report
│   └── fcc_defect_catalogue.py #   Composition catalogue of FCC defect cores
├── spectral_mass/              # Hadron masses from the cluster elastic spectrum
│   ├── spectral_classifier.py  #   FCC cluster builders + O_h irrep decomposition
│   ├── cosserat_classifier.py  #   Full Cosserat dynamical matrix
│   ├── hadron_spectral_mass.py #   Mass law m = N m_0 − N(4−λ) m_e
│   ├── bare_shell_irreps.py    #   Bare 13-shell spectrum by O_h irrep (reference)
│   ├── parity_flip_rule.py     #   J^P → parent irrep map
│   ├── retrospective_adiabatic.py # Excited-baryon adiabatic audit
│   └── delta1600_dual_orbit.py #   Δ(1600) dual-orbit closure (+ supporting builders)
├── decays/                     # Decay rates and lifetimes
│   ├── cosserat_decay_engine.py#   Full decay rate calculator
│   └── test_decay_engine.py    #   Regression tests
├── neutrinos/                  # Neutrino masses, mixing, error budget
│   ├── neutrino_predictions.py #   ν mass and mixing predictions
│   └── neutrino_error_budget.py#   Error propagation
├── gravity/                    # Black holes and gravity
│   ├── bh_entropy_derivation.py#   Bekenstein–Hawking entropy
│   └── bh_statistical_mechanics.py # BH statistical mechanics
└── d4/                         # Four-dimensional interpretation
    ├── d4_lattice_sums.py      #   D3 vs D4 elastic lattice sums
    ├── d4_mu_prime.py          #   μ' = 2 self-consistency: KK reduction
    └── d4_scales.py            #   Temperature hierarchy, G × K_sf = c⁴/ℓ²
```

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

Please go read it: [doi.org/10.5281/zenodo.18636501](https://doi.org/10.5281/zenodo.18636501)

---

## Code Licence

MIT

---

Warp drive doesn't jump out of this framework (unfortunately). But if the vacuum really is a crystal, then the universe is (now) an engineering problem — and that's a start.

---

## A note on the commit history

The git log carries about fifteen different author names, and all of them are
me. Most are old machines and placeholder emails I never got round to fixing.
The rest come from commits fired off from a random chat on my phone, where
Claude, asked to write the commit, would confidently supply an author it had
invented on the spot. My favourite is *Mitch Kovari*, which has never been my
surname. A few commits are signed by Claude itself, which at least has the
virtue of honesty.

So: one person, several laptops, and a language model with rather more
enthusiasm than access to my passport. The physics was checked more carefully
than the metadata.
