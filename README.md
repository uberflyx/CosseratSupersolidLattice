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
│   └── three_mechanisms.py     #   Vacuum energy mechanisms (Ch. 13)
├── decays/                     # Ch. 12: decay rates and lifetimes
│   ├── cosserat_decay_engine.py#   Full decay rate calculator
│   └── test_decay_engine.py    #   Regression tests
├── neutrinos/                  # Ch. 8: neutrino sector
│   ├── neutrino_predictions.py #   ν mass and mixing predictions
│   └── neutrino_error_budget.py#   Error propagation
└── gravity/                    # Ch. 15–16: black holes, gravity
    ├── bh_entropy_derivation.py#   Bekenstein-Hawking entropy
    └── bh_statistical_mechanics.py # BH statistical mechanics
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

The fine structure constant is the Boltzmann factor for tunnelling through the Peierls–Nabarro barrier. Newton's constant requires coherent tunnelling of all 19 nodes in the FCC coordination cluster, giving *G* ∝ α¹⁹ — the 10³⁸ hierarchy between electromagnetism and gravity is a counting problem (1 node vs 19). A universal mass formula *m* = *m*ₑ(*N*/α*ᵏ* + *Q*) yields the full particle spectrum from five geometric building blocks and a four-term electromagnetic correction, each derived from the defect geometry with no reference to experiment.

The obvious objection — that a crystal breaks Lorentz invariance — is addressed quantitatively in the monograph and in a [dedicated paper](https://doi.org/10.5281/zenodo.18739953) (submitted, under review). Two suppression mechanisms intrinsic to the lattice (the Peierls–Nabarro form factor and Debye–Waller smearing from quantum zero-point delocalisation) combine to push Bragg scattering five orders of magnitude below the Fermi-LAT gamma-ray bound. The lattice is there, but it hides well.

---

## Code Licence

MIT

---

## `cosserat_graph.py` — Constructive Graph Calculator

The centrepiece of this repository. A 1099-line Python script that computes hadron masses from edge counting on the FCC cuboctahedral graph. No lookup tables, no fitted parameters — every integer in every mass formula is computed from loops over actual 3D lattice coordinates.

**Three inputs:** *c*, *ℏ*, *mₑ* (speed of light, Planck's constant, electron mass).

**One equation:** *m = N × m₀ + Q × mₑ*, where *m₀ = mₑ/α* and α is derived self-consistently from the Peierls–Nabarro tunnelling amplitude.

**One graph:** the 12-vertex cuboctahedron (first coordination shell of the FCC lattice), built from three primitive vectors and analysed with networkx.

### How it works

1. **Build the lattice.** The code generates ~340 FCC lattice sites from primitive vectors a₁=(1,1,0), a₂=(1,0,1), a₃=(0,1,1), identifies the 12-vertex cuboctahedral shell, and computes all graph invariants: coordination number Z₁=12, chromatic number N_c=3, edge count E=24, hex-cap size N_hex=7, etc.

2. **Assemble the defect.** Quantum numbers (B, S, I, J, P, n_charm) determine which FCC building blocks to activate: coordination shell (baryons), hex-cap stacking fault (kaons), crossed faults (vector mesons), etc. Each node gets actual 3D FCC coordinates.

3. **Count edges.** The electromagnetic correction Q is computed by looping over actual edges of the embedded defect: core bonds (centre → shell), boundary nodes (exposed surface), antipodal pairs (dislocation character), common-NN bonds (void modification), and colour-spatial cross-edges (charm coupling). The counting method is determined by the spin: J=0 counts vertices, J=1 counts edges, J≥2 counts faces — the simplicial decomposition of the Cosserat field on the graph.

4. **Return the mass.** *m = N × m₀ + Q × mₑ*, both integers from the graph.

### Quick start

```python
from cosserat_graph import predict, predict_molecular, QN

# Proton
r = predict(QN(B=1, S=0, I=0.5, I3=0.5, J=0.5, P=+1))
print(f"Proton: {r.mass:.1f} MeV (PDG: 938.3)")   # 939.2 MeV

# J/ψ meson
r = predict(QN(B=0, I=0, I3=0, J=1, P=-1, n_charm=2))
print(f"J/ψ: {r.mass:.1f} MeV (PDG: 3096.9)")     # 3097.0 MeV

# Ξ_cc⁺⁺ (doubly charmed baryon, discovered LHCb 2017)
r = predict(QN(B=1, S=0, I=0.5, I3=0.5, J=0.5, P=+1, n_charm=2))
print(f"Ξ_cc⁺⁺: {r.mass:.1f} MeV (PDG: 3621.2)")  # 3621.3 MeV

# Ξ_cc⁺ (discovered LHCb Moriond 2026)
r = predict(QN(B=1, S=0, I=0.5, I3=-0.5, J=0.5, P=+1, n_charm=2))
print(f"Ξ_cc⁺: {r.mass:.1f} MeV (LHCb: 3619.97)")  # 3620.3 MeV

# X(3872) molecular exotic — docking mode selected automatically
D0  = QN(B=0, S=0, I=0.5, I3=0.5, J=0, P=-1, n_charm=1)
D0s = QN(B=0, S=0, I=0.5, I3=0.5, J=1, P=-1, n_charm=1)
r = predict_molecular(D0, D0s)
print(f"X(3872): {r.mass:.1f} MeV (PDG: 3871.7)")  # 3871.8 MeV
```

### Key graph-theoretic results

| Result | Formula | Value | What it determines |
|--------|---------|-------|-------------------|
| Coordination number | Z₁ = len(shell) | 12 | Proton Q, base coupling |
| Chromatic number | χ(G) = greedy colouring | 3 | Number of QCD colours |
| Hyperfine coupling | E×χ(G) + Z₁ | 84 | J/ψ – η_c splitting (113 MeV) |
| Surplus-edge theorem | Remove edge, count induced | 27 = N_c³ | Dibaryon binding (83.8 MeV) |
| Charm coupling target | N_charm × (N_coord − N_c×|S|) | 235, 180, 98, 8 | All charm baryon Q values |
| Pauli void activation | Count identical light quarks + spin | bool | Σ vs Λ mass difference |

### Master table: 52/52 within 1%

The code reproduces every entry in the monograph's three-part master prediction table:
- 21 light hadrons (π through Ω⁻, d*(2380))
- 8 charmonium states (η_c, J/ψ, ψ(2S), η_c(2S), χ_c0, χ_c1, h_c, χ_c2)
- 6 open charm mesons (D⁰, D⁺, D*⁰, D_s, D_s*, B_c)
- 8 charm baryons (Λ_c, Σ_c, Σ_c*, Ξ_c, Ω_c, Ξ_cc⁺⁺, Ξ_cc⁺, Ω_ccc)
- 4 molecular exotics (X(3872), T_cc⁺, P_c(4457), P_c(4440))
- 4 bottom splittings (Υ−η_b, B*−B, Bs*−Bs, Υ(2S)−Υ(1S))
- 1 cage tetraquark (T_4c(6600) from stacking-fault tetrahedron)

Worst residual: proton at 0.10%. Best: Ξ_cc⁺ at 0.009%.

### Dependencies

Python 3.8+, numpy, networkx. No other packages required.

```bash
pip install numpy networkx
python cosserat_graph.py
```

---

Warp drive doesn't jump out of this framework (unfortunately). But if the vacuum really is a crystal, then the universe is (now) an engineering problem — and that's a start.
