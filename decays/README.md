# Decays

Everything that computes a decay width or a stability classification lives here.
If you are looking for a lifetime, a branching pattern, or "why is this particle
stable", start in this folder, not in `spectral_mass/` and not at the repo root.

## Start here (authoritative)

| File | What it does | Accuracy |
|------|--------------|----------|
| `cosserat_decay_engine.py` | Universal engine. Parses parent+daughters to a defect graph, takes literal master-formula mode sums, selects one of five structural combination rules by graph topology, multiplies by relativistic phase space. Derives `g_piNN = N_H = 13` and the decuplet widths with **no fitted coupling**. | Δ −3.7%, Σ* −2.4%, ρ +1.1%, π⁰→γγ +0.5%; 31/38 channels within 15% |
| `test_decay_engine.py` | Regression table, 38 channels against PDG. | run it to check |
| `cosserat_decay.py` | Tier 1/2 graph engine: Laplacian/Fiedler stability classification and the factorisation theorem Γ = 2(\|∂E\|mₑ + \|∂V\|m₀/π), plus decuplet void-deactivation and crossed-fault P-wave healing. | the structural backbone |

The **derivations** behind these numbers are in `monograph_decays.tex` (carved out
of `monograph_dynamics.tex`): `g_piNN = N_H = 13` at `sec:gA_derivation`, the Δ
width at `eq:delta_ab_initio`, the strange extension at `eq:decuplet_ab_initio_NLO`,
and the closed-form cubic vertices at `sec:vertex_spectral_synthesis`.

## exploratory/ (session work, superseded — kept for the record)

These are explorations from the synthesis sessions. They are **not** the
production path; their absolute scale is less accurate than the engine above.
Each carries a SUPERSEDED banner. Kept so the reasoning is not lost.

| File | What it explored | Why superseded |
|------|------------------|----------------|
| `synth_test.py` | Decay as parent→daughter spectral-mode overlap across the Fiedler cut; proton stability by absence of a lighter daughter. | Structural test only; the engine does this and gives absolute widths. |
| `cosserat_vertex.py` | Rotation (φφu) cubic vertex for the baryon sector; on-shell pion projection (ℓ=r_e), P-wave selection rule. | Analytic vertex; engine uses the literal node sum, which is calibration-free. |
| `gamma3_vertex.py` | Displacement cubic vertex; ρ→ππ P-wave selection rule. | Same — analytic counterpart of the engine's sum. |
| `decay_width.py` | Outgoing-pion projection at p* via ℓ = r_e. | Folded into the engine's phase space. |
| `tier2_complex_mass.py` | Mass = Re, width = Im of one complex eigenvalue. | Channel-blind widths; conceptual only. |
| `strong_widths.py` | n_free² coherence model for the decuplet. | Inferior to the engine's (λ₂/λ₂^Δ)(1−\|S\|σ) form. |
| `mode_normalised_width.py` | Absolute scale from mode zero-point amplitudes; finds kappa3 ≈ N mₑ but lands Δ at 62 MeV (−47%). | A factor ~1.9 below the engine's −3.7%. Documents the normalisation question. |
| `fiedler_check.py` | Scratch check of void-swap λ₂. | Numbers now live in the engine's `_void_swap_lambda2`. |

## One-line map

```
decays/
  cosserat_decay_engine.py   <- production: widths + g_piNN=13, no fit
  test_decay_engine.py       <- 38-channel regression
  cosserat_decay.py          <- Tier 1/2: stability + factorisation
  exploratory/               <- superseded session scripts (see table)
```
