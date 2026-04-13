# D4 Compact Direction Analysis

Scripts for the D4 (four-dimensional) interpretation of the FCC Cosserat lattice.

## Scripts

### `layer_entropy_catalogue.py`

Computes the D4 stacking-layer decomposition for every hadron in the constructive
Cosserat mass calculator. Assigns each node of each defect graph to one of three
stacking layers (A, B, C), classifies bonds as spatial or temporal, computes the
layer entropy and variance, and tests the prediction that the mass-formula residual
correlates with temporal asymmetry.

**Key results:**
- Every baryon graph fragments when temporal bonds are removed (confinement = temporal binding)
- The decuplet baryons (Δ, Σ*, Ξ*, Ω⁻) show perfect rank correlation (ρ = 1.000)
  between layer entropy and prediction accuracy
- The Ω⁻ (H = 1.000, perfectly balanced across all three layers) has the smallest
  mass residual of any light baryon (0.001%)

**Usage:**
```
python layer_entropy_catalogue.py          # full analysis
python layer_entropy_catalogue.py --plot   # also generate scatter plot
```

**Dependencies:** numpy, scipy, matplotlib (optional), cosserat_graph

## Reference

M. A. Cox, *The Cosserat Supersolid*, Sec. layer_decomposition (decay chapter).
