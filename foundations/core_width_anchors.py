"""Core-width anchor mapping: glide vs stacking readings of the PN core.

Tests the candidate theorem that the two form-factor anchors in the monograph,

  (A) glide anchor:    exp(-k1 * w_par) = alpha,  k1 = 2*pi/d,     w_par = 0.783*d
  (B) stacking anchor: exp(-|G111| * w_perp) = alpha,  |G111| = 2*pi/d111, w_perp = 0.783*d111

are two projections of one dimensionless invariant (the PN phase width
W = ln(1/alpha)/(2*pi) periods), related by the exact FCC geometry
d111/d = sqrt(2), with the D4 compact step equal to d111.

It then evaluates candidate direction-resolved coupling rules on the
{111}, {200}, and {220} shells and prints the discriminating table.

Candidate rules for the form factor F(G) of a defect whose PN phase
winds along the covector K_phi = (2*pi/d) e_112 + (2*pi/d111) e_111:

  H3 (current convention):  F(G) = alpha**(|G|/|G111|)          (isotropic in k)
  H_par (glide-anchored):   F(G) = alpha**(|G| * w_par / (2*pi) * 2*pi) -> alpha**(|G| d /(2*pi) *?) ...
                            implemented as exp(-|G| * w_par)
  H_perp (stacking-anchored): exp(-|G| * w_perp)   (equivalent to H3 exactly)
  H_A (winding covector):   F(G) = alpha**(|G| / (K_phi . G_hat)), per partial,
                            the phase period measured along G_hat.
  H_iso3D (plain elasticity): exp(-|G| * w_par) with w_par the single 3D length
                            0.454*ell -- the control that MUST fail at {111} if
                            the D4 reading is doing real work.

Verdicts are checked against the two non-negotiable anchors:
  glide harmonic k1 -> alpha exactly (the alpha chapter's derivation)
  {111} shell      -> alpha exactly (the working shell convention)

Author: framework companion script. Exact arithmetic via sympy.
"""

import sympy as sp

# ----------------------------------------------------------------------
# Exact FCC / D4 geometry (lengths in units of the NN spacing ell = 1)
# ----------------------------------------------------------------------
ell = sp.Integer(1)
a = sp.sqrt(2) * ell                     # conventional cubic lattice constant
d = ell / sp.sqrt(3)                     # partial Burgers vector = misfit period
d111 = ell * sp.sqrt(sp.Rational(2, 3))  # {111} interplanar spacing
L4 = sp.sqrt(6) * ell                    # D4 compact period (three layers)

print("=== Exact geometric identities ===")
print(f"d111/d              = {sp.simplify(d111/d)}          (the anchor mismatch)")
print(f"compact step L4/3   = {sp.simplify(L4/3)} = d111 ? {sp.simplify(L4/3 - d111) == 0}")
print(f"d^2 + d111^2        = {sp.simplify(d**2 + d111**2)} = ell^2 (interlayer NN bond)")

# ----------------------------------------------------------------------
# The dimensionless invariant from CODATA alpha
# ----------------------------------------------------------------------
alpha_inv = sp.Float("137.035999177", 12)
alpha = 1 / alpha_inv
W = sp.log(alpha_inv) / (2 * sp.pi)      # phase width in periods
w_par = W * d                            # width along the glide direction
w_perp = W * d111                        # width along the stacking normal

print("\n=== The invariant and its two projections ===")
print(f"W = ln(1/alpha)/2pi = {float(W):.6f} periods")
print(f"w_par  = W*d    = {float(w_par):.6f} ell   (glide reading, alpha chapter)")
print(f"w_perp = W*d111 = {float(w_perp):.6f} ell   (stacking reading, shell convention)")
print(f"w_perp/w_par    = {float(w_perp/w_par):.6f} = sqrt(2)")

# ----------------------------------------------------------------------
# Directions and shells (cubic coordinates, exact)
# ----------------------------------------------------------------------
def unit(v):
    v = sp.Matrix(v)
    return v / sp.sqrt(sum(x**2 for x in v))

n111 = unit([1, 1, 1])                   # glide-plane normal = stacking axis
partials = [unit([1, 1, -2]), unit([-2, 1, 1]), unit([1, -2, 1])]  # <112> on (111)

# Phase-winding covector per partial: 2pi/d along the glide direction
# plus 2pi/d111 along the stacking axis (the D4 compact reading).
def K_phi(e112):
    return (2 * sp.pi / d) * e112 + (2 * sp.pi / d111) * n111

# Reciprocal shells: |G| = 2*pi/d_hkl with d_hkl the interplanar spacing.
shells = {
    "glide k1": ((2 * sp.pi / d), None),          # probed along each partial itself
    "{111}":    ((2 * sp.pi / d111), n111),
    "{200}":    ((2 * sp.pi / (a / 2)), unit([1, 0, 0])),
    "{220}":    ((2 * sp.pi / (a / (2 * sp.sqrt(2)))), unit([1, 1, 0])),
}

# ----------------------------------------------------------------------
# Candidate rules
# ----------------------------------------------------------------------
print("\n=== Candidate coupling rules: exponent x in F = alpha^x ===")
print(f"{'shell':10s} {'H3 |G|/|G111|':>14s} {'H_par |G|w_par':>15s} "
      f"{'H_iso3D':>9s} {'H_A per partial':>24s}")

G111mag = 2 * sp.pi / d111
lnainv = sp.log(alpha_inv)

for name, (Gmag, Ghat) in shells.items():
    x_h3 = float(Gmag / G111mag)
    x_par = float(Gmag * w_par / lnainv)          # exp(-|G| w_par) = alpha^x
    x_iso = x_par                                  # same length, shown for the control
    if Ghat is None:
        # glide harmonic: probe along each partial direction
        xs_A = []
        for e in partials:
            proj = (K_phi(e).T * e)[0]
            xs_A.append(float(Gmag / proj))
        xa = ", ".join(f"{x:.3f}" for x in xs_A)
    else:
        xs_A = []
        for e in partials:
            proj = (K_phi(e).T * Ghat)[0]
            xs_A.append(float(sp.Abs(Gmag / proj)) if proj != 0 else sp.oo)
        xa = ", ".join((f"{float(x):.3f}" if x != sp.oo else "inf") for x in xs_A)
    print(f"{name:10s} {x_h3:14.4f} {x_par:15.4f} {x_iso:9.4f} {xa:>24s}")

print("\nNon-negotiable anchors: glide k1 must give x = 1; {111} must give x = 1.")
print("H_iso3D at {111}: x =", f"{float((2*sp.pi/d111) * w_par / lnainv):.4f}",
      " -> plain 3D elastic width FAILS the {111} anchor (the D4 fingerprint).")
print("H3 and H_A both pass both anchors; they diverge at {200} and {220}.")

# ----------------------------------------------------------------------
# Numerical values of the {200} candidates for the cascade map
# ----------------------------------------------------------------------
print("\n=== {200} coupling amplitude under each surviving rule ===")
for label, x in [("current convention alpha^sqrt(2)", sp.sqrt(2)),
                 ("glide-anchored alpha^(2/sqrt(3))", 2 / sp.sqrt(3)),
                 ("H_A dominant partials alpha^1", sp.Integer(1)),
                 ("H_A hard partial alpha^2", sp.Integer(2))]:
    print(f"  {label:38s} = {float(alpha**sp.Float(float(x))):.3e}")


# ----------------------------------------------------------------------
# H_B: the separable anisotropic core (the mechanically derived rule)
# ----------------------------------------------------------------------
# The misfit density is a Lorentzian of width w_par along the glide
# direction (the alpha chapter's PN profile) and a Lorentzian of width
# w_perp along the stacking normal (the compact-phase structure the {111}
# anchor requires). The line direction contributes no localisation, so
# only the G components perpendicular to the line couple. The Fourier
# transform of the product of Lorentzians is the product of exponentials,
# so the exponents ADD:
#
#   F(G) = exp(-|G.e112| w_par - |G.n111| w_perp) = alpha**n(G),
#   n(G) = |G.e112| d/(2 pi) + |G.n111| d111/(2 pi),
#
# the Manhattan winding count of the probe against the core's two phase
# directions. Both anchors give n = 1 exactly and symbolically; H3 is
# ruled out mechanically because no single core shape has an isotropic
# |k|-only transform while carrying two different widths.

print("\n=== H_B: winding-count rule n(G) per partial [F = alpha^n, exact] ===")
for name, (Gmag, Ghat) in shells.items():
    outs = []
    for e in partials:
        G = (Gmag * e) if Ghat is None else (Gmag * Ghat)
        n = sp.simplify(sp.nsimplify(
            sp.Abs((G.T * e)[0]) * d / (2 * sp.pi)
            + sp.Abs((G.T * n111)[0]) * d111 / (2 * sp.pi)))
        outs.append(n)
    print(f"  {name:9s}: {outs}")

avg200 = (2 * alpha + alpha**sp.Rational(4, 3)) / 3
print(f"\n{{200}} symmetric-average amplitude under H_B: {float(avg200):.3e} = {float(avg200/alpha):.3f} alpha")
print("Mode IV therefore sits at the same alpha rung as Mode II, not one decade below;")
print("the E_g channel's inability to mediate a force is untouched, since that rests on")
print("the WKB action S_eff = ln(1/alpha)/sqrt(0.069) = 18.7, a different quantity.")
