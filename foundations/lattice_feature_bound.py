"""The dressed defect propagator, and the lattice feature it permits in alpha(q^2).

THE QUESTION.  The paper assumes the dressed defect propagator keeps its
continuum Dirac form past the zone scale and flags the assumption as underived
(sec:status).  The monograph has in fact already computed the propagator on the
periodic background at nearly-free order (ch05, eq-nearly-free-band): the
defect's dispersion is a Bloch band with an avoided crossing at the glide-zone
boundary, and the barrier that opens it is dressed by the medium's zero-point
motion, V_eff = sqrt(S) V_0 (eq-suppressed-barrier), the same convention the
antiprotonic-helium check validated for the proton (374 MeV -> ~11 keV,
five orders inside the Hori bound).  What no one has done is push that band
structure through the vacuum-polarisation loop to get the FEATURE SIZE in
alpha(q^2).  That is this calculation.  Everything upstream is an established
monograph result; the loop propagation is the one new step.

THE SETUP.  The PN potential is periodic along each glide direction with the
misfit period d = l/sqrt(3), so the defect sees Bragg planes at

    p_G = hbar G_P / 2,   G_P = 2 pi / d   ->   p_G c = pi sqrt(3) m_0 c^2.

FCC has 4 {111} planes x 3 <112> partial directions = 12 misfit wavevectors.
At each Bragg plane an avoided crossing opens of size

    Delta = 2 |V_G^eff| = sqrt(S) V_0        (V(x) = V_0 cos(G_P x), V_G = V_0/2)

with V_0 = m_0 c^2 the electron's ceiling barrier (ch05) and S = 8.0e-10 the
condensate-free floor.  The floor is a lower bound on the suppression, so
Delta is an UPPER bound on the gap; the full supersolid suppression makes it
smaller still.  The coherent (gap-opening) part of the potential is the
Debye-Waller-dressed one whatever the probe's speed, because gap formation
needs the periodic average; the fluctuating remainder scatters diffusively and
is the friction question sec:proton settles separately.

THE LOOP.  A gap Delta at |p.n| = p_G modifies pair production gamma* -> e+e-
at sqrt(s) ~ 2 E_G.  Phase space: on the sphere |p| = p(s) the measure in p_par
is uniform, d(cos) = dp_par/p, and the crossing disturbs the band
|E(p) - E(p - 2 p_G n)| < Delta, whose half-width in p_par is
Delta E/(4 c^2 p_G).  Per wavevector the sphere cuts both +/- planes, so the
disturbed solid-angle fraction at s = 4 E_G^2 (where p ~ p_G ~ E_G/c) is

    f_1 = 2 x [2 x Delta E/(4 c^2 p_G)] / (2p) = Delta/(2 E_G),
    f   = 12 f_1 = 6 Delta / E_G.

Within that fraction the pair energies shift by O(Delta), so spectral weight
R(s) f is redistributed over a window delta(s) = 4 E_G Delta; the zeroth moment
of delta-R vanishes (states are shifted, not created), which is why every
smooth observable responds at second order.  The magnitude of the local
feature in the running is bounded by the total displaced weight:

    |delta Dalpha| <= (alpha/3pi) x R(s_G) x f x (delta s / s_G)
                    = (alpha/3pi) x R x 24 (Delta/2E_G)^2 x ... -> see below,

computed exactly with the numbers rather than in the head.
"""
import math

# ---------------------------------------------------------------- constants
ALPHA_INV = 137.035999177
ALPHA = 1.0 / ALPHA_INV
M_E = 0.51099895069            # MeV
M_0 = M_E * ALPHA_INV          # node mass, MeV
S_FLOOR = 8.0e-10              # condensate-free suppression (intensity)
N_GLIDE = 12                   # 4 {111} planes x 3 <112> partials

line = lambda t: print("\n" + t + "\n" + "=" * len(t))

line("SCALES, ALL FROM ESTABLISHED MONOGRAPH RESULTS")
V0 = M_0                                        # electron ceiling barrier, ch05
E_G = math.pi * math.sqrt(3.0) * M_0            # glide zone edge, ch09: 381 MeV
sqrtS = math.sqrt(S_FLOOR)
Delta = sqrtS * V0                              # dressed gap (upper bound)
print(f"  m_0 c^2              = {M_0:.3f} MeV")
print(f"  glide zone edge E_G  = pi sqrt(3) m_0 = {E_G:.1f} MeV")
print(f"  pair threshold        sqrt(s_G) = 2 E_G = {2*E_G:.0f} MeV")
print(f"  bare barrier V_0     = {V0:.2f} MeV")
print(f"  sqrt(S_floor)        = {sqrtS:.3e}")
print(f"  dressed gap Delta    = sqrt(S) V_0 = {Delta*1e3:.2f} keV   (upper bound)")
print(f"  proton check (374 MeV -> {sqrtS*374.0*1e3:.0f} keV) reproduces ch05's 11 keV")

line("THE DISTURBED PHASE-SPACE FRACTION AND WINDOW")
f1 = Delta / (2.0 * E_G)
f = N_GLIDE * f1
ds_over_s = Delta / E_G          # (4 E_G Delta)/(4 E_G^2)
print(f"  per wavevector  f_1 = Delta/2E_G = {f1:.3e}")
print(f"  all {N_GLIDE} wavevectors f = {f:.3e}")
print(f"  window width    ds/s = Delta/E_G = {ds_over_s:.3e}")

line("THE FEATURE IN THE RUNNING (upper bound, weight-displacement)")
# R for a lepton pair near s_G: beta(3-beta^2)/2 with beta = sqrt(1-4m^2/s) ~ 1
s_G = (2.0 * E_G) ** 2
beta = math.sqrt(1.0 - 4.0 * M_E ** 2 / s_G)
R_lep = 0.5 * beta * (3.0 - beta ** 2)
dDalpha = (ALPHA / (3.0 * math.pi)) * R_lep * f * ds_over_s
print(f"  R_e(s_G) = {R_lep:.6f}")
print(f"  |delta Dalpha|_max <= (alpha/3pi) R f (ds/s) = {dDalpha:.2e}")
print(f"  relative to Dalpha_had = 0.0275     : {dDalpha/0.0275:.1e}")
print(f"  relative to MUonE design (0.3% of 8.33e-4 = 2.5e-6): "
      f"{dDalpha/2.5e-6:.1e}")

line("THE SAME BOUND UNDER WEAKER AND STRONGER DRESSINGS")
for nm, Sval in (("bare barrier, S = 1 (excluded by antiprotonic He)", 1.0),
                 ("form-factor only, S = |F(G)|^2 = 9.5e-4", 9.5e-4),
                 ("condensate-free floor, S = 8.0e-10  <- adopted", 8.0e-10)):
    D = math.sqrt(Sval) * V0
    ff = N_GLIDE * D / (2.0 * E_G)
    dd = (ALPHA / (3.0 * math.pi)) * R_lep * ff * (D / E_G)
    perturbative = "  [non-perturbative, order-of-magnitude only]" if D / E_G > 0.05 else ""
    print(f"  {nm}")
    print(f"      Delta = {D:.3e} MeV   |dDalpha| <= {dd:.1e}{perturbative}")

line("WHERE THE FEATURE SITS, AND THE ON-SHELL COMPANION")
print(f"  timelike: pair threshold feature at sqrt(s) = {2*E_G:.0f} MeV,")
print(f"            width in sqrt(s) ~ Delta = {Delta*1e3:.1f} keV")
print(f"  spacelike (MUonE): the dispersion integral smears the window,")
print(f"            so |delta Dalpha(t)| < the bound above at every t")
print(f"  on-shell: the electron's own dispersion carries the avoided")
print(f"            crossing, a {Delta*1e3:.1f} keV anomaly at p = {E_G:.0f} MeV/c")
print(f"            in 12 lattice-fixed directions: fractional {Delta/E_G:.1e},")
print(f"            far below any beam-energy scan's resolution")

line("REFINEMENTS DELIBERATELY NOT TAKEN")
mask = M_E / E_G
print(f"  1. vector-coupled Bragg backscattering of an ultrarelativistic Dirac")
print(f"     particle is chirality-suppressed by m_e/E_G = {mask:.1e}")
print(f"     (the graphene/Klein mechanism); taking it would shrink Delta to")
print(f"     {Delta*mask*1e6:.1f} eV and the feature to {dDalpha*mask**2:.0e}.")
print(f"     Not used: the defect's Dirac embedding is exactly the part not")
print(f"     yet derived, so the scalar-coupling (unsuppressed) case is kept.")
print(f"  2. the full supersolid suppression (NCRI-enhanced) is stronger than")
print(f"     the condensate-free floor; not used, floor kept.")
print(f"  3. higher harmonics of the PN potential carry further form-factor")
print(f"     penalties e^-2(G-G_P)w and are negligible against shell one.")
