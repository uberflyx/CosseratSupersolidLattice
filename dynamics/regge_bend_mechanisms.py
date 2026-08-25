"""What bends the trajectory?  Three candidate mechanisms, tested against the data.

THE OBSERVATION.  The chiral ladder M_N^2 = (2N-1) pi sigma places the first two rungs
to 0.6% and 0.14% and misses the third by 3.3%.  Equivalently the measured rung
spacings run 1.229 then 1.021 GeV^2 against a constant 2 pi sigma = 1.216, so the
trajectory bends after the first step.

Three mechanisms have been proposed for that bend, and this script tests each against
the spectrum rather than arguing about it.

  (1) Massive endpoints.  Excluded already by regge_massive_endpoints.py: loading the
      rotating string raises every step ABOVE 2 pi sigma while the data sits below.

  (2) Unitarity self-energy, the corpus's current claim: wider states are pulled
      further below their bare mass as more channels open, so the deficit should grow
      with the width.  Tested here.

  (3) 1/N_c.  The ladder is a narrow-resonance construction, exact only at large N_c,
      so the departure should scale with Gamma/M.  Tested here.

Data: PDG.  The rung level is read two ways where both are available, the leading
J = N state and the centre of the chirally split J = 1 pair, which the level structure
requires to be degenerate.
"""
import numpy as np

alpha = 1 / 137.035999177
m_0 = 0.51099895069 / alpha
PS = np.pi * (2 * np.pi * m_0)**2 / 1e6          # GeV^2


def centre(a, v):
    return np.sqrt((a * a + v * v) / 2)


# rung, level mass [MeV], mean width of the states on it [MeV], how the level is read
RUNGS = [
    (1, 778.96, (3 * 147.4 + 8.68) / 4, "rho(770) and omega(782), isospin-weighted"),
    (2, centre(1230., 1465.), (400. + 400.) / 2, "a1(1260)/rho(1450) centre"),
    (3, centre(1655., 1720.), (250. + 250.) / 2, "a1(1640)/rho(1700) centre"),
]
LEAD = [(1, 775.26, 147.4, "rho(770)"), (2, 1318.2, 107., "a2(1320)"),
        (3, 1688.8, 161., "rho3(1690)"), (4, 1967., 324., "a4(1970)")]

print("=== the rung levels, and mechanism (2) ===")
print("If wider states are pulled further down, the deficit must grow with the width.\n")
print("  N   level     width    Gamma/M   ladder    dev(M)   deficit(M^2)  deficit/(Gamma M)")
for N, M, G, how in RUNGS:
    lad = 1e3 * np.sqrt((2 * N - 1) * PS)
    d = (lad**2 - M**2) / 1e6
    print("  %d  %7.1f  %7.1f   %6.3f  %7.1f  %+6.2f%%   %+8.4f    %+8.3f"
          % (N, M, G, G / M, lad, 100 * (lad / M - 1), d, d / (G * M / 1e6)))
print("\n  The widths run 148, 400, 250: rung 2 carries the BROADEST states on the")
print("  trajectory and is the rung the ladder places best, at -0.14%.  Rung 3 is")
print("  narrower and is missed by 3.3%.  The deficit therefore ANTI-correlates with")
print("  the width, and mechanism (2) predicts the opposite.  Falsified.")

print("\n=== the leading trajectory, same test ===")
print("  N   state        width   ladder     dev      deficit   deficit/(Gamma M)")
for N, M, G, nm in LEAD:
    lad = 1e3 * np.sqrt((2 * N - 1) * PS)
    d = (lad**2 - M**2) / 1e6
    print("  %d  %-11s %6.1f  %7.1f  %+6.2f%%  %+8.4f    %+8.3f"
          % (N, nm, G, lad, 100 * (lad / M - 1), d, d / (G * M / 1e6)))
print("\n  Here the widths do grow, and for N = 2, 3, 4 the ratio is nearly constant at")
print("  0.61 to 0.69.  The rho breaks it by a factor of ten, and the a2 carries the")
print("  SMALLEST width on the trajectory while showing twelve times the rho's deficit.")
print("  So the regularity holds over three points and fails on the fourth, which is a")
print("  regularity and not a mechanism.")

print("\n=== mechanism (3): does the departure scale with Gamma/M? ===")
print("  N   Gamma/M   |dev| in M^2   ratio")
for N, M, G, how in RUNGS:
    lad2 = (2 * N - 1) * PS * 1e6
    rel = abs(lad2 - M**2) / M**2
    print("  %d   %6.3f      %6.3f      %6.3f" % (N, G / M, rel, rel / (G / M)))
print("\n  No scaling: rung 2 has the largest Gamma/M and the smallest departure.")
print("  What survives is a BOUND rather than an explanation.  A narrow-resonance")
print("  construction carries an accuracy of order (Gamma/M)^2 in M^2 when no threshold")
print("  sits nearby, which is 2% to 9% across these rungs; the observed departures are")
print("  0.3%, 0.6% and 6.7%, so the ladder performs at or inside its own accuracy on")
print("  the first two rungs and at the edge of it on the third.")

print("\n=== what is actually established ===")
print("  Endpoint masses are excluded, width-driven self-energy is falsified on the")
print("  rungs, and 1/N_c does not scale.  The bend is not explained.  It is also not")
print("  a defect peculiar to this framework: imposing the Adler zero on the")
print("  Lovelace-Shapiro amplitude gives an intercept of 1/2 where fits to the")
print("  physical spectrum return about 0.42, a tension the dual-resonance literature")
print("  has carried since 1968 and attributes to the model being tree-level.")


# ----------------------------------------------------------------------
# 4. Is the "-7% Regge slope" a defect of the tension at all?
# ----------------------------------------------------------------------
print("\n=== the slope, with the intercept held at the Adler value ===")
print("A linear fit has two parameters.  Fitting both to a bending trajectory spreads")
print("the bend across them, so part of the apparent slope error is the intercept")
print("absorbing curvature.  The chiral ladder fixes the intercept independently at")
print("exactly 1/2, which lets each rung state a tension of its own.\n")
TPS_GEV = 2 * np.pi * (2 * np.pi * 0.51099895069e-3 * 137.035999177)**2
PS = TPS_GEV / 2.0
ISO = {1: [775.26, 782.66], 2: [1230., 1229.5, 1318.2, 1439., 1465.],
       3: [1655., 1688.8, 1706., 1720.]}
M2 = {N: np.mean([(v * 1e-3)**2 for v in st]) for N, st in ISO.items()}
b_fit, c_fit = np.polyfit(np.array([1, 2, 3.]),
                          np.array([M2[1], M2[2], M2[3]]), 1)
print("  free two-parameter fit: slope %.4f GeV^2 (%+.1f%% from the derived %.4f),"
      % (b_fit, 100 * (TPS_GEV / b_fit - 1), TPS_GEV))
print("     intercept a = %.4f against the Adler value of 0.5." % (-c_fit / b_fit))
print("  intercept held at 1/2, each rung on its own:")
for n in (1, 2, 3):
    print("     rung %d  ->  2 pi sigma = %.6f GeV^2  (derived %+.2f%%)"
          % (n, M2[n] / (n - 0.5), 100 * (TPS_GEV / (M2[n] / (n - 0.5)) - 1)))
print("\n  The first rung confirms the derived tension to 0.23%.  The seven per cent is")
print("  the bend counted as a slope error, not a separate defect of sigma.")
d = [TPS_GEV * (n - 0.5) - M2[n] for n in (1, 2, 3)]
print("  And no single tension can absorb it: the M^2 deficits divided by (N-1/2) run")
print("  %s, which is not constant."
      % ", ".join("%.4f" % (x / (n - 0.5)) for x, n in zip(d, (1, 2, 3))))


# ----------------------------------------------------------------------
# 5. The standard analyticity form, which is the fourth candidate
# ----------------------------------------------------------------------
print("\n=== square-root trajectory: alpha(s) = 1/2 + a1(sqrt(s0) - sqrt(s0-s)) ===")
print("This is the form analyticity and unitarity motivate: linear at low s, bending to")
print("sqrt(s) asymptotically, with one new parameter s0 setting where the bend begins.")
print("Hold the low-s slope at the framework's own value and each rung then solves for")
print("s0 on its own, via s0 = u^2 / (2u - M^2) with u = (N - 1/2) pi sigma.\n")
for n in (1, 2, 3):
    u = (n - 0.5) * PS
    s0 = u * u / (2 * u - M2[n])
    print("  rung %d  ->  s0 = %7.2f GeV^2   (sqrt(s0) = %.1f GeV)" % (n, s0, np.sqrt(s0)))
print("\n  Not one scale but three, falling by a factor of five, so the form is excluded")
print("  as soon as the framework's slope is imposed.  Letting the slope float fits, as")
print("  two parameters on three points must, and buys nothing: the slope it then wants")
print("  is low by the same margin a straight line wanted.  Fourth candidate, excluded.")
