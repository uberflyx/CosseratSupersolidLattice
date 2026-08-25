"""Where the chiral ladder stops working, and what the ribbon is doing there.

THE OBSERVATION.  The ladder M_N^2 = (2N-1) pi sigma places the first rung to 0.11% and
the third to 3.0%, and four candidate mechanisms for that drift have been excluded by
calculation: massive endpoints, width-driven self-energy, 1/N_c scaling, and the
square-root trajectory analyticity motivates (regge_bend_mechanisms.py,
regge_massive_endpoints.py).  All four treat the meson as a string.

THE OBJECT.  It is not a string.  A meson here is a patch of stacking fault with a
partial dislocation at each end, so its energy is the fault energy times an area, and
sigma = gamma w is that area energy folded with the ribbon's own equilibrium width.  A
ribbon has something a string does not: a length at which it stops being one object.
The framework already states the condition, that the ribbon cannot break by itself and
can only break by nucleating a fresh partial-antipartial pair from the lattice.  That
pair costs 2 m_q, so the ribbon reaches breaking at

    L* = 2 m_q / sigma,

with m_q = N_c^2 m_0 / 2 the constituent partial mass and sigma = (2 pi m_0)^2, both of
which the framework derives.  Nothing here is fitted.

THE TEST.  A classical rotating ribbon with light ends carries M = pi sigma l over a
half-length l, so its full length is L = 2M/(pi sigma).  Comparing that with L* asks
whether a given rung is a single ribbon at all.

THE REDUCTION.  The tension cancels between the two lengths:

    L / L* = M / (pi m_Q),

so whether a rotating ribbon has passed breaking depends only on its mass and the
constituent mass at its ends.  That turns the light-sector observation into a statement
about every sector, and the heavy ones test it.
"""
import numpy as np

alpha = 1 / 137.035999177
m_0 = 0.51099895069e-3 / alpha          # GeV
sigma = (2 * np.pi * m_0)**2            # GeV^2
m_q = 3.0**2 * m_0 / 2                  # constituent partial, N_c^2 m_0 / 2
HBARC = 0.1973269804                    # GeV fm

print("sigma = %.5f GeV^2 (sqrt = %.1f MeV)" % (sigma, 1e3 * np.sqrt(sigma)))
print("m_q = %.1f MeV, so a fresh pair costs 2 m_q = %.0f MeV" % (1e3 * m_q, 2e3 * m_q))
print("L* = 2 m_q / sigma = %.3f GeV^-1 = %.3f fm\n" % (2 * m_q / sigma,
                                                        2 * m_q / sigma * HBARC))

# isovector rung levels, quadratic mean over the states on each rung
ISO = {1: [775.26, 782.66], 2: [1230., 1229.5, 1318.2, 1439., 1465.],
       3: [1655., 1688.8, 1706., 1720.]}
print("=== the light rungs ===")
print("  N   level     ladder    drift     ribbon L    L/L*")
for N, st in ISO.items():
    M = np.sqrt(np.mean([(v * 1e-3)**2 for v in st]))
    lad = np.sqrt((2 * N - 1) * np.pi * sigma)
    L = 2 * M / (np.pi * sigma)
    print("  %d  %7.1f  %7.1f   %+5.2f%%    %.3f fm    %.2f"
          % (N, 1e3 * M, 1e3 * lad, 100 * (lad / M - 1), L * HBARC, L / (2 * m_q / sigma)))
print("\n  The ladder holds while the ribbon is shorter than breaking and degrades once")
print("  it is longer.  It also says why the width tests failed: what matters is the")
print("  ribbon's length, not how fast the state decays.\n")

print("=== the reduction, and what it predicts elsewhere ===")
print("  L/L* = [2M/(pi sigma)] / [2 m_Q/sigma] = M / (pi m_Q).  The tension cancels.\n")
HEAVY = [("light q qbar", m_q, [("rung 1", 0.7790), ("rung 2", 1.3401), ("rung 3", 1.6926)]),
         ("charmonium", 1.5, [("J/psi", 3.0969), ("chi_c1(3510)", 3.5107),
                              ("chi_c2(3556)", 3.5562), ("psi(3770)", 3.7737)]),
         ("bottomonium", 4.7, [("Upsilon(1S)", 9.4604), ("chi_b2(1P)", 9.9122),
                               ("Upsilon(2S)", 10.0234), ("Upsilon(3S)", 10.3551)])]
for nm, mQ, states in HEAVY:
    print("  %s, m_Q = %.3f GeV:" % (nm, mQ))
    for s, M in states:
        flag = "past breaking" if M / (np.pi * mQ) > 1 else ""
        print("     %-13s M = %7.4f GeV   L/L* = %.2f   %s" % (s, M, M / (np.pi * mQ), flag))
print("\n  Every heavy state sits below breaking where the light rungs 2 and 3 sit above.")
print("  So the picture predicts heavy-quark trajectories are straighter than light ones,")
print("  which is what the spectroscopy has always found and what any string picture")
print("  alone leaves unexplained.")
print("\nSTATUS.  This identifies the scale and gets the ordering right; it does not")
print("compute the size of the drift.  The relative M^2 deficits run 0.23%, 1.59% and")
print("6.14% against L/L* of 0.79, 1.35 and 1.71, and no simple power of the excess")
print("length reproduces them, so the mixing with the broken-ribbon configurations is")
print("named and not calculated.")
