"""The chiral vector ladder M_L^2 = (2L+1) pi sigma, and what anchoring it costs.

BACKGROUND.  regge_massive_endpoints.py excluded endpoint masses as the cause of the
Regge bend and, along the way, recorded an observation: with the Adler-zero intercept
alpha_rho(0) = 1/2, the framework's string alone fixes the bottom of the vector
trajectory.  This script chases that observation on fixed inputs.

THE ANCHOR.  The Lovelace-Shapiro amplitude for pi pi scattering,
A = Gamma(1-a(s)) Gamma(1-a(t)) / Gamma(1-a(s)-a(t)) with a(s) = a_0 + a' s, is the
minimal dual amplitude with an Adler zero: it vanishes at zero pion momentum, which is
what a Goldstone boson requires.  Imposing the zero forces a_0 = 1/2, and the rho at
the first pole gives m_rho^2 = 1/(2 a').  The framework supplies a' = 1/(2 pi sigma)
from its own string tension and supplies the Adler zero from its own pion, the
Goldstone mode of stacking slip, whose uniform displacement costs nothing.  Combining
the two: m_rho^2 = pi sigma, and the chiral trajectory is the odd ladder

    M_L^2 = (2L + 1) pi sigma,   pi sigma = 4 pi^3 m_0^2,   m_rho = 2 pi^{3/2} m_0,

with m_0 = m_e / alpha and nothing else.

Everything below is computed on that one line.  Nothing is fitted.
"""
import numpy as np

alpha = 1 / 137.035999177
m_e = 0.51099895069e-3                       # GeV
m_0 = m_e / alpha
sigma = (2 * np.pi * m_0)**2
PS = np.pi * sigma                           # the ladder unit, pi sigma
m_pi = 0.138007                              # isospin average, 2 m_0 - 4 m_e
m_rho_meas = 0.77526                         # PDG, e+e-
m_rho_cluster = 11 * m_0 * 1e3 - 11 * (4 - 4.891) * m_e * 1e3   # MeV, spectral route

print("pi sigma = %.5f GeV^2      sqrt = %.2f MeV" % (PS, 1e3 * np.sqrt(PS)))
print("m_rho(ladder) = 2 pi^(3/2) m_0 = %.2f MeV = %.4f m_0"
      % (1e3 * 2 * np.pi**1.5 * m_0, 2 * np.pi**1.5))
print("m_rho(cluster) = %.2f MeV = %.4f m_0   [the spectral route]"
      % (m_rho_cluster, m_rho_cluster / (1e3 * m_0)))
print("m_rho(PDG)     = %.2f MeV" % (1e3 * m_rho_meas))

# ----------------------------------------------------------------------
# 1. The ladder against the chiral-limit objects it actually predicts
# ----------------------------------------------------------------------
print("\n=== 1. the odd ladder against the chiral centres ===")
print("The doublet centre sqrt((m_V'^2 + m_A^2)/2) is the corpus's own chiral-limit")
print("construction, so it is the like-for-like comparison; the physical rho is not.")
rows = [("L=0  (physical rho, category-mismatched)", 0, m_rho_meas, None),
        ("L=0  (physical omega, same caveat)", 0, 0.78266, None),
        ("L=1  centre of a1(1260)/rho(1450)", 1, 1.35261, 0.023),
        ("L=2  centre of a1(1640)/rho(1700)", 2, 1.68781, 0.013)]
for nm, L, M, err in rows:
    pred = np.sqrt((2 * L + 1) * PS)
    dev = 100 * (pred / M - 1)
    sig_txt = ""
    if err:
        sig_txt = "  (%.1f sigma)" % abs((pred - M) / err)
    print("  %-42s pred %7.1f  obs %7.1f  %+5.2f%%%s"
          % (nm, 1e3 * pred, 1e3 * M, dev, sig_txt))

# ----------------------------------------------------------------------
# 2. What the physical pion does to the anchor, stated honestly
# ----------------------------------------------------------------------
print("\n=== 2. the anchor with the physical pion ===")
print("The Adler zero sits at the pion mass shell, not at zero, once m_pi > 0:")
print("alpha(m_pi^2) = 1/2 gives m_rho^2 = pi sigma + m_pi^2.")
mr2 = PS + m_pi**2
print("  m_rho(physical-pion anchor) = %.1f MeV  (+%.1f%% on the PDG rho)"
      % (1e3 * np.sqrt(mr2), 100 * (np.sqrt(mr2) / m_rho_meas - 1)))
print("  Lattice-QCD chiral extrapolations put the chiral-limit rho near 760 MeV,")
print("  below the physical 775, so a chiral-limit formula SHOULD overshoot the")
print("  chiral rho and undershoot nothing; the ladder's %.1f against ~760 is +2.6%%."
      % (1e3 * np.sqrt(PS)))
print("  The clean statements are therefore the centres of section 1, not the rho.")

# ----------------------------------------------------------------------
# 3. The near-identity the corpus carries becomes an identity
# ----------------------------------------------------------------------
print("\n=== 3. eq-wsr-coincidence closed ===")
print("The corpus records 8 pi^3 =? 2 (m_rho/m_0)^2, 248.05 against 245.16, 1.2%,")
print("as an unexplained near-identity.  On the ladder it is exact:")
print("  2 (m_rho^2/m_0^2) = 2 pi sigma / m_0^2 = 2 pi (2 pi)^2 = 8 pi^3 identically.")
gap = 8 * np.pi**3 / (2 * (m_rho_cluster / (1e3 * m_0))**2) - 1
print("  The measured 1.2%% is then the cluster rho's own offset from sqrt(pi sigma):")
print("  (m_rho_ladder/m_rho_cluster)^2 - 1 = %+.2f%%, i.e. the same number." % (100 * gap))

# ----------------------------------------------------------------------
# 4. m_a1^2 = 3 m_rho^2 exact: the two Gamma_ee routes merge
# ----------------------------------------------------------------------
print("\n=== 4. the Weinberg factor becomes exactly 3/2 ===")
f_pi = 3**0.25 * m_0 * (1 + alpha / np.pi)
G_tree = 1e3 * 8 * np.pi * alpha**2 * f_pi**2 / (3 * m_rho_meas) * 1e3   # keV, on GeV inputs
G_closed = 1e3 * 4 * np.pi * alpha**2 * f_pi**2 / m_rho_meas * 1e3
m_a1_regge = np.sqrt(m_rho_cluster**2 * 1e-6 + 2 * PS)
wein_now = 1 / (1 - (m_rho_cluster * 1e-3)**2 / m_a1_regge**2)
fin = 1 - m_pi**2 / (3 * PS)
print("  On the ladder m_a1^2/m_rho^2 = 3 exactly, so the sum-rule route is")
print("  KSRF x 3/2 = closed form identically; the corpus's 0.4%% inter-route gap")
print("  (7.34 against 7.36 keV) is the same cluster offset as section 3.")
print("  current Weinberg factor on corpus masses: %.4f   ladder: 1.5 exact" % wein_now)
print("  Gamma_ee on the ladder: %.3f x %.5f = %.3f keV  (corpus closed-form 7.288)"
      % (G_closed, fin, G_closed * fin))

# ----------------------------------------------------------------------
# 5. Side predictions and their honest status
# ----------------------------------------------------------------------
print("\n=== 5. side checks ===")
print("  pion trajectory (LS daughter, intercept 0): M^2 = 2 pi sigma L")
print("    L=1: %.0f MeV against the b1(1235): %+.1f%%   [weak; recorded, not claimed]"
      % (1e3 * np.sqrt(2 * PS * 1), 100 * (np.sqrt(2 * PS) / 1.2295 - 1)))
print("  bottom of the vector ladder has no axial partner: the a1 first appears at")
print("    L = 1, which is the parity-doubling-up-the-spectrum pattern the chiral")
print("    literature reports; the ladder builds it in structurally.")
print("  the cluster integer: the spectral route counts N = 11 nodes; the ladder gives")
print("    the pure number 2 pi^(3/2) = %.4f, within 1.2%% of 11." % (2 * np.pi**1.5))
print("    The two routes agree on the rho to 0.6%% but their fine structure differs:")
print("    the ladder sits 18.7 m_e above 11 m_0, the cluster 9.8 m_e above it, so")
print("    the alpha-scale corrections are where the routes genuinely part company.")
print("    Recorded as a curiosity, not a closure.")
print("  Lovelace-Shapiro is reported ghost-free with critical dimension FOUR")
print("    (Bianchi et al.), the one hadronic dual amplitude distinguished in d = 4;")
print("    the framework's lattice is four-dimensional.  Recorded, not built on.")


# ----------------------------------------------------------------------
# 6. The level structure: leading state and J=1 daughter must be degenerate
# ----------------------------------------------------------------------
print("\n=== 6. the LS level structure, and the degeneracy it predicts ===")
print("At level N the amplitude's residue carries spins 0..N, all at (2N-1) pi sigma.")
print("So the leading J = N state and the J = 1 daughter sit together.  The J = 1")
print("daughter is an axial/vector pair split by chiral symmetry breaking, and its")
print("centre is exactly the doublet construction the corpus already uses.")
print("Nothing below imposes the degeneracy; it is a prediction of the level structure.\n")


def centre(a, ea, v, ev):
    c = np.sqrt((a * a + v * v) / 2)
    dc = np.hypot((a / (2 * c)) * ea, (v / (2 * c)) * ev)
    return c, dc


print("  N  leading J=N          J=1 daughter pair        centre       degeneracy")
for N, ln, lm, le, dn, a, ea, v, ev in (
        (2, "a2(1320)", 1318.2, 0.6, "a1(1260)/rho(1450)", 1230., 40., 1465., 25.),
        (3, "rho3(1690)", 1688.8, 2.1, "a1(1640)/rho(1700)", 1655., 16., 1720., 20.)):
    ct, dct = centre(a, ea, v, ev)
    d = ct - lm
    print("  %d  %-11s %7.1f  %-19s %7.1f+-%4.1f   %+6.1f MeV (%.2f sigma)"
          % (N, ln, lm, dn, ct, dct, d, abs(d) / np.hypot(dct, le)))
print("\n  N = 3 is the sharp one: two multiplets sharing no state land on top of")
print("  each other, 1688.8 against 1687.8, a tenth of a standard deviation.")
print("  N = 2 is broken at 1.5 sigma, the a2 sitting low.")

print("\n=== 7. what the ladder's centre forces on the doublet ===")
M_A1 = 1227.37     # Born-cluster eigenvalue, an independent derivation
for nm, ct in (("corpus centre, sqrt(m_rho,cluster^2 + 2 pi sigma)", 1348.11),
               ("ladder centre, sqrt(3 pi sigma)", 1e3 * np.sqrt(3 * PS))):
    rp = np.sqrt(2 * ct**2 - M_A1**2)
    half = (ct**2 - M_A1**2) / 1e6
    print("  %-46s centre %7.2f MeV" % (nm, ct))
    print("     rho' = %7.1f against 1465 +- 25  (%+.2f%%, %.2f sigma)"
          % (rp, 100 * (rp / 1465 - 1), abs(rp - 1465) / 25))
    print("     half-splitting %.4f against 0.3167 GeV^2  (%+.2f%%)"
          % (half, 100 * (half / 0.31666 - 1)))

# ----------------------------------------------------------------------
# 8. The degeneracy tested across the WHOLE rung, not a chosen pair
# ----------------------------------------------------------------------
print("\n=== 8. the full-rung test ===")
print("Section 6 compared the leading J = N state with the J = 1 daughter centre and")
print("found 1 MeV at rung 3 against 34 MeV at rung 2.  Neither figure should be read")
print("alone: the amplitude puts EVERY spin from 0 to N on rung N, so the test is")
print("whether all of them sit together.  Isovector states, rung by rung:\n")
RUNG_STATES = {
    2: [("a1(1260)", 1, 1230.0), ("b1(1235)", 1, 1229.5), ("a2(1320)", 2, 1318.2),
        ("a0(1450)", 0, 1439.0), ("rho(1450)", 1, 1465.0)],
    3: [("a1(1640)", 1, 1655.0), ("rho3(1690)", 3, 1688.8), ("a2(1700)", 2, 1706.0),
        ("rho(1700)", 1, 1720.0)],
}
spreads = {}
for N, st in RUNG_STATES.items():
    m = np.array([x[2] for x in st])
    spreads[N] = (m.max() - m.min()) / m.mean()
    print("  rung %d: %s" % (N, ", ".join("%s J=%d %.0f" % (n, j, v) for n, j, v in st)))
    print("     spread %.0f MeV = %.1f%% of the mean %.0f, ladder %.0f MeV (%+.1f%%)\n"
          % (m.max() - m.min(), 100 * spreads[N], m.mean(),
             1e3 * np.sqrt((2 * N - 1) * PS),
             100 * (1e3 * np.sqrt((2 * N - 1) * PS) / m.mean() - 1)))
print("  The degeneracy is broken by %.1f%% at rung 2 and %.1f%% at rung 3, improving by"
      % (100 * spreads[2], 100 * spreads[3]))
print("  a factor of %.1f in one step.  Two consequences, pulling opposite ways:"
      % (spreads[2] / spreads[3]))
print("    - the a2's 34 MeV is not a puzzle about the a2; it is one slice of a rung")
print("      split by 236 MeV, and any pair from that rung would disagree similarly.")
print("    - rung 3's 1 MeV is better than its own rung warrants, the four states")
print("      there scattering over 65 MeV, so that pair agrees more closely than the")
print("      prediction can support.")
print("  What the ladder gets right is the trend.  It is spin-blind, and the departures")
print("  are spin-dependent forces falling away with excitation; the framework does not")
print("  compute them, so the rate of the fall is observed rather than predicted.")
