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
