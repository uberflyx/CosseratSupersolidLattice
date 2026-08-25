"""The framework's rotating string against the light-meson Regge data, no free parameters.

FINDING.  The monograph previously said the massive-endpoint correction to the Regge
slope "runs in the right direction".  It does not, and this script is the computation
that replaced that sentence: massive ends raise every M^2 step above 2 pi sigma, while
the measured steps sit below it and shrink with excitation, so no endpoint mass can
close the -7% slope gap or the second-step overshoot.  The shrinking steps instead
track the growth of the resonance widths along the trajectory, which sizes a unitarity
self-energy shift the framework does not yet compute.

OBSERVATION, recorded for follow-up and not yet in the monograph.  With the standard
chiral-limit Adler-zero intercept alpha_rho(0) = 1/2 (Lovelace-Shapiro), the string
alone puts the chiral trajectory at M_L^2 = (2L+1) pi sigma: 779.8 MeV at L = 0 against
the rho's 775.3 (+0.6%), and 1350.7 at L = 1 against the chiral doublet centre
1352.6 (-0.14%), with nothing fitted.

The question.  The monograph records a -7% miss on the Regge slope (1/(2 pi sigma) =
0.822 GeV^-2 against the fitted 0.88-0.95) and a +19% overshoot on the second Regge
step, and says the massive-endpoint correction "runs in the right direction".  This
script checks that sentence and then asks what, on the framework's own fixed inputs,
actually reproduces the measured states.

The classical rotating Nambu-Goto string with tension sigma and equal endpoint masses m,
each end moving at speed beta:

    boundary condition   sigma / gamma = gamma m omega beta
    E(beta) = 2 m gamma (1 + gamma beta arcsin beta)
    J(beta) = (m^2 gamma^3 beta^2 / sigma) (beta + gamma arcsin beta)

with gamma = 1/sqrt(1-beta^2).  Massless limit: J = E^2/(2 pi sigma).

Fixed inputs (nothing scanned, nothing fitted except where labelled as a diagnostic):
    sigma = (2 pi m_0)^2,  m_0 = m_e/alpha       [the framework's string tension]
    m_q in {0, 315 MeV = N_c^2 m_0 / 2}          [current vs constituent ends]

Data: the leading light trajectory rho(775.26) J=1, a2(1318.2) J=2, rho3(1688.8) J=3,
a4(1967) J=4 [PDG], and the chiral doublet centres 1352.6 (L=1) and 1687.8 (L=2).
"""
import numpy as np
from scipy import optimize

alpha = 1 / 137.035999177
m_e = 0.51099895069e-3            # GeV
m_0 = m_e / alpha
sigma = (2 * np.pi * m_0)**2      # GeV^2
TPS = 2 * np.pi * sigma           # the Regge step 2 pi sigma
m_con = 4.5 * m_0                 # N_c^2 m_0 / 2 = 315 MeV

print("sigma = %.5f GeV^2   sqrt(sigma) = %.1f MeV   2 pi sigma = %.4f GeV^2"
      % (sigma, 1e3 * np.sqrt(sigma), TPS))
print("1/(2 pi sigma) = %.4f GeV^-2\n" % (1 / TPS))


def EJ(beta, m, sig=sigma):
    """Energy and orbital angular momentum of the classical rotating loaded string."""
    g = 1 / np.sqrt(1 - beta**2)
    E = 2 * m * g * (1 + g * beta * np.arcsin(beta))
    J = (m**2 * g**3 * beta**2 / sig) * (beta + g * np.arcsin(beta))
    return E, J


def E_of_J(J, m, sig=sigma):
    """Invert J(beta) for the loaded string; massless limit if m = 0."""
    if m == 0.0:
        return np.sqrt(2 * np.pi * sig * J)
    f = lambda b: EJ(b, m, sig)[1] - J
    b = optimize.brentq(f, 1e-9, 1 - 1e-12, xtol=1e-15)
    return EJ(b, m, sig)[0]


# ----------------------------------------------------------------------
# 1. The direction of the classical massive-endpoint correction
# ----------------------------------------------------------------------
print("=== 1. direction check: chord slope of the loaded string, m = 315 MeV ===")
print("   L    E [MeV]    M^2 [GeV^2]   step in M^2   local dJ/dM^2")
prev = None
for L in (1, 2, 3, 4):
    E = E_of_J(L, m_con)
    step = E**2 - prev if prev is not None else float('nan')
    dE = (E_of_J(L + 1e-4, m_con) - E_of_J(L - 1e-4, m_con)) / 2e-4
    print("  %2d   %7.1f     %7.4f      %7.4f       %7.4f"
          % (L, 1e3 * E, E**2, step, 1 / (2 * E * dE)))
    prev = E**2
print("  massless slope for comparison: dJ/dM^2 = %.4f everywhere" % (1 / TPS))
print("  -> steps EXCEED 2 pi sigma and the chord slope sits BELOW 1/(2 pi sigma):")
print("     massive ends push the fitted alpha' DOWN, away from the observed 0.88-0.95.\n")

# ----------------------------------------------------------------------
# 2. The observed pattern, laid out
# ----------------------------------------------------------------------
data_J = [("rho(770)", 1, 0.77526), ("a2(1320)", 2, 1.3182),
          ("rho3(1690)", 3, 1.6888), ("a4(1970)", 4, 1.967)]
print("=== 2. the measured leading trajectory ===")
print("   state        J    M [MeV]   M^2      step     needed slope")
prev = None
for nm, J, M in data_J:
    step = M**2 - prev if prev is not None else float('nan')
    print("  %-12s %d   %7.1f  %6.4f   %6.4f      %6.4f"
          % (nm, J, 1e3 * M, M**2, step, step / 1.0))
    prev = M**2
M2 = np.array([M**2 for _, _, M in data_J])
Jv = np.array([J for _, J, _ in data_J])
A = np.vstack([M2, np.ones_like(M2)]).T
slope, icept = np.linalg.lstsq(A, Jv, rcond=None)[0]
print("  linear fit: alpha' = %.4f GeV^-2, intercept = %.4f" % (slope, icept))
print("  steps run %.3f, %.3f, %.3f: all BELOW 2 pi sigma = %.3f and shrinking.\n"
      % (M2[1] - M2[0], M2[2] - M2[1], M2[3] - M2[2], TPS))

# ----------------------------------------------------------------------
# 3. Fixed-input predictions, three hypotheses, no dials
# ----------------------------------------------------------------------
print("=== 3. state-by-state predictions on fixed inputs ===")
print("hypothesis A: massless ends, string alone, M^2 = 2 pi sigma (J - 1/2)")
print("   (the 1/2 is the Adler-zero intercept of the Lovelace-Shapiro amplitude,")
print("    alpha_rho(0) = 1/2 in the chiral limit: standard, parameter-free)")
print("hypothesis B: massless ends, M^2 = m_rho^2 + 2 pi sigma (J - 1)")
print("hypothesis C: constituent ends m = 315 MeV, classical loaded string")
m_rho = 0.775286
print("   state        J    obs      A [dev]        B [dev]        C [dev]")
for nm, J, M in data_J:
    a = np.sqrt(TPS * (J - 0.5))
    b = np.sqrt(m_rho**2 + TPS * (J - 1))
    c = E_of_J(J, m_con)
    print("  %-12s %d  %7.1f  %7.1f %+5.1f%%  %7.1f %+5.1f%%  %7.1f %+5.1f%%"
          % (nm, J, 1e3 * M, 1e3 * a, 100 * (a / M - 1),
             1e3 * b, 100 * (b / M - 1), 1e3 * c, 100 * (c / M - 1)))

# the doublet centres, read as L = 1 and L = 2 on the same string
print("\n   doublet centres against the same three hypotheses (L = J here):")
for nm, L, M in (("centre L=1", 1, 1.35261), ("centre L=2", 2, 1.68781)):
    a = np.sqrt(TPS * (L + 0.5))          # J = L+1 stretched? no: centre carries L; use L + 1/2
    b = np.sqrt(m_rho**2 + TPS * L)
    c = E_of_J(L, m_con)
    print("  %-12s %d  %7.1f  %7.1f %+5.1f%%  %7.1f %+5.1f%%  %7.1f %+5.1f%%"
          % (nm, L, 1e3 * M, 1e3 * a, 100 * (a / M - 1),
             1e3 * b, 100 * (b / M - 1), 1e3 * c, 100 * (c / M - 1)))

# ----------------------------------------------------------------------
# 4. Diagnostic only: what endpoint mass would the data demand, sigma held fixed?
# ----------------------------------------------------------------------
print("\n=== 4. diagnostic inversion (one parameter, labelled as such) ===")


def chi2(m):
    return sum((E_of_J(J, m) - M)**2 for _, J, M in data_J)


res = optimize.minimize_scalar(chi2, bounds=(1e-6, 0.6), method='bounded')
m_fit = res.x
print("  best-fit endpoint mass at fixed sigma: m = %.1f MeV, rms = %.1f MeV"
      % (1e3 * m_fit, 1e3 * np.sqrt(res.fun / 4)))
print("  (against the constituent 315 and the current-quark few MeV)")

# and the reverse: massless ends, what sigma would the data demand?
res2 = optimize.minimize_scalar(lambda s: sum((np.sqrt(2 * np.pi * s * J) - M)**2
                                              for _, J, M in data_J),
                                bounds=(0.1, 0.3), method='bounded')
print("  massless ends, best-fit sigma: %.4f GeV^2, sqrt = %.1f MeV (framework %.1f)"
      % (res2.x, 1e3 * np.sqrt(res2.x), 1e3 * np.sqrt(sigma)))
