"""Independent recomputation of the pion-form-factor truncation systematic on a_mu.

Companion check on formfactor_truncation.py, written so the two share no code.

Written from the closed forms in the monograph rather than by reusing the companion
module, so that agreement between the two is evidence rather than a copy.  Nothing is
imported from the public repository.

The question.  The pion form factor is a coupling-weighted sum over the vector tower,
F_pi(s) = sum_n c_n A_n(s) with sum_n c_n = 1 (charge conservation, F_pi(0) = 1).  The
framework truncates at the first recurrence, so c_2 is an unknown bounded but not
derived.  How far does a_mu move as c_2 sweeps its bound?

The answer depends entirely on WHICH normalisation condition fixes c_0:

  AREA rule  - c_0 is solved so the channel's total spectral area returns Gamma_ee.
               Adding c_2 then rescales the whole channel, so a_mu moves a lot.
  PEAK rule  - c_0 is solved so |F_pi(m_rho^2)| returns the value the leptonic width
               demands at the pole.  A third resonance then displaces the second rather
               than rescaling the channel, and a_mu barely moves.

Both are computed below.  The amplitudes are unit-normalised relativistic Breit-Wigners
with P-wave energy-dependent widths, deliberately NOT the Gounaris-Sakurai form the
production calculation uses: the sensitivity being tested is a property of how the
constraint propagates, so it must survive a change of line shape to mean anything.

Checks performed before any sensitivity is quoted:
  1. F_pi(0) = 1 identically, to machine precision, for every (c_0, c_1, c_2).
  2. The peak condition reproduces the leptonic width through the narrow-resonance
     relation R(m_rho^2) = 9 Gamma_ee / (alpha^2 Gamma_rho).
  3. The muon kernel K(r) reproduces its analytic limits: K -> 1/(3r) for r >> 1 and
     K(0) = 1/2.
"""

import numpy as np
from scipy import integrate, optimize

# ----------------------------------------------------------------------
# Inputs.  alpha is the Peierls-Nabarro tunnelling amplitude; m_e is the single
# scale-setting input.  Everything hadronic below is built from the two.
# ----------------------------------------------------------------------
alpha = 1.0 / 137.035999177
m_e = 0.51099895069            # MeV, CODATA 2022
m_mu = 105.6583755             # MeV, PDG 2026

m_0 = m_e / alpha              # node mass, 70.02 MeV
N_c = 3.0

# Pion: spectral mass formula gives the isospin average; the charged mass carries
# an imported splitting.  The two-pion threshold uses the charged mass, the
# Weinberg sum rules use the isospin average.
m_pi_iso = 2.0 * m_0 - 4.0 * m_e
m_pi = m_pi_iso + (139.57039 - 134.9768) / 3.0

f_pi = N_c**0.25 * m_0 * (1.0 + alpha / np.pi)

# Rho: spectral crossed-fault closure, eleven-node cluster, graph eigenvalue 4.891.
RHO_EIG = 4.891
m_rho = 11.0 * m_0 - 11.0 * (4.0 - RHO_EIG) * m_e


def beta_pi(s, mpi=m_pi):
    """Pion velocity in the pair centre-of-mass frame."""
    return np.sqrt(np.maximum(0.0, 1.0 - 4.0 * mpi**2 / np.asarray(s, float)))


# KSRF coupling and the width it implies.
g2 = m_rho**2 / (2.0 * f_pi**2)
Gamma_rho = g2 * m_rho * beta_pi(m_rho**2)**3 / (48.0 * np.pi)

# Leptonic width: chiral-limit Weinberg value times the derived finite-pion factor.
sigma_string = (2.0 * np.pi * m_0)**2
m_a1_sq = m_rho**2 + 2.0 * np.pi * sigma_string
Gee_rho = (4.0 * np.pi * alpha**2 * f_pi**2 / m_rho) * (1.0 - m_pi_iso**2 / m_a1_sq)

# Tower: one state per Regge interval of width 2 pi sigma.
TPS = 2.0 * np.pi * sigma_string
m1 = np.sqrt(m_rho**2 + 1.0 * TPS)
m2 = np.sqrt(m_rho**2 + 2.0 * TPS)
G1, G2 = 400.0, 300.0          # MeV, stated widths for the two recurrences

print("inputs")
print("  m_0        = %8.4f MeV" % m_0)
print("  m_pi       = %8.4f MeV   (charged)" % m_pi)
print("  f_pi       = %8.4f MeV" % f_pi)
print("  m_rho      = %8.4f MeV" % m_rho)
print("  Gamma_rho  = %8.4f MeV" % Gamma_rho)
print("  Gee_rho    = %8.5f keV" % (Gee_rho * 1e3))
print("  m_rho(1),(2) = %.1f, %.1f MeV" % (m1, m2))


# ----------------------------------------------------------------------
# Amplitudes.  Relativistic Breit-Wigner with a P-wave energy-dependent width,
# divided through by its own value at s = 0 so that A_n(0) = 1 exactly.  The
# constraint sum_n c_n = 1 is then F_pi(0) = 1 identically.
# ----------------------------------------------------------------------
def bw_unit(m, G):
    def raw(s):
        s = np.asarray(s, dtype=complex)
        k = np.sqrt(np.maximum(0.0, s.real - 4.0 * m_pi**2)) / 2.0
        k_m = np.sqrt(m**2 - 4.0 * m_pi**2) / 2.0
        Gs = G * (m / np.sqrt(np.maximum(s.real, 1e-9))) * (k / k_m)**3
        Gs = np.where(s.real > 4.0 * m_pi**2, Gs, 0.0)
        return m**2 / (m**2 - s - 1j * np.sqrt(np.maximum(s.real, 0.0)) * Gs)
    n0 = raw(0.0)
    return lambda s: raw(s) / n0


A0, A1, A2 = bw_unit(m_rho, Gamma_rho), bw_unit(m1, G1), bw_unit(m2, G2)


def F_pi(s, c0, c1, c2):
    return c0 * A0(s) + c1 * A1(s) + c2 * A2(s)


def R_pipi(s, c0, c1, c2):
    """Two-pion spectral function.  The beta^3 is the P-wave threshold."""
    return 0.25 * beta_pi(s)**3 * np.abs(F_pi(s, c0, c1, c2))**2


# ----------------------------------------------------------------------
# The muon kernel.  K(r) = int_0^1 dx x^2 (1-x) / (x^2 + r(1-x)) with r = s/m_mu^2,
# and a_mu = (1/3)(alpha/pi)^2 int ds/s K(s/m_mu^2) R(s).
# ----------------------------------------------------------------------
_XG, _WG = np.polynomial.legendre.leggauss(400)
_X = 0.5 * (_XG + 1.0)
_W = 0.5 * _WG


def K_kernel(r):
    r = np.atleast_1d(np.asarray(r, float))[:, None]
    return np.sum(_W * _X**2 * (1.0 - _X) / (_X**2 + r * (1.0 - _X)), axis=1)


def amu_channel(c0, c1, c2, s_lo, s_hi):
    """Integrate in t = ln s so the resonance is resolved without a fine grid."""
    def integrand(t):
        s = np.exp(t)
        return K_kernel(s / m_mu**2)[0] * R_pipi(s, c0, c1, c2)
    val, _ = integrate.quad(integrand, np.log(s_lo), np.log(s_hi), limit=400)
    return (1.0 / 3.0) * (alpha / np.pi)**2 * val


# ----------------------------------------------------------------------
# Check 3: kernel limits.
# ----------------------------------------------------------------------
assert abs(K_kernel(1e-12)[0] - 0.5) < 1e-3, "K(0) != 1/2"
big = 1e6
assert abs(K_kernel(big)[0] * 3.0 * big - 1.0) < 1e-3, "K -> 1/(3r) fails"
print("\nkernel checks passed: K(0) = %.6f (1/2), 3 r K(r) -> %.5f (1)"
      % (K_kernel(1e-12)[0], K_kernel(big)[0] * 3.0 * big))

# ----------------------------------------------------------------------
# The two normalisation rules.
# ----------------------------------------------------------------------
S_R = m_rho**2
PEAK_TARGET = np.sqrt(36.0 * Gee_rho / (alpha**2 * beta_pi(S_R)**3 * Gamma_rho))
S_LO, S_HI = 4.0 * m_pi**2, 1909.0**2   # duality matching point, N = 2


def solve_peak(c2):
    """c_0 such that |F_pi(m_rho^2)| meets the value the leptonic width demands."""
    def resid(c0):
        return abs(F_pi(S_R, c0, 1.0 - c0 - c2, c2)) - PEAK_TARGET
    c0 = optimize.brentq(resid, 0.3, 3.0, xtol=1e-14)
    return c0, 1.0 - c0 - c2


def area_target():
    """Spectral area a single narrow rho of this leptonic width would carry."""
    return 9.0 * np.pi * Gee_rho * m_rho / alpha**2


def solve_area(c2):
    """c_0 such that the channel's total spectral area returns Gamma_ee."""
    tgt = area_target()

    def area(c0):
        def integrand(t):
            s = np.exp(t)
            return s * R_pipi(s, c0, 1.0 - c0 - c2, c2)
        val, _ = integrate.quad(integrand, np.log(S_LO), np.log(S_HI), limit=400)
        return val
    c0 = optimize.brentq(lambda c: area(c) - tgt, 0.3, 3.0, xtol=1e-10)
    return c0, 1.0 - c0 - c2


# ----------------------------------------------------------------------
# Check 1 and 2.
# ----------------------------------------------------------------------
c0p, c1p = solve_peak(0.0)
assert abs(F_pi(0.0, c0p, c1p, 0.0).real - 1.0) < 1e-12, "F_pi(0) != 1"
R_at_pole = R_pipi(S_R, c0p, c1p, 0.0)
R_expected = 9.0 * Gee_rho / (alpha**2 * Gamma_rho)
print("F_pi(0) = %.15f   (must be 1)" % F_pi(0.0, c0p, c1p, 0.0).real)
print("R(m_rho^2) = %.3f against 9 Gee/(alpha^2 Gamma) = %.3f"
      % (R_at_pole, R_expected))
assert abs(R_at_pole / R_expected - 1.0) < 1e-8

# ----------------------------------------------------------------------
# The bound on c_2.  The photon couplings come from the Cornell tower ratios,
# f_n proportional to sqrt(Gamma_ee,n m_n); the ratio f_2/f_1 = 0.8515 is the
# framework's and is carried here as a stated input, since re-deriving the
# Cornell solve is a separate calculation.
# ----------------------------------------------------------------------
F2_OVER_F1 = 0.8515
c2max = abs(c1p) * F2_OVER_F1

print("\nc_1 = %+.4f at c_2 = 0, so |c_2| <= %.4f" % (c1p, c2max))

for label, solver in (("peak", solve_peak), ("area", solve_area)):
    print("\n%s normalisation" % label)
    print("   c_2        c_0       c_1        a_mu(pipi) x 1e10")
    vals = []
    for c2 in (-c2max, -0.5 * c2max, 0.0, 0.5 * c2max, c2max):
        c0, c1 = solver(c2)
        a = amu_channel(c0, c1, c2, S_LO, S_HI) * 1e10
        vals.append(a)
        print("  %+.4f   %+.4f   %+.4f      %8.2f" % (c2, c0, c1, a))
    print("   spread across the bound: %.1f units" % (max(vals) - min(vals)))
