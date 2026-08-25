"""Numeric audit of the hadronic-vacuum-polarisation and running-coupling chapter.

Every framework-derived decimal the chapter prints is recomputed here from the two
inputs, alpha and m_e, through the closed forms the chapter itself states, and compared
with the printed value at the printed precision.  Nothing is imported from the
dispersion machinery, so a PASS is independent evidence rather than a copy.

Comparison values (PDG, CODATA, data-driven determinations) are listed as inputs and
are not recomputed.  They are checked against their sources separately.

Two conventions matter and are easy to get wrong; both are set deliberately here.

  SPACELIKE.  The running is a function of momentum transfer, so every evaluation of
  the one-loop shift is at spacelike q^2 = -Q^2.  Evaluated timelike instead, the real
  part near a pair threshold goes negative and the charge-census bound at 45 GeV comes
  out meaningless, which is an artefact of the wrong analytic continuation and not a
  defect in the chapter.

  THE PION MASS.  The Weinberg sum rules are isospin statements about the pion pole, so
  they take the isospin average 2 m_0 - 4 m_e.  The two-pion threshold takes the charged
  mass.  The bare 2 m_0 is the pair before the -4 m_e correction and is not a pion mass.

Run:  python hvp_chapter_audit.py
Exits non-zero if any literal fails to reproduce.
"""

import sys
import numpy as np
from scipy import integrate

FAILS = []


def chk(label, printed, computed, tol=None, unit=""):
    """Compare a printed literal with its recomputation at the printed precision."""
    if tol is None:
        s = repr(float(printed))
        nd = len(s.split(".")[1].rstrip("0")) if "." in s else 0
        tol = 0.5 * 10**(-nd) * 1.001
    ok = abs(printed - computed) <= tol
    if not ok:
        FAILS.append(label)
    print("  %-50s %-13s %12.6g  %s"
          % (label, ("%g" % printed) + unit, computed, "PASS" if ok else "**FAIL**"))


# ======================================================================
# 1.  Inputs
# ======================================================================
alpha = 1.0 / 137.035999177        # derived: the Peierls-Nabarro tunnelling amplitude
m_e = 0.51099895069                # MeV, CODATA 2022: the scale-setting input
N_c = 3.0

m_mu = 105.6583755                 # MeV, PDG, for comparison only
m_tau = 1776.93
MZ = 91187.9

m_0 = m_e / alpha                  # node mass


def dalpha_inv(Q2, m, Nc=1.0, Q=1.0):
    """One-loop shift in the inverse coupling at SPACELIKE momentum transfer Q.

    Delta-alpha^-1_f(-Q^2) = (Nc Q_f^2 / pi) Int_0^1 dx 2x(1-x) ln[1 + x(1-x) Q^2/m^2].

    Positive by construction: a charged loop screens, so the inverse coupling falls.
    """
    v = Q2 / m**2
    f = lambda x: 2 * x * (1 - x) * np.log(1.0 + x * (1 - x) * v)
    val, _ = integrate.quad(f, 0.0, 1.0, limit=400)
    return (Nc * Q**2 / np.pi) * val


print("\n=== the node mass and the scales built on it ===")
chk("m_0 = m_e/alpha", 70.0253, m_0, unit=" MeV")
chk("Lambda_QCD = pi m_0", 220.0, np.pi * m_0, tol=0.5, unit=" MeV")
chk("sqrt(sigma) = 2 pi m_0, the string tension", 440.0, 2 * np.pi * m_0, tol=0.5,
    unit=" MeV")
chk("E_BZ = pi sqrt2 m_0, the D4 zone face", 311.1, np.pi * np.sqrt(2) * m_0,
    unit=" MeV")

# ======================================================================
# 2.  The light hadron sector
# ======================================================================
m_pi_iso = 2.0 * m_0 - 4.0 * m_e                       # spectral mass formula
m_pi = m_pi_iso + (139.57039 - 134.9768) / 3.0         # charged, splitting imported
f_pi = N_c**0.25 * m_0 * (1.0 + alpha / np.pi)
RHO_EIG = 4.891                                        # cluster graph eigenvalue
m_rho = 11.0 * m_0 - 11.0 * (4.0 - RHO_EIG) * m_e
sigma = (2.0 * np.pi * m_0)**2
m_a1 = np.sqrt(m_rho**2 + 2.0 * np.pi * sigma)         # Regge, L = 1


def beta_pi(s, mpi=m_pi):
    """Pion velocity in the pair centre-of-mass frame."""
    return np.sqrt(np.maximum(0.0, 1.0 - 4.0 * mpi**2 / np.asarray(s, float)))


g2 = m_rho**2 / (2.0 * f_pi**2)                        # KSRF
Gamma_rho = g2 * m_rho * beta_pi(m_rho**2)**3 / (48.0 * np.pi)

print("\n=== the light hadrons ===")
chk("m_pi, isospin average", 138.007, m_pi_iso, unit=" MeV")
chk("f_pi = 3^(1/4) m_0 (1 + alpha/pi)", 92.4, f_pi, unit=" MeV")
chk("m_rho, spectral crossed-fault closure", 775.3, m_rho, unit=" MeV")
chk("m_rho, the Z^2/(Z+1) charge-screening estimate", 775.7,
    (12.0**2 / 13.0) * m_0, unit=" MeV")
chk("m_a1 on the Regge trajectory", 1348.0, m_a1, tol=0.6, unit=" MeV")
chk("  m_a1^2 as displayed, in 1e6 MeV^2", 1.817, m_a1**2 / 1e6)
chk("Gamma_rho from KSRF", 147.1, Gamma_rho, unit=" MeV")
chk("  as a fraction of the mass", 0.190, Gamma_rho / m_rho)

# ======================================================================
# 3.  The rho leptonic width, both routes
# ======================================================================
print("\n=== the rho leptonic width ===")
Gee_tree = 1e3 * 8.0 * np.pi * alpha**2 * f_pi**2 / (3.0 * m_rho)      # keV
chk("Gamma_ee tree level, VMD universality", 4.91, Gee_tree)
chk("  shortfall against the measured 7.04 keV", 30.0,
    100 * (1 - Gee_tree / 7.04), tol=0.6, unit="%")

wein = 1.0 / (1.0 - m_rho**2 / m_a1**2)
chk("m_rho^2/m_a1^2", 0.331, m_rho**2 / m_a1**2)
chk("Weinberg factor, on the printed 0.331", 1.495, 1.0 / (1.0 - 0.331))
chk("Gamma_ee, sum-rule route in the chiral limit", 7.34, Gee_tree * wein, tol=0.006)

Gee_closed = 1e3 * 4.0 * np.pi * alpha**2 * f_pi**2 / m_rho
chk("Gamma_ee, closed form in the chiral limit", 7.36, Gee_closed, tol=0.006)
chk("  the two routes agree iff m_a1^2 = 3 m_rho^2, i.e. 8 pi^3 = 2(m_rho/m_0)^2",
    248.05, 8 * np.pi**3)
chk("    against", 245.16, 2 * (m_rho / m_0)**2)
chk("    so the near-identity holds to", 1.2,
    100 * (8 * np.pi**3 / (2 * (m_rho / m_0)**2) - 1), tol=0.06, unit="%")

fin = 1.0 - m_pi_iso**2 / m_a1**2                      # derived finite-pion correction
chk("finite-pion correction m_pi^2/m_a1^2", 0.0105, 1.0 - fin, tol=5e-5)
chk("  as a percentage", 1.05, 100 * (1.0 - fin), tol=0.006, unit="%")
chk("  the m_pi^2/m_rho^2 a dimensional estimate would give", 0.032,
    m_pi_iso**2 / m_rho**2)
chk("Gamma_ee sum-rule route, corrected", 7.26, 4.91 * 1.495 * fin, tol=0.006)
chk("  its excess over the pre-2023 average", 3.2,
    100 * (4.91 * 1.495 * fin / 7.040 - 1), tol=0.06, unit="%")
Gee_rho = Gee_closed * fin
chk("Gamma_ee closed-form route, corrected", 7.288, Gee_rho, tol=0.0011)
chk("  its excess over the pre-2023 average", 3.5,
    100 * (Gee_rho / 7.040 - 1), tol=0.06, unit="%")
chk("pion mass that would close the residual entirely", 270.0,
    m_a1 * np.sqrt(1.0 - 7.040 / (Gee_tree * wein)), tol=3.0, unit=" MeV")

chk("Gamma_ee(phi), charge ratio 2/9", 1.24,
    1e3 * (2 / 9) * 4 * np.pi * alpha**2 * f_pi**2 / 1019.460, tol=0.006)
chk("omega-rho split, charge scaling 12 m_e", 6.1, 12 * m_e)
chk("omega-rho split, spectral 11 x 1.350 m_e", 7.59, 11 * 1.350 * m_e)

# ======================================================================
# 4.  The a_1 / rho' doublet, and where the trajectory stops working
# ======================================================================
print("\n=== the chirally split doublet ===")
M_A1_SPECTRAL = 1227.37     # Born-cluster eigenvalue route, an independent derivation
centre1 = np.sqrt((1465.0**2 + 1230.0**2) / 2)
chk("doublet centre from the PDG pair", 1353.0, centre1, tol=0.6, unit=" MeV")
chk("  trajectory residual against it", 0.33, abs(100 * (m_a1 / centre1 - 1)),
    tol=0.006, unit="%")
chk("  in standard deviations", 0.2, abs(m_a1 - centre1) / 23.0, tol=0.06)
rhoprime = np.sqrt(2 * m_a1**2 - M_A1_SPECTRAL**2)
chk("rho' predicted from the splitting", 1459.0, rhoprime, tol=0.6, unit=" MeV")
chk("  residual against the observed 1465", -0.42, 100 * (rhoprime / 1465 - 1))
chk("half-splitting implied", 0.311, (m_a1**2 - M_A1_SPECTRAL**2) / 1e6,
    unit=" GeV^2")
chk("  measured", 0.317, ((1465.0**2 - 1230.0**2) / 2) / 1e6, unit=" GeV^2")

centre2 = np.sqrt((1720.0**2 + 1655.0**2) / 2)
m2 = np.sqrt(m_rho**2 + 2 * 2 * np.pi * sigma)
chk("second doublet centre", 1688.0, centre2, tol=0.6, unit=" MeV")
chk("  trajectory there", 1742.0, m2, tol=0.6, unit=" MeV")
chk("  overshoot", 3.2, 100 * (m2 / centre2 - 1), tol=0.06, unit="%")
chk("2 pi sigma", 1.216, 2 * np.pi * sigma / 1e6, unit=" GeV^2")
chk("  first step, rho to the first centre", 1.228, (centre1**2 - m_rho**2) / 1e6,
    unit=" GeV^2")
chk("  second step", 1.019, (centre2**2 - centre1**2) / 1e6, unit=" GeV^2")
chk("  slope overshoot on the second step", 19.0,
    100 * (2 * np.pi * sigma / 1e6 / 1.019 - 1), tol=0.6, unit="%")

# ======================================================================
# 5.  The tower, and where the bound-state treatment stops
# ======================================================================
print("\n=== the vector tower and the duality matching ===")
TPS = 2.0 * np.pi * sigma
m1 = np.sqrt(m_rho**2 + TPS)
chk("first recurrence against the observed rho(1450)", -8.0, 100 * (m1 / 1465 - 1),
    tol=0.6, unit="%")
chk("second recurrence against the observed rho(1700)", 1.3, 100 * (m2 / 1720 - 1),
    tol=0.06, unit="%")
chk("constituent mass N_c^2 m_0 / 2", 315.0, N_c**2 * m_0 / 2, tol=0.6, unit=" MeV")

s_over = np.pi * sigma * m_rho / Gamma_rho
chk("overlap bound pi sigma m_rho / Gamma_rho", 3.21, s_over / 1e6, tol=0.006,
    unit=" GeV^2")
chk("  its closed form in f_pi and m_rho", 3.21,
    96 * np.pi**2 * sigma * f_pi**2 / (m_rho**2 * beta_pi(m_rho**2)**3) / 1e6,
    tol=0.006, unit=" GeV^2")
chk("  crossing", 1.79, np.sqrt(s_over) / 1e3, unit=" GeV")
# level by level, starting from the rho itself: Gamma_n / Delta m_n with
# Gamma_n = (Gamma_rho/m_rho) m_n and Delta m_n = pi sigma / m_n
for n, printed in enumerate((0.19, 0.57, 0.95, 1.33, 1.70)):
    mn = np.sqrt(m_rho**2 + n * TPS)
    chk("  overlap ratio at level n = %d" % n, printed,
        (Gamma_rho / m_rho) * mn**2 / (np.pi * sigma))
chk("Gamma/m that reaching N = 3 would require", 0.143,
    (np.pi * sigma) / (m_rho**2 + 3 * TPS) * 1.0, tol=0.0011)
for N, printed in ((1, 1558.0), (2, 1909.0)):
    chk("matching point at N = %d" % N, printed,
        np.sqrt(m_rho**2 + (N + 0.5) * TPS), tol=1.0, unit=" MeV")

# the pure-linear ceiling, tested against measured radial pairs
print("\n  the 1/m^2 ceiling against measurement (Gamma_ee = B_ee x Gamma_tot, PDG):")
for nm, gn, g0, mn_, m0_, printed in (
        ("psi(2S)/J-psi", 7.95e-3 * 293e-3, 5.971e-2 * 92.6e-3, 3686.097, 3096.900, 0.60),
        ("Upsilon(2S)/(1S)", 1.91e-2 * 31.98e-3, 2.39e-2 * 54.02e-3, 10023.4, 9460.40, 0.53),
        ("Upsilon(3S)/(1S)", 2.18e-2 * 20.32e-3, 2.39e-2 * 54.02e-3, 10355.1, 9460.40, 0.41)):
    chk("    measured / (1/m^2 law), " + nm, printed, (gn / g0) / ((m0_ / mn_)**2),
        tol=0.006)

# ======================================================================
# 6.  The running
# ======================================================================
print("\n=== the running ===")
lep = sum(dalpha_inv(MZ**2, m) for m in (m_e, m_mu, m_tau))
chk("sum_f Delta-alpha^-1_f at M_Z, three leptons", 4.306, lep, tol=0.002)
chk("  as a fraction, Delta-alpha_lep at one loop", 0.03142, lep * alpha, tol=6e-6)

E_BZ = np.pi * np.sqrt(2) * m_0
chk("ln(pi sqrt2 / alpha)", 6.412, np.log(np.pi * np.sqrt(2) / alpha))
ll_e = (2.0 / (3 * np.pi)) * np.log(np.pi * np.sqrt(2) / alpha)
ll_mu = (1.0 / (3 * np.pi)) * np.log(E_BZ**2 / m_mu**2)
chk("electron shift at the zone, leading log", -1.36, -ll_e)
chk("muon shift at the zone, leading log", -0.23, -ll_mu)
chk("  their sum", -1.59, -(ll_e + ll_mu), tol=0.006)
chk("alpha^-1(E_BZ), leading log", 135.4, 137.036 - ll_e - ll_mu, tol=0.06)
chk("alpha^-1(E_BZ), on-shell thresholds restored", 135.7,
    1.0 / alpha - dalpha_inv(E_BZ**2, m_e) - dalpha_inv(E_BZ**2, m_mu), tol=0.06)
chk("two-pion threshold 2 m_pi", 279.0, 2 * m_pi, tol=0.6, unit=" MeV")
chk("  its distance below the zone face", 32.0, E_BZ - 2 * m_pi, tol=0.6, unit=" MeV")

chk("Landau exponent 3 pi / alpha, in log10", 561.0, (3 * np.pi / alpha) / np.log(10),
    tol=0.5)
chk("  q_L^2, in log10 GeV^2", 554.0,
    np.log10(m_e**2 * 1e-6) + (3 * np.pi / alpha) / np.log(10), tol=0.5)

# ======================================================================
# 7.  How little of the spectral function the running resolves
# ======================================================================
print("\n=== the spectral-sensitivity table ===")


def da_disp(Rfun, s0, Q2=MZ**2, n=400001):
    """Once-subtracted dispersion integral: Delta-alpha^-1 at spacelike -Q^2."""
    t = np.linspace(np.log(s0 * (1 + 1e-12)), np.log(1e18), n)
    s = np.exp(t)
    return (Q2 / (3 * np.pi)) * np.trapezoid(Rfun(s) / (s + Q2), t)


s0 = 4 * m_e**2
b = lambda s: np.sqrt(np.maximum(0.0, 1.0 - 4 * m_e**2 / s))
ref = None
for lbl, f, printed in (
        ("beta(3-beta^2)/2, the point-charge form", lambda s: b(s) * (3 - b(s)**2) / 2, 2.38918),
        ("beta, the threshold power alone", lambda s: b(s), 2.35381),
        ("1, a bare step at threshold", lambda s: np.ones_like(s), 2.41892),
        ("beta^3, a P-wave threshold", lambda s: b(s)**3, 2.28307),
        ("beta^3/4, a scalar pair", lambda s: b(s)**3 / 4, 0.57077)):
    v = da_disp(f, s0)
    if ref is None:
        ref = v
    chk("  " + lbl, printed, v, tol=0.002)
chk("dispersive form against the Feynman-parameter form", 2.38918,
    dalpha_inv(MZ**2, m_e), tol=0.002)
chk("  step function moves alpha^-1(M_Z) by", 0.030, 2.41892 - 2.38918, tol=0.0011)
chk("  a pure beta^3 moves it by", 4.4, 100 * (1 - 2.28307 / 2.38918), tol=0.06,
    unit="%")

# ======================================================================
# 8.  The defect's spin, read backwards from the measurement
# ======================================================================
print("\n=== the spin read backwards ===")
AINV0 = 1.0 / alpha
DA_LEP3, DA_HAD, DA_TOP = 0.031498, 0.027526, -0.7e-4
chk("alpha^-1(M_Z), the framework's value", 128.957,
    AINV0 * (1 - DA_LEP3 - DA_HAD - DA_TOP), tol=0.0011)
d_ferm = AINV0 - 128.946
chk("Delta-alpha^-1 a Dirac fermion must supply", 8.090, d_ferm, tol=0.0011)
chk("  a complex scalar, exactly one quarter", 2.022, d_ferm / 4, tol=0.0011)
chk("  the coupling it would leave", 135.013, AINV0 - d_ferm / 4, tol=0.0011)
chk("  its departure", 404.0, (AINV0 - d_ferm / 4 - 128.946) / 0.015, tol=1.5,
    unit=" sigma")

# ======================================================================
# 9.  The node form factor
# ======================================================================
print("\n=== the node form factor ===")
d_mis = 1.0 / np.sqrt(3.0)                      # misfit period, in units of ell
w = (np.log(1.0 / alpha) / (2 * np.pi)) * d_mis  # dressed core half-width
chk("dressed w/d = ln(1/alpha)/(2 pi)", 0.783, np.log(1 / alpha) / (2 * np.pi))
chk("  dressed w", 0.452, w, unit=" ell")
chk("  bare Cosserat w = (pi/4) d", 0.453, (np.pi / 4) * d_mis, unit=" ell")
chk("  the gap between them", 0.3, 100 * ((np.pi / 4) * d_mis / w - 1), tol=0.06,
    unit="%")
chk("sech half-width 2 ln(1+sqrt2)/(pi w)", 1.24,
    2 * np.log(1 + np.sqrt(2)) / (np.pi * w), tol=0.006, unit="/ell")
# FCC {111}: |G| = 2 pi sqrt3 / a with cubic constant a = sqrt2 ell, so |G| = pi sqrt6 / ell
G111 = np.pi * np.sqrt(6.0)
chk("2 |G_111| d / (2 pi), an identity", np.sqrt(2), 2 * G111 * d_mis / (2 * np.pi),
    tol=1e-9)
chk("|F(G_111)|^2 = alpha^sqrt2, at the dressed width", 9.5e-4, alpha**np.sqrt(2),
    tol=5e-6)
chk("  at the bare pi/4 width instead, no power of alpha", 9.3e-4,
    np.exp(-2 * G111 * (np.pi / 4) * d_mis), tol=5e-6)
chk("1/w in energy", 155.0, m_0 / w, tol=3.0, unit=" MeV")

# ======================================================================
# 10.  The closure census
# ======================================================================
print("\n=== the closure census ===")
resid = AINV0 * (1 - DA_LEP3 - DA_HAD - DA_TOP) - 128.946
sig = np.hypot(0.009, 0.015)
chk("residual against the measurement", 0.011, resid, tol=0.0011)
chk("  its uncertainty, in quadrature", 0.017, sig, tol=0.0011)
bound = 0.011 + 2 * sig            # two-sigma ceiling on an unaccounted contribution
chk("  two-sigma ceiling on a missed species", 0.046, bound, tol=0.0011)
for m, pr_nq, pr_q in ((131.0, 0.038, 0.19), (148.0, 0.039, 0.20),
                       (1e3, 0.059, 0.24), (1e4, 0.154, 0.39), (45e3, 0.733, 0.86)):
    nq = bound / dalpha_inv(MZ**2, m)
    chk("  N_c Q^2 bound at %6.0f MeV" % m, pr_nq, nq, tol=0.0011)
    chk("    as a charge, colour singlet", pr_q, np.sqrt(nq), tol=0.006, unit=" e")

# ======================================================================
print("\n" + "=" * 82)
if FAILS:
    print("%d literal(s) did not reproduce:" % len(FAILS))
    for f in FAILS:
        print("   ", f)
    sys.exit(1)
print("all literals reproduce")
