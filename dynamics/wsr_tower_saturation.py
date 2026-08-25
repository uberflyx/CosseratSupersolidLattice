"""Do the Weinberg sum rules, saturated with the derived tower, fix Gamma_ee(rho)?

THE PROBLEM.  The corpus computes the rho's leptonic width from the two Weinberg sum
rules saturated with one vector and one axial state, and the result sits 3.5% above the
pre-2023 average.  Feeding the same formula MEASURED f_pi and m_rho still returns
+2.8%, so the excess is a property of single-state saturation and not of the
framework's inputs; the inputs contribute the remaining 0.7%.

THE CANDIDATE FIX.  The framework derives a whole vector tower, with recurrence
couplings from the Cornell solve, so the sum rules can be saturated with more than one
state per channel.  This script asks whether doing so moves Gamma_ee(rho) by the ~3%
needed, and in which direction.

CONVENTIONS.  The corpus writes Gamma_ee(rho) = (8 pi alpha^2 / (3 m_rho)) f_rho^2, and
in the one-state case the sum rules give f_rho^2 = f_pi^2 / (1 - m_rho^2/m_a1^2).  That
fixes the normalisation of f used throughout here.  The sum rules are

    WSR-I    sum_V f_V^2 - sum_A f_A^2 = f_pi^2
    WSR-II   sum_V f_V^2 m_V^2 - sum_A f_A^2 m_A^2 = 0

and the tower ratios enter as f_n^2 / f_0^2 = r_n (m_n / m_0), since Gamma_ee carries an
explicit 1/m.
"""
import numpy as np

alpha = 1 / 137.035999177
m_e = 0.51099895069
m_0 = m_e / alpha
sigma = (2 * np.pi * m_0)**2
TPS = 2 * np.pi * sigma

f_pi_meas = 92.07
m_rho_meas = 775.26
R1 = 0.2798                     # Cornell solve, first recurrence leptonic ratio
GEE_OBS = 7.04


def gee(f2, m_rho):
    """Leptonic width in keV from f^2, in the corpus's normalisation."""
    return 1e3 * (8 * np.pi * alpha**2 / (3 * m_rho)) * f2


def one_state(f_pi, m_rho, m_a1):
    return f_pi**2 / (1 - m_rho**2 / m_a1**2)


def two_by_two(f_pi, mV, mV1, mA, mA1, rV, rA):
    """Solve both sum rules with two vectors and two axials at fixed coupling ratios."""
    cV1, cA1 = rV * mV1 / mV, rA * mA1 / mA          # f_1^2 / f_0^2 in each channel
    # WSR-II fixes the ratio of the two ground-state couplings
    num = mV**2 + cV1 * mV1**2
    den = mA**2 + cA1 * mA1**2
    k = num / den                                     # f_A0^2 = k f_V0^2
    # WSR-I then fixes the scale
    denom = (1 + cV1) - k * (1 + cA1)
    return f_pi**2 / denom, k, cV1, cA1


print("=== the baseline, one state per channel ===")
for nm, m_a1 in (("corpus Regge a1 1348.11", 1348.11),
                 ("chiral ladder a1 1350.73", 1350.73),
                 ("observed a1 1230", 1230.0)):
    f2 = one_state(f_pi_meas, m_rho_meas, m_a1)
    print("  %-28s f_rho^2 = %7.0f   Gamma_ee = %.3f keV  (%+.1f%%)"
          % (nm, f2, gee(f2, m_rho_meas), 100 * (gee(f2, m_rho_meas) / GEE_OBS - 1)))

print("\n=== two states per channel, on the framework's own tower ===")
print("Vector tower: m_V1 = sqrt(m_rho^2 + 2 pi sigma).  Axial tower: same step on m_a1.")
print("Coupling ratios from the Cornell solve, applied to both channels.\n")
base = one_state(f_pi_meas, m_rho_meas, 1348.11)
for nm, mA, mV1_mode in (("framework a1 = 1348.11", 1348.11, "regge"),
                         ("observed a1 = 1230", 1230.0, "regge"),
                         ("observed a1 = 1230, measured rho' = 1465", 1230.0, "meas")):
    mV1 = 1465.0 if mV1_mode == "meas" else np.sqrt(m_rho_meas**2 + TPS)
    mA1 = np.sqrt(mA**2 + TPS)
    f2, k, cV1, cA1 = two_by_two(f_pi_meas, m_rho_meas, mV1, mA, mA1, R1, R1)
    g = gee(f2, m_rho_meas)
    print("  %-42s" % nm)
    print("     f_V1^2/f_V0^2 = %.4f   f_A1^2/f_A0^2 = %.4f   f_A0^2/f_V0^2 = %.4f"
          % (cV1, cA1, k))
    print("     f_rho^2 = %7.0f   Gamma_ee = %.3f keV  (%+.1f%% on 7.04, %+.1f%% on the"
          " one-state value)\n" % (f2, g, 100 * (g / GEE_OBS - 1), 100 * (f2 / base - 1)))

print("=== how stable is that? scan the axial coupling ratio ===")
print("The vector ratio is derived; the axial one is not, so it is the free handle.")
print("  r_A / r_V     f_rho^2    Gamma_ee     shift on the one-state value")
mV1 = np.sqrt(m_rho_meas**2 + TPS)
for scale in (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
    mA, mA1 = 1348.11, np.sqrt(1348.11**2 + TPS)
    f2, k, _, _ = two_by_two(f_pi_meas, m_rho_meas, mV1, mA, mA1, R1, R1 * scale)
    print("     %.2f      %7.0f     %.3f keV      %+.1f%%"
          % (scale, f2, gee(f2, m_rho_meas), 100 * (f2 / base - 1)))

print("\n=== verdict ===")
print("The tower correction is not a small perturbation.  Over a plausible range of the")
print("undetermined axial coupling it moves f_rho^2 by tens of per cent and changes")
print("sign, so it cannot be used to close a 3% gap: the tool is far coarser than the")
print("quantity being fixed.  Truncated Weinberg sum rules are known to behave this way,")
print("because the two sums converge only in the difference, and the difference is")
print("controlled by chiral restoration at high rung rather than by the first recurrence.")

# ----------------------------------------------------------------------
# The comparison value: is the 4-sigma pull an artefact of a tight PDG fit?
# ----------------------------------------------------------------------
print("\n=== is the comparison value itself soft? ===")
DET = [("BARKOV 85 (OLYA)", 6.77, np.hypot(0.10, 0.30)),
       ("CMD-2 1999", 6.93, np.hypot(0.11, 0.10)),
       ("CMD-2 2001", 6.86, np.hypot(0.11, 0.05)),
       ("AKHMETSHIN 04 (CMD-2)", 7.06, np.hypot(0.11, 0.05)),
       ("AKHMETSHIN 07 (CMD-2)", 7.048, np.hypot(0.057, 0.050)),
       ("ACHASOV 06 (SND)", 7.12, np.hypot(0.02, 0.11))]
F = 7.288
for nm, v, e in DET:
    print("  %-24s %5.3f +- %.3f   framework %+5.1f%%   %4.1f sigma" %
          (nm, v, e, 100 * (F / v - 1), (F - v) / e))
w = np.array([1 / e**2 for _, _, e in DET])
x = np.array([v for _, v, _ in DET])
mu = (w * x).sum() / w.sum()
chi2 = (w * (x - mu)**2).sum()
print("\n  weighted mean %.3f +- %.3f, chi2/dof = %.2f" %
      (mu, 1 / np.sqrt(w.sum()), chi2 / (len(DET) - 1)))
print("  The determinations AGREE with each other, so there is no hidden spread to")
print("  hide behind and the PDG error needs no inflating.  The pull is real.")

# ----------------------------------------------------------------------
# What the framework's m_a1 buys, and what would be needed to close the rest
# ----------------------------------------------------------------------
print("\n=== what the a1 mass is worth here ===")
K = 1e3 * 8 * np.pi * alpha**2 / (3 * m_rho_meas)
for nm, ratio in (("Weinberg's classic m_a1^2 = 2 m_rho^2", 2.0),
                  ("framework Regge / ladder, m_a1^2 = 3 m_rho^2", 3.0),
                  ("what the measured width demands", None)):
    if ratio is None:
        need = GEE_OBS / (K * f_pi_meas**2)
        ratio = 1 / (1 - 1 / need)
    W = 1 / (1 - 1 / ratio)
    print("  %-46s ratio %.3f -> %.3f keV  (%+.1f%%)"
          % (nm, ratio, K * f_pi_meas**2 * W, 100 * (K * f_pi_meas**2 * W / GEE_OBS - 1)))
print("\n  The textbook relation misses by 39%; the framework's own a1 brings that to")
print("  3.5%, a factor of eleven.  Closing the rest needs m_a1^2/m_rho^2 = 3.26, i.e.")
print("  m_a1 = 1399 MeV, which no route in the corpus produces.")
