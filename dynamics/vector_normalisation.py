"""Is the Gamma_ee excess specific to the rho, or common to all three vector families?

THE QUESTION.  The corpus computes the rho's leptonic width from the Weinberg sum rules
and lands 3.5% high.  It computes the omega and phi from the same VMD expression scaled
by quark-charge weights, and a companion pass fixed the omega with an omega-phi mixing
angle.  That fix works when the two isoscalars are scaled off the MEASURED rho and
fails when they are scaled off the framework's own f_pi, and the difference between the
two routes is exactly the rho's own 3.5%.  So the isoscalars carry information about
whether that 3.5% is a rho defect or a normalisation common to the whole sector.

THE TEST.  Fit one overall normalisation k and one mixing angle delta to all three
measured widths at once.  The rho is isovector and does not mix, so it constrains k
alone; the two isoscalars constrain both.  Three measurements, two parameters, one
degree of freedom.

  Gamma_ee(rho) = k * B_rho
  Gamma_ee(omega) = k * (1/9) * B_omega * (cos d - sqrt2 sin d)^2
  Gamma_ee(phi)   = k * (2/9) * B_phi   * (cos d + sin d / sqrt2)^2

with B_V = 4 pi alpha^2 f_pi^2 / m_V evaluated on the framework's own derived f_pi.

If k comes out near 1, the rho excess is a rho defect and the isoscalars are healthy.
If k comes out near 0.96, the excess is a single normalisation shared by all three, and
one correction removes three residuals rather than one.
"""
import numpy as np
from scipy import optimize

alpha = 1 / 137.035999177
m_e = 0.51099895069
m_0 = m_e / alpha
f_pi = 3**0.25 * m_0 * (1 + alpha / np.pi)          # derived
m_pi_iso = 2 * m_0 - 4 * m_e

M_RHO, M_OMEGA, M_PHI = 775.286, 782.66, 1019.461
FIN = 1 - m_pi_iso**2 / 1348.11**2                   # finite-pion factor, rho channel

# measurements: PDG, isoscalars as B_ee x Gamma_tot
OBS = {"rho": (7.04, 0.06), "omega": (7.38e-5 * 8.68e3, 0.22e-5 * 8.68e3),
       "phi": (2.979e-4 * 4.249e3, 0.033e-4 * 4.249e3)}


def base(mV, w=1.0):
    """Ideal-mixing VMD leptonic width in keV, on the framework's own f_pi."""
    return 1e3 * w * 4 * np.pi * alpha**2 * f_pi**2 / mV


B = {"rho": base(M_RHO) * FIN, "omega": base(M_OMEGA, 1 / 9), "phi": base(M_PHI, 2 / 9)}

print("=== ideal mixing, no normalisation adjustment ===")
for k in ("rho", "omega", "phi"):
    o, e = OBS[k]
    print("  %-6s framework %6.3f   measured %6.3f +- %.3f   %+6.1f%%  %5.1f sigma"
          % (k, B[k], o, e, 100 * (B[k] / o - 1), (B[k] - o) / e))
print("  Opposite signs of very different size, which is what a single mixing angle")
print("  alone cannot produce.\n")


def model(p):
    k, d = p
    Rw = np.cos(d) - np.sqrt(2) * np.sin(d)
    Rp = np.cos(d) + np.sin(d) / np.sqrt(2)
    return {"rho": k * B["rho"], "omega": k * B["omega"] * Rw**2,
            "phi": k * B["phi"] * Rp**2}


def chi2(p):
    m = model(p)
    return sum(((m[k] - OBS[k][0]) / OBS[k][1])**2 for k in OBS)


res = optimize.minimize(chi2, [0.95, 0.05], method="Nelder-Mead",
                        options=dict(xatol=1e-10, fatol=1e-12))
k_fit, d_fit = res.x
print("=== joint fit: one normalisation, one angle, three widths ===")
print("  k     = %.4f   (1 would mean the rho excess is a rho defect)" % k_fit)
print("  delta = %.3f degrees" % np.degrees(d_fit))
print("  chi2  = %.2f on 1 degree of freedom\n" % res.fun)
m = model(res.x)
for kk in ("rho", "omega", "phi"):
    o, e = OBS[kk]
    print("  %-6s fitted %6.3f   measured %6.3f +- %.3f   %+6.2f%%  %5.2f sigma"
          % (kk, m[kk], o, e, 100 * (m[kk] / o - 1), (m[kk] - o) / e))

# what each channel wants on its own
print("\n=== what each channel asks for separately ===")
print("  rho alone wants k = %.4f, i.e. the framework is %.1f%% high" %
      (OBS["rho"][0] / B["rho"], 100 * (B["rho"] / OBS["rho"][0] - 1)))
rat = np.sqrt((OBS["phi"][0] / B["phi"]) / (OBS["omega"][0] / B["omega"]))
t = optimize.brentq(lambda d: (np.cos(d) + np.sin(d) / np.sqrt(2)) /
                    (np.cos(d) - np.sqrt(2) * np.sin(d)) - rat, 0.0, 0.5)
print("  the two isoscalars ALONE fix delta = %.3f degrees from their ratio, with the"
      % np.degrees(t))
Rw = np.cos(t) - np.sqrt(2) * np.sin(t)
print("  normalisation then k = %.4f from the omega." % (OBS["omega"][0] / (B["omega"] * Rw**2)))
print("  So the isoscalar pair and the rho ask for the same normalisation to within")
print("  %.1f per cent, having been given no chance to agree." %
      abs(100 * (OBS["omega"][0] / (B["omega"] * Rw**2) / (OBS["rho"][0] / B["rho"]) - 1)))

print("\n=== three candidate sources for k, all excluded ===")
print("  (a) the chiral-limit decay constant.  eq-gee-final is a chiral-limit relation")
print("      fed the physical f_pi, so F_0 is the natural suspect.  Direction right,")
print("      size badly wrong: k wants F_0/f_pi = %.4f where the lattice sits near"
      % np.sqrt(k_fit))
for r in (0.93, 0.94, 0.95):
    print("      0.94.  At %.2f, k = %.4f and the rho lands %+.1f%%."
          % (r, r * r, 100 * (B["rho"] * r * r / OBS["rho"][0] - 1)))
x = np.log(np.sqrt(k_fit) * 3**0.25) / np.log(3)
print("  (b) the N_c exponent of eq-fpi-exponent.  k would need N_c^%.4f where the" % x)
print("      topological susceptibility confirms N_c^0.249.  Excluded.")
print("  (c) each channel with its OWN chiral partner instead of a common Weinberg")
print("      factor.  This fails differently and the failure is informative:")
for nm, w, mV, mA in (("rho", 1.0, M_RHO, 1348.11), ("omega", 1 / 9, M_OMEGA, 1281.9),
                      ("phi", 2 / 9, M_PHI, 1426.3)):
    W = 1 / (1 - mV**2 / mA**2)
    G = w * (2 / 3) * base(mV) * W
    print("      %-6s m_A = %7.1f  W = %.4f  ->  %6.3f keV vs %6.3f  (%+.1f%%)"
          % (nm, mA, W, G, OBS[nm][0], 100 * (G / OBS[nm][0] - 1)))
print("      The two isoscalars overshoot by the same 34%, agreeing to half a per")
print("      cent, where the rho sits at +4.2%.  That says the isoscalar pair shares")
print("      a factor the rho does not, which is the opposite of a sector-wide k, so")
print("      the two readings cannot both be the whole story.")
print("\n  What survives is the measurement and not its explanation: a single factor")
print("  near 0.96 is asked for independently by the rho and by the isoscalar pair,")
print("  and no candidate yet accounts for it.")
