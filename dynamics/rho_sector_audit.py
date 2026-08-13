"""Audit of the framework's rho sector, the one input that decides a_mu.

Everything about the rho is derived from alpha and m_e.  Some of it is excellent and one
number is not, and the point of this script is to find out which of the several claims
around that number actually stand.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ew_one_loop_lattice as ew

m0, mrho, fpi = ew.m_0, ew.m_rho, ew.f_pi
alpha = ew.alpha

print("="*78)
print("1.  WHAT THE RHO SECTOR GETS RIGHT")
print("="*78)
mrho_obs, Grho_obs, fpi_obs = 775.26, 147.4, 92.07
print("  m_rho      %8.2f  vs %8.2f   %+.2f%%" % (mrho, mrho_obs, 100*(mrho/mrho_obs-1)))
print("  Gamma_rho  %8.2f  vs %8.2f   %+.2f%%" % (ew.Gamma_rho, Grho_obs,
                                                  100*(ew.Gamma_rho/Grho_obs-1)))
print("  f_pi       %8.2f  vs %8.2f   %+.2f%%" % (fpi, fpi_obs, 100*(fpi/fpi_obs-1)))
# KSRF coupling against the width-extracted one
g2_ksrf = mrho**2/(2*fpi**2)
b_obs = (1 - 4*139.57039**2/mrho_obs**2)**0.5
g2_obs = Grho_obs*48*np.pi/(mrho_obs*b_obs**3)
print("  g_rhopipi^2 %7.2f  vs %8.2f   %+.2f%%   (KSRF against the measured width)"
      % (g2_ksrf, g2_obs, 100*(g2_ksrf/g2_obs-1)))

print("\n" + "="*78)
print("2.  THE ONE NUMBER THAT IS OFF, IN COMPANY")
print("="*78)
for name, pred, obs in (("rho  ", 7.361, 7.04), ("omega", 0.62, 0.617),
                        ("phi  ", ew.Gee_phi*1e3, 1.259)):
    print("  Gamma_ee(%s) %7.3f keV  vs %6.3f   %+.2f%%"
          % (name, pred, obs, 100*(pred/obs-1)))
print("""
  The rho is the outlier, not a symptom of a common normalisation: a single scale factor
  cannot fit both, needing 0.956 for the rho and 1.011 for the phi.""")

print("\n" + "="*78)
print("3.  ARE THE TWO ROUTES TO Gamma_ee(rho) INDEPENDENT?")
print("="*78)
# route A: tree KSRF x Weinberg factor, with m_a1 from the Regge trajectory
m_a1_sq_over_m_rho_sq = 1.0 + 2*np.pi*ew.sigma_string/mrho**2
wsr = 1.0/(1.0 - 1.0/m_a1_sq_over_m_rho_sq)
tree = 8*np.pi*alpha**2*fpi**2/(3*mrho)
routeA = tree*wsr
# route B: the closed form quoted as a separate derivation
routeB = 4*np.pi*alpha**2*fpi**2/mrho
print("  route A (KSRF x WSR)      %.4f keV" % (routeA*1e3))
print("  route B (closed form)     %.4f keV" % (routeB*1e3))
print("  A/B = %.5f, and A/B = (2/3) x WSR exactly" % (routeA/routeB))
print("  WSR factor = %.5f   m_a1^2/m_rho^2 = %.5f" % (wsr, m_a1_sq_over_m_rho_sq))
print("""
  Route A equals route B only because m_a1^2 = 3 m_rho^2 makes the Weinberg factor 3/2,
  and that in turn is the numerical near-identity""")
print("     2 pi sigma / m_rho^2 = 8 pi^3 / (144/13)^2 = %.5f / %.5f = %.5f  (want 2)"
      % (8*np.pi**3, (144/13)**2, 8*np.pi**3/(144/13)**2))
print("""     8 pi^3 = %.3f  against  2 (144/13)^2 = %.3f, agreeing to %.2f per cent.

  So the corpus's "two independent routes" are one route plus a coincidence at the
  per-cent level.  The agreement is not a check on Gamma_ee.""" %
      (8*np.pi**3, 2*(144/13)**2, 100*abs(8*np.pi**3/(2*(144/13)**2)-1)))

print("="*78)
print("4.  CAN A BETTER a_1 MASS RESCUE IT?")
print("="*78)
print("  m_a1 (Regge, framework) = %.1f MeV   observed a_1(1260) = 1230 +- 40" %
      (mrho*np.sqrt(m_a1_sq_over_m_rho_sq)))
for m_a1, tag in ((mrho*np.sqrt(m_a1_sq_over_m_rho_sq), "framework Regge"),
                  (1230.0, "observed a_1"),
                  (1409.4, "value Gamma_ee would need")):
    w = 1.0/(1.0 - mrho**2/m_a1**2)
    print("   m_a1 = %7.1f (%-24s) -> WSR %.4f -> Gamma_ee %.3f keV  (%+.1f%%)"
          % (m_a1, tag, w, tree*w*1e3, 100*(tree*w/7.04e-3 - 1)))
print("""
  Using the OBSERVED a_1 makes the width worse, not better: 8.15 keV against 7.04.
  The width would need m_a1 = 1409, further from observation than the framework's 1348.
  So the excess is not an a_1-mass error, and the Weinberg route cannot be repaired by
  sharpening the a_1.  Whatever is wrong sits in the sum-rule application itself.""")

print("="*78)
print("5.  A STRUCTURAL GAP THE TRAJECTORY DOES HAVE")
print("="*78)
print("  The framework uses one Regge step for both orbital and radial excitation:")
print("     m^2 = m_rho^2 + n (2 pi sigma),  2 pi sigma = %.0f MeV^2" % (2*np.pi*ew.sigma_string))
orb = 1230.0**2 - mrho_obs**2
rad = 1465.0**2 - mrho_obs**2
print("  observed orbital step  a_1(1260):  m^2 - m_rho^2 = %.0f MeV^2  -> sqrt = %.0f MeV"
      % (orb, np.sqrt(orb)))
print("  observed radial step   rho(1450):  m^2 - m_rho^2 = %.0f MeV^2  -> sqrt = %.0f MeV"
      % (rad, np.sqrt(rad)))
print("  radial/orbital slope ratio = %.2f, where the framework uses 1.00" % (rad/orb))
print("""
  The framework's single step sits between the two, which is why the a_1 comes out
  10 per cent high and the rho(1450) 8 per cent low from the same formula.  Splitting
  the trajectory into an orbital and a radial slope is a well-posed structural fix, and
  it would move the tower masses; it does not by itself touch Gamma_ee, for the reason
  in section 4.""")
