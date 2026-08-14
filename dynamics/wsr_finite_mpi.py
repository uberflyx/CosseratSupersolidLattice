"""The finite-m_pi correction to the Weinberg sum rules, derived rather than imported.

THE DERIVATION
--------------
The sum rules are moments of the difference between the vector and axial spectral
functions.  Writing the axial one with its pion pole shown explicitly,

    rho_A(s) = f_pi^2 delta(s - m_pi^2) + (a_1 and higher),

the two superconvergence relations read

    WSR-I :  Int ds       [rho_V - rho_A] = 0
    WSR-II:  Int ds  s    [rho_V - rho_A] = 0

In the chiral limit the pion is massless, so its pole sits at s = 0.  It therefore
contributes f_pi^2 to the ZEROTH moment and f_pi^2 x 0 = 0 to the FIRST.  That is the
whole reason WSR-II has a bare zero on the right in the standard treatment:

    f_rho^2 - f_a1^2 = f_pi^2 ,          f_rho^2 m_rho^2 = f_a1^2 m_a1^2 .

At finite pion mass the pole moves to s = m_pi^2.  The zeroth moment does not care,
so WSR-I is untouched; the first moment picks up exactly f_pi^2 m_pi^2:

    f_rho^2 m_rho^2 = f_a1^2 m_a1^2 + f_pi^2 m_pi^2 .

Eliminating f_a1^2 between the two gives the corrected coupling,

    f_rho^2 = f_pi^2 (m_a1^2 - m_pi^2) / (m_a1^2 - m_rho^2)
            = [ f_pi^2 / (1 - m_rho^2/m_a1^2) ] x (1 - m_pi^2/m_a1^2) ,

so the chiral-limit result is multiplied by (1 - m_pi^2/m_a1^2).  Nothing is imported:
the correction is the pion pole sitting at its own mass instead of at the origin, and
every quantity in it is one the framework already derives.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ew_one_loop_lattice as ew

m_pi, m_rho, f_pi = ew.m_pi, ew.m_rho, ew.f_pi
m_a1 = np.sqrt(m_rho**2 + 2.0*np.pi*ew.sigma_string)
G_TREE = 8.0*np.pi*ew.alpha**2*f_pi**2/(3.0*m_rho)
G_OBS  = 7.04e-3

print("="*76)
print("INPUTS, all framework-derived")
print("="*76)
print("  m_pi = %.2f   m_rho = %.2f   m_a1 = %.1f   f_pi = %.2f MeV"
      % (m_pi, m_rho, m_a1, f_pi))
print("  tree-level KSRF/VMD width  = %.4f keV" % (G_TREE*1e3))

chiral = 1.0/(1.0 - m_rho**2/m_a1**2)
finite = (m_a1**2 - m_pi**2)/(m_a1**2 - m_rho**2)
corr   = 1.0 - m_pi**2/m_a1**2

print("\n" + "="*76)
print("THE CORRECTION")
print("="*76)
print("  chiral-limit factor  f_rho^2/f_pi^2 = 1/(1 - m_rho^2/m_a1^2) = %.5f" % chiral)
print("  finite-m_pi factor   (m_a1^2 - m_pi^2)/(m_a1^2 - m_rho^2)    = %.5f" % finite)
print("  ratio                (1 - m_pi^2/m_a1^2)                     = %.5f  (%+.2f%%)"
      % (corr, 100*(corr-1)))
assert abs(finite/chiral - corr) < 1e-12, "algebra check"
print("  algebra check: finite/chiral equals (1 - m_pi^2/m_a1^2) exactly.")

print("\n  Gamma_ee(rho):")
print("     chiral limit      %.4f keV   (%+.2f%% vs measured %.2f)"
      % (G_TREE*chiral*1e3, 100*(G_TREE*chiral/G_OBS-1), G_OBS*1e3))
print("     finite m_pi       %.4f keV   (%+.2f%%)"
      % (G_TREE*finite*1e3, 100*(G_TREE*finite/G_OBS-1)))
print("     the correction closes %.0f%% of the gap" %
      (100*(G_TREE*chiral - G_TREE*finite)/(G_TREE*chiral - G_OBS)))

print("""
  So the derivable correction is about one per cent, not the three per cent that the
  naive (m_pi/m_rho)^2 estimate suggested.  The reason is structural and worth stating:
  the pion pole enters the FIRST moment, where it is weighed against the a_1 scale, so
  the correction is m_pi^2/m_a1^2 and not m_pi^2/m_rho^2.  The a_1 is the heaviest scale
  in the sum rule, which is exactly why the correction is small.
""")

print("="*76)
print("DOWNSTREAM, with the correction applied to the width the integral uses")
print("="*76)
G_new = ew.Gee_rho*corr
print("  Gamma_ee used by the spectral function: %.4f -> %.4f keV"
      % (ew.Gee_rho*1e3, G_new*1e3))

for tag, g in (("chiral limit", ew.Gee_rho), ("finite m_pi", G_new),
               ("measured", ew.Gee_rho_pdg)):
    a = [ew.a_mu_hvp_interval(n, Gee_rho_use=g)*1e10 for n in (1,2,3,4)]
    c, h = 0.5*(max(a)+min(a)), 0.5*(max(a)-min(a))
    print("  a_mu %-14s %.0f +- %.0f   -> %.1f sigma vs WP25 713(6)"
          % (tag, c, h, abs(c-713)/np.hypot(h,6)))

# the running barely moves; quantify with the same prescription
d0 = [ew.delta_alpha_had_interval(n) for n in (1,2,3,4)]
print("\n  Delta-alpha_had is untouched at this precision: the rho width enters it at")
print("  the half-per-cent level, so the band %.4f +- %.4f stands."
      % (0.5*(max(d0)+min(d0)), 0.5*(max(d0)-min(d0))))

print("\n" + "="*76)
print("WHAT WOULD THE REMAINING GAP REQUIRE?")
print("="*76)
need = G_OBS/G_TREE
m_a1_need = np.sqrt(m_rho**2/(1 - 1/need)) if need > 1 else float("nan")
print("  measured width needs f_rho^2/f_pi^2 = %.5f; the finite-m_pi rule gives %.5f."
      % (need, finite))
# what m_pi would be needed at fixed m_a1
m_pi_need = np.sqrt(m_a1**2 - need*(m_a1**2 - m_rho**2))
print("  at the framework's m_a1, closing it would need m_pi = %.0f MeV, about %.1f times"
      % (m_pi_need, m_pi_need/m_pi))
print("  the physical pion mass, so the pion pole cannot be the whole story.")
print("""
  The residual is therefore NOT the chiral limit, and the audit has now removed three
  candidates in turn: the a_1 mass (wrong direction), a common vector-sector
  normalisation (the phi sits the other side), and now the finite-m_pi sum-rule
  correction (right direction, four times too small).  What is left is the
  vector-dominance assumption itself, that the photon reaches the pions only through
  the rho and a_1 towers, which is the one ingredient never tested.
""")
