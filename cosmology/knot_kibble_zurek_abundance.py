"""
Kibble-Zurek relic abundance of the condensate knot sector.

The superfluid chapter leaves the knot sector's present abundance as an honest
unknown: Kibble-Zurek makes the objects at the crystallisation transition, and
reconnection removes some fraction, and neither number was derived. This script
supplies the first one from the framework's own scales, and shows that it
brackets the dark-matter density rather than sitting far above or below it.

The Zurek freeze-out length. A quench through a continuous transition freezes
the order parameter when the correlation length can no longer keep up with the
cooling. The frozen length is

    xi_hat = xi_0 (tau_Q / tau_0)^sigma,   sigma = nu / (1 + z nu),

with xi_0 = l (the healing length ~ lattice spacing), tau_0 = l/c the
microscopic time, and tau_Q the quench time. For a cosmological quench
tau_Q ~ t_c, the age at crystallisation. The framework puts crystallisation near
t_c ~ 1e-5 s, so tau_Q/tau_0 ~ 1e18, and the mean-field exponents give
sigma = 1/4 (z=2) or 1/3 (z=1). The frozen length is then ~1e4 to 1e6 lattice
spacings: microscopic on cosmic scales, but large compared with the core.

From the frozen length to a relic. One closed loop forms per correlation volume,
so n ~ 1/xi_hat^3. Each loop is a length of vortex line, so its mass is the
line tension (about one node energy m0 c^2 = 70 MeV per lattice spacing) times
its length ~ 2 pi xi_hat. The comoving yield Y = n/s (s the entropy density at
T_c) is frozen in, and the relic density follows from the standard

    Omega h^2 = m_knot * Y * (s0 / (rho_c/h^2)),   s0/(rho_c/h^2) = 2.74e8 /GeV.

The result depends on the frozen length through Y ~ 1/xi_hat^3 and m ~ xi_hat,
so Omega ~ 1/xi_hat^2: a longer freeze-out length means fewer, heavier loops and
a lower density. The correlation length that reproduces Omega_DM h^2 = 0.12 is
computed and compared with the Zurek range above.

Honest uncertainties, all flagged in the output: the quench exponent (a factor
~100 in xi_hat, hence ~1e4 in Omega), the survival fraction against reconnection
(taken as 1 here; any decay lowers Omega), the geometric O(1) in the defect
density, and g_* at the transition. The claim is not that Omega is pinned, but
that Omega_DM sits inside the plausible band, so the sector is a viable
dark-matter-scale relic, not a negligible trace.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Fundamentals from the framework.
# ---------------------------------------------------------------------------
alpha = 1.0 / 137.035999177
m_e_MeV = 0.51099895069
m0_MeV = m_e_MeV / alpha                 # node mass 70.03 MeV
Tc_MeV = m0_MeV / np.sqrt(6.0)           # crystallisation temperature T_geom
l_over_c_s = 2.8179403205e-15 / 2.99792458e8   # microscopic time l/c [s]
t_c_s = 1.0e-5                           # crystallisation age [s] (first 1e-5 s)

# Line-tension coefficient: T = c_T * (m0 c^2 / l), c_T = (4pi/5) ln(R/xi).
def c_T(lnRxi=1.0):
    return (4.0 * np.pi / 5.0) * lnRxi

# Loop mass in MeV for a ring of size xi_hat = b * l (length 2 pi xi_hat).
def m_knot_MeV(b, lnRxi=1.0):
    return c_T(lnRxi) * (2.0 * np.pi * b) * m0_MeV     # T * length, in units m0

# Standard relic conversion: Omega h^2 = m[GeV] * Y * K, K = s0/(rho_c/h^2).
K_relic = 2.742e8            # /GeV  (s0 = 2891/cm^3, rho_c/h^2 = 1.054e-5 GeV/cm^3)

def yield_and_omega(b, g_star=10.0, f_geo=1.0, lnRxi=1.0):
    """Y = n/s and Omega h^2 for frozen length xi_hat = b*l."""
    # xi_hat * T_c is dimensionless in natural units: l = 1/m0, T_c = m0/sqrt6,
    # so xi_hat * T_c = b / sqrt6. Hence (xi_hat T_c)^3 = b^3 / 6^1.5.
    xiT3 = (b / np.sqrt(6.0)) ** 3
    s_over_Tc3 = (2.0 * np.pi ** 2 / 45.0) * g_star     # s = this * T_c^3
    Y = f_geo / (xiT3 * s_over_Tc3)                     # n/s, dimensionless
    m_GeV = m_knot_MeV(b, lnRxi) / 1e3
    Omega_h2 = m_GeV * Y * K_relic
    return Y, m_GeV, Omega_h2

def zurek_length(sigma):
    return (t_c_s / l_over_c_s) ** sigma                # xi_hat / l = b

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def line():
    print("-" * 76)

print(__doc__.strip().splitlines()[0])
line()
print("Scales from the framework:")
print(f"  node mass m0 = {m0_MeV:.2f} MeV   T_c = m0/sqrt6 = {Tc_MeV:.2f} MeV")
print(f"  microscopic time l/c = {l_over_c_s:.2e} s   quench time t_c ~ {t_c_s:.0e} s")
print(f"  tau_Q/tau_0 = {t_c_s/l_over_c_s:.2e}")
line()
print("Zurek freeze-out length b = xi_hat/l for the mean-field exponents:")
for name, sig in [("z=2 (sigma=1/4)", 0.25), ("z=1 (sigma=1/3)", 1.0/3.0)]:
    b = zurek_length(sig)
    print(f"  {name:<18s}  b = {b:.2e}   xi_hat = {b*2.818e-15*1e9:.2g} nm")
line()
print("Relic density for frozen length b (g_* = 10, f_geo = 1, ln = 1):")
print(f"  {'b = xi_hat/l':>14s} {'m_knot [GeV]':>13s} {'Y':>11s} {'Omega h^2':>12s}")
for b in [3.0e4, 1.0e5, 3.0e5, 1.0e6]:
    Y, mG, Om = yield_and_omega(b)
    print(f"  {b:>14.1e} {mG:>13.1f} {Y:>11.2e} {Om:>12.2e}")
line()

# Correlation length that reproduces Omega_DM = 0.12.
Omega_DM = 0.12
# Omega ~ C / b^2  =>  b* = sqrt(C / Omega_DM). Find C from a reference point.
Y_ref, m_ref, Om_ref = yield_and_omega(1.0e5)
C = Om_ref * (1.0e5) ** 2
b_star = np.sqrt(C / Omega_DM)
Y_s, m_s, Om_s = yield_and_omega(b_star)
print("Matching the dark-matter density:")
print(f"  Omega_DM h^2 = {Omega_DM}")
print(f"  required frozen length  b* = {b_star:.2e}  (xi_hat = {b_star*2.818e-15*1e9:.2g} nm)")
print(f"  loop mass there         m_knot = {m_s:.1f} GeV")
print(f"  check: Omega(b*) h^2 = {Om_s:.3f}")
b_lo, b_hi = zurek_length(0.25), zurek_length(1.0/3.0)
inside = min(b_lo, b_hi) <= b_star <= max(b_lo, b_hi)
print(f"  Zurek band is b = {min(b_lo,b_hi):.1e} to {max(b_lo,b_hi):.1e};"
      f"  b* inside band: {inside}")
line()
# sigma/m of the frozen loop itself (size xi_hat = b* l).
R_star_cm = b_star * 2.8179403205e-13
sigma_star = np.pi * R_star_cm ** 2
som_star = sigma_star / (m_s * 1.782662e-24)     # cm^2/g, m_s in GeV
print("A caveat the abundance forces into the open:")
print(f"  at b*, the frozen loop is heavy ({m_s:.0e} GeV) and its self-interaction")
print(f"  is sigma/m ~ {som_star:.0e} cm^2/g, far ABOVE the SIDM window (~1 cm^2/g).")
print("  So the abundance and the SIDM cross-section come from different points")
print("  in the size spectrum. The frozen loops carry the right mass density but")
print("  are collisional; the SIDM window belongs to the small (~GeV) knots at")
print("  the bottom of the reconnection cascade. Grinding the tangle down to")
print("  those also radiates line length into second sound, so it lowers Omega.")
line()
print("Reading it:")
print("  Two things are robust. The Kibble-Zurek abundance brackets the dark-")
print("  matter density: the frozen length that matches Omega_DM sits inside the")
print("  Zurek band, so the sector is a dark-matter-scale relic, not a negligible")
print("  trace and not a gross overclosure. And the sector is strongly self-")
print("  interacting at every scale. What is not fixed is the endpoint of the")
print("  reconnection cascade, which sets both the surviving fraction and the")
print("  effective sigma/m -- collisional if loops stay large, SIDM-window if")
print("  they grind to ~GeV knots.")
print()
print("  With the framework's bulk dark matter already the collisionless")
print("  vacancy, the natural reading is a MIXED dark sector: vacancy CDM for")
print("  the smooth bulk, plus a knot component that is self-interacting. If the")
print("  cascade lands near the GeV knots, that component is SIDM at the tens-of-")
print("  percent level the small-scale problems call for; if it stalls at large")
print("  loops, it is a heavier, more collisional minority. Either way the bulk")
print("  stays with the vacancy and the knots supply the self-interaction.")
line()
print("Still uncertain (unchanged in kind, narrowed in range):")
print("  - quench exponent sigma: a factor ~100 in b, ~1e4 in Omega;")
print("  - survival fraction against reconnection: taken 1, any decay lowers Omega;")
print("  - geometric O(1) in the defect density and g_* at T_c.")
print("  The abundance is bracketed, not pinned. Pinning it needs the")
print("  reconnection dynamics of the crystallisation front.")
