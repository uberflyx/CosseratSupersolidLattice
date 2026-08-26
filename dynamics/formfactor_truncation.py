"""Bound the truncation systematic from the framework's own tower couplings.

c_n is a product of the state's coupling to the PHOTON and its coupling to PIONS.  The
framework derives the first through the Cornell solve and does not derive the second, so
c_2 cannot be predicted.  It CAN be bounded.  Taking the pion coupling to be the same for
every tower member, which is the most generous assumption because the physical
recurrences decay dominantly to four pions and so couple to two pions LESS as n rises,
gives c_2/c_1 = f_2/f_1 with f_n proportional to sqrt(Gamma_ee,n m_n).  That is an upper
bound on |c_2|, and the truncation systematic follows.

NORMALISATION.  This script used to bisect c_0 until the channel's total spectral AREA
returned Gamma_ee, which is the prescription R_pipi_factory rejects.  Gamma_ee is a
Breit-Wigner quantity and fixes the PEAK of the whole form factor, so c_0 is solved from

    | c_0 A_0(m^2) + c_1 A_1(m^2) + c_2 A_2(m^2) | = sqrt(36 Gee/(alpha^2 beta^3 Gamma)),

the same condition c0_normalisation imposes, generalised to carry c_2.  The systematic
comes out much smaller under that rule, because c_0 is re-solved whenever c_2 moves and
so the peak is pinned: a third resonance displaces the second rather than rescaling the
channel.  The admissible matching points are also N_ADMISSIBLE, not 1..4.
"""
import sys, os, numpy as np
from scipy import integrate, optimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ew_one_loop_lattice as ew

m_pi, m_rho, G_rho = ew.m_pi, ew.m_rho, ew.Gamma_rho
TPS = 2.0*np.pi*ew.sigma_string
SC, SB = 3739.0**2, 10560.0**2
masses, ratios = ew.tower_ratios()
m1, m2 = masses[1], masses[2]
f1 = np.sqrt(ratios[1]*m1/m_rho)
f2 = np.sqrt(ratios[2]*m2/m_rho)
print("  Cornell tower leptonic ratios : n=1 %.4f   n=2 %.4f" % (ratios[1], ratios[2]))
print("  photon couplings f_n/f_0      : n=1 %.4f   n=2 %.4f" % (f1, f2))
print("  upper bound on |c_2/c_1|      : %.4f" % (f2/f1))

# each amplitude normalised to 1 at s = 0, so that sum c_n = 1 is F_pi(0) = 1
def _unit(m, G):
    A = ew.gs_amplitude(m, G, m_pi)
    n = A(0.0)
    return lambda s: A(s)/n

A0 = _unit(m_rho, G_rho)
A1 = _unit(m1, ew.GAMMA_RHO1)
A2 = _unit(m2, 300.0)

_sR = m_rho**2
_TARGET = np.sqrt(36.0*ew.Gee_rho
                  / (ew.alpha**2 * ew.beta_pi(_sR, m_pi)**3 * G_rho))

def build(c2):
    """Solve the peak condition on the full form factor at this c_2."""
    def R_of(c0):
        c1 = 1.0 - c0 - c2
        def R(s):
            F = c0*A0(s) + c1*A1(s) + c2*A2(s)
            b = np.sqrt(np.maximum(1.0 - 4.0*m_pi*m_pi/np.asarray(s, float), 0.0))
            return 0.25*b**3*abs(F)**2
        return R
    c0 = optimize.brentq(
        lambda c: abs(c*A0(_sR) + (1.0 - c - c2)*A1(_sR) + c2*A2(_sR)) - _TARGET,
        0.3, 2.5, xtol=1e-14)
    return R_of(c0), c0, 1.0 - c0 - c2

def observables(Rp, c1=None):
    """Integrate the anomaly and the running for one member of the c_2 family.

    c_1 is passed in because the recurrence's two-pion branching is not a free
    input alongside c_2: ew.gamma_pipi_recurrence reads it off c_1 itself, so
    holding B fixed while c_2 moves would count one freedom twice and would
    make the sweep disagree with the chapter's own table.
    """
    B = ew.B_PIPI_RHO1 if c1 is None else ew.gamma_pipi_recurrence(c1=c1)/ew.GAMMA_RHO1
    out = {}
    for N in ew.N_ADMISSIBLE:
        sm = min(m_rho**2 + (N+0.5)*TPS, SC)
        a  = ew.amu_disp(Rp, 4*m_pi**2, sm)
        d  = ew.disp(Rp, 4*m_pi**2, sm, ew.MZ2)
        a += ew.amu_narrow(ew.m_omega_pdg, ew.Gee_omega) + ew.amu_narrow(ew.m_phi_pdg, ew.Gee_phi)
        d += ew.narrow(ew.m_omega_pdg, ew.Gee_omega) + ew.narrow(ew.m_phi_pdg, ew.Gee_phi)
        _, Rt = ew.tower_ratios()
        # gee_tower, not G*r: above the first recurrence the couplings come from
        # duality (Gamma_ee(V_n) m_n constant), not from the Cornell ratio, and
        # using the Cornell law up the whole tower reopens the matching band from
        # 0.000042 to 0.000178 and misstates every number below.
        for k, (m0V, G) in enumerate(((m_rho, ew.Gee_rho), (ew.m_omega_pdg, ew.Gee_omega), (ew.m_phi_pdg, ew.Gee_phi))):
            for n, r in enumerate(Rt[1:], 1):
                mn = np.sqrt(m0V**2 + n*TPS)
                if mn**2 >= sm: continue
                w = (1.0-B) if (n==1 and abs(m0V-m_rho)<1.0) else 1.0
                gn = ew.gee_tower(k, n, mn, G, r)
                a += w*ew.amu_narrow(mn, gn); d += w*ew.narrow(mn, gn)
        if sm < SC:
            a += ew.amu_disp(lambda s: 2.0*(1+ew.delta_qcd(s,3)), sm, SC)
            d += ew.disp(lambda s: 2.0*(1+ew.delta_qcd(s,3)), sm, SC, ew.MZ2)
        a += ew.amu_disp(lambda s: (10/3)*(1+ew.delta_qcd(s,4)), SC, SB)
        d += ew.disp(lambda s: (10/3)*(1+ew.delta_qcd(s,4)), SC, SB, ew.MZ2)
        a += ew.amu_disp(lambda s: (11/3)*(1+ew.delta_qcd(s,5)), SB, 4e12)
        d += ew.disp(lambda s: (11/3)*(1+ew.delta_qcd(s,5)), SB, 4e12, ew.MZ2)
        for nm, mV, Gd, Gm, Bee, tau, isd in ew.quarkonia:
            a += ew.amu_narrow(mV, Gd if isd else Gm, 1-(3 if tau else 2)*Bee)
            d += ew.narrow(mV, Gd if isd else Gm, 1-(3 if tau else 2)*Bee)
        out[N] = (a*1e10, d)
    return out

print("\n  c_2       a_mu (N band)          Delta-alpha_had")
allA, allD = [], []
R0, c0_0, c1_0 = build(0.0)
c2flat = abs(c1_0)*f2/f1                       # flat pion coupling up the tower: the ceiling
# c_1 measures the pion coupling at the first recurrence rather than assuming it:
# g_1/g_0 = (c_1/c_0)/(f_1/f_0).  Extrapolating one step and padding by two gives the
# adopted bound; the flat-tower figure is kept as the ceiling if that step fails.
g1g0   = abs(c1_0/c0_0)/f1
c2max  = 2.0*abs(c1_0)*(f2/f1)*g1g0
print("  (c_1 = %.4f at c_2 = 0; g_1/g_0 = %.3f, so |c_2| <= %.4f, ceiling %.4f)"
      % (c1_0, g1g0, c2max, c2flat))
for c2 in (0.0, -0.5*c2max, -c2max, +0.5*c2max, +c2max):
    Rp, c0, c1 = build(c2)
    o = observables(Rp, c1=c1)
    av = [o[n][0] for n in ew.N_ADMISSIBLE]; dv = [o[n][1] for n in ew.N_ADMISSIBLE]
    allA += av; allD += dv
    print("  %+.3f   %.0f to %.0f            %.5f to %.5f"
          % (c2, min(av), max(av), min(dv), max(dv)))

aL, aH = min(allA), max(allA)
dL, dH = min(allD), max(allD)
ac, ah = 0.5*(aL+aH), 0.5*(aH-aL)
dc, dh = 0.5*(dL+dH), 0.5*(dH-dL)
print("\n" + "="*74)
print("  a_mu            = %.0f +- %.0f   -> %.2f sigma vs WP25 713(6)"
      % (ac, ah, abs(ac-713)/np.hypot(ah,6)))
# the top-quark term is carried so the comparison is like for like: without it the
# coupling comes out a hundredth of a unit low and the reference cannot be reproduced
lep, top = 0.031498, -0.7e-4
ai = (1/ew.alpha)*(1-lep-dc-top); au = (1/ew.alpha)*dh
print("  Delta-alpha_had = %.5f +- %.5f" % (dc, dh))
print("  alpha^-1(M_Z)   = %.3f +- %.3f  -> %.2f sigma vs KNT 128.946(15)"
      % (ai, au, abs(ai-128.946)/np.hypot(au,0.015)))
