"""Bound the truncation systematic from the framework's own tower couplings.

c_n is a product of the state's coupling to the PHOTON and its coupling to PIONS.  The
framework derives the first through the Cornell solve and does not derive the second, so
c_2 cannot be predicted.  It CAN be bounded.  Taking the pion coupling to be the same for
every tower member, which is the most generous assumption because the physical
recurrences decay dominantly to four pions and so couple to two pions LESS as n rises,
gives c_2/c_1 = f_2/f_1 with f_n proportional to sqrt(Gamma_ee,n m_n).  That is an upper
bound on |c_2|, and the truncation systematic follows.
"""
import sys, os, numpy as np
from scipy import integrate
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

A0 = ew.gs_amplitude(m_rho, G_rho, m_pi)
A1 = ew.gs_amplitude(m1, ew.GAMMA_RHO1, m_pi)
A2 = ew.gs_amplitude(m2, 300.0, m_pi)

def build(c2):
    def R_of(c0):
        c1 = 1.0 - c0 - c2
        def R(s):
            F = c0*A0(s) + c1*A1(s) + c2*A2(s)
            b = np.sqrt(np.maximum(1.0 - 4.0*m_pi*m_pi/s, 0.0))
            return 0.25*b**3*abs(F)**2
        return R
    def gee(c0):
        ar, _ = integrate.quad(R_of(c0), 4*m_pi**2, (m_rho+12*G_rho)**2, limit=250)
        return ar*ew.alpha**2/(9*np.pi*m_rho)
    lo, hi = 0.5, 3.0
    for _ in range(45):
        mid = 0.5*(lo+hi)
        if gee(mid) < ew.Gee_rho: lo = mid
        else: hi = mid
    c0 = 0.5*(lo+hi)
    return R_of(c0), c0, 1.0 - c0 - c2

def observables(Rp):
    out = {}
    for N in (1,2,3,4):
        sm = min(m_rho**2 + (N+0.5)*TPS, SC)
        a  = ew.amu_disp(Rp, 4*m_pi**2, sm)
        d  = ew.disp(Rp, 4*m_pi**2, sm, ew.MZ2)
        a += ew.amu_narrow(ew.m_omega_pdg, ew.Gee_omega) + ew.amu_narrow(ew.m_phi_pdg, ew.Gee_phi)
        d += ew.narrow(ew.m_omega_pdg, ew.Gee_omega) + ew.narrow(ew.m_phi_pdg, ew.Gee_phi)
        _, Rt = ew.tower_ratios()
        for m0V, G in ((m_rho, ew.Gee_rho), (ew.m_omega_pdg, ew.Gee_omega), (ew.m_phi_pdg, ew.Gee_phi)):
            for n, r in enumerate(Rt[1:], 1):
                mn = np.sqrt(m0V**2 + n*TPS)
                if mn**2 >= sm: continue
                w = (1.0-ew.B_PIPI_RHO1) if (n==1 and abs(m0V-m_rho)<1.0) else 1.0
                a += w*ew.amu_narrow(mn, G*r); d += w*ew.narrow(mn, G*r)
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
c2max = abs(c1_0)*f2/f1
print("  (c_1 = %.4f at c_2 = 0, so the bound is |c_2| <= %.4f)" % (c1_0, c2max))
for c2 in (0.0, -0.5*c2max, -c2max, +0.5*c2max, +c2max):
    Rp, c0, c1 = build(c2)
    o = observables(Rp)
    av = [o[n][0] for n in (1,2,3,4)]; dv = [o[n][1] for n in (1,2,3,4)]
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
lep = 0.031498
ai = (1/ew.alpha)*(1-lep-dc); au = (1/ew.alpha)*dh
print("  Delta-alpha_had = %.4f +- %.4f" % (dc, dh))
print("  alpha^-1(M_Z)   = %.2f +- %.2f  -> %.2f sigma vs 128.943(14)"
      % (ai, au, abs(ai-128.943)/np.hypot(au,0.014)))
