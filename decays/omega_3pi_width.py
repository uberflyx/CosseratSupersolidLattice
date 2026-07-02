#!/usr/bin/env python3
"""
omega_3pi_width.py
===================

The omega width from the winding vertex, with every input lattice-derived.

The omega shares the crossed-fault geometry of the rho but pairs the two
{111} planes in the isoscalar combination, and for that combination the
direct two-pion door is shut (the Standard Model's G-parity).  Its
dominant channel is three pions, and the vertex that opens it is the
same Z_3 winding trace that fixes pi0 -> gamma gamma: the anomalous
(epsilon-tensor) coupling of one isoscalar vector, one isovector vector,
and one pion,

    g_{omega rho pi} = N_c g^2 / (8 pi^2 f_pi),

with every factor already derived in the chapter: N_c = 3 stacking
positions, g = 5.938 the crossed-fault healing coupling (KSRF, derived
from the strain overlap), f_pi = 92.37 MeV the stacking-edge decay
constant.  Substituting the quark-level slip coupling g/2 (universality)
for the two charges e Q_u, e Q_d of the pi0 trace and summing the
isoscalar flavour trace reproduces the coefficient; see the chapter text
for the three-line substitution.

The decay then runs omega -> (rho pi)_virtual -> 3 pi, one rho per pion
pair, all three charge channels adding coherently.  The rho propagators
carry the chapter's own derived width Gamma_rho = 147.3 MeV, made
energy-dependent by the derived P-wave law Gamma(s) ~ p(s)^3.

This script:
  1. verifies the two polarisation-sum contractions symbolically with
     sympy (rank-4 epsilon contractions are where sign errors live);
  2. integrates the 3-pion Dalitz plot numerically;
  3. computes omega -> pi0 gamma through the same vertex with one rho
     leg converted to a photon by the rho-gamma mixing e/g (the same
     mixing the chapter already used for omega -> e+ e-);
  4. prints the partial and total widths against PDG.

All masses PDG 2024; all couplings lattice-derived.
"""

import numpy as np

# ---------------------------------------------------------------- inputs
# Lattice-derived (dynamics chapter):
ALPHA   = 1.0 / 137.035999177
G_RHO   = 5.938                    # crossed-fault healing coupling (KSRF)
F_PI    = 92.37                    # MeV, stacking-edge decay constant
GAMMA_RHO_0 = 147.3                # MeV, derived rho width at the pole
N_C     = 3                        # stacking positions A, B, C
E_CHARGE = np.sqrt(4.0 * np.pi * ALPHA)

# PDG 2024 masses [MeV]:
M_OMEGA = 782.66
M_RHO   = 775.26
M_PIC   = 139.5704                 # pi+/-
M_PI0   = 134.9768

# PDG 2024 observed widths [MeV]:
GAMMA_OMEGA_TOTAL_OBS = 8.68       # +/- 0.13
BR_3PI  = 0.892                    # +/- 0.007
BR_PI0G = 0.0834                   # +/- 0.0026
GAMMA_3PI_OBS  = GAMMA_OMEGA_TOTAL_OBS * BR_3PI
GAMMA_PI0G_OBS = GAMMA_OMEGA_TOTAL_OBS * BR_PI0G

HBAR_MEV_S = 6.582119569e-22

# The winding-trace vertex [MeV^-1]:
G_ORP = N_C * G_RHO**2 / (8.0 * np.pi**2 * F_PI)


# ---------------------------------------------------- symbolic checks
def check_polarisation_sums():
    """Verify the two epsilon-tensor contractions with sympy.

    (a) omega -> 3pi:  T^mu = eps^{mu nu rho sigma} p+_nu p-_rho p0_sigma.
        In the omega rest frame T^0 = 0 (the three momenta are coplanar
        in the sense p0 = -p+ - p-), and the massive-vector polarisation
        sum gives  sum_pol |eps.T|^2 = |T_spatial|^2
                 = m_omega^2 |p+ x p-|^2.

    (b) omega -> pi0 gamma:  M = C eps_{mu nu rho sigma}
        eps_omega^mu eps_gamma^{nu*} k_omega^rho k_gamma^sigma.
        Summing photon polarisations (-g^{nu nu'}) and averaging the
        omega gives  |Mbar|^2 = (2/3) C^2 (k_omega . k_gamma)^2
                              = (2/3) C^2 m_omega^2 E_gamma^2,
        hence Gamma = C^2 E_gamma^3 / (12 pi).
    """
    import sympy as sp
    from sympy import LeviCivita

    # metric signature (+,-,-,-)
    g = sp.diag(1, -1, -1, -1)

    # (a) generic non-collinear momenta in the rest frame
    m, e1, e2 = sp.symbols('m e1 e2', positive=True)
    p1x, p2x, p2y = sp.symbols('p1x p2x p2y', real=True)
    Pp = sp.Matrix([e1, p1x, 0, 0])
    Pm = sp.Matrix([e2, p2x, p2y, 0])
    P0 = sp.Matrix([m - e1 - e2, -p1x - p2x, -p2y, 0])

    T = [sum(LeviCivita(mu, nu, rho, sig) * Pp[nu] * Pm[rho] * P0[sig]
             for nu in range(4) for rho in range(4) for sig in range(4))
         for mu in range(4)]
    T = [sp.simplify(t) for t in T]
    assert sp.simplify(T[0]) == 0, "T^0 should vanish in the rest frame"

    # sum over massive polarisations: -g + kk/m^2 acting on T with k=(m,0):
    # the kk/m^2 piece needs k.T = m T^0 = 0, so the sum is just |T_spatial|^2
    Tsq = sp.simplify(sum((-g[mu, mu]) * T[mu]**2 for mu in range(1, 4)))
    cross_sq = sp.expand((p1x * 0 - 0 * p2x)**2 + (0 * p2x - p1x * 0)**2
                         + (p1x * p2y - 0 * p2x)**2)  # |p+ x p-|^2 for these vectors
    assert sp.simplify(Tsq - m**2 * cross_sq) == 0, "|T|^2 != m^2 |p+ x p-|^2"

    # (b) V -> P gamma contraction
    C, k = sp.symbols('C k', positive=True)   # k = photon energy
    mV = sp.symbols('mV', positive=True)
    kV = sp.Matrix([mV, 0, 0, 0])
    kG = sp.Matrix([k, 0, 0, k])
    # sum over photon pol: -g^{nu nu'}; sum over V pol: -g^{mu mu'} + kV kV/mV^2
    total = 0
    for mu in range(4):
        for mup in range(4):
            for nu in range(4):
                for nup in range(4):
                    e1_ = sum(LeviCivita(mu, nu, r, s) * kV[r] * kG[s]
                              for r in range(4) for s in range(4))
                    e2_ = sum(LeviCivita(mup, nup, r, s) * kV[r] * kG[s]
                              for r in range(4) for s in range(4))
                    polV = -g[mu, mup] + kV[mu] * kV[mup] / mV**2
                    polG = -g[nu, nup]
                    total += e1_ * e2_ * polV * polG
    total = sp.simplify(total)
    expected = 2 * (mV * k)**2          # 2 (kV . kG)^2 with kV.kG = mV k
    assert sp.simplify(total - expected) == 0, "V->Pgamma pol sum failed"

    # (c) the channel factor 2: the rho -> pi pi current is (p_a - p_b),
    # and contracting the vertex's (p_a + p_b) with it inside the epsilon
    # gives eps(..., p_a+p_b, p_a-p_b, ...) = -2 eps(..., p_a, p_b, ...),
    # one factor of 2 in amplitude for every channel.
    a0, a1, a2, a3, b0, b1, b2, b3 = sp.symbols('a0 a1 a2 a3 b0 b1 b2 b3')
    c0, c1, c2, c3 = sp.symbols('c0 c1 c2 c3')
    A = [a0, a1, a2, a3]; B = [b0, b1, b2, b3]; Cv = [c0, c1, c2, c3]
    lhs = sum(LeviCivita(0, nu, rho, sig) * (A[nu] + B[nu]) * (A[rho] - B[rho]) * Cv[sig]
              for nu in range(4) for rho in range(4) for sig in range(4))
    rhs = -2 * sum(LeviCivita(0, nu, rho, sig) * A[nu] * B[rho] * Cv[sig]
                   for nu in range(4) for rho in range(4) for sig in range(4))
    assert sp.simplify(lhs - rhs) == 0, "channel factor 2 identity failed"
    return True


# ---------------------------------------------------- rho propagator
def p_cm(s, m1, m2):
    """Daughter momentum for a pair of invariant mass sqrt(s) [MeV]."""
    s = np.asarray(s, dtype=float)
    val = (s - (m1 + m2)**2) * (s - (m1 - m2)**2)
    return np.sqrt(np.clip(val, 0.0, None)) / (2.0 * np.sqrt(s))


def bw_rho(s, m1, m2):
    """Relativistic rho Breit-Wigner with the derived P-wave running
    width Gamma(s) = Gamma_0 (p/p0)^3 (m_rho/sqrt(s))  [MeV^-2]."""
    p  = p_cm(s, m1, m2)
    p0 = p_cm(M_RHO**2, m1, m2)
    gam = GAMMA_RHO_0 * (p / p0)**3 * (M_RHO / np.sqrt(s))
    return 1.0 / (M_RHO**2 - s - 1j * np.sqrt(s) * gam)


# ---------------------------------------------------- Dalitz integral
def gamma_omega_3pi(n=1200):
    """Integrate dGamma = |Mbar|^2 / (32 (2 pi)^3 M^3) ds dt over the
    Dalitz region, with
      |Mbar|^2 = (1/3) g_orp^2 g_rho^2 |BW(s)+BW(t)+BW(u)|^2
                 * M^2 |p+ x p-|^2 .
    Labels: 1 = pi+, 2 = pi-, 3 = pi0;  s = m12^2, t = m13^2, u = m23^2,
    with the rho0 in s and the charged rhos in t and u."""
    M = M_OMEGA
    m1 = m2 = M_PIC
    m3 = M_PI0
    s_lo, s_hi = (m1 + m2)**2, (M - m3)**2
    s_vals = np.linspace(s_lo, s_hi, n)
    ds = s_vals[1] - s_vals[0]
    total = 0.0
    for s in s_vals[:-1] + ds / 2.0:
        # t boundaries at fixed s (standard Dalitz limits)
        E2 = (s - m1**2 + m2**2) / (2.0 * np.sqrt(s))
        E3 = (M**2 - s - m3**2) / (2.0 * np.sqrt(s))
        p2 = np.sqrt(max(E2**2 - m2**2, 0.0))
        p3 = np.sqrt(max(E3**2 - m3**2, 0.0))
        t_lo = (E2 + E3)**2 - (p2 + p3)**2
        t_hi = (E2 + E3)**2 - (p2 - p3)**2
        if t_hi <= t_lo:
            continue
        t_vals = np.linspace(t_lo, t_hi, n)
        dt = t_vals[1] - t_vals[0]
        t = t_vals[:-1] + dt / 2.0
        u = M**2 + m1**2 + m2**2 + m3**2 - s - t

        # rest-frame energies from pair invariants
        E1 = (M**2 + m1**2 - u) / (2.0 * M)     # u pairs (2,3), recoils 1
        E2r = (M**2 + m2**2 - t) / (2.0 * M)    # t pairs (1,3), recoils 2
        q1 = np.sqrt(np.clip(E1**2 - m1**2, 0, None))
        q2 = np.sqrt(np.clip(E2r**2 - m2**2, 0, None))
        # p1 . p2 from s = (p1+p2)^2
        dot = E1 * E2r - (s - m1**2 - m2**2) / 2.0
        cross_sq = np.clip(q1**2 * q2**2 - dot**2, 0.0, None)

        F = bw_rho(s, m1, m2) + bw_rho(t, m1, m3) + bw_rho(u, m2, m3)
        # The factor 2 per channel comes from the rho -> pi pi current:
        # eps(..., p_a+p_b, p_a-p_b, ...) = -2 eps(..., p_a, p_b, ...),
        # verified symbolically in check_polarisation_sums().
        Msq = (1.0 / 3.0) * (2.0 * G_ORP * G_RHO)**2 * np.abs(F)**2 * M**2 * cross_sq
        total += np.sum(Msq) * dt
    total *= ds
    return total / (32.0 * (2.0 * np.pi)**3 * M**3)


def gamma_omega_pi0_gamma():
    """omega -> pi0 gamma through the same vertex with one rho leg
    converted to a photon by the rho-gamma mixing e/g (the identity the
    chapter's omega -> e+ e- already used):
        C = (e / g) g_orp,  Gamma = C^2 E_gamma^3 / (12 pi)."""
    C = (E_CHARGE / G_RHO) * G_ORP
    k = (M_OMEGA**2 - M_PI0**2) / (2.0 * M_OMEGA)
    return C**2 * k**3 / (12.0 * np.pi), k


# ---------------------------------------------------------------- main
if __name__ == '__main__':
    print("symbolic polarisation-sum checks:", end=" ")
    check_polarisation_sums()
    print("both pass\n")

    print(f"g_(omega rho pi) = N_c g^2/(8 pi^2 f_pi) "
          f"= {G_ORP*1e3:.2f} GeV^-1   (all factors lattice-derived)\n")

    g3 = gamma_omega_3pi()
    r3 = (g3 - GAMMA_3PI_OBS) / GAMMA_3PI_OBS
    print(f"Gamma(omega -> 3pi)      = {g3:.2f} MeV   "
          f"obs {GAMMA_3PI_OBS:.2f}   residual {100*r3:+.1f}%")

    gg, k = gamma_omega_pi0_gamma()
    rg = (gg - GAMMA_PI0G_OBS) / GAMMA_PI0G_OBS
    print(f"Gamma(omega -> pi0 g)    = {gg:.3f} MeV  "
          f"obs {GAMMA_PI0G_OBS:.3f}  residual {100*rg:+.1f}%   (E_gamma = {k:.1f} MeV)")

    ge_obs = 0.000617   # MeV, already closed in the chapter at +0.9%
    tot = g3 + gg + ge_obs
    rt = (tot - GAMMA_OMEGA_TOTAL_OBS) / GAMMA_OMEGA_TOTAL_OBS
    print(f"Gamma(omega) total       = {tot:.2f} MeV   "
          f"obs {GAMMA_OMEGA_TOTAL_OBS:.2f}   residual {100*rt:+.1f}%")

    tau = HBAR_MEV_S / tot
    print(f"\nomega lifetime           = {tau:.3e} s")
