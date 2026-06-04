#!/usr/bin/env python3
"""
Critical-collapse verification for the Cosserat lattice black-hole chapter.

The black-hole chapter argues that the horizon is a second-order critical
point (Sec. scrambling_time) and that this predicts Choptuik echoing at the
threshold of formation (Sec. critical_collapse).  This script does not derive
the lattice echoing period; that requires the full nonlinear collapse
simulation.  It establishes the three things the chapter actually leans on:

  (A) The Ecker-Ecker-Grumiller leading-order large-D critical solution
      (PRL 136, 191401, 2026) genuinely solves the reduced Einstein-massless-
      Klein-Gordon field equations.  Verified symbolically: the leading-order
      fields satisfy d_x Pi = -x Pi^3 and d_x ln f = x Pi^2 with zero residual,
      and the self-similar-horizon condition f(x=1) = 1 holds.

  (B) Their next-to-leading-order admissibility condition,
          Delta = |beta''| / (3 |beta'|)   at the zero of beta,
      reproduces the published amplitude A = 15.9476 that normalises the
      echoing period to one, for the worked two-mode example
          beta(tau) = cos(2 pi tau) + sin(6 pi tau)/A.
      The subtlety is that the condition must be evaluated at the true zero of
      beta, which the second Fourier mode shifts away from tau = 1/4; the zero
      is therefore located by root-finding rather than assumed.

  (C) The Koike-Hara-Adachi relation gamma = 1/Re(kappa) is checked as a
      consistency statement, and the self-similar growth rate Re(kappa) is
      contrasted with the framework's boundary Lyapunov exponent
      lambda_L = c/(2 r_s).  These are different objects (one dimensionless in
      self-similar time, one a real-time rate in inverse seconds), so the
      chapter must not claim that fixing lambda_L fixes gamma.

Symbolic where possible (sympy) for exact verification, with a numba-accelerated
numerical cross-check on the period condition.  All quantities in SI unless
noted.  The lattice spacing is l = r_e (the classical electron radius), as
established in the monograph.

Mitchell A. Cox / Cosserat supersolid monograph
"""

import numpy as np
import sympy as sp

print("=" * 72)
print("CRITICAL COLLAPSE PROBE  --  large-D anatomy + period + exponent")
print("=" * 72)

# ----------------------------------------------------------------------
# (A) Verify the leading-order solution solves the reduced EOM
# ----------------------------------------------------------------------
# EEG coordinates: tau = log-time, x = self-similar radius, ' = d/dtau.
# Fields at LO: Pi_LO(tau,x), Omega_LO = Pi_LO^2, Psi_LO = D1 Pi_LO / f_LO,
# with D1 = d/dtau + x d/dx + 1, and f_LO from (3b).
#
# Their closed forms (Eqs. 10, 13):
#   Pi_LO = beta / sqrt(1 + beta^2 x^2)
#   f_LO  = sqrt( (1 + beta^2 x^2) / (1 + beta^2) )
#
# We check the two LO PDEs that must hold:
#   (i)   d_x Pi_LO = -x Pi_LO^3                                  [their (9)]
#   (ii)  d_x ln f_LO = x Pi_LO^2                                 [their (11)]
# and the SSH boundary condition f_LO(tau, x=1) = 1.
print("\n(A) Leading-order solution check (symbolic, exact)")

tau, x = sp.symbols('tau x', real=True)
beta = sp.Function('beta')(tau)

Pi_LO = beta / sp.sqrt(1 + beta**2 * x**2)
f_LO = sp.sqrt((1 + beta**2 * x**2) / (1 + beta**2))

# (i) radial matter equation  d_x Pi = -x Pi^3
lhs_i = sp.diff(Pi_LO, x)
rhs_i = -x * Pi_LO**3
res_i = sp.simplify(lhs_i - rhs_i)
print(f"   (i)  d_x Pi + x Pi^3 = 0 ?   residual = {res_i}")

# (ii) metric equation  d_x ln f = x Pi^2
lhs_ii = sp.diff(sp.log(f_LO), x)
rhs_ii = x * Pi_LO**2
res_ii = sp.simplify(lhs_ii - rhs_ii)
print(f"   (ii) d_x ln f - x Pi^2 = 0 ? residual = {res_ii}")

# SSH boundary condition
ssh = sp.simplify(f_LO.subs(x, 1))
print(f"   SSH: f_LO(x=1) = {ssh}   (must be 1)")

# Omega_LO = Pi_LO^2 consistency (their Eq. 7): the 1/epsilon terms in (3a)
# force Omega_LO = Pi_LO^2.  We just record the identity used downstream.
Omega_LO = Pi_LO**2
print(f"   Omega_LO = Pi_LO^2 imposed (their Eq. 7).")

allA = (res_i == 0) and (res_ii == 0) and (ssh == 1)
print(f"   => (A) leading-order solution verified: {allA}")

# ----------------------------------------------------------------------
# (B) NLO echoing-period admissibility condition
# ----------------------------------------------------------------------
# EEG show (their Eq. 19) that near a zero of beta and near the SSH,
#   f_LO + eps f_NLO = 1 - eps (1-x) beta (beta'' + 3 beta') + ...
# Admissibility f <= 1 then requires (beta'' + 3 beta') to vanish where beta=0,
# which for a rescaled argument beta(tau/Delta) fixes the period:
#       Delta = |beta''| / (3 |beta'|)   evaluated at beta = 0.            (20)
#
# Worked 2-mode example (their Eq. 21):
#   beta(tau) = cos(2 pi tau) + sin(6 pi tau)/A,  with Delta = 1 forced.
# We solve for A that makes Eq. (20) give Delta = 1.
print("\n(B) NLO echoing-period condition (their Eqs. 19-21)")

from scipy.optimize import brentq

twopi, sixpi = 2.0*np.pi, 6.0*np.pi
def beta_f(t, A):  return np.cos(twopi*t) + np.sin(sixpi*t)/A
def b1_f(t, A):    return -twopi*np.sin(twopi*t) + (sixpi/A)*np.cos(sixpi*t)
def b2_f(t, A):    return -(twopi**2)*np.cos(twopi*t) - ((sixpi**2)/A)*np.sin(sixpi*t)

def zero_of_beta(A):
    # beta(0)=+1, beta(1/2)=-1: a single sign change in (0, 1/2).
    # The second mode shifts the zero away from t=1/4, so it must be found,
    # not assumed -- that was the subtlety.
    return brentq(lambda t: beta_f(t, A), 1e-9, 0.5 - 1e-9)

def Delta_at_zero(A):
    t0 = zero_of_beta(A)
    return abs(b2_f(t0, A)) / (3.0*abs(b1_f(t0, A))), t0

# Solve Delta_at_zero(A) = 1 (their condition with the period normalised to 1).
def cond(A):
    d, _ = Delta_at_zero(A)
    return d - 1.0

grid = np.linspace(2.0, 40.0, 4000)
vals = [cond(A) for A in grid]
roots = []
for i in range(len(grid)-1):
    if vals[i] == 0 or vals[i]*vals[i+1] < 0:
        roots.append(brentq(cond, grid[i], grid[i+1]))
bestA = roots[0] if roots else float('nan')
t0 = zero_of_beta(bestA)
print(f"   true zero of beta lies at t0 = {t0:.5f}  (NOT 1/4; second mode shifts it)")
print(f"   beta(t0)   = {beta_f(t0, bestA):.2e}")
print(f"   beta'(t0)  = {b1_f(t0, bestA):.5f}")
print(f"   beta''(t0) = {b2_f(t0, bestA):.5f}")
print(f"   |beta''|/(3|beta'|) at t0 = {abs(b2_f(t0,bestA))/(3*abs(b1_f(t0,bestA))):.6f}")
print(f"   => A that fixes Delta = 1:  A = {bestA:.4f}   (EEG quote A = 15.9476)")

# ----------------------------------------------------------------------
# (C) Koike-Hara-Adachi exponent relation and the lattice Lyapunov number
# ----------------------------------------------------------------------
# gamma = 1 / Re(kappa), kappa = growth rate of the single unstable self-similar
# mode.  Reference DSS scalar value (4D): gamma ~ 0.374, Delta ~ 3.4453.
# => Re(kappa) = 1/gamma.
print("\n(C) Exponent relation gamma = 1/Re(kappa)  (Koike-Hara-Adachi)")
gamma_4D_scalar = 0.374    # Choptuik / Gundlach value, massless scalar, D=4
Re_kappa = 1.0/gamma_4D_scalar
print(f"   D=4 massless scalar: gamma ~ {gamma_4D_scalar}")
print(f"   => Re(kappa) = 1/gamma ~ {Re_kappa:.4f}  (self-similar growth rate)")
print(f"   Delta(4D, scalar) ~ 3.4453 (echoing period, for reference)")

# Framework's own boundary Lyapunov exponent (already in monograph / code):
#   lambda_L = c/(2 r_s) = c^3/(4 G M)  -- a *real-time* rate, units 1/s.
# NOTE: this is NOT the same object as Re(kappa) (dimensionless, self-similar
# time).  We print both to make the distinction explicit and to show the
# overreach we are about to fix in the tex.
c = 2.99792458e8
G = 6.67430e-11
Msun = 1.98892e30
for Msol in (1.0, 30.0):
    M = Msol*Msun
    rs = 2*G*M/c**2
    lamL = c/(2*rs)
    print(f"   M={Msol:>4.0f} Msun:  r_s={rs:.3e} m,  lambda_L=c/(2 r_s)={lamL:.3e} 1/s "
          f"(real-time boundary rate)")
print("   -> lambda_L [1/s] and Re(kappa) [dimensionless] are different objects;")
print("      the tex must not assert lambda_L 'fixes gamma' directly.")

print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)
print(f"  (A) LO solution solves reduced EOM + SSH BC : {allA}")
print(f"  (B) period condition reproduces A ~ {bestA:.2f} (EEG: 15.95), Delta=1 : "
      f"{abs(bestA-15.9476) < 0.1}")
print(f"  (C) gamma=1/Re(kappa): Re(kappa)~{Re_kappa:.3f} for gamma~{gamma_4D_scalar};")
print(f"      lambda_L (framework) is a distinct real-time rate.")
