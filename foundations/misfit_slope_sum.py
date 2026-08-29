#!/usr/bin/env python3
"""pi/4 reduces to one number: the misfit slope sum S, and alpha demands 1.001067.

THE EXACT REDUCTION (sympy-verified, all m)
For the Peierls-Nabarro arctangent profile spanning one period d, every
misfit harmonic integrates in closed form:
    Int [1 - cos(2 pi m u/d)] dx = 2 pi m w        (exactly, every m)
so for a general periodic misfit gamma(u) = sum_m V_m [1-cos(2 pi m u/d)]
the misfit force dE/dw = 2 pi sum_m m V_m is width-independent and the
variational equilibrium (alphaG eq:pn_equilibrium) generalises to the
Frenkel condition with the right-hand side scaled by ONE number:
    S = sum_m m V_m / V_1.
The entire shape of the gamma-surface enters the core width through S alone.

WHAT ALPHA SAYS S MUST BE
Validation: S = 1 (pure Frenkel) reproduces the paper's w/d = 0.786160 to
all printed digits. Inverting at w/d = pi/4 exactly:
    S_alpha = 1.001067,
with sensitivity d(w/d)/dS = -0.713, so the +0.107% slope surplus is the
same bookkeeping as the known +0.097% width residual. alpha's 2.4e-11 band
on w/d pins S to 3.4e-11 in principle; the meaningful near-term target is
1e-3, the level of the residual itself. Deriving pi/4 from first
principles now means: show the FCC contact gamma-surface has S = 1.0011.

THE FIRST ATTEMPT, FALSIFIED BY ITS OWN ARTEFACT (recorded)
A normal-springs-only surface with hard nearest-three contact switching at
the bridge gives alternating slowly-decaying harmonics
(V_m/V_1 = 1, -0.236, +0.106, -0.059, ...) whose slope sum climbs with
truncation (0.609, 0.641, 0.658 at m = 4, 6, 8): the signature of a
derivative kink at the switch, which is a prescription artefact, not
lattice physics. The model also omits the tangential-rotational channel,
which carries kappa_c/2 of the slide stiffness mu-bar = mu + kappa_c/2,
i.e. 46.7 per cent of it, and is smooth through the bridge. Conclusion:
the honest first-principles surface needs (a) the framework's actual
cohesion mechanism for tangent-sphere contacts, looked up not invented,
and (b) the finite-amplitude tangential energy with node rotations
relaxed, exactly the two ingredients the alpha paper's open-calculation
flag names.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

N2 = 1/np.pi
P2 = 4*N2*(1-N2)/(1-2*N2); Q2 = 4*N2/(1-2*N2)
MUT = 1/(1-2*N2); MUBAR = (1-N2)/(1-2*N2)
D111 = np.sqrt(2/3); D = 1/np.sqrt(3)

def kernel_int(w):
    f = lambda k: 0.5*MUT*(k*k+P2)/(k*k+Q2)*np.exp(-2*k*w)
    return quad(f, 0, np.inf, limit=400)[0]

def w_of_S(S):
    return brentq(lambda w: kernel_int(w) - S*MUBAR/(2*D111), 0.2, 0.9, xtol=1e-13)

def main():
    w1 = w_of_S(1.0)
    S_alpha = kernel_int(np.pi/4*D)/(MUBAR/(2*D111))
    h = 1e-6
    slope = (w_of_S(S_alpha+h) - w_of_S(S_alpha-h))/(2*h)/D
    print(f"validation, S=1:      w/d = {w1/D:.6f}   (paper: 0.786160)")
    print(f"inversion at pi/4:    S_alpha = {S_alpha:.6f}")
    print(f"sensitivity:          d(w/d)/dS = {slope:+.4f}")
    print(f"alpha band on S:      {2.4e-11/abs(slope):.1e}   (near-term target 1e-3)")
    print(f"kinked v1 surface:    S ~ 0.658 and truncation-unstable; rejected")

if __name__ == "__main__":
    main()
