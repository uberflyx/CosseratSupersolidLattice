#!/usr/bin/env python3
"""The neutrino magnetic moment from the chirally dressed edge: the Cosserat
angular overlap, computed rather than power-counted.

THE QUESTION
The edge dislocation is neutral and, in the achiral lattice, moment-free: the
photon block is anti-plane, the edge's fields are in-plane, and the angular
integral of its in-plane rotation vanishes. The chirality W_ch = 2 C_ch e_(ij)
k_(ij) dresses the edge with rotational content, and the paper brackets the
resulting moment between 2 alpha^2 theta_ch (m_nu/m_e) mu_B ~ 1e-17 and
2 alpha theta_ch^2 mu_B ~ 1e-12, the discriminator being "the full Cosserat
angular overlap of the defect with the photon mode in the presence of the
chirality". This script computes that overlap.

THE MECHANISM, ORDER BY ORDER
Order 0: the edge's in-plane strain e0(r, theta), harmonics cos/sin(theta).
Order 1 (in theta_ch): W_ch is linear in the microrotation curvature, so the
  edge's strain acts as a static source and the rotational sector responds
  coherently: phi1 = -G_rot * S, with G_rot the screened rotational Green's
  function (screening length = the kernel crossover, mass m_D) and S the
  source built from e0. phi1 inherits the source's odd harmonics, so every
  moment integral at this order carries int cos/sin d(theta) = 0. The first
  order vanishes by the angular integral, which is the same statement as the
  paper's parity argument, now at the level of components.
Order 2: phi1 re-enters W_ch and the kappa_c rotational lock, generating an
  anti-plane displacement u_z2 and an in-plane microrotation drag. Products
  of two first-order harmonics contain the s-wave, so the angular integral
  survives. The moment is the s-wave flow dipole of the induced curl,
      D = int d^2r (curl u)_z-projections, converted to mu_B by the
  electron's own normalisation (one circulation quantum b over the ring area
  pi a_star^2 == mu_B, monograph, magnetic-moment-as-geometry).

WHAT IS COMPUTED
1. The angular algebra: which harmonic products reach the s-wave, exactly.
2. The radial overlap: the s-wave projection integral over the Cosserat-
   screened profiles, with the two physical cutoffs the framework supplies,
   the core (ell_c = 1/2) and the rotational screening length (1/Q, Q = 1.87).
3. Assembly: mu_nu / mu_B = O * theta_ch^2 * (flavour projection w_i^2),
   with O the computed overlap, and the cross-check against the admixture
   route |c|^2 = K_w / m_D.
4. Confrontation: XENONnT (lab) and the tip of the red-giant branch
   (plasmon decay), the latter the tighter and previously uncited.

Units: ell = 1, mu_shear = 1, so m0 = mu ell^3 = 1 and energies are in m0.

HONESTY BOX, read before quoting the number.
This is the first computation of the overlap and it collapses the tensor
chain (which strain component sources which phi component, and how the
kappa_c lock drags u_z) to its scalar s-wave skeleton: source-squared times
the screened rotational propagator, with the exact angular algebra but
scalar radial structure, and the rotational stiffness absorbed into the
screening scale Q rather than carried as an explicit 1/gamma. The angular
conclusions are exact: the first order vanishes, the induced moment is
axial, the flavour pattern is w_i^2. The absolute scale carries a model
spread we estimate at a factor of a few either way. What the computation
excludes at any credible spread is an order-unity overlap: the profiles
and the screening cost three orders of magnitude, so the bracket's upper
branch 2 alpha theta_ch^2 (which treats the overlap as ~2 alpha) and the
patch estimate (overlap ~0.1) both overstate the moment. A full tensor
solve is the remaining refinement and would sharpen the coefficient, not
the order.
"""

import numpy as np
from scipy.integrate import quad
from scipy.special import k0, k1

# ---------------------------------------------------------------- constants
ALPHA   = 1 / 137.035999177
THETA   = ALPHA**2 / (2 * np.pi)          # constitutive chirality
N2      = 1 / np.pi                       # Cosserat coupling number
ELL_C   = 0.5                             # Cosserat characteristic length [ell]
Q_ROT   = 1.87                            # rotational screening wavevector [1/ell]
A_STAR  = 1.073                           # electron mu_B ring radius [ell]
POISSON = 0.25                            # FCC Cosserat effective Poisson ratio
KW_MEV  = 6.628e-3 * 1e-6                 # screw-channel energy [MeV] -> for cross-check
M_D_MEV = 130.9                           # rotational gap seen as mass [MeV]
M0_MEV  = 0.51099895 / ALPHA              # node mass [MeV]

# flavour screw projections w_i^2 on (perp, engaged, engaged)
W2 = np.array([1.0, 0.25, 0.25])

# ------------------------------------------------- the edge strain harmonics
# Classical edge, b along x, line along z (Cosserat corrections to the far
# field are O(ell_c^2/r^2) and enter only through the regularised radials
# below). In polar coordinates the independent in-plane strain components are
#   e_rr = -A sin(theta)/r,  e_tt = ... , e_rt = A cos(theta)/r,
# with A = b (1 - 2 nu)/(4 pi (1 - nu)) for the dilatational parts and
# A' = b /(4 pi (1 - nu)) for the shear part. Every component is a single
# first harmonic: the edge is a pure l = 1 source.
A_DIL = (1 - 2 * POISSON) / (4 * np.pi * (1 - POISSON))   # x b, dilatational
A_SH  = 1.0 / (4 * np.pi * (1 - POISSON))                  # x b, shear

# ------------------------------------------------------------- angular algebra
def angular_overlap_table():
    """Products of two l = 1 harmonics reaching the s-wave.

    <cos^2> = <sin^2> = 1/2, <sin cos> = 0 over the circle. The second-order
    source for the z-directed moment is the s-wave part of
    (phi1 x source) products; each term is (first harmonic)^2, so the s-wave
    projection is exactly 1/2 per matched pair and 0 for mixed pairs. The
    moment along the in-plane directions needs an l = 1 part of a product of
    two l = 1 fields, which is absent (products give l = 0 and l = 2 only),
    so the induced moment is purely axial: along the dislocation line.
    """
    return 0.5

# ------------------------------------------------------ radial overlap integral
def radial_overlap(r_max=200.0):
    """The s-wave radial integral of the second-order induced curl.

    Chain: source e0 ~ A/r (regularised inside ell_c by the Cosserat core,
    factor [1 - (r/l_c) K1(r/l_c)] as in the screw solution); the screened
    rotational response convolves with exp screening K0/K1 at 1/Q; the
    second-order curl density then scales as
        rho(r) = e0_reg(r) * G_screen(r) * e0_reg(r) * r-measure.
    The moment integral D = int rho(r) 2 pi r dr converges at both ends:
    the core regularisation removes the r -> 0 divergence and the screening
    removes the log at infinity, so no arbitrary cutoff enters. That closure
    is itself the result: the Cosserat medium supplies both lengths.
    """
    def e0_reg(r):
        # 1/r far field with the Cosserat core closure (screw-form regulator)
        x = r / ELL_C
        return (1.0 / r) * (1.0 - x * k1(x)) if r > 1e-12 else 0.0

    def g_screen(r):
        # screened rotational propagator weight at separation r
        return k0(Q_ROT * r) * Q_ROT**2 / (2 * np.pi)

    def integrand(r):
        return e0_reg(r)**2 * g_screen(r) * r * 2 * np.pi

    val, err = quad(integrand, 1e-8, r_max, limit=400)
    return val, err

# --------------------------------------------------------------- assembly
def main():
    ang = angular_overlap_table()
    rad, rad_err = radial_overlap()
    # geometric strain weights: the two independent edge amplitudes both
    # contribute matched-harmonic squares; their s-wave weight is the sum of
    # squares (mixed terms vanish by <sin cos> = 0)
    geom = A_DIL**2 + A_SH**2
    # the second-order induced circulation, in units of b * theta_ch^2:
    D_ind = ang * geom * rad
    # electron normalisation: circulation b over ring area pi a_star^2 == mu_B
    O_total = D_ind / (np.pi * A_STAR**2)
    mu_perp = O_total * THETA**2            # w^2 = 1 flavour, units of mu_B

    print("=" * 68)
    print("Cosserat angular overlap for the edge moment")
    print("=" * 68)
    print(f"  angular s-wave weight            : {ang:.3f}")
    print(f"  strain harmonic weight (b^2)     : {geom:.5f}")
    print(f"  radial overlap integral          : {rad:.5f}  (quad err {rad_err:.1e})")
    print(f"  induced dipole D [b theta^2]     : {D_ind:.5f}")
    print(f"  overlap O = D/(pi a*^2)          : {O_total:.5f}")
    print()
    print(f"  mu(nu_2)  [w^2 = 1  ]            : {mu_perp:.3e} mu_B")
    print(f"  mu(nu_1) = mu(nu_3) [w^2 = 1/4]  : {0.25*mu_perp:.3e} mu_B")
    mu_vec = mu_perp * W2
    mu_quad = np.sqrt(np.sum(mu_vec**2))
    print(f"  quadrature sum (plasmon decay)   : {mu_quad:.3e} mu_B")
    print()
    # cross-check: admixture route
    c2 = KW_MEV / M_D_MEV
    print("cross-check, admixture route:")
    print(f"  |c|^2 = K_w/m_D                  : {c2:.3e}")
    print(f"  |c|^2 / theta_ch^2               : {c2/THETA**2:.3f}   (paper: 0.70)")
    print(f"  mu via |c|^2 * O_patch(0.092)    : {c2*0.092:.3e} mu_B")
    print()
    print("confrontation:")
    print(f"  XENONnT (2022, lab)              : < 6.4e-12 mu_B")
    print(f"  TRGB plasmon (Capozzi-Raffelt)   : < 1.5e-12 mu_B  [quadrature]")
    for name, mu in [("field route, quadrature", mu_quad),
                     ("admixture route, nu_2", c2 * 0.092)]:
        for bname, b in [("XENONnT", 6.4e-12), ("TRGB", 1.5e-12)]:
            print(f"    {name:28s} vs {bname:8s}: ratio {mu/b:6.2f}")

if __name__ == "__main__":
    main()
