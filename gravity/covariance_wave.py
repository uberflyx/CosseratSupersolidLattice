"""Spin-2 without a spin-2 quantum: the covariance wave, proved four ways.

The medium has no linear helicity-2 branch (theorem, part 1). The helicity-2
object is the second central moment of a gapless transverse field: the
covariance Q_ij = <f_i f_j> - (1/3) delta_ij <f^2>, where f is the transverse
microrotation phi (the algebra is identical for any transverse vector field).
This script proves, symbolically and numerically:

  1. NO-GO. For a plane wave along z, any rank-2 tensor built as a derivative
     of a vector field has helicities 0 and +-1 only; its TT projection
     vanishes identically. No linear branch of any elastic medium radiates
     the LIGO strain pattern.
  2. HELICITY 2. The traceless covariance transforms under rotations about
     the propagation axis with e^{-+2i theta}: exactly helicity +-2, because
     the rotation acts twice on a quadratic form.
  3. MASSLESSNESS. A travelling covariance modulation is carried by pairs of
     collinear gapless quanta: omega = c(k1+k2) = c|k| identically, so the
     mode is gapless exactly when its constituents are, and it moves at c.
  4. MOMENT ORTHOGONALITY. In a Gaussian state, displacing the mean
     (a coherent wave: electromagnetism / gravitomagnetism, the FIRST moment)
     leaves the central covariance untouched, and squeezing (a covariance
     wave: the tensor mode, the SECOND central moment) leaves the mean
     untouched. The two radiations are orthogonal observables of one medium,
     so no "two-photon graviton" arises: a tensor wave has <f> = 0.

Consequence for the black-hole log correction: the medium's fundamental
fields are the two gapless gauge vectors (transverse displacement = photon,
transverse microrotation). The metric is a composite moment of them, so it
contributes no independent loop; adding one would double-count the phonon
determinant. C_local = -26/90 (photon) - 26/90 (microrotation vector)
= -26/45, against +199/45 (GR: fundamental tensor graviton + photon) and
approximately -3/2 (LQG). The discriminator survives with a repaired
justification: not "the graviton is a vector" but "there is no graviton
field to run in the loop".

No observation requires the quantum: Dyson (2013) and Rothman & Boughn
(2006) on the infeasibility of single-graviton detection; Carney, Domcke &
Rodd (2024) on why even a single-graviton-sensitive click would not by
itself establish quantisation.
"""

import numpy as np
import sympy as sp

print("=" * 72)
print("[1] NO-GO: linear branches carry no helicity 2")
print("=" * 72)
k, phix, phiy, phiz = sp.symbols('k phi_x phi_y phi_z')
# plane wave along z: every gradient is (0, 0, ik); field amplitude arbitrary
grad = sp.Matrix([0, 0, sp.I * k])
f = sp.Matrix([phix, phiy, phiz])
S = (grad * f.T + f * grad.T) / 2          # symmetric derivative tensor
print("  S_xx - S_yy =", sp.simplify(S[0, 0] - S[1, 1]),
      ";  S_xy =", sp.simplify(S[0, 1]), "   (helicity-2 slots empty)")
# TT projector for k along z kills every k-carrying component:
P = sp.diag(1, 1, 0)
TT = P * S * P - sp.Rational(1, 2) * P * sp.trace(P * S)
print("  TT projection of S =", sp.simplify(TT), " (zero tensor)")

print()
print("=" * 72)
print("[2] HELICITY 2 of the covariance")
print("=" * 72)
th = sp.symbols('theta', real=True)
Qxx, Qyy, Qxy = sp.symbols('Q_xx Q_yy Q_xy', real=True)
Q = sp.Matrix([[Qxx, Qxy], [Qxy, Qyy]])
R = sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]])
Qp = R * Q * R.T
lhs = (Qp[0, 0] - Qp[1, 1]) + 2 * sp.I * Qp[0, 1]
rhs = sp.exp(2 * sp.I * th) * ((Qxx - Qyy) + 2 * sp.I * Qxy)   # active rotation; passive convention flips the sign
resid = sp.simplify((lhs - rhs).rewrite(sp.exp).expand())
print("  (Q'_xx - Q'_yy) + 2i Q'_xy  -  e^{2i theta} x original =", resid)
assert resid == 0
print("  -> (h_+, h_x) = (Q_xx - Q_yy, 2 Q_xy) is exactly helicity -+2")
psi = np.linspace(0, 2 * np.pi, 9)
print("  antenna pattern of the differential arm, h_+ cos2psi + h_x sin2psi:")
print("  pi-periodic (quadrupolar):",
      np.allclose(np.cos(2 * psi), np.cos(2 * (psi + np.pi))),
      "| a helicity-1 pattern cos(psi) is not:",
      np.allclose(np.cos(psi), np.cos(psi + np.pi)))

print()
print("=" * 72)
print("[3] MASSLESSNESS: collinear pair kinematics")
print("=" * 72)
k1, k2, c = sp.symbols('k1 k2 c', positive=True)
omega = c * (k1 + k2)          # two gapless quanta, omega_i = c k_i
ktot = k1 + k2                 # collinear momenta add
print("  omega - c|k_total| =", sp.simplify(omega - c * ktot),
      "  -> the covariance mode is gapless and moves at c, identically")

print()
print("=" * 72)
print("[4] MOMENT ORTHOGONALITY: coherent vs squeezed")
print("=" * 72)
rng = np.random.default_rng(3)
N = 400_000
base = rng.normal(size=(N, 2))                      # isotropic zero-point cloud
coh = base + np.array([1.7, 0.0])                   # coherent displacement (EM)
sq = base @ np.diag([1.25, 0.8])                    # squeeze x, anti-squeeze y (GW)
for name, X in (("vacuum   ", base), ("coherent ", coh), ("squeezed ", sq)):
    m = X.mean(0)
    C = np.cov(X.T)
    print(f"  {name}: <f> = ({m[0]:+.3f},{m[1]:+.3f})   "
          f"Q_xx - Q_yy = {C[0,0]-C[1,1]:+.3f}   Q_xy = {C[0,1]:+.3f}")
print("  -> displacement moves the mean and not the covariance;")
print("     squeezing moves the covariance and not the mean. Orthogonal channels.")

print()
print("=" * 72)
print("[5] LOG-CORRECTION COUNT (no graviton loop)")
print("=" * 72)
photon, microrot = sp.Rational(-26, 90), sp.Rational(-26, 90)
lattice = photon + microrot
gr = sp.Rational(424 - 26, 90)
print(f"  lattice: {photon} + {microrot} = {lattice} = {sp.nsimplify(lattice)}")
print(f"  GR     : (424 - 26)/90 = {gr}")
print(f"  LQG    : approximately -3/2")
print("  Three theories, three numbers; the lattice value now rests on the")
print("  absence of a fundamental metric field, not on a vector graviton.")
