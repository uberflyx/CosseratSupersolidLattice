#!/usr/bin/env python3
"""The planar chiral channel is exactly zero, and that retires my own 0.036.

WHAT THIS SETTLES
Solving the localised core-libration problem properly (bound state plus full
Dalgarno-Lewis mode sum, no truncation) on the settled edge width
w = 0.590 ell returned a mass coefficient of 0.036 to 0.071 depending on
which core strain sourced it. That number is an ARTEFACT and is retracted
here, because the channel it describes vanishes identically.

The index argument, checked term by term:
  edge with b along x, line along z, glide normal y
  nonzero strain    : eps_(xx), eps_(xy), eps_(yy)   -> both indices in {x,y}
  microrotation     : only phi_z (it tracks the macrorotation, cos = 0.994)
  nonzero curvature : kappa_(zx), kappa_(zy)         -> one index is z
  matching contractions in eps_(ij) kappa_(ij): NONE
So the planar vertex is zero term by term, not merely small. My 1D reduction
broke that structure and leaked a spurious 0.036. Any calculation that
collapses the tensor indices to one dimension will do the same.

WHAT THE MECHANISM MUST THEREFORE BE
The monograph is explicit and this computation confirms the need for it: the
chirality generates out-of-plane corrections u_z, phi_x, phi_y at O(theta),
and only those produce matching index pairs. The chiral energy evaluated on
the relaxed configuration is then O(theta) from the vertex prefactor times
O(theta) from the induced field, i.e. O(theta^2), which is the power counting
m_1 = theta_ch^2 m_0 without any planar overlap ever appearing. Equivalently,
it is the relaxation energy -S^2/(2K) of the chirally sourced out-of-plane
sector, which is why the second-order and first-order-on-relaxed-state
descriptions coincide.

WHAT I COULD NOT DO, STATED PLAINLY
Two attempts at that out-of-plane channel today, both by guessing which
gradient sources which field, returned 1e-4 and are not reported as results.
Guessing the tensor chain does not converge on the right answer; the chain
has to come out of the hemitropic constitutive equations. The correct
calculation is the coupled system

    (mu + kappa_c) grad^2 u + kappa_c curl phi + C_ch (chiral cross terms) = 0
    gamma grad^2 phi - 2 kappa_c phi + kappa_c curl u + C_ch (chiral cross) = 0

solved as a boundary-value problem on the PN disregistry with the O(theta)
out-of-plane components retained, which is exactly the "effective Cosserat
boundary-value problem" the monograph flags as open. It is a substantial
symbolic derivation, not an evening's numerics.

STATUS OF THE COEFFICIENT: still open, and every value I have produced this
week (0.01, 0.036, 0.06, 0.30, 1.29) is retracted with a stated reason. The
one thing now certain is that the planar channel contributes exactly nothing,
which removes an entire class of wrong attempts.
"""

def matching_contractions():
    eps = {('x', 'x'), ('x', 'y'), ('y', 'y')}
    kappa = {('z', 'x'), ('z', 'y')}
    sym = lambda p: frozenset([p, p[::-1]])
    ksym = {sym(q) for q in kappa}
    return {p for p in eps if sym(p) in ksym}

def main():
    m = matching_contractions()
    print("planar edge, eps_(ij) kappa_(ij):")
    print(f"  matching index pairs = {m or 'NONE'}")
    print("  => the planar chiral channel vanishes term by term.")
    print("  => the 0.036 from the 1D Dalgarno-Lewis run is index leakage,")
    print("     and is retracted. The mechanism is the O(theta) out-of-plane")
    print("     sector, which needs the full hemitropic BVP.")
    assert not m

if __name__ == "__main__":
    main()
