"""Slip-patch geometry of the neutrino's edge dislocation, from FCC periods.

The neutrino is a full edge dislocation on a close-packed {111} plane: Burgers
vector b along <110> with |b| = ell, line direction xi along the in-plane
<112>.  Two lengths fix the patch.  The width L_b is the plane-strain
Peierls-Nabarro core half-width, solved in edge_core_width.py.  The length
L_xi is what this script establishes, two independent ways, and then
propagates into the coherent form factor.

Route 1, lattice commensurability.  A segment periodic along its own line must
be periodic under a lattice translation along that line, so L_xi is a whole
primitive period of <112>.  The primitive periods are found exactly here by
scanning multiples of a/6 along each direction and testing the face-centred
cubic membership condition: written in units of a/2, a translation must have
all-integer components whose sum is even.

Route 2, the compact wrap.  A step along the compact D_4 direction is one ABC
stacking step, which shifts the layer's registry in-plane by the Shockley
vector (a/6)<112>.  Three steps close the compact circle (L_4 = 3 d_111), so
the brane slides 3 x (a/6)<112> = (a/2)<112> along the line as it wraps once.

Both give L_xi = sqrt(3) ell, and they agree because three Shockley steps along
<112> are the primitive period of that direction.

Units: lengths in ell (nearest-neighbour spacing), with the conventional cubic
cell edge a = sqrt(2) ell.  Numerical output is in fm at ell = r_e.
"""

import numpy as np
import sympy as sp
from scipy.special import j1

# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------
ELL_FM = 2.8179403205          # classical electron radius r_e [fm], CODATA 2022
HBARC = 197.3269804            # hbar c [MeV fm]
W_EDGE = 0.5902556948955939    # plane-strain PN core half-width [ell]
                               # (edge_core_width.solve_width(4/3))

A = sp.sqrt(2)                 # cubic cell edge, in units of ell


# --------------------------------------------------------------------------
# Route 1: primitive periods of the FCC lattice
# --------------------------------------------------------------------------
def is_fcc_translation(components_in_half_a):
    """FCC membership test.

    In units of a/2 a face-centred cubic translation has all-integer
    components whose sum is even; the odd-sum integer vectors belong to the
    complementary simple-cubic sublattice and are not lattice translations.
    """
    comps = [sp.nsimplify(c) for c in components_in_half_a]
    if not all(c.is_integer for c in comps):
        return False
    return sp.Integer(sum(comps)) % 2 == 0


def primitive_period(direction, max_sixths=60):
    """Shortest FCC translation along `direction`, in units of ell.

    Scans multiples of a/6 so that Shockley partials, which are a/6<112> and
    deliberately not translations, are inside the search grid and are seen to
    fail the test rather than being excluded by a coarse step.
    """
    d = sp.Matrix(direction)
    for n in range(1, max_sixths + 1):
        # n*(a/6)*direction, expressed in units of a/2
        comps = [sp.Rational(n, 3) * c for c in d]
        if is_fcc_translation(comps):
            return sp.simplify(sp.Rational(n, 3) * d.norm() * A / 2)
    raise RuntimeError(f"no translation found along {direction}")


# --------------------------------------------------------------------------
# Route 2: the slide accumulated over one wrap of the compact direction
# --------------------------------------------------------------------------
def compact_wrap_slide():
    """In-plane displacement carried by one wrap of the compact direction.

    Returns (slide_vector_in_ell, magnitude_in_ell, n_steps).  One step along
    the compact direction advances the stacking label by one, which in the
    close-packed plane is a rigid registry shift by (a/6)<112>.  The compact
    circumference is L_4 = 3 d_111, so one wrap is three such steps.
    """
    shockley = A / 6 * sp.Matrix([1, 1, -2])
    d111 = sp.simplify(A / sp.sqrt(3))
    n_steps = sp.simplify(3 * d111 / d111)          # L_4 / d_111 = 3, kept explicit
    slide = n_steps * shockley
    return sp.simplify(slide), sp.simplify(slide.norm()), int(n_steps)


# --------------------------------------------------------------------------
# Downstream: the uniform-disk form factor of the patch
# --------------------------------------------------------------------------
def effective_radius(L_xi_fm, L_b_fm):
    """Radius of the disk of equal area, which is what scattering resolves."""
    return np.sqrt(L_xi_fm * L_b_fm / np.pi)


def coherent_deficit(R_fm, q_MeV):
    """1 - |F(q)|^2 for a uniform disk, F(q) = 2 J_1(qR)/(qR)."""
    x = q_MeV * R_fm / HBARC
    return 1.0 - (2.0 * j1(x) / x) ** 2


def main():
    print("FCC primitive periods (units of ell, a = sqrt(2) ell)")
    for label, direction in [("<110>  Burgers vector, screw line", (1, -1, 0)),
                             ("<112>  edge line direction      ", (1, 1, -2)),
                             ("<111>  stacking normal          ", (1, 1, 1))]:
        p = primitive_period(direction)
        print(f"  {label}: {p}  = {float(p):.5f} ell")

    shockley_mag = sp.simplify(A / 6 * sp.Matrix([1, 1, -2]).norm())
    print(f"\n  Shockley partial (a/6)<112>       : {shockley_mag}"
          f"  = {float(shockley_mag):.5f} ell   [not a translation]")

    slide, mag, n = compact_wrap_slide()
    print(f"\nCompact wrap: {n} stacking steps, slide = {list(slide.T)} ell")
    print(f"  |slide| = {mag} = {float(mag):.5f} ell")
    assert sp.simplify(mag - primitive_period((1, 1, -2))) == 0, \
        "the two routes must give the same length"
    print("  agrees with the <112> primitive period: OK")

    L_xi = float(mag) * ELL_FM
    L_b = W_EDGE * ELL_FM
    d = ELL_FM / np.sqrt(3)                     # partial period
    R = effective_radius(L_xi, L_b)
    print(f"\nSlip patch at ell = {ELL_FM:.4f} fm")
    print(f"  L_xi  = {L_xi:.3f} fm = {L_xi/d:.3f} d")
    print(f"  L_b   = {L_b:.3f} fm = {L_b/d:.3f} d")
    print(f"  R_eff = {R:.3f} fm")

    print("\nCoherent cross-section deficit, 1 - |F(q)|^2")
    for q in (10, 30, 50, 100):
        print(f"  q = {q:3d} MeV : {100*coherent_deficit(R, q):5.2f}%")

    # The width convention is the remaining softness: the PN slip density is a
    # Lorentzian of half-width w, so 2w is an equally defensible patch width.
    R2 = effective_radius(L_xi, 2 * L_b)
    print(f"\nWidth convention L_b = 2w: R_eff = {R2:.3f} fm,"
          f" deficit at 100 MeV = {100*coherent_deficit(R2, 100):.1f}%")


if __name__ == "__main__":
    main()
