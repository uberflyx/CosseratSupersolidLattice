#!/usr/bin/env python3
"""
birefringence_sign_audit.py
===========================
The load-bearing sign check for the cosmic-birefringence mechanism
(Confrontation chapter, "Cosmic birefringence").

Question: for a photon crossing a stack of alternating chirality domains
(+,-,+,-,... along the line of sight), do the screw-dislocation twist walls
carry a SINGLE handedness (so an odd, twist-driven rotation would add), or
does the handedness ALTERNATE (so any odd rotation telescopes)?

This settles it with explicit Frank-Bilby bookkeeping in one fixed lab frame,
using no orientation convention that could pre-decide the answer. The result
is that successive walls carry OPPOSITE screw handedness: they are a
dislocation dipole, forced by the fact that the net twist across (+ -> - -> +)
is zero. Any per-crossing rotation odd in the wall's twist therefore
telescopes, exactly like the axion pseudoscalar, and cannot produce the
observed ~0.3 deg isotropic angle. The accumulating mechanism must instead
live in the EVEN part of the coupling (see birefringence_ribbon_coupling.py).

Units: theta in units of theta_ch; lengths in the lattice spacing ell.
"""
import numpy as np

np.set_printoptions(precision=3, suppress=True)


def wall_screw_pseudoscalar(phi_below, phi_above):
    """Frank-Bilby content of a pure twist wall between two grains whose
    helical phases (about the fixed global axis +z) are phi_below (behind)
    and phi_above (ahead of) the wall, for a photon travelling +z.

    Relative rotation of the grain ahead w.r.t. the grain behind:
        Omega = (phi_above - phi_below) * z_hat.
    Net Burgers vector crossed by a unit in-plane segment V = x_hat:
        B = Omega x V   (Frank's formula, small angle).
    The screw handedness pseudoscalar is b . xi with line sense xi = +y;
    its SIGN is what tells left- from right-handed screw grids, and it is
    invariant under the choice of xi and V (both flip together)."""
    Omega = np.array([0.0, 0.0, phi_above - phi_below])
    V = np.array([1.0, 0.0, 0.0])
    B = np.cross(Omega, V)
    xi = np.array([0.0, 1.0, 0.0])
    return float(np.dot(B, xi))


def run_audit(n_show=6):
    phis = [(+1.0 if i % 2 == 0 else -1.0) for i in range(n_show + 1)]
    print("wall   grains        b.xi     handedness")
    handed = []
    for i in range(n_show):
        h = wall_screw_pseudoscalar(phis[i], phis[i + 1])
        handed.append(np.sign(h))
        print(f"  {i+1:2d}   {phis[i]:+.0f} -> {phis[i+1]:+.0f}    {h:+.2f}    "
              f"{'RIGHT' if h > 0 else 'LEFT'}-handed")
    alt = all(handed[i] == -handed[i + 1] for i in range(len(handed) - 1))
    print(f"\nhandedness alternates wall to wall: {alt}")
    print("forced, because net twist across (+ -> - -> +) = 0: the two walls")
    print("are a dislocation dipole and their screw content cancels.\n")
    return alt


def telescoping_consequence(N=940, n_trials=4000, seed=0):
    """Show that an ODD (twist-driven) per-wall rotation, with any incidence
    weighting, telescopes to an endpoint-bounded value: dead against 0.3 deg."""
    rng = np.random.default_rng(seed)
    alpha = 1 / 137.035999177
    th_deg = np.degrees(alpha**2 / (2 * np.pi))
    signs = np.array([(-1) ** i for i in range(N)])       # alternating handedness
    # flux-weighted incidence cos: p(c) = 2c on [0,1]  ->  c = sqrt(U)
    iso = []
    for _ in range(n_trials):
        c = np.sqrt(rng.random(N))
        iso.append((signs * c).sum())                     # odd coupling * weight
    iso = np.array(iso)
    print("ODD coupling (rotation proportional to signed wall twist):")
    print(f"  unweighted alternating sum        = {signs.sum():d} theta_ch (exactly 0)")
    print(f"  |cos|-weighted: mean {iso.mean():+.2f}, rms {iso.std():.1f} theta_ch")
    print(f"  => isotropic |beta| <~ 1 theta_ch = {th_deg:.1e} deg")
    print(f"     (measured ~0.3 deg: the odd mechanism is excluded by ~600x)")


if __name__ == "__main__":
    run_audit()
    telescoping_consequence()
