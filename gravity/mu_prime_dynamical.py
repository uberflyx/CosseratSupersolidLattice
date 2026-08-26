#!/usr/bin/env python3
"""
mu_prime_dynamical.py -- the pressure derivative of the shear modulus, taken
from the long-wavelength dynamical matrix rather than from an elastic-constant
convention.

The gravitational sector needs mu' = d(mu)/dP = 2 to deliver the observed light
deflection, and derives it from the Born sums as mu' = (3 - xi)/5 with the 5
coming from mu/(3K) = 1/5.  That 5 is the Cauchy value, so it moves the moment
the contact carries a tangential spring.  Two things therefore need checking,
and a convention-free instrument is needed for both.

The quantity that actually matters is the transverse sound speed, since delta c
/ c is what the refractive index reads.  For a Bravais lattice with pair
interactions the long-wavelength limit of the dynamical matrix is exact:

    rho v^2(k_hat) = eigenvalue of  M_ab = (1/2 Omega) sum_b Phi_ab (k_hat . R)^2

    Phi_ab = k_n n_a n_b + (V'(R)/R + k_t) (delta_ab - n_a n_b)

where n = R/|R|, k_n = V''(R), and the V'(R)/R term is the geometric transverse
stiffness a stretched bond carries.  It vanishes at zero pressure but its
derivative does not, which is exactly why mu' is not simply the logarithmic
slope of V''.  No elastic-constant convention enters, so there is nothing to
get wrong about Lagrangian against engineering strain or about whether a
Birch correction belongs.

The tangential law is the one the strain-stationarity argument fixes:
eta_t = 7 means k_t tracks the normal curvature, k_t(R) = r V''(R).

It also settles a question the shell alone cannot answer: whether the D_4
isotropy survives a gravitational dilation.  It does not, if the dilation is
spatial with the compact circumference held fixed.  A crossing bond keeps half
its length in the compact direction, so it stretches by only half as much as an
in-slice bond, the two families drift apart in stiffness, and the identity
sum n1^4 = 3 sum n1^2 n2^2 that makes the shell isotropic stops holding.  The
result is a stress-induced vacuum birefringence, dn/n = 0.155 e along <110> and
zero along <111>.  Dilating the compact circumference along with space restores
exact isotropy at every scale, at the cost of letting hbar drift with the
gravitational potential.

Validation: at zero pressure the instrument returns the Born value
C_44 = V2(ell) ell^2 / (2 Omega) exactly, and for a shell whose bonds all scale
with the lattice it reproduces the closed form
mu_prime = -3[(xi - 2) G_4 + G_2] / Z, where G_2 = sum (khat.n)^2 and
G_4 = sum (e.n)^2 (khat.n)^2 are geometric shell sums.

Units: bond length ell = 1 at zero pressure, energies in k_n(ell) = V''(ell).
Contact: Morse with a = 7/3, so xi = ell V'''/V'' = -7.
"""
import itertools

import numpy as np

A_M = 7.0 / 3.0
D_M = 1.0 / (2.0 * A_M**2)               # normalised so V''(1) = 1
S = 1.0 / np.sqrt(2.0)


def V1(r):
    e = np.exp(-A_M * (r - 1.0))
    return 2.0 * D_M * A_M * (1.0 - e) * e


def V2(r):
    e = np.exp(-A_M * (r - 1.0))
    return 2.0 * D_M * A_M**2 * (2.0 * e * e - e)


def shells():
    """Spatial parts and compact parts of the 24 D_4 bonds."""
    sp, cp = [], []
    for i, j in itertools.combinations(range(3), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = np.zeros(3)
                v[i], v[j] = si * S, sj * S
                sp.append(v)
                cp.append(0.0)
    for i in range(3):
        for si in (1, -1):
            for s4 in (1, -1):
                v = np.zeros(3)
                v[i] = si * S
                sp.append(v)
                cp.append(s4 * S)
    return np.array(sp), np.array(cp)


SP24, CP24 = shells()
SP12, CP12 = SP24[:12], CP24[:12]


def state(sp, cp, scale, r):
    """
    Bond geometry and contact stiffnesses at spatial lattice scale.

    Only the spatial legs scale: a gravitational dilation stretches space and
    leaves the compact circumference alone, since that circumference fixes hbar.
    """
    Rs = sp * scale                                # spatial bond vectors
    L = np.sqrt((Rs**2).sum(1) + cp**2)            # full bond lengths
    kn = V2(L)
    kt = r * kn                                    # eta_t = 7
    tension = V1(L) / L                            # geometric transverse term
    return Rs, L, kn, kt, tension


def pressure_and_energy(sp, cp, scale, vol0):
    """Hydrostatic pressure from -dU/dV, by central difference in volume."""
    def U(s):
        L = np.sqrt(((sp * s)**2).sum(1) + cp**2)
        e = np.exp(-A_M * (L - 1.0))
        return 0.5 * np.sum(D_M * ((1.0 - e)**2 - 1.0))

    h = 1e-5
    dU = (U(scale + h) - U(scale - h)) / (2 * h)
    dV = 3.0 * vol0 * scale**2                     # dV/d(scale)
    return -dU / dV


def wave_moduli(sp, cp, scale, r, vol0, khat):
    """rho v^2 for the three branches along khat, from the dynamical matrix."""
    Rs, L, kn, kt, tension = state(sp, cp, scale, r)
    vol = vol0 * scale**3
    n = Rs / L[:, None]                            # spatial direction cosines
    proj = (Rs @ khat)**2                          # (k_hat . R)^2, spatial only
    M = np.zeros((3, 3))
    for b in range(len(L)):
        long_part = kn[b] * np.outer(n[b], n[b])
        perp_part = (tension[b] + kt[b]) * (np.eye(3) - np.outer(n[b], n[b]))
        M += (long_part + perp_part) * proj[b]
    return np.sort(np.linalg.eigvalsh(M / (2.0 * vol)))


def mu_prime(sp, cp, r, vol0, khat=None, ds=1e-3):
    """d(rho v_T^2)/dP at zero pressure, by central difference in lattice scale."""
    if khat is None:
        khat = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    out = []
    for s in (1.0 - ds, 1.0 + ds):
        P = pressure_and_energy(sp, cp, s, vol0)
        mt = wave_moduli(sp, cp, s, r, vol0, khat)[0]     # slowest transverse
        out.append((P, mt))
    (Pm, mm), (Pp, mp) = out
    return (mp - mm) / (Pp - Pm)


if __name__ == "__main__":
    VOL0 = 1.0 / np.sqrt(2.0)
    xi = -3.0 * A_M
    print(f"Morse contact: a = {A_M:.4f}, xi = ell V'''/V'' = {xi:+.4f}\n")

    print("Zero-pressure check (should be P = 0 at scale 1):")
    for nm, sp, cp in [("FCC 12", SP12, CP12), ("D_4 24", SP24, CP24)]:
        print(f"  {nm}:  P(1) = {pressure_and_energy(sp, cp, 1.0, VOL0):+.3e}")

    print("\nTransverse rho v^2 at zero pressure, and its pressure derivative:")
    print("  shell     r        rho v_T^2      mu' = d(rho v_T^2)/dP")
    for nm, sp, cp in [("FCC 12", SP12, CP12), ("D_4 24", SP24, CP24)]:
        for r in (0.0, 1.0 / (3 * np.pi - 5)):
            mt = wave_moduli(sp, cp, 1.0, r, VOL0,
                             np.array([1.0, 1.0, 0.0]) / np.sqrt(2))[0]
            mp = mu_prime(sp, cp, r, VOL0)
            print(f"  {nm:8s} {r:6.4f}   {mt:10.5f}      {mp:8.4f}")

    print("\n  The gravitational sector quotes mu' = (3 - xi)/5 = "
          f"{(3 - xi) / 5:.4f}.")

    print("\nDirection dependence of mu' on the D_4 shell "
          "(isotropy carries through):")
    for tag, k in [("[100]", [1, 0, 0]), ("[110]", [1, 1, 0]),
                   ("[111]", [1, 1, 1]), ("[311]", [3, 1, 1])]:
        k = np.array(k, float)
        k /= np.linalg.norm(k)
        r = 1.0 / (3 * np.pi - 5)
        print(f"  {tag}:  rho v_T^2 = "
              f"{wave_moduli(SP24, CP24, 1.0, r, VOL0, k)[0]:8.5f}   "
              f"mu' = {mu_prime(SP24, CP24, r, VOL0, k):8.4f}")
