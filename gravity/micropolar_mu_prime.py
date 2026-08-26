#!/usr/bin/env python3
"""
micropolar_mu_prime.py -- the photon's modulus under pressure, with the
microrotation allowed to do its work.

The plain dynamical matrix freezes the node rotors, so its transverse branch
returns mu_tot = mu + kappa_c, the short-wavelength shear.  The photon lives at
long wavelength, where the rotors relax adiabatically and the governing modulus
drops to mu_bar.  This script builds the six-degree-of-freedom (three
translations, three rotor angles) dynamical matrix for the D_4 spatial slice,
eliminates the rotor block by its Schur complement in the k -> 0 limit, and
differentiates the surviving transverse eigenvalue against pressure.

Contact model, per bond of 4-vector R = (R_s, R_4), L = |R|, n = R/L:
  normal spring   k_n = V''(L), loaded by the component of the relative centre
                  displacement along n (translations carry no compact part, so
                  this is Delta . R_s / L);
  tangential      k_t = r V''(L)  (eta_t = 7: the tangential stiffness tracks
                  the normal curvature), loaded by the relative displacement of
                  the contact surfaces perpendicular to n, which is where the
                  rotors enter with lever arm R_s / 2;
  bond tension    V'(L)/L acting on the perpendicular centre displacement, the
                  geometric term a stretched bond carries; zero at equilibrium
                  but its pressure derivative is not.

Pressure comes from the same total contact energy under uniform spatial
dilation (rotors vanish there by cubic symmetry), so the numerator and the
denominator of mu' = d(rho v_T^2)/dP share one energy function.

Validation: at r = 0 the rotors decouple, and the relaxed and frozen values
must agree with the plain instrument; at r > 0 and k along [100], the frozen
value must reproduce Lambda (1 + 5r) and the relaxed value Lambda (1 + 2r),
the corpus's mu_tot and mu_bar.
"""
import itertools

import numpy as np

A_M = 7.0 / 3.0
D_M = 1.0 / (2.0 * A_M**2)
S = 1.0 / np.sqrt(2.0)
VOL0 = 1.0 / np.sqrt(2.0)


def V0f(L):
    e = np.exp(-A_M * (L - 1.0))
    return D_M * ((1.0 - e)**2 - 1.0)


def V1f(L):
    e = np.exp(-A_M * (L - 1.0))
    return 2.0 * D_M * A_M * (1.0 - e) * e


def V2f(L):
    e = np.exp(-A_M * (L - 1.0))
    return 2.0 * D_M * A_M**2 * (2.0 * e * e - e)


def shell():
    """(spatial part, compact part) of the 24 D_4 bonds at unit scale."""
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


SP, CP = shell()


def cross_mat(v):
    """Matrix form of the spatial cross product v x (.)."""
    return np.array([[0.0, -v[2], v[1]],
                     [v[2], 0.0, -v[0]],
                     [-v[1], v[0], 0.0]])


def energy_uniform(scale, r):
    """Total contact energy per node under uniform spatial dilation.

    The centre displacement across a bond is (scale - 1) R_s, purely spatial,
    while the bond axis n is the 4-vector: crossing bonds therefore load their
    tangential springs even under a pure dilation, and that storage belongs in
    the pressure."""
    E = 0.0
    for b in range(24):
        Rs0, c4 = SP[b], CP[b]
        Rs = Rs0 * scale
        L = np.sqrt(Rs @ Rs + c4 * c4)
        E += V0f(L)
        # tangential storage: perpendicular part of the affine displacement,
        # measured against the *reference* bond frame the spring was set in
        L0 = 1.0
        n_s, n_4 = Rs0 / L0, c4 / L0
        d = (scale - 1.0) * Rs0                       # spatial displacement
        d_par = d @ n_s                               # component along n (4D)
        d_perp2 = d @ d - d_par**2
        E += 0.5 * (r * V2f(L0)) * d_perp2
    return 0.5 * E                                    # half: bonds shared


def pressure(scale, r, h=1e-5):
    dE = (energy_uniform(scale + h, r) - energy_uniform(scale - h, r)) / (2 * h)
    return -dE / (3.0 * VOL0 * scale**2)


def blocks(kvec, scale, r):
    """H_UU, H_UW, H_WW of the quadratic energy for a plane wave at kvec."""
    HUU = np.zeros((3, 3), complex)
    HUW = np.zeros((3, 3), complex)
    HWW = np.zeros((3, 3), complex)
    for b in range(24):
        Rs = SP[b] * scale
        c4 = CP[b]
        L = np.sqrt(Rs @ Rs + c4 * c4)
        kn, kt = V2f(L), r * V2f(L)
        tension = V1f(L) / L
        m = Rs / L                                    # spatial direction cosines
        ph = np.exp(1j * (kvec @ Rs))
        a = ph - 1.0                                  # Delta = a U
        s2 = 1.0 + ph                                 # omega_i + omega_j = s2 W
        X = cross_mat(Rs / 2.0)                       # (W x R_s/2) = -X^T W = X W? define below
        # surface relative displacement g = a U - s2 (W x R_s/2)
        # with (W x v) = -cross_mat(v) W  =>  g = a U + s2 X W, X = cross_mat(Rs/2)
        # normal component: n . g uses only the spatial part of n, i.e. m . g
        # E_b = 1/2 kn |m.g_n|^2 ... careful: normal loading is centre-to-centre
        #   (rotations do not stretch the bond at linear order): delta_n = m . (aU)
        # tangential: delta_t = g - n (n.g); |delta_t|^2 = |g|^2 - (m.g)^2
        #   (g has no compact part, n.g = m.g)
        # tension: acts on centre displacement only: |aU|^2 - (m.aU)^2
        aU = a * np.eye(3)
        # normal spring
        v = m[:, None] * a                            # row vector effect: (m.aU) = a m^T U
        HUU += kn * np.outer(m, m) * (a.conjugate() * a)
        # tension on centre displacement
        HUU += tension * (np.eye(3) - np.outer(m, m)) * (a.conjugate() * a)
        # tangential on surface displacement g = a U + s2 X W
        # |delta_t|^2 = g^H (I - m m^T) g
        Pp = np.eye(3) - np.outer(m, m)
        HUU += kt * Pp * (a.conjugate() * a)
        HUW += kt * (a.conjugate() * s2) * (Pp @ X)
        HWW += kt * (s2.conjugate() * s2) * (X.T @ Pp @ X)
    # halves: energy per node counts each bond once with factor 1/2 twice over
    return 0.5 * HUU, 0.5 * HUW, 0.5 * HWW


def relaxed_speeds(khat, scale, r, kmag=1e-4):
    """rho v^2 for the three acoustic branches with rotors eliminated."""
    kvec = khat * kmag
    HUU, HUW, HWW = blocks(kvec, scale, r)
    if r == 0.0:                              # no tangential springs: rotors
        Heff = HUU                            # decouple and nothing relaxes
    else:
        Heff = HUU - HUW @ np.linalg.solve(HWW, HUW.conjugate().T)
    vol = VOL0 * scale**3
    w = np.linalg.eigvalsh(Heff) / (kmag**2 * vol)
    return np.sort(w.real)


def frozen_speeds(khat, scale, r, kmag=1e-4):
    kvec = khat * kmag
    HUU, _, _ = blocks(kvec, scale, r)
    vol = VOL0 * scale**3
    return np.sort(np.linalg.eigvalsh(HUU).real / (kmag**2 * vol))


def mu_prime(khat, r, ds=1e-3, relaxed=True):
    f = relaxed_speeds if relaxed else frozen_speeds
    out = []
    for s in (1.0 - ds, 1.0 + ds):
        out.append((pressure(s, r), f(khat, s, r)[0]))
    (Pm, mm), (Pp, mp) = out
    return (mp - mm) / (Pp - Pm)


if __name__ == "__main__":
    LAM = 1.0 / (2.0 * VOL0)                 # Born prefactor V'' l^2 / (2 Omega)
    r_star = 1.0 / (3.0 * np.pi - 5.0)
    k100 = np.array([1.0, 0.0, 0.0])
    k110 = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    k111 = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)

    print("Validation at zero pressure (expect mu_tot = L(1+5r), "
          "mu_bar = L(1+2r)):\n")
    print("     r      frozen [100]   L(1+5r)    relaxed [100]   L(1+2r)")
    for r in (0.0, 0.1, r_star, 0.3):
        fz = frozen_speeds(k100, 1.0, r)[0]
        rx = relaxed_speeds(k100, 1.0, r)[0]
        print(f"  {r:6.4f}   {fz:11.6f}  {LAM*(1+5*r):9.6f}"
              f"   {rx:12.6f}  {LAM*(1+2*r):9.6f}")

    print("\nIsotropy of the relaxed transverse branch at the rolling point:")
    for tag, k in (("[100]", k100), ("[110]", k110), ("[111]", k111)):
        print(f"  {tag}: rho v_T^2 = {relaxed_speeds(k, 1.0, r_star)[0]:.6f}")

    print("\nPressure derivative of the photon's modulus, "
          "d(rho v_T^2)/dP at P = 0:\n")
    print("     r      frozen      relaxed")
    for r in (0.0, 0.1, r_star, 0.3):
        fz = mu_prime(k110, r, relaxed=False)
        rx = mu_prime(k110, r, relaxed=True)
        print(f"  {r:6.4f}   {fz:8.4f}    {rx:8.4f}")
    print("\n  Direction dependence of the relaxed mu' at the rolling point:")
    for tag, k in (("[100]", k100), ("[110]", k110), ("[111]", k111)):
        print(f"    {tag}: mu' = {mu_prime(k, r_star):.4f}")

    print("\nRelaxed birefringence under spatial dilation, L_4 fixed:")
    for s in (1.000, 1.005, 1.010):
        w = relaxed_speeds(k110, s, r_star)
        split = np.sqrt(w[1] / w[0]) - 1.0
        print(f"  e = {3*(s-1):+6.3f}  [110]: vT1^2 = {w[0]:.6f}  "
              f"vT2^2 = {w[1]:.6f}   dn/n = {split:.3e}")
