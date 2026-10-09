#!/usr/bin/env python3
"""
contact_vertex_lever.py
=======================

Where a contact-level reading of the Peierls misfit puts the vertex's lever arm,
compared with the lever arm whose static links are wavelength-independent.

Background
----------
The rank-one Peierls vertex reads the medium through a covector (1, lambda) in
the node's displacement u and rotation phi.  In the isotropic Cosserat continuum
the only real covector that the static crystal answers like a plain solid of the
relaxed modulus mubar, at every wavelength, has

    lambda^2 = N^2 gamma / mubar,

and with the tangent-sphere curvature modulus gamma = 4 mu R^2 (R = l/2) and
mu/mubar = (1 - 2N^2)/(1 - N^2) this is

    lambda^2 / R^2 = 4 N^2 (1 - 2N^2) / (1 - N^2).

Its inertial (f-sum) weight is 1 + lambda^2 rho/J, which with J = gamma/c^2 is
1 + N^2 = 1 + 1/pi at the rolling point.

What the script checks
----------------------
1. A discrete FCC contact lattice (nodes with displacement and rotation, normal
   and tangential contact springs acting on the relative displacement of the
   touching surfaces, delta = u_j - u_i - R (phi_i + phi_j) x n) reproduces the
   framework's homogenised moduli mu_tot and kappa_c, and N^2 = 1/pi at
   r = k_t/k_n = 1/(2 pi - 3), along [100] and [110].  Along [111] the FCC slice
   is cubic-anisotropic; the comparison below is made in the isotropic continuum.
2. A misfit acting at the contact points across a {111} glide plane reads the
   rotation through the lever R n_perp = h/2 = l/sqrt(6), with h = d_111, for
   both in-plane slip directions.
3. The selected lever equals the contact lever only at N^2 = 1/4 and 1/3; in
   between it is larger, by at most a factor (12 - 8 sqrt 2)/(2/3) in lambda^2.
4. Within the identification alpha_G = C_G alpha^19, Newton's constant fixes
   C_0 - 1 = lambda^2 mubar/gamma and hence the lever arm (for the adopted
   gamma = mu l^2).  The script reports the measured lever, the distance of the
   contact value in standard deviations, C_0 for the centre-registry,
   contact-point and selected readings, the contact endpoint when its static
   links sit at the misfit wavevector, and the selected lever for the
   gamma/kappa_c that the density modulation gives instead.
   The contact lever is a one-sided, face-centred-cubic-slice result for a wave
   crossing the plane; for a wave in the plane the levers cancel.
5. The static link weight of the contact-point covector in the isotropic
   continuum, against wavevector.

Usage:    python3 contact_vertex_lever.py
Requires: numpy, sympy
Runtime:  a few seconds.

Repository: https://github.com/uberflyx/CosseratSupersolidLattice
Authors:    Mitchell Cox, Warren Carlson (University of the Witwatersrand)
"""
import itertools

import numpy as np
import sympy as sp

ELL = 1.0                     # nearest-neighbour distance
R = ELL / 2                   # sphere radius, tangent spheres
V_CELL = ELL**3 / np.sqrt(2)  # FCC primitive cell volume
N2_ROLL = 1.0 / np.pi         # rolling-constraint coupling number
C0_G = (0.318313, 0.000030)   # C_0 - 1 from CODATA 2022 G within the identification


def fcc_bonds():
    """The twelve FCC nearest-neighbour unit vectors."""
    vecs = set()
    for p in itertools.permutations((1, 1, 0)):
        for s1, s2 in itertools.product((1, -1), repeat=2):
            v = list(p)
            nz = [i for i in range(3) if v[i] != 0]
            v[nz[0]] *= s1
            v[nz[1]] *= s2
            vecs.add(tuple(v))
    out = [np.array(v, float) / np.sqrt(2) for v in sorted(vecs)]
    assert len(out) == 12
    return out


def cross_matrix(n):
    """Matrix of n x (.)."""
    return np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])


def contact_stiffness(k, kn, kt):
    """6x6 Fourier stiffness per node, each bond counted from both of its nodes.

    The relative displacement of the touching surfaces is
    delta = u_j - u_i - R (phi_i + phi_j) x n = u_j - u_i + R [n]x (phi_i + phi_j),
    with energy (1/2) delta.C.delta, C = kn nn + kt (1 - nn).  Counting every bond
    from each node matches the homogenisation convention C = (l^2/V) sum K n n.
    """
    M = np.zeros((6, 6), complex)
    for n in fcc_bonds():
        C = kn * np.outer(n, n) + kt * (np.eye(3) - np.outer(n, n))
        phase = np.exp(1j * ELL * (k @ n))
        Pj = np.hstack([np.eye(3), R * cross_matrix(n)])
        Pi = np.hstack([-np.eye(3), R * cross_matrix(n)])
        D = Pj * phase + Pi
        M += D.conj().T @ C @ D
    return M / V_CELL


def lattice_check():
    """Effective transverse moduli of the contact lattice along three axes."""
    kn = 1.0
    kt = kn / (2 * np.pi - 3)
    mu = np.sqrt(2) * (kn - kt) / ELL
    kappa = 4 * np.sqrt(2) * kt / ELL
    mu_tot = mu + kappa
    z = np.ones(3) / np.sqrt(3)
    x = np.array([1.0, -1.0, 0.0]) / np.sqrt(2)
    rows = []
    for name, kd, ud, pd in (("[100]", np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0])),
                             ("[110]", np.array([1.0, 1, 0]) / np.sqrt(2), np.array([0, 0, 1.0]),
                              np.array([1.0, -1, 0]) / np.sqrt(2)),
                             ("[111]", z, x, np.cross(z, x))):
        kk = 1e-3
        M = contact_stiffness(kk * kd, kn, kt)
        eu = np.r_[ud, 0, 0, 0]
        ep = np.r_[0, 0, 0, pd]
        k_uu = (eu @ M @ eu).real / kk**2
        k_up = abs(eu @ M @ ep) / kk
        k_pp = (ep @ M @ ep).real
        rows.append((name, k_uu, k_up, k_pp / 2, k_up / (2 * k_uu)))
    for name, k_uu, k_up, k_pp2, n2 in rows[:2]:
        assert abs(k_uu / mu_tot - 1) < 1e-5 and abs(k_up / kappa - 1) < 1e-5
        assert abs(k_pp2 / kappa - 1) < 1e-5 and abs(n2 - N2_ROLL) < 1e-5
    return mu_tot, kappa, rows


def contact_levers():
    """Slip-direction displacement of each contact point across (111) per unit rotation."""
    z = np.ones(3) / np.sqrt(3)
    up = [n for n in fcc_bonds() if n @ z > 1e-9]
    out = {}
    for name, slip in (("<110>", np.array([1.0, -1.0, 0.0]) / np.sqrt(2)),
                       ("<112>", np.array([1.0, 1.0, -2.0]) / np.sqrt(6))):
        axis = np.cross(z, slip)
        out[name] = [float(np.cross(axis, R * n) @ slip) for n in up]
    return len(up), out


def in_plane_levers():
    """Levers for a wave travelling in the (111) plane, perpendicular to a <110> slip.

    The transverse rotation axis is then the plane normal, and the three contact
    levers sum to zero.
    """
    z = np.ones(3) / np.sqrt(3)
    slip = np.array([1.0, -1.0, 0.0]) / np.sqrt(2)
    up = [n for n in fcc_bonds() if n @ z > 1e-9]
    levers = [float(np.cross(z, R * n) @ slip) for n in up]
    assert abs(sum(levers)) < 1e-12
    return levers


def lever_of_coupling():
    """Selected lever^2/R^2 against N^2, its crossings with 2/3 and its maximum."""
    x = sp.symbols('x', positive=True)
    f = 4 * x * (1 - 2 * x) / (1 - x)
    crossings = sorted(sp.solve(sp.Eq(f, sp.Rational(2, 3)), x))
    x_max = [s for s in sp.solve(sp.diff(f, x), x) if s.is_real and 0 < float(s) < 0.5][0]
    f_max = sp.simplify(f.subs(x, x_max))
    return crossings, x_max, f_max, float(f.subs(x, N2_ROLL))


def g_measures_lever():
    """Within the identification, C_0 - 1 = lambda^2 mubar/gamma with J = gamma/c^2."""
    j_over_rho = (np.pi - 2) / (np.pi - 1)      # gamma/mubar in units of l^2
    c0m, dc0 = C0_G
    lam = np.sqrt(c0m * j_over_rho)
    dlam = 0.5 * lam * dc0 / c0m
    selected = np.sqrt(N2_ROLL * j_over_rho)
    contact = 1 / np.sqrt(6)
    readings = {"centre registry (slip only)": 1.0,
                "contact points across {111}": 1 + contact**2 / j_over_rho,
                "selected (static links scale-free)": 1 + selected**2 / j_over_rho,
                "contact radius l/2": 1 + R**2 / j_over_rho}
    assert abs(readings["contact points across {111}"] - (1 + (np.pi - 1) / (6 * (np.pi - 2)))) < 1e-14
    return lam, dlam, selected, contact, (lam - contact) / dlam, readings


def static_weight(lever, x):
    """Stripped static link weight of covector (1, lever) in the isotropic continuum."""
    a = np.pi / (np.pi - 1)
    b = 4 / (np.pi - 2)
    n2 = lever**2 / ((np.pi - 2) / (np.pi - 1))
    return (x * x * (1 + n2 * a) + b) / (a * x * x + b)


def main():
    mu_tot, kappa, rows = lattice_check()
    print("1. Contact lattice, transverse block (framework units, k_n = 1, r = 1/(2pi - 3))")
    print(f"   framework: mu_tot = {mu_tot:.6f}, kappa_c = {kappa:.6f}, N^2 = 1/pi = {N2_ROLL:.6f}")
    for name, k_uu, k_up, k_pp2, n2 in rows:
        print(f"   {name}: mu_tot {k_uu:.6f}  kappa_c {k_up:.6f} / {k_pp2:.6f}  N^2 {n2:.6f}")
    n_up, levers = contact_levers()
    print(f"\n2. Contact-point lever across (111): {n_up} crossing bonds")
    for name, ls in levers.items():
        print(f"   slip {name}: {np.round(ls, 6)}  (h/2 = l/sqrt6 = {1/np.sqrt(6):.6f})")
    crossings, x_max, f_max, f_roll = lever_of_coupling()
    print(f"   wave in the plane, perpendicular to the slip: levers {np.round(in_plane_levers(), 4)}")
    print("\n3. Selected lever^2/R^2 = 4N^2(1-2N^2)/(1-N^2) against the contact value 2/3")
    print(f"   equal at N^2 = {crossings}; maximum {f_max} = {float(f_max):.6f} at N^2 = {x_max}")
    print(f"   at N^2 = 1/pi: {f_roll:.6f}, {100*(f_roll/(2/3)-1):.2f} per cent above 2/3")
    lam, dlam, selected, contact, nsig, readings = g_measures_lever()
    print("\n4. Within the identification, G measures the lever")
    print(f"   lambda_G = {lam:.6f} +- {dlam:.6f} l;  selected {selected:.6f} l;  contact {contact:.6f} l "
          f"({nsig:.0f} standard deviations)")
    for name, c0 in readings.items():
        print(f"   C_0 for {name:36s} {c0:.6f}  ({1e6*(c0/(1+N2_ROLL)-1):+.0f} ppm)")
    w_misfit = static_weight(contact, 2 * np.pi * np.sqrt(3))
    print(f"   contact endpoint with all 18 links at the misfit wavevector: "
          f"{readings['contact points across {111}']:.5f} x {w_misfit:.5f}^18 = "
          f"{readings['contact points across {111}'] * w_misfit**18:.4f}")
    j_over_rho_derived = 0.509 * 2 / (np.pi - 2) * (np.pi - 2) / (np.pi - 1)
    print(f"   selected lever with the derived gamma/kappa_c = 0.509 l^2: "
          f"{np.sqrt(N2_ROLL * j_over_rho_derived):.4f} l")
    print("\n5. Static link weight of the contact covector against k l")
    for x in (0.01, 1.0, 2 * np.pi, 2 * np.pi * np.sqrt(3), 1e3):
        print(f"   k l = {x:8.3f}: {static_weight(contact, x):.6f}  (selected {static_weight(selected, x):.6f})")


if __name__ == "__main__":
    main()
