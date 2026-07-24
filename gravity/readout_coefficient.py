"""The readout coefficient from first principles.

The gravitational-wave section leaves one number to matching: how strongly an
anisotropic displacement covariance delta<u_i u_j> shifts the speed of a
transverse light wave. This script derives it from the pair potential and the
D4 shell, with no matching.

Mechanism. A probe shear wave (propagation zhat, polarisation ehat) stretches
each bond by delta_1 = eps d (n.e)(n.z) at first order in the strain eps, and
by the transverse geometric term at second order. The zero-point cloud dresses
every bond: with bond-relative covariance q_ij = <Du_i Du_j> [m^2], the
fluctuation ensemble shifts the mean bond length and its stiffness,

    <dr>   = (tr q - n.q.n) / (2d)          (transverse fluctuations lengthen)
    <dr^2> = n.q.n                          (longitudinal fluctuations, leading)

    V''_eff = V'' + V''' <dr> + (1/2) V'''' <dr^2>     (dressed stiffness)
    T_eff   = <V'(r)> = V'' <dr> + (1/2) V''' <dr^2>   (fluctuation tension)

and the probe's effective modulus per bond is

    d^2E/deps^2 = V''_eff d^2 (n.e)^2 (n.z)^2 + T_eff d (n.z)^2 [1 - (n.e)^2],

the second term being the induced tension acting through the probe's own
transverse geometry (a string under tension resists sideways displacement).
Summing the D4 shell and inserting the framework's anharmonicity constants,
gamma_3 = V''' d / V'' = -7 (from gravity, mu' = 2) and
gamma_4 = V'''' d^2 / V'' = 343/9 (the Morse value consistent with gamma_3),
gives the dimensionless response kernel

    delta mu(ehat) / mu = Lambda_read * (qtilde_ij e_i e_j - trace part)
                        + Lambda_iso  * tr qtilde,     qtilde = q / d^2,

and the differential light speed between the two transverse polarisations of
an h_+ covariance pattern (qtilde ~ diag(+1,-1,0) q0) is

    [c(x) - c(y)] / c = Lambda_read * q0.

The number, its sign, and the 12-bond-slice comparison are printed below.
Everything uses the corrected 24-bond D4 sums.
"""

import numpy as np
from itertools import permutations, product

G3 = -7.0            # gamma_3 = V''' d / V''   (gravity chapter, mu' = 2)
G4 = 343.0 / 9.0     # gamma_4 = V'''' d^2/V''  (Morse, consistent with gamma_3)

def d4_bonds():
    v = set()
    for p in permutations(range(4), 2):
        for s1, s2 in product((1, -1), repeat=2):
            w = [0, 0, 0, 0]
            w[p[0]], w[p[1]] = s1, s2
            v.add(tuple(w))
    return np.array(sorted(v), float) * np.sqrt(0.5)   # |n| d = d, here d = 1

BONDS = d4_bonds()                 # 24 bonds, unit length
SLICE = BONDS[BONDS[:, 3] == 0]    # the 12 in-slice bonds, for contrast

def mu_eff(bonds, e4, z4, q3):
    """Probe modulus (arbitrary units) for polarisation e4, propagation z4,
    with spatial bond-relative covariance q3 (3x3, units of d^2).
    Fluctuations are spatial; bond projections use full 4D geometry."""
    tot = 0.0
    for n in bonds:
        ns = n[:3]                                  # spatial leg
        nqn = ns @ q3 @ ns
        perp = np.trace(q3) - nqn                   # perturbation of <|Du_perp|^2>:
        #      may be negative for a traceless squeeze; clipping it would be wrong
        dr = perp / 2.0                             # d = 1
        dr2 = nqn
        v2eff = 1.0 + G3 * dr + 0.5 * G4 * dr2      # in units of V''
        teff = dr + 0.5 * G3 * dr2                  # <V'>/V'' d
        ne, nz = n @ e4, n @ z4
        tot += v2eff * ne**2 * nz**2 + teff * nz**2 * (1.0 - ne**2)
    return tot

ex = np.array([1.0, 0, 0, 0]); ey = np.array([0, 1.0, 0, 0]); ez = np.array([0, 0, 1.0, 0])
Q0 = 1e-6                                           # linear-response probe amplitude

for label, bonds in (("D4, 24 bonds", BONDS), ("FCC slice, 12 bonds", SLICE)):
    mu0 = mu_eff(bonds, ex, ez, np.zeros((3, 3)))
    qplus = Q0 * np.diag([1.0, -1.0, 0.0])          # h_+ covariance pattern
    dmx = mu_eff(bonds, ex, ez, qplus) - mu0        # x-polarised probe
    dmy = mu_eff(bonds, ey, ez, qplus) - mu0        # y-polarised probe
    lam_read = (dmx - dmy) / (2.0 * mu0 * Q0)       # delta(c_x - c_y)/c per q0
    qiso = Q0 * np.eye(3)
    lam_iso = (mu_eff(bonds, ex, ez, qiso) - mu0) / (mu0 * Q0)
    print(f"{label}:  mu0 = {mu0:.3f} V'' d^2")
    print(f"    Lambda_read = {lam_read:+.4f}   [ (c_x - c_y)/c = Lambda_read * q0/d^2 ]")
    print(f"    Lambda_iso  = {lam_iso:+.4f}   (isotropic softening per unit tr q/d^2)")

print()
print("Exact value: Lambda_read = 239/18 (rational because gamma_3, gamma_4 are).")
print("Sign: Lambda_read > 0: fuzz fattened along x STIFFENS x-polarised light")
print("(the positive quartic dominates), so the two polarisations split with a")
print("derived gain, and the differential arm reads the h_+ pattern directly.")
print("Structural: the crossing bonds contribute nothing to the anisotropic kernel")
print("(their single spatial leg cannot supply both the polarisation and the")
print("propagation projection), while they do enter the isotropic softening:")
print("the same one-leg rule that sorted the elastic-constant sweep.")
print()
print("Decomposition at the D4 value (per unit q0/d^2):")
mu0 = mu_eff(BONDS, ex, ez, np.zeros((3, 3)))
for name, g3, g4 in (("geometric (V'' tension)", 0.0, 0.0),
                     ("Murnaghan (V''')", G3, 0.0),
                     ("quartic (V'''')", 0.0, G4)):
    globals()['G3'], globals()['G4'] = g3, g4
    d = (mu_eff(BONDS, ex, ez, 1e-6*np.diag([1.,-1.,0.])) -
         mu_eff(BONDS, ey, ez, 1e-6*np.diag([1.,-1.,0.]))) / (2*mu0*1e-6)
    print(f"    {name:<26} {d:+.4f}")
G3, G4 = -7.0, 343.0/9.0
