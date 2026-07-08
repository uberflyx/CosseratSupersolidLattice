"""
The Hagedorn skin of a black hole: where the tangle ensemble is prepared.

Physics context (monograph, black-holes chapter; paper "Black holes as
superfluid droplets"): outside the horizon the cost of adding one unit of
mutual winding between two defect lines is the elastic winding gap

    epsilon(r) = epsilon_0 * psi(r),    psi = mu_eff/mu = 1 - r_s/r = xi^2,

because every elastic energy is linear in the shear modulus, while the
temperature a static thermometer reads climbs by the Tolman-Ehrenfest law

    T_loc(r) = T_H / xi(r).

The free energy of one winding, Delta F = epsilon_0 xi^2 - k_B T_loc s_w,
with s_w ~ 2 ln N the configurational entropy per winding (the choice of
pair among ~N^2/2), therefore changes sign at a unique depth

    xi_*^3 = 2 ln N * k_B T_H / epsilon_0 ,

the local Hagedorn condition of the gauged permutation-invariant matrix
model (O'Connor-Ramgoolam, JHEP 07 (2024) 152) realised as a physical
skin.  Proper depth of the skin below r at redshift xi: rho ~ 2 r_s xi.

The script evaluates the skin depth for a range of hole masses, with and
without the entropy factor, and the healing-length hierarchy at the skin
edge.  epsilon_0 is taken at the node scale m_0 c^2 = m_e c^2 / alpha;
the cube root makes the depth insensitive to its exact value.
"""

from math import log, pi

# CODATA / framework constants
hbar = 1.054571817e-34          # J s
c = 2.99792458e8                # m / s
G = 6.67430e-11                 # m^3 / (kg s^2)
kB = 1.380649e-23               # J / K
M_sun = 1.98841e30              # kg
M_P = 2.176434e-8               # kg (Planck mass)
m_e = 9.1093837015e-31          # kg
alpha = 7.2973525693e-3
m0c2 = m_e * c**2 / alpha       # node rest energy [J], ~70.02 MeV
ell = 2.8179403262e-15          # lattice spacing = r_e [m]


def skin(M, eps0=m0c2, entropy_factor=True):
    """Return (xi_*, proper depth rho_* [m], r_s [m], T_H [K], ln N)."""
    r_s = 2 * G * M / c**2
    T_H = hbar * c**3 / (8 * pi * G * M * kB)
    N = (8 * pi / log(2)) ** 0.5 * M / M_P          # active line budget, k = 2
    s_w = 2 * log(N) if entropy_factor else 1.0
    xi3 = s_w * kB * T_H / eps0
    xi = xi3 ** (1 / 3)
    rho = 2 * r_s * xi
    return xi, rho, r_s, T_H, log(N)


if __name__ == "__main__":
    print(f"epsilon_0 = m0 c^2 = {m0c2/1.602176634e-13:.1f} MeV")
    print(f"{'M/M_sun':>10} {'r_s':>10} {'T_H [K]':>10} {'ln N':>7} "
          f"{'xi_*':>10} {'skin rho_* ':>12} {'zeta(skin) ':>12}")
    for Ms in (1, 30, 4.3e6, 6.5e9):
        M = Ms * M_sun
        xi, rho, r_s, T_H, lnN = skin(M)
        zeta = ell / xi                      # healing length at the skin edge
        def fmt(x):
            for u, s in ((1e3, 'km'), (1.0, 'm'), (1e-2, 'cm'), (1e-3, 'mm'),
                         (1e-9, 'nm'), (1e-15, 'fm')):
                if abs(x) >= u:
                    return f"{x/u:8.2f} {s}"
            return f"{x:.2e} m"
        print(f"{Ms:>10.3g} {fmt(r_s):>10} {T_H:>10.2e} {lnN:>7.1f} "
              f"{xi:>10.2e} {fmt(rho):>12} {fmt(zeta):>12}")
    # robustness: cube-root dependence on epsilon_0 and the entropy factor
    xi1, rho1, *_ = skin(M_sun, entropy_factor=False)
    xi2, rho2, *_ = skin(M_sun, eps0=10 * m0c2)
    _, rho0, *_ = skin(M_sun)
    print(f"\nRobustness at 1 M_sun: rho_* = {rho0*100:.2f} cm (fiducial); "
          f"{rho1*100:.2f} cm without the entropy factor; "
          f"{rho2*100:.2f} cm with epsilon_0 x 10.")
