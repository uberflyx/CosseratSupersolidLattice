"""
The Hagedorn skin of a black hole: where the record of the interior tangle
would thaw, and where the crystal's own nodes sit against it.

Physics context (monograph, black-holes chapter, Hagedorn-skin subsection).
Both quantities that cross in the skin are measured in one frame, that of a
static observer.  Such an observer measures the crystal's own shear modulus
at every radius where there is crystal, so the elastic cost of adding one
unit of mutual winding between two defect lines is a fixed gap

    epsilon(r) = epsilon_0,    at and above the crystal's lowest nodes,

while the temperature a static thermometer reads climbs by the
Tolman-Ehrenfest law

    T_loc(r) = T_H / xi(r),    xi = sqrt(1 - r_s/r)  (the lapse).

The free energy of one winding, Delta F = epsilon_0 - k_B T_loc s_w, with
s_w ~ 2 ln N the configurational entropy per winding (the choice of pair
among ~N^2/2), therefore changes sign at a unique depth

    xi_* = 2 ln N * k_B T_H / epsilon_0 ,

the local Hagedorn condition of the gauged permutation-invariant matrix
model (O'Connor-Ramgoolam, JHEP 07 (2024) 152).  Near the horizon the
proper height is rho = 2 r_s xi, and with k_B T_H = hbar c/(4 pi r_s) the
Schwarzschild radius cancels:

    rho_* = (ln N / pi) (hbar c / epsilon_0),

about 80 fm at a solar mass, depending on the mass only through ln N.

Taking the crystal's coordinates as the Schwarzschild areal chart, a node
one lattice spacing above the horizon (r = r_s + l) sits at proper height
rho_1 = 2 sqrt(r_s l), 5.8 micrometres for a solar mass.  A static
thermometer there reads k_B T_loc = hbar c/(2 pi rho_1).  A winding there
costs at least epsilon_0, so it carries a Boltzmann exponent of at least

    epsilon_0 / (k_B T_loc) = 4 pi sqrt(r_s l) epsilon_0 / (hbar c),

about 1.3e10 for a solar mass, against the 2 ln N ~ 178 at which the first
winding would thaw.  The lower bound on the cost is clearest in the
crystal's own time, where the heat is k_B T_H at every node: a winding
about a line running out from the horizon strains the bonds along the
horizon, whose stiffness in crystal time is mu_0 with no factor of the
lapse, so no power of the lapse lowers its cost, while the bound needs
only xi_1 epsilon_0 there, with xi_1 = sqrt(l/r_s) ~ 1e-9 (solar mass);
xi_1 epsilon_0 / (k_B T_H) is the same exponent.  (The radial Debye energy,
by contrast, falls as xi^2 in crystal time, which is why its margin over
the heat a full spacing up is only 4 pi^2.)  The whole skin lies within
rho_*^2/(4 r_s) of the horizon in areal height, about 2e-16 l for a solar
mass, so on every node from a full spacing up the record of the tangle
stays frozen; what windings would do at any node the horizon leaves that
low is not derived.

The harmonic chart puts the nodes at the same heights and only widens
their spacing along the horizon, which can only raise the gap.  The
isotropic chart puts a node n spacings up at proper height 4 n l at every
mass.  With the gap epsilon_0 the exponent there is
8 pi n l epsilon_0/(hbar c) ~ 25 n, so the skin as computed reaches the
seventh node; the static observer reads every spacing as 4 l, which raises
a node-scale winding to ~ 4^3 epsilon_0 and the exponent at the lowest
node to ~ 1600, only nine times the 178 needed to thaw, a margin inside
the uncertainty of epsilon_0.  Even on that gap, a node within a ninth of
a spacing of the horizon lies inside the skin.

The Tolman heating itself needs the near-horizon fields in the state that
is regular for freely falling observers.  The crystal's own shear field has
a largest wavenumber and ends, after its formation transient, in its
Killing ground state (horizon-cutoff section), which heats no node at all;
on that reading the record is frozen a fortiori.

The script evaluates the depth for a range of hole masses, with and without
the entropy factor, its (linear) sensitivity to epsilon_0, which is taken
at the node scale m_0 c^2 = m_e c^2 / alpha, and the comparison with a
node one lattice spacing above the horizon in the areal and isotropic
charts.
"""

from math import log, pi, sqrt

# CODATA / framework constants
hbar = 1.054571817e-34          # J s
c = 2.99792458e8                # m / s
G = 6.67430e-11                 # m^3 / (kg s^2)
kB = 1.380649e-23               # J / K
eV = 1.602176634e-19            # J
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
    xi = s_w * kB * T_H / eps0                    # linear, not a cube root
    rho = 2 * r_s * xi                            # = (ln N / pi) hbar c / eps0
    return xi, rho, r_s, T_H, log(N)


def node_one_spacing_up(M, eps0=m0c2):
    """Node one spacing above the horizon (r = r_s + l, areal radius):
    proper height rho_1 [m],
    local temperature k_B T_loc [J], lower bound on the Boltzmann exponent
    of one winding (gap >= eps0), the same bound in crystal time, and the
    gap below which the first winding would thaw there [J]."""
    _, _, r_s, T_H, lnN = skin(M, eps0)
    rho1 = 2 * sqrt(r_s * ell)                    # leading order in l/r_s
    kT_loc = hbar * c / (2 * pi * rho1)
    exponent = eps0 / kT_loc                      # = 4 pi sqrt(r_s l) eps0/(hbar c)
    xi1 = sqrt(ell / r_s)
    exponent_crystal = xi1 * eps0 / (kB * T_H)    # the same bound in crystal time
    gap_thaw = 2 * lnN * kT_loc
    return rho1, kT_loc, exponent, exponent_crystal, gap_thaw


if __name__ == "__main__":
    print(f"epsilon_0 = m0 c^2 = {m0c2/1.602176634e-13:.1f} MeV")
    print(f"{'M/M_sun':>10} {'r_s':>10} {'T_H [K]':>10} {'ln N':>7} "
          f"{'xi_*':>10} {'skin rho_* ':>12}")
    fm = 1e-15
    for Ms in (1, 30, 4.0e6, 6.5e9):
        M = Ms * M_sun
        xi, rho, r_s, T_H, lnN = skin(M)
        print(f"{Ms:>10.3g} {r_s:>10.3e} {T_H:>10.2e} {lnN:>7.1f} "
              f"{xi:>10.2e} {rho/fm:>9.2f} fm")
    # sensitivity: the depth is linear in the entropy factor and in 1/epsilon_0
    xi1, rho1, *_ = skin(M_sun, entropy_factor=False)
    xi2, rho2, *_ = skin(M_sun, eps0=10 * m0c2)
    _, rho0, *_ = skin(M_sun)
    print(f"\nSensitivity at 1 M_sun: rho_* = {rho0/fm:.1f} fm (fiducial); "
          f"{rho1/fm:.2f} fm without the entropy factor; "
          f"{rho2/fm:.1f} fm with epsilon_0 x 10.")
    print("Closed form check: (ln N/pi)(hbar c/eps0) =",
          f"{log((8*pi/log(2))**0.5*M_sun/M_P)/pi*hbar*c/m0c2/fm:.1f} fm")

    print("\nA node one spacing above the horizon (r = r_s + l) against the skin:")
    print(f"{'M/M_sun':>10} {'rho_1':>12} {'rho_1/rho_*':>12} {'k_B T_loc':>12} "
          f"{'T_loc [K]':>10} {'eps0/kT_loc':>12} {'(crystal)':>12} {'2 ln N':>7} "
          f"{'thaw gap':>12}")
    for Ms in (1, 30, 4.0e6, 6.5e9):
        M = Ms * M_sun
        _, rho, _, _, lnN = skin(M)
        r1, kT, E, Ec, gth = node_one_spacing_up(M)
        print(f"{Ms:>10.3g} {r1:>10.3e} m {r1/rho:>12.3e} {kT/eV:>9.3e} eV "
              f"{kT/kB:>10.3g} {E:>12.4e} {Ec:>12.4e} {2*lnN:>7.1f} "
              f"{gth/eV:>9.3e} eV")
    r_s = 2 * G * M_sun / c**2
    _, rho_sun, *_ = skin(M_sun)
    print(f"Every node higher up is colder still.  The skin lies within "
          f"rho_*^2/(4 r_s) = {rho_sun**2/(4*r_s)/ell:.2e} l of the horizon "
          f"in areal height at 1 M_sun.")

    # isotropic chart: node n at proper height 4 n l, static spacing 4 l
    lnN = skin(M_sun)[4]
    e1 = 8 * pi * ell * m0c2 / (hbar * c)          # exponent per n with gap eps0
    print("\nIsotropic chart (node n at proper height 4 n l), 1 M_sun:")
    print(f"  gap eps0:      exponent = {e1:.2f} n; nodes inside the skin: "
          f"n <= {int(2 * lnN // e1)}")
    print(f"  gap 4^3 eps0:  exponent = {64 * e1:.0f} n; margin at n = 1: "
          f"{64 * e1 / (2 * lnN):.2f}; nodes within "
          f"{2 * lnN / (64 * e1):.3f} of a spacing lie inside the skin")
