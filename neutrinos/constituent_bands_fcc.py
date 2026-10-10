"""A first pass at the constituent bands on the FCC slice, for the reading of the
antineutrino as a hole in the filled band of the vacuum's constituent fermions.

Three questions, each answered by a short calculation:

1. Where does a four-component constituent with nearest-neighbour hopping cross,
   and how fast does it move there?  With bond vectors n (12 of them, units
   s = l/sqrt 2 = 1) the odd hopping gives
       d_i(k) = 2 t sin(k_i) (cos k_j + cos k_k),   H = alpha . d(k),
   an isotropic crossing at Gamma with speed 4 t s/hbar, anisotropic crossings
   at the four L points (speeds in the ratio 2:1:1, opposite chirality) and
   nodal lines through X and W.  Setting the Gamma speed to c fixes
   t = hbar c/(4 s) = 24.8 MeV and the bandwidth at 198 MeV, so the group
   speed falls below c by 5e-5 at 1 MeV and 5e-3 at 10 MeV.

2. What a uniform condensate in each fermion bilinear does to a crossing:
   scalar and pseudoscalar gap it; a polar vector shifts it; an axial vector
   (spin density) splits it into two isotropic Weyl cones at +-b, which a mass
   m < b cannot close; a tensor spreads it into a ring.

3. A registry term of Wilson type, r sum (1 - cos k.n), gaps every doubler
   (16 r on the X-W lines, 12 r at L) and leaves Gamma massless at tree level,
   but it removes the chiral symmetry [gamma5, H] = 0, so an on-site scalar
   moves the Gamma mass linearly.  A scalar node texture with zero average,
   m_t cos(Q.r) with Q = X, is (-1)^{n_1} on the sites, since 2Q is a
   reciprocal vector; it couples Gamma to X with the full amplitude m_t and
   gives the Gamma fermion a mass m_t^2/M_X through the gapped doubler.
   Keeping that below m_1 c^2 = 5 meV with M_X = m_0 c^2 bounds m_t at about
   0.6 keV.

All numbers quoted in the monograph are asserted here.
"""

import itertools

import numpy as np

HBARC = 197.3269804          # MeV fm
L_FM = 2.8179403             # lattice spacing l = r_e, fm
ALPHA = 1 / 137.035999177
M0C2 = 0.51099895 / ALPHA    # node rest energy, MeV
M1C2 = 5.03e-9               # lightest neutrino rest energy, MeV
S_FM = L_FM / np.sqrt(2)     # D4 integer unit

I2, Z = np.eye(2), np.zeros((2, 2))
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.diag([1, -1]).astype(complex)
BETA = np.block([[I2, Z], [Z, -I2]]).astype(complex)
ALPHAS = [np.block([[Z, s], [s, Z]]) for s in (SX, SY, SZ)]
GAMMA5 = np.block([[Z, I2], [I2, Z]]).astype(complex)

NN = np.array([v for v in itertools.product((-1, 0, 1), repeat=3)
               if sorted(map(abs, v)) == [0, 1, 1]])


def d_vec(k, t=1.0):
    """Odd nearest-neighbour hopping on FCC: (t/2) sum_n sin(k.n) n."""
    return (t / 2) * np.array([np.sum(np.sin(NN @ k) * NN[:, i]) for i in range(3)])


def wilson(k, r):
    """Registry (Wilson-type) term r sum_n (1 - cos k.n), zero at Gamma."""
    return r * np.sum(1 - np.cos(NN @ k))


def hamiltonian(k, t=1.0, r=0.0, m0=0.0):
    return sum(a * di for a, di in zip(ALPHAS, d_vec(k, t))) + (wilson(k, r) + m0) * BETA


def min_abs_energy(hmat):
    return np.min(np.abs(np.linalg.eigvalsh(hmat)))


def crossings():
    """Question 1: the closed form, the crossings and their speeds."""
    k = np.array([0.4, -1.1, 0.7])
    closed = 2 * np.array([np.sin(k[0]) * (np.cos(k[1]) + np.cos(k[2])),
                           np.sin(k[1]) * (np.cos(k[0]) + np.cos(k[2])),
                           np.sin(k[2]) * (np.cos(k[0]) + np.cos(k[1]))])
    assert np.allclose(d_vec(k), closed)

    def jacobian(z, eps=1e-5):
        return np.array([(d_vec(z + eps * np.eye(3)[j]) - d_vec(z - eps * np.eye(3)[j])) / (2 * eps)
                         for j in range(3)]).T

    gamma, ell, x, w = np.zeros(3), np.full(3, np.pi / 2), np.array([np.pi, 0, 0]), np.array([np.pi, np.pi / 2, 0])
    jg, jl = jacobian(gamma), jacobian(ell)
    assert np.allclose(np.linalg.svd(jg, compute_uv=False), [4, 4, 4]) and np.linalg.det(jg) > 0
    assert np.allclose(np.linalg.svd(jl, compute_uv=False), [4, 2, 2]) and np.linalg.det(jl) < 0
    # nodal line (pi, q, 0) through X and W
    for q in np.linspace(-np.pi, np.pi, 17):
        assert np.linalg.norm(d_vec(np.array([np.pi, q, 0]))) < 1e-12
    print("1. Gamma: isotropic, speed 4 t s/hbar, chirality +; L: speeds 4:2:2, chirality -;"
          " d = 0 along the lines (pi, q, 0) through X and W")
    t_c = HBARC / (4 * S_FM)
    assert abs(t_c - 24.76) < 0.01
    ks = np.linspace(-np.pi, np.pi, 41)
    emax = max(np.linalg.norm(d_vec(np.array(kk))) for kk in itertools.product(ks, repeat=3))
    assert abs(emax - 4.0) < 1e-9
    print(f"   v = c at Gamma needs t = hbar c/(4 s) = {t_c:.2f} MeV; bandwidth 2|d|max = {2 * emax * t_c:.0f} MeV")
    for energy, expect in ((1.0, -5.1e-5), (10.0, -5.1e-3)):
        u = np.array([1.0, 0, 0])
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo + hi) / 2
            if np.linalg.norm(d_vec(mid * u)) * t_c < energy:
                lo = mid
            else:
                hi = mid
        q, eps = (lo + hi) / 2, 1e-6
        v = (np.linalg.norm(d_vec((q + eps) * u)) - np.linalg.norm(d_vec((q - eps) * u))) / (2 * eps) * t_c / (HBARC / S_FM)
        assert abs((v - 1) / expect - 1) < 0.05
        print(f"   group speed along [100] at {energy:4.0f} MeV: v/c - 1 = {v - 1:+.1e}")


def channels():
    """Question 2: uniform condensates at a continuum Dirac crossing."""
    k, b, m = np.array([0.3, -0.2, 0.5]), 0.7, 0.4

    def spec(kk, mm, v):
        return np.sort(np.linalg.eigvalsh(sum(a * q for a, q in zip(ALPHAS, kk)) + mm * BETA + v))

    g3g5 = ALPHAS[2] @ GAMMA5               # beta gamma^3 gamma5 = alpha_3 gamma5
    axial = b * g3g5
    e = spec(k, 0.0, axial)
    pred = np.sort([s1 * np.sqrt(k[0]**2 + k[1]**2 + (k[2] + s2 * b)**2) for s1 in (1, -1) for s2 in (1, -1)])
    assert np.allclose(e, pred)
    e = spec(k, m, axial)
    pred = np.sort([s1 * np.sqrt(k[0]**2 + k[1]**2 + (np.sqrt(k[2]**2 + m**2) + s2 * b)**2) for s1 in (1, -1) for s2 in (1, -1)])
    assert np.allclose(e, pred)
    sigma12 = np.block([[SZ, Z], [Z, -SZ]])   # beta sigma^{12} = beta Sigma_z
    e = spec(k, m, b * sigma12)
    pred = np.sort([s1 * np.sqrt(k[2]**2 + (np.sqrt(k[0]**2 + k[1]**2 + m**2) + s2 * b)**2) for s1 in (1, -1) for s2 in (1, -1)])
    assert np.allclose(e, pred)
    for name, v in (("scalar", b * BETA), ("pseudoscalar", b * BETA @ (1j * GAMMA5))):
        assert np.allclose(np.abs(spec(k, 0.0, v)), np.sqrt(k @ k + b**2))
    print("2. scalar and pseudoscalar: gap sqrt(k^2 + b^2); axial (spin density): two Weyl cones at k_z = -+b,"
          " E^2 = k_perp^2 + (sqrt(k_z^2 + m^2) -+ b)^2, open while b > m; tensor: ring at k_perp^2 = b^2 - m^2")


def registry_and_texture():
    """Question 3: Wilson term, lost chiral symmetry, and the texture bound."""
    r = 0.1
    for z, expect in ((np.zeros(3), 0.0), (np.array([np.pi, 0, 0]), 1.6), (np.array([np.pi, np.pi / 2, 0]), 1.6), (np.full(3, np.pi / 2), 1.2)):
        assert abs(min_abs_energy(hamiltonian(z, 1.0, r)) - expect) < 1e-9
    k = np.array([0.3, 0.2, -0.4])
    assert np.allclose(GAMMA5 @ hamiltonian(k) - hamiltonian(k) @ GAMMA5, 0)
    assert not np.allclose(GAMMA5 @ hamiltonian(k, 1, r) - hamiltonian(k, 1, r) @ GAMMA5, 0)
    assert abs(min_abs_energy(hamiltonian(np.zeros(3), 1.0, r, 0.05)) - 0.05) < 1e-9
    print("3. Wilson r: gaps 16 r on the X-W lines and 12 r at L, Gamma massless; [gamma5, H] = 0 lost;"
          " an on-site scalar m0 gives Gamma the mass m0 (unprotected)")
    # Folded two-site cell for m(r) = m_t (-1)^{n_1}: the sublattices n_1 even
    # and odd see +-m_t, and the Bloch Hamiltonian at k in the folded zone is
    # built from the bonds that stay on, or change, the sublattice.
    q_x = np.array([np.pi, 0, 0])
    same = NN[NN[:, 0] == 0]          # bonds with n_1 = 0 keep the sublattice
    cross = NN[NN[:, 0] != 0]         # bonds with n_1 = +-1 change it

    def folded(k, t, r, mt):
        """Per-bond hopping T_n = (i t/2) alpha.n - r beta, on-site 12 r beta."""
        def hop(bonds, kk):
            return sum(((1j * t / 2) * sum(a * ni for a, ni in zip(ALPHAS, n)) - r * BETA) * np.exp(-1j * kk @ n)
                       for n in bonds)
        onsite = 12 * r * BETA
        h_aa = onsite + hop(same, k) + mt * BETA
        h_bb = onsite + hop(same, k) - mt * BETA
        h_ab = hop(cross, k)
        return np.block([[h_aa, h_ab], [h_ab.conj().T, h_bb]])

    # check: at m_t = 0 the folded spectrum at k is that of H(k) and H(k + X) together
    k = np.array([0.3, -0.2, 0.4])
    ref = np.sort(np.concatenate([np.linalg.eigvalsh(hamiltonian(k, 1.0, r)),
                                  np.linalg.eigvalsh(hamiltonian(k + q_x, 1.0, r))]))
    assert np.allclose(np.sort(np.linalg.eigvalsh(folded(k, 1.0, r, 0.0))), ref)
    for mt in (0.02, 0.04, 0.08):
        gap = min_abs_energy(folded(np.zeros(3), 1.0, r, mt))
        exact = np.sqrt((16 * r)**2 / 4 + mt**2) - 16 * r / 2      # two-level result
        assert abs(gap / exact - 1) < 1e-6 and abs(gap / (mt**2 / (16 * r)) - 1) < 0.03
    mt_max = np.sqrt(M1C2 * M0C2)
    assert abs(mt_max * 1e6 - 593) < 3
    print(f"   staggered scalar texture m_t (-1)^n1: Gamma mass m_t^2/M_X; with M_X = m0 c^2 and"
          f" m_eff <= m1 c^2 the texture is bounded at m_t <= {mt_max * 1e6:.0f} eV")


if __name__ == "__main__":
    crossings()
    channels()
    registry_and_texture()
