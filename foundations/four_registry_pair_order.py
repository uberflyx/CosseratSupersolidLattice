"""Four-registry pair order on the D4 lattice: geometry, finite contact and Hartree states.

Companion to the monograph chapter "Internal order and sublattices". The model is
conditional: each node carries a rigid pair of two constituents at +/- a n from its
centre, neighbouring nodes interact through four equal central Morse contacts
between their constituents, and the pair rotor has rotational energy B l(l+1).
Nothing here fixes the vacuum's constituent forces; the script checks the stated
model's numbers.

It computes, in order:

1. The registry lattice of the tetrahedral texture. Translations n in D4 change
   the three carrier signs by the charge c(n) = (n_x+n_y, n_x+n_z, n_y+n_z) mod 2.
   The kernel is H D4 = 2 D4*, of index 4, and the 24 shortest registry-neutral
   vectors (length sqrt2 ell) each split into two orthogonal parent roots in
   exactly three ways, one for each nonzero charge class.
2. The stiffness calibration. A node of mass m0 on the 24-bond shell has transverse
   sound speed c^2 = (1+2r) k_n ell^2 / (2 m0) once the tangential contact
   (r = k_t/k_n = 1/(3 pi - 5)) is engaged, so k_n ell^2 / (m0 c^2) = 2/(1+2r).
   With two half-masses m0/2 and hbar = m0 c ell, the small-pair crossing budget
   Lambda_pair (a/ell)^6 = 0.07904488148 extrapolates to a/ell = 0.78249.
3. The quadrupole coefficient C22 of the unexpanded contact against its leading
   small-pair value -284 (a/ell)^4 / 27, in units of k_n ell^2.
4. Unrestricted even-rotor Hartree states in the four-registry cell (fixed centres,
   real even spherical harmonics through degree L at every site), with the
   constrained four-site wavefunction Hessian as a local stability test.

Units: ell = hbar = m0 = c = 1. Energies are per node unless stated. The Hessian
is that of the total four-site energy. Runtime is about a minute.
"""
import itertools
import os

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh
from scipy.special import eval_legendre, sph_harm_y

BETA = 7/3                      # Morse inverse range times ell (anharmonicity -7)
R_ROLL = 1/(3*np.pi - 5)        # k_t/k_n at the D4 rolling point
KN = 2/(1 + 2*R_ROLL)           # k_n ell^2 / (m0 c^2)
PAIRS = [(0, 1), (0, 2), (1, 2)]
SIGNS = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])
DIRECTORS = SIGNS[:, [0, 2, 1]]/np.sqrt(3)   # tetrahedral axes of the four registries


# ---------------------------------------------------------------- geometry

def roots():
    """The 24 shortest D4 vectors in integer coordinates (lengths sqrt2)."""
    out = []
    for i, j in itertools.combinations(range(4), 2):
        for si, sj in itertools.product((-1, 1), repeat=2):
            v = np.zeros(4, dtype=int)
            v[i], v[j] = si, sj
            out.append(v)
    return np.array(out)


def charge(n):
    """Registry charge: which of the three carrier signs a translation flips."""
    return tuple(int((n[i] + n[j]) % 2) for i, j in PAIRS)


def registry_geometry():
    rs = roots()
    H = np.array([[1, 1, 0, 0], [1, -1, 0, 0], [0, 0, 1, 1], [0, 0, 1, -1]])
    assert np.array_equal(H.T @ H, 2*np.eye(4)) and round(np.linalg.det(H)) == 4
    classes = {}
    for r in rs:
        classes.setdefault(charge(r), []).append(r)
    assert sorted(len(v) for v in classes.values()) == [8, 8, 8]
    assert (0, 0, 0) not in classes
    neutral = np.array([v for v in itertools.product(range(-2, 3), repeat=4)
                        if sum(x*x for x in v) == 4 and charge(v) == (0, 0, 0)])
    assert len(neutral) == 24
    assert {tuple(v) for v in neutral} == {tuple(H @ r) for r in rs}
    for v in neutral:
        routes = [(a, b) for k, a in enumerate(rs) for b in rs[k:]
                  if np.array_equal(a + b, v)]
        assert len(routes) == 3
        assert all(a @ b == 0 and charge(a) == charge(b) for a, b in routes)
        assert len({charge(a) for a, _ in routes}) == 3
    return {'index': 4, 'root_classes': {k: len(v) for k, v in classes.items()},
            'neutral_vectors': len(neutral), 'routes_each': 3}


# ---------------------------------------------------------------- contact

def angular_grid(n):
    """Gauss-Legendre in cos(theta) times 2n uniform azimuths; weights sum to one."""
    z, w = leggauss(n)
    phi = np.arange(2*n)*np.pi/n
    zz, pp = np.repeat(z, 2*n), np.tile(phi, n)
    s = np.sqrt(1 - zz*zz)
    pts = np.array([s*np.cos(pp), s*np.sin(pp), zz]).T
    return pts, np.repeat(w, 2*n)/(4*n)


def morse(dist):
    """Morse contact with W(ell) = -1/(2 beta^2), W'(ell) = 0 and W''(ell) = 1."""
    e = np.exp(-BETA*(dist - 1))
    return (e*e - 2*e)/(2*BETA**2)


def bond_kernels(a, n):
    """Constituent-averaged contact kernels W(|R + a(n' - n)|), summed by registry charge.

    Even angular probabilities make the four constituent sign choices equal, so one
    displacement a(n' - n) represents all four contacts of a bond.
    """
    pts, w = angular_grid(n)
    delta = pts[None, :, :] - pts[:, None, :]
    kernels = {}
    for root in roots():
        R = root/np.sqrt(2)
        dist = np.sqrt(np.sum((R[:3] + a*delta)**2, axis=2) + R[3]**2)
        kernels[charge(root)] = kernels.get(charge(root), 0) + morse(dist)
    return pts, w, kernels


def neighbour_registry(c):
    """For a bond of charge c, the registry index reached from each of the four sites."""
    flip = np.array([(-1)**x for x in c])
    return [int(np.flatnonzero(np.all(SIGNS == s*flip, axis=1))[0]) for s in SIGNS]


def c22(a, n=32):
    """Coefficient of s^2 in the tetrahedral contact energy, in units of k_n ell^2.

    Each site's orientation probability is written in Legendre moments about its own
    director; C22 couples the P2 moments of neighbouring registries.
    """
    pts, w, kernels = bond_kernels(a, n)
    p2 = np.array([5*eval_legendre(2, pts @ d)*w for d in DIRECTORS])
    total = 0.
    for c, K in kernels.items():
        for site, other in enumerate(neighbour_registry(c)):
            total += p2[site] @ K @ p2[other]/8      # half bonds, four sites
    return total


# ---------------------------------------------------------------- Hartree

def real_harmonics(pts, L):
    """Real even spherical harmonics, normalised to one under the unit sphere measure."""
    phi = np.arctan2(pts[:, 1], pts[:, 0])
    th = np.arccos(np.clip(pts[:, 2], -1, 1))
    Y, ll = [], []
    for l in range(0, L + 1, 2):
        for m in range(-l, l + 1):
            y = sph_harm_y(l, abs(m), th, phi)
            if m == 0:
                y = y.real
            elif m > 0:
                y = np.sqrt(2)*(-1)**m*y.real
            else:
                y = np.sqrt(2)*(-1)**m*y.imag
            Y.append(np.sqrt(4*np.pi)*y)
            ll.append(l*(l + 1))
    return np.array(Y).T, np.array(ll, dtype=float)


def hartree(a, n, L, eta, seed, hessian=False, maxiter=600):
    """Self-consistent four-registry product state at pair size a and inertia multiplier eta.

    B = m0 c^2 / (2 eta (a/ell)^2): rotor energy of two half-masses m0/2 at +/- a,
    with the inertia multiplied by eta. Spin-orbit is held at lambda = 2B, where the
    spin-one ground manifold reduces exactly to this scalar rotor.
    """
    pts, w, kernels = bond_kernels(a, n)
    Y, T = real_harmonics(pts, L)
    assert np.max(abs(Y.T @ (w[:, None]*Y) - np.eye(len(T)))) < 1e-12
    B = 1/(2*a*a*eta)
    target = {c: neighbour_registry(c) for c in kernels}

    def fields(g):
        p = w[None, :]*(g @ Y.T)**2
        v = np.zeros_like(p)
        for c, K in kernels.items():
            v += p[target[c]] @ K.T
        return p, v, float(np.mean((g*g) @ (B*T)) + KN*np.sum(p*v)/8)

    g = np.zeros((4, len(T)))
    g[:, 0] = 1
    if seed in ('tetra', 'aligned'):
        for i in range(4):
            d = DIRECTORS[i] if seed == 'tetra' else np.array([0, 0, 1.])
            V = -20*eval_legendre(2, pts @ d)
            g[i] = eigh(np.diag(T) + Y.T @ ((w*V)[:, None]*Y), subset_by_index=[0, 0])[1][:, 0]
            g[i] *= np.sign(g[i, 0])
    for _ in range(maxiter):
        p, v, energy = fields(g)
        new = np.array([eigh(np.diag(B*T) + KN*Y.T @ ((w*v[i])[:, None]*Y),
                             subset_by_index=[0, 0])[1][:, 0] for i in range(4)])
        new *= np.sign(np.sum(new*g, axis=1))[:, None]
        if np.max(abs(new - g)) < 2e-10:
            break
        step = .5                                   # damped map with backtracking
        while True:
            trial = (1 - step)*g + step*new
            trial /= np.linalg.norm(trial, axis=1)[:, None]
            if fields(trial)[2] <= energy + 2e-14 or step < 1e-5:
                break
            step *= .5
        g = trial
    p, v, energy = fields(g)
    Q = np.einsum('an,ni,nj->aij', p, pts, pts) - np.eye(3)/3
    staggered = max(abs(Q[i, x, y]) for i in range(4) for x in range(3) for y in range(3) if x != y)
    p0 = np.tile(w, (4, 1))
    v0 = sum(p0[target[c]] @ K.T for c, K in kernels.items())
    out = {'energy': energy, 'uniform_sphere_energy': float(KN*np.sum(p0*v0)/8),
           'staggered_quadrupole': float(staggered), 's': float(3*staggered)}
    if hessian:
        ng = len(T) - 1
        hh = np.zeros((4*ng, 4*ng))
        tangent = []
        for i in range(4):
            eig, vec = eigh(np.diag(B*T) + KN*Y.T @ ((w*v[i])[:, None]*Y))
            hh[i*ng:(i + 1)*ng, i*ng:(i + 1)*ng] = np.diag(2*(eig[1:] - eig[0]))
            tangent.append(2*(w*(Y @ g[i]))[:, None]*(Y @ vec[:, 1:]))
        for c, K in kernels.items():
            for i, j in enumerate(target[c]):
                if i < j:
                    block = KN*tangent[i].T @ K @ tangent[j]
                    hh[i*ng:(i + 1)*ng, j*ng:(j + 1)*ng] = block
                    hh[j*ng:(j + 1)*ng, i*ng:(i + 1)*ng] = block.T
        out['hessian_min'] = float(np.linalg.eigvalsh(hh)[0])
    return out


if __name__ == '__main__':
    print('Registry lattice:', registry_geometry())

    budget = 0.07904488148                 # Lambda_pair (a/ell)^6 at the small-pair crossing
    lam = KN/4                             # Lambda_pair with mu_red = m0/4 and hbar = m0 c ell
    a_ext = (budget/lam)**(1/6)
    print(f'\nk_n ell^2/(m0 c^2) = 2/(1+2r) = {KN:.8f}; Lambda_pair = {lam:.8f}; '
          f'extrapolated a/ell = {a_ext:.8f}')

    print('\n a/ell    leading C22    unexpanded C22    ratio')
    for a in (0.3, 0.5, 0.73, 0.7825, 0.83):
        lead = -284*a**4/27
        full = c22(a)
        print(f'{a:6.4f} {lead:14.6f} {full:16.6f} {full/lead:10.5f}')

    print('\nEqual half-masses (eta = 1) at the extrapolated radius:')
    for seed in ('tetra', 'aligned'):
        r = hartree(0.7825, 28, 10, 1., seed, hessian=(seed == 'tetra'))
        print(f'  {seed:8s} E - E_uniform = {r["energy"] - r["uniform_sphere_energy"]:.7f}'
              f'  staggered Q = {r["staggered_quadrupole"]:.1e}'
              + (f'  Hessian min = {r["hessian_min"]:.5f}' if 'hessian_min' in r else ''))

    print('\nHigh-inertia diagnostic (a/ell = 0.3, eta = 600):')
    t = hartree(0.3, 20, 10, 600., 'tetra', hessian=True)
    u = hartree(0.3, 20, 10, 600., 'aligned')
    print(f'  tetrahedral s = {t["s"]:.6f}; below the aligned-start state by '
          f'{u["energy"] - t["energy"]:.8f}; Hessian min = {t["hessian_min"]:.7f}')

    print('\nCompeting order (a/ell = 0.5, eta = 50):')
    t = hartree(0.5, 20, 10, 50., 'tetra')
    u = hartree(0.5, 20, 10, 50., 'aligned')
    print(f'  aligned state below the tetrahedral stationary state by {t["energy"] - u["energy"]:.5f}')
