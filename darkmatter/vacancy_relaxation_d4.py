#!/usr/bin/env python3
"""
Vacancy relaxation and the D_4 conversion of the dark-matter chapter.

Part A: DISCRETE RELAXATION.  Does the shell around a missing node move inward,
and by how much?  This decides two things at once: the vacancy formation energy
(hence the dark-matter mass) and the vacancy's relaxation volume (hence its
elastic dipole, hence whether it couples to a nucleon's pressure field at all).
Two interaction models are run on the same clusters (NOTE: the harmonic
result below is exact; for the quantitative anharmonic bands use the
equilibrated, control-subtracted protocol in darkmatter/d4_conversion.py,
which supersedes the illustrative LJ run here):
  (i)  harmonic nearest-neighbour springs, which is the model the framework
       actually specifies through the contact stiffnesses k_n and k_t;
  (ii) an anharmonic pair potential with second neighbours, the generic case.

Part B: D_4 versus FCC.  The dark-matter chapter was written when the lattice
was FCC (Z = 12).  The framework now derives D_4 (Z = 24) as the foundation.
Every coordination-counted number in the chapter therefore has to be recomputed:
the bond quantum, the tight-binding band edges, the bandwidth and the band mass.

Part C: consequences for the elastic vacancy-nucleon channel.
"""

import numpy as np
from itertools import product
from scipy.optimize import minimize

ME = 0.510_998_950
AINV = 137.035_999_177
ALPHA = 1.0 / AINV
M0 = ME * AINV                    # node mass energy [MeV]


# ---------------------------------------------------------------- lattices
def lattice_points(dim, rmax):
    """Checkerboard lattice D_n: integer points with even coordinate sum."""
    rng = range(-rmax, rmax + 1)
    pts = [p for p in product(rng, repeat=dim)
           if sum(p) % 2 == 0 and np.dot(p, p) <= rmax**2]
    return np.array(pts, dtype=float)


def shells(pts, nmax=3):
    """Distinct neighbour distances and their multiplicities about the origin."""
    d = np.linalg.norm(pts, axis=1)
    d = d[d > 1e-9]
    vals, counts = np.unique(np.round(d, 6), return_counts=True)
    return list(zip(vals[:nmax], counts[:nmax]))


# ------------------------------------------------------- interaction models
def lj_energy(positions, sigma, eps, cutoff):
    """Lennard-Jones total energy, pair minimum at r = sigma."""
    diff = positions[:, None, :] - positions[None, :, :]
    r = np.linalg.norm(diff, axis=-1)
    iu = np.triu_indices(len(positions), k=1)
    rr = r[iu]
    rr = rr[rr < cutoff]
    s6 = (sigma / rr) ** 6
    return eps * np.sum(s6**2 - 2.0 * s6)


def harmonic_energy(positions, pairs, r0, k):
    """Harmonic springs on a fixed bond list, natural length r0."""
    d = np.linalg.norm(positions[pairs[:, 0]] - positions[pairs[:, 1]], axis=1)
    return 0.5 * k * np.sum((d - r0) ** 2)


def relax_vacancy(dim, rmax, model, r_free, verbose=True):
    """Remove the central node, relax the inner region, report the first-shell move."""
    pts = lattice_points(dim, rmax)
    origin = np.argmin(np.linalg.norm(pts, axis=1))
    pts = np.delete(pts, origin, axis=0)          # the vacancy
    radii = np.linalg.norm(pts, axis=1)
    free = radii <= r_free
    d1 = np.sqrt(2.0)                              # nearest-neighbour distance

    x0 = pts[free].ravel()
    fixed = pts[~free]

    def total(x):
        p = np.vstack([x.reshape(-1, dim), fixed])
        return model(p)

    res = minimize(total, x0, method="L-BFGS-B",
                   options={"maxiter": 4000, "ftol": 1e-14, "gtol": 1e-12})
    relaxed = res.x.reshape(-1, dim)

    r_new = np.linalg.norm(relaxed, axis=1)
    r_old = radii[free]
    first = np.abs(r_old - d1) < 1e-6
    u1 = np.mean(r_new[first] - r_old[first])      # radial move, negative = inward
    return u1, u1 / d1, res.fun


# ------------------------------------------------------------- D_4 band
def band_sum(k, dim):
    """Sum of exp(i k.delta) over the checkerboard nearest neighbours.

    The neighbours are (+-1, +-1, 0, ...) and permutations, so the sum
    factorises into 4 cos(k_i) cos(k_j) over all coordinate pairs i < j.
    """
    c = np.cos(k)
    tot = 0.0
    for i in range(dim):
        for j in range(i + 1, dim):
            tot += 4.0 * c[i] * c[j]
    return tot


def band_properties(dim):
    """Band edges, width and curvature for the checkerboard tight-binding band.

    Energy is E(k) = E_f - t * S(k), with S the neighbour sum above.  Expanding
    S about k = 0 gives S = Z - (Z/dim) * (k a / 2)^2 * ... ; the curvature is
    extracted numerically to avoid an algebra slip.
    """
    Z = 4 * dim * (dim - 1) / 2
    s0 = band_sum(np.zeros(dim), dim)
    # zone scan for the minimum of S (which is the band maximum)
    grid = np.linspace(0, np.pi, 41)
    smin = np.inf
    for combo in product(grid, repeat=dim):
        s = band_sum(np.array(combo), dim)
        smin = min(smin, s)
    # curvature: S ~ s0 - C k^2 for small k, in units where the neighbour
    # offset is (1,1,0,..)/1, i.e. the integer lattice
    h = 1e-4
    kk = np.zeros(dim)
    kk[0] = h
    curv = (s0 - band_sum(kk, dim)) / h**2
    return Z, s0, smin, curv


if __name__ == "__main__":
    print("=" * 74)
    print("PART A.  DOES THE SHELL RELAX INWARD?")
    print("=" * 74)
    for dim, name in ((3, "FCC (D_3)"), (4, "D_4")):
        pts = lattice_points(dim, 4)
        print(f"\n{name}: {len(pts)} sites within r = 4;  shells "
              f"{[(round(v,4), int(c)) for v, c in shells(pts)]}")

        # (i) harmonic nearest-neighbour springs on the surviving bond list
        full = lattice_points(dim, 3)
        origin = np.argmin(np.linalg.norm(full, axis=1))
        keep = np.delete(full, origin, axis=0)
        d = np.linalg.norm(keep[:, None, :] - keep[None, :, :], axis=-1)
        pairs = np.array(np.where(np.triu(np.abs(d - np.sqrt(2)) < 1e-6, 1))).T
        model_h = lambda p: harmonic_energy(p, pairs, np.sqrt(2.0), 1.0)
        u1, frac, _ = relax_vacancy(dim, 3, model_h, r_free=2.9)
        print(f"  harmonic NN springs      : first-shell move {u1:+.3e} "
              f"({100*frac:+.4f} % of the bond)")

        # (ii) anharmonic pair potential, second neighbours included
        model_lj = lambda p: lj_energy(p, sigma=np.sqrt(2.0), eps=1.0, cutoff=2.1)
        u1, frac, _ = relax_vacancy(dim, 3, model_lj, r_free=2.9)
        print(f"  Lennard-Jones, shells 1-2: first-shell move {u1:+.3e} "
              f"({100*frac:+.4f} % of the bond)")

    print("\n" + "=" * 74)
    print("PART B.  D_4 VERSUS FCC: EVERY COORDINATION-COUNTED NUMBER")
    print("=" * 74)
    rows = []
    for dim, name in ((3, "FCC (chapter as written)"), (4, "D_4 (framework foundation)")):
        Z, s0, smin, curv = band_properties(dim)
        eps_bond = M0 / (Z / 2)
        e_bottom = M0 * (1 - Z * ALPHA)
        width = (s0 - smin)
        # E = E_f - t S ;  E - E_bottom = t (S0 - S) = t * curv * k^2
        # hbar^2 k^2 / 2 m* = t curv k^2  ->  m* = hbar^2 / (2 t curv a_int^2)
        # with hbar = m0 c ell and the integer spacing a_int = ell/sqrt(2)*sqrt(2)=ell... 
        # neighbour offset (1,1,0..) has length sqrt(2) = ell  ->  integer unit = ell/sqrt(2)
        # k in integer units means k_phys = k / (ell/sqrt(2))
        m_star = M0 / (2 * ALPHA * curv * 0.5)     # = m0 /(alpha * curv)
        rows.append((name, Z, eps_bond, e_bottom, width, m_star))
        print(f"\n{name}")
        print(f"  nearest neighbours Z              : {int(Z)}")
        print(f"  bond quantum  eps = m_0 c^2/(Z/2) : {eps_bond:6.3f} MeV")
        print(f"  di-vacancy binding (one bond)     : {eps_bond:6.3f} MeV")
        print(f"  band bottom  m_0 c^2 (1 - Z alpha): {e_bottom:6.2f} MeV")
        print(f"  band top     E_f + {-smin:.0f} t")
        print(f"  bandwidth    {width:.0f} t = {width*ALPHA*M0:6.2f} MeV")
        print(f"  band (inertial) mass m*           : {m_star:8.1f} MeV "
              f"= m_0/({M0/m_star:.2f} alpha)")
        print(f"  EP ratio m*/m_g                   : {m_star/e_bottom:6.1f}")

    print("\n" + "=" * 74)
    print("PART C.  WHAT THIS DOES TO THE ELASTIC CHANNEL")
    print("=" * 74)
    print("The elastic dipole of a point defect is set by its RELAXATION volume.")
    print("Part A returns that volume as zero for harmonic nearest-neighbour")
    print("springs, in any dimension, because a bond removed from a lattice whose")
    print("bonds all sit at their natural length exerts no force on anything.")
    print("The coupling to a nucleon's pressure field is therefore zero at the")
    print("order the framework specifies its interactions, and the whole")
    print("vacancy-nucleon elastic channel is an anharmonic effect.")
