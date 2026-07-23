#!/usr/bin/env python3
"""
The D_4 conversion of the dark-matter chapter, computed.

A. THE LIFT QUESTION (numbers for the resolution).  A dislocation must wrap the
   compact circle because a line cannot end; a vacancy carries no topological
   charge, so nothing forces its lift to extend along the compact direction.
   The single missing node is the particle; the 'missing column' (one vacancy
   per compact layer) is just a tri-vacancy cluster, and its mass is computed
   here to show it is a heavier bound state, not the fundamental object.

B. THE D_4 DEFECTON BAND.  Tight binding on the 24-neighbour shell with one
   amplitude t = alpha m0 c^2 (the 24 bonds are a single orbit of the lattice
   symmetry).  Band bottom, exact-isotropy effective mass, bandwidth, ceiling
   group velocity, and the compact (Kaluza-Klein) gap from the 6+6 crossing
   bonds at Bloch phases 0, +-2pi/3.

C. CLUSTER LADDER ON THE D_4 GRAPH.  Exact maximal internal bond counts
   B_max(n) for n <= 6 by exhaustive search, checked against the Turan bound
   (no K5 exists in the D_4 contact graph, proved en route).

D. RELAXATION FROM THE FRAMEWORK'S OWN BOND.  The bond potential's depth is
   the bond quantum eps = m0 c^2/12 and its curvature is the contact stiffness
   k_n = mu ell/(sqrt2 (1-r)) from eq:fcc_moduli; a Morse profile matched to
   (depth, curvature) then fixes the anharmonicity with no free choice except
   the profile family, which is varied as a check.  The vacancy relaxation,
   exactly zero at harmonic order, is recomputed with the matched potentials.

E. THE ELASTIC CHANNEL, RE-PRICED.  The derived relaxation volume replaces the
   guessed one; the Born validity parameter is checked; the cross-section band
   is compared with the Milky Way satellite bound.

F. REMNANT STOPPING.  Geometric contact cross-section of the stair-rod cage
   against the overburden threshold for reaching a deep laboratory.
"""

import numpy as np
from itertools import product, combinations
from scipy.optimize import minimize

ME = 0.510_998_950              # electron mass [MeV]
AINV = 137.035_999_177
ALPHA = 1.0 / AINV
M0 = ME * AINV                  # node mass [MeV] = 70.025
T_HOP = ALPHA * M0              # hop amplitude [MeV] = 0.511
ELL_FM = 2.817_940_3205         # lattice spacing [fm]
HBARC = 197.326_980_4           # [MeV fm]
C_KMS = 299_792.458
NA = 6.022_140_76e23

# ------------------------------------------------------------------ lattice
def d4_points(r2max):
    rng = range(-3, 4)
    return np.array([p for p in product(rng, repeat=4)
                     if sum(p) % 2 == 0 and 0 < np.dot(p, p) <= r2max])


NN = d4_points(2)                                  # the 24-cell shell
assert len(NN) == 24


def band_S(k):
    """Structure factor S(k) = sum over the 24 neighbours of exp(i k.delta)."""
    return np.cos(NN @ k).sum(axis=-1) if k.ndim == 1 else None


def part_B():
    print("=" * 72)
    print("B. THE D_4 DEFECTON BAND")
    print("=" * 72)
    # exact curvature from the 4-design property: sum_delta delta_i delta_j
    M2 = NN.T @ NN                                  # integer units, NN dist sqrt2
    iso = np.allclose(M2, M2[0, 0] * np.eye(4))
    print(f"second-moment tensor isotropic: {iso};  sum d_i d_j = {M2[0,0]} delta_ij")
    # physical: sum |d|^2 per axis = 12 (integer) -> sum (k.d)^2 = 12 k^2 a^2,
    # a = ell/sqrt2  ->  E = E0 + (t/2) * 6 ell^2 k^2  ->  m* = m0/(6 alpha)
    mstar = M0 / (6 * ALPHA)
    print(f"band bottom   E_f - 24 t  ->  m_DM = m0(1 - 24a) = "
          f"{M0*(1-24*ALPHA):.2f} MeV")
    print(f"effective mass m* = m0/(6 alpha) = {mstar:.1f} MeV (exactly isotropic)")
    print(f"EP ratio m*/m_g = 1/[6a(1-24a)] = {1/(6*ALPHA*(1-24*ALPHA)):.2f}")
    # bandwidth: S extremes on random + structured sample (S in [-8, 24] exact)
    rng = np.random.default_rng(2)
    K = rng.uniform(-np.pi, np.pi, size=(400000, 4))
    K = np.vstack([K, np.array([[np.pi/2, -np.pi/2, np.pi/2, -np.pi/2]]),
                   np.zeros((1, 4))])
    S = np.cos(K @ NN.T).sum(axis=1)
    print(f"S range sampled: [{S.min():.2f}, {S.max():.2f}]  ->  bandwidth = "
          f"{(S.max()-S.min()):.0f} t = {(S.max()-S.min())*T_HOP:.2f} MeV")
    sn = np.sin(K @ NN.T)
    gradS = -(sn @ NN)
    vmax = np.linalg.norm(gradS, axis=1).max()
    v_kms = vmax * T_HOP / HBARC * (ELL_FM / np.sqrt(2)) * C_KMS
    print(f"ceiling group velocity ~ {v_kms:,.0f} km/s (= {v_kms/C_KMS:.4f} c)")
    # KK gap from the 6+6 crossing bonds, phases 2*pi*n/3
    gap = 12 * T_HOP * (1 - np.cos(2 * np.pi / 3))
    print(f"compact (KK) partners at +18 t = {gap:.2f} MeV above the zero mode")
    print(f"compactification correction to crossing hops ~ e^(-L4/w) = "
          f"{np.exp(-np.sqrt(6)/0.452):.1e}")
    return mstar


def part_C():
    print("\n" + "=" * 72)
    print("C. CLUSTER LADDER: EXACT B_max(n) ON THE D_4 GRAPH")
    print("=" * 72)
    pts = np.vstack([[0, 0, 0, 0], d4_points(4)])   # 49 sites: origin + 2 shells
    n_sites = len(pts)
    d2 = ((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1)
    adj = d2 == 2
    adj_bits = [int("".join("1" if adj[i, j] else "0"
                for j in range(n_sites))[::-1], 2) for i in range(n_sites)]

    def bonds(subset):
        return sum(adj[i, j] for a, i in enumerate(subset)
                   for j in subset[a + 1:])

    # no K5: try to extend every K4 by brute force
    k4 = [(0, i, j, k) for i, j, k in combinations(range(1, n_sites), 3)
          if adj[0, i] and adj[0, j] and adj[0, k]
          and adj[i, j] and adj[i, k] and adj[j, k]]
    k5 = any(all(adj[m, x] for x in q) for q in k4 for m in range(n_sites)
             if m not in q)
    print(f"K4 cliques through origin: {len(k4)};  extendable to K5: {k5}")

    best = {2: 1, 3: 3, 4: 6}
    # n = 5, 6 exhaustive with connectivity implicit via bond count
    for n in (5, 6):
        turan = n * (n - 1) // 2 - {5: 1, 6: 2}[n]  # max edges, clique <= 4
        bmax, arg = 0, None
        for combo in combinations(range(n_sites), n - 1):
            subset = (0,) + combo                    # anchor at origin, WLOG
            b = bonds(subset)
            if b > bmax:
                bmax, arg = b, subset
                if bmax == turan:
                    break
        best[n] = bmax
        print(f"n = {n}: B_max = {bmax}  (Turan bound {turan})  "
              f"example {[tuple(pts[i]) for i in arg]}")
    eps = M0 / 12
    print(f"\nbond quantum eps = m0 c^2/12 = {eps:.3f} MeV")
    print("ladder (growth-step releases m*eps):",
          ", ".join(f"{m*eps:.2f}" for m in (1, 2, 3, 4)), "MeV")
    print(f"tetrahedron binding 6 eps = {6*eps:.2f} MeV = m0 c^2/2;"
          f"  bound tetra-vacancy weighs 4 m0 - m0/2 = {3.5*M0:.1f} MeV")
    tri = 3 * M0 - 3 * eps
    print(f"the 'compact column' = tri-vacancy ring: mass {tri:.1f} MeV "
          f"against the single node's {M0:.1f} MeV: a cluster, not the particle")
    return eps


def relax_with_control(V, s_eq, remove):
    """Relax a D_4 ball at equilibrated scale s_eq; optionally remove the origin.

    Returns (energy recovered by relaxation, mean first-shell radial move).
    The caller subtracts the remove=False run as a boundary control, which
    isolates the vacancy's own contribution from the finite-cluster surface.
    """
    pts = d4_points(9).astype(float) * s_eq
    if not remove:
        pts = np.vstack([[0.0, 0.0, 0.0, 0.0], pts])
    rad = np.sqrt((pts ** 2).sum(1))
    free = rad <= 2.9 * s_eq
    fixed = pts[~free]
    cut = 2.25 * s_eq                      # between second and third shells

    def E(x):
        p = np.vstack([x.reshape(-1, 4), fixed])
        d = p[:, None, :] - p[None, :, :]
        rr = np.sqrt((d ** 2).sum(-1))[np.triu_indices(len(p), 1)]
        return np.sum(V(rr[rr < cut]))

    x0 = pts[free].ravel()
    e0 = E(x0)
    res = minimize(E, x0, method="L-BFGS-B",
                   options={"maxiter": 4000, "ftol": 1e-14})
    rel = res.x.reshape(-1, 4)
    rn = np.sqrt((rel ** 2).sum(1))
    first = np.abs(rad[free] - s_eq * np.sqrt(2)) < 1e-9
    u1 = (rn[first] - rad[free][first]).mean() if first.any() else 0.0
    return e0 - res.fun, u1


def part_D(eps):
    print("\n" + "=" * 72)
    print("D. RELAXATION FROM THE FRAMEWORK'S OWN BOND (equilibrated, controlled)")
    print("=" * 72)
    from scipy.optimize import minimize_scalar
    r1, r2 = np.sqrt(2.0), 2.0
    NN2 = 24                               # second-shell count in D_4

    results = []
    for lam, tag in ((17.3, "slice dictionary r = 0.305"),
                     (15.5, "D_4 dictionary r = 0.226")):
        aM = np.sqrt(lam / 2) / r1         # Morse: k r1^2 = 2 (a r1)^2 = Lambda
        nL = np.sqrt(lam / 2)              # gen. LJ: k r1^2 = 2 n^2 = Lambda
        for V, pname in ((lambda r: (lambda e: e * e - 2 * e)(np.exp(-aM * (r - r1))), "Morse "),
                         (lambda r: (lambda x: x * x - 2 * x)((r1 / r) ** nL), "genLJ ")):
            bulk = lambda sc: 0.5 * (24 * V(sc * r1) + NN2 * V(sc * r2))
            s_eq = minimize_scalar(bulk, bounds=(0.90, 1.02), method="bounded").x
            dEv, u1v = relax_with_control(V, s_eq, remove=True)
            dE0, u10 = relax_with_control(V, s_eq, remove=False)
            u1 = u1v - u10                 # vacancy-only first-shell move
            dE = dEv - dE0                 # vacancy-only recovered energy [eps]
            Ecoh = -bulk(s_eq)             # formation energy in this model [eps]
            dV4 = 2 * np.pi ** 2 * (s_eq * r1) ** 3 * u1 / 2.0   # / Omega_4
            dV3 = 0.289 * dV4              # Omega_4/(L_4 Omega_3) = 1/(2 sqrt 3)
            results.append(dV3)
            print(f"{tag} {pname}: s* = {s_eq:.4f}  u1 = {100*u1/(s_eq*r1):+.2f}%  "
                  f"delta_relax = {100*dE/Ecoh:.1f}% of E_f  dV3/O3 = {dV3:+.3f}")
    print("Within the framework's strict nearest-neighbour contacts all of the")
    print("above is exactly zero; these bands are the model uncertainty from")
    print("giving the bond a generic second-shell tail at the framework's own")
    print("depth (eps = m0 c^2/12) and curvature (k_n l^2 = 90-100 MeV).")
    return max(abs(v) for v in results)


def part_EF(mstar_d4, eps, dv_max):
    print("\n" + "=" * 72)
    print("E. THE ELASTIC CHANNEL, RE-PRICED")
    print("=" * 72)
    MN = 938.918
    c0, epsd = (0.05, 0.20), (0.03, 0.10)
    for mstar, tag in ((M0 * (1 - 24 * ALPHA), "condensate transport, m* ~ m_g"),
                      (mstar_d4, "crystalline transport, m* = m0/6a")):
        mu = mstar * MN / (mstar + MN)
        for hi, name in ((0, "low "), (1, "high")):
            ceff = c0[hi] * epsd[hi] * (0.0 if hi == 0 else dv_max)
            ceff = max(ceff, c0[0] * epsd[0] * 0.05)   # nominal floor for print
            A2 = ceff * M0 * ELL_FM ** 3
            sig = (64 * np.pi / 45) * (mu * A2 / HBARC ** 2) ** 2 * 1e-26
            born = 2 * mu * (ceff * M0) * ELL_FM ** 2 / HBARC ** 2
            print(f"  {tag:36s} {name}: Ceff = {ceff:.1e}, "
                  f"sigma = {sig:.1e} cm^2, Born parameter {born:.3f}")
    print("  Strict contact model: sigma = 0 identically.")
    print("  Milky Way satellite bound ~ 1e-28 cm^2 (Nadler 2019).")

    print("\n" + "=" * 72)
    print("F. REMNANT STOPPING")
    print("=" * 72)
    sig_geo = 6 * ELL_FM * 1.27 * 1e-26
    thresh = 1.0 / (2.7 * NA * 1.4e5)
    print(f"  stair-rod cage contact ~ {sig_geo:.1e} cm^2 "
          f"({sig_geo/1e-24:.2f} barn, geometric)")
    print(f"  ceiling to reach a deep laboratory: {thresh:.1e} cm^2")
    print(f"  margin: stopped by a factor ~ {sig_geo/thresh:.0e}")


if __name__ == "__main__":
    ms = part_B()
    ep = part_C()
    dv = part_D(ep)
    part_EF(ms, ep, dv)
