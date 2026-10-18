#!/usr/bin/env python3
"""
d4_local_order.py - local order around a core: triangle counts, cluster energetics, and what the
Koide lattice sums need of the first shell
=================================================================================================
Companion script to "The Cosserat Supersolid" monograph (M. A. Cox).

Three checks, all in lattice units (nearest-neighbour distance 1):

1. Resonant triangles. For the first-shell wavevector sets of several lattices (and of the
   icosahedral 12-shell, which is not a lattice shell) count the triples k1 + k2 + k3 = 0 with
   all three on the shell. The cubic term of a Landau (Alexander-McTague) free energy pays only
   on such triangles, so 6T/N ranks the candidates near a weak crystallisation.

2. Cluster energetics (Frank's argument) with the Morse contact
       V(r) = D [exp(-2a(r-1)) - 2 exp(-a(r-1))],  a = 7/3, D = 9/98  (V''(1) = 1).
   A centre plus its Z shell nodes, all pairs interacting, relaxed by L-BFGS from the lattice
   shell and from random starts. 3D: cuboctahedron (FCC), anticuboctahedron (HCP), icosahedron.
   4D: the 24-cell against whatever random starts find. The contact range is scanned.

3. The Koide quartic lattice sums (sum a1^2, W4^(V''), C3, C4, with a1 = pq, a2 = p^2(1-q^2)/2,
   a3 = -p^3 q (1-q^2)/2, p = n.b, q = n.xi) on the cuboctahedral shell, on a shell of twelve
   neighbours with no angular order (sphere average), on an orientation-averaged icosahedral
   shell and on a randomly rotated cuboctahedron.

Run: python3 d4_local_order.py   (about a minute)
"""
import itertools
import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

rng = np.random.default_rng(0)


def shell_vectors(dim):
    """Minimal vectors of the root lattice D_dim (D3 = FCC, D4 = 24-cell), length sqrt2."""
    return np.array([v for v in itertools.product([-1, 0, 1], repeat=dim) if sum(abs(x) for x in v) == 2], float)


def icosahedron():
    phi = (1 + 5**0.5) / 2
    v = np.array([[0, s1, s2 * phi] for s1 in (-1, 1) for s2 in (-1, 1)], float)
    v = np.vstack([v, v[:, [2, 0, 1]], v[:, [1, 2, 0]]])
    return v / np.linalg.norm(v, axis=1)[:, None]


def triangles(vectors):
    S = {tuple(np.round(v, 6)) for v in vectors}
    n = 0
    for a, b in itertools.combinations(list(S), 2):
        c = tuple(np.round(-(np.array(a) + np.array(b)), 6))
        if c in S and c != a and c != b:
            n += 1
    return n // 3


# ---------------------------------------------------------------- 1. triangles
print("1. resonant triangles on the first shell")
cands = {
    "Z^4 hypercubic": np.array([v for v in itertools.product([-1, 0, 1], repeat=4) if sum(abs(x) for x in v) == 1], float),
    "A4 (FCC analogue)": None,
    "D4 (24-cell)": shell_vectors(4),
    "3D BCC crystal (its shell = 12 FCC-type vectors)": shell_vectors(3),
    "3D FCC crystal (its shell = 8 BCC-type vectors)": np.array(list(itertools.product([-1, 1], repeat=3)), float),
    "icosahedral 12-shell": icosahedron(),
}
# A4 minimal vectors: e_i - e_j in R^5 restricted to the sum-zero hyperplane (20 vectors)
e = np.eye(5)
cands["A4 (FCC analogue)"] = np.array([e[i] - e[j] for i in range(5) for j in range(5) if i != j])
for name, V in cands.items():
    N = len(V); T = triangles(V)
    print(f"   {name:48s} N = {N:2d}  T = {T:2d}  6T/N = {6*T/N:.2f}")

# ---------------------------------------------------------------- 2. clusters
print("\n2. cluster energetics with the Morse contact (centre + shell, all pairs)")
D = 9 / 98


def make(a):
    def V(r):
        ex = np.exp(-a * (r - 1)); return D * (ex * ex - 2 * ex)

    def energy(x, dim):
        P = np.vstack([np.zeros(dim), x.reshape(-1, dim)])
        d = np.linalg.norm(P[:, None] - P[None], axis=-1)
        return V(d[np.triu_indices(len(P), 1)]).sum()

    def relax(P0, dim):
        r = minimize(energy, P0.ravel(), args=(dim,), method="L-BFGS-B",
                     options=dict(maxiter=20000, ftol=1e-14, gtol=1e-10))
        return r.fun
    return relax


cubo = shell_vectors(3) / np.sqrt(2)
hexa = np.array([[np.cos(t), np.sin(t), 0] for t in np.arange(6) * np.pi / 3])
tri = np.array([[np.cos(t), np.sin(t), 0] for t in np.arange(3) * 2 * np.pi / 3 + np.pi / 6]) / np.sqrt(3)
hcp = np.vstack([hexa, tri + [0, 0, np.sqrt(2 / 3)], tri - [0, 0, np.sqrt(2 / 3)]])
ico = icosahedron()
cell24 = shell_vectors(4) / np.sqrt(2)


def random_shell(Z, dim):
    P = rng.normal(size=(Z, dim)); P /= np.linalg.norm(P, axis=1)[:, None]
    return P * rng.uniform(0.9, 1.1, (Z, 1))


print("   a*ell |  3D: cubo      hcp      ico   best random |  4D: 24-cell  best random  24-cell lowest?")
for a in (1.5, 7 / 3, 3.5, 5.0, 8.0):
    relax = make(a)
    e3 = [relax(P, 3) for P in (cubo, hcp, ico)]
    r3 = min(relax(random_shell(12, 3), 3) for _ in range(15))
    e4 = relax(cell24, 4)
    r4 = min(relax(random_shell(24, 4), 4) for _ in range(25))
    print(f"   {a:5.2f} | {e3[0]:9.4f} {e3[1]:8.4f} {e3[2]:8.4f} {r3:9.4f}    | {e4:10.4f} {r4:11.4f}   {'yes' if e4 <= r4 + 1e-6 else 'NO'}")

# ---------------------------------------------------------------- 3. Koide sums
print("\n3. Koide quartic lattice sums (sum a1^2, W4^(V''), C3, C4)")
b = np.array([-1, 1, 0]) / np.sqrt(2); xi = np.array([1, 1, -2]) / np.sqrt(6)


def koide_sums(N):
    p = N @ b; q = N @ xi
    a1 = p * q; a2 = 0.5 * p**2 * (1 - q**2); a3 = -0.5 * p**3 * q * (1 - q**2)
    S1 = (a1**2).sum()
    return S1, ((2 * a1 * a3 + a2**2).sum()) / S1, (a1**2 * a2).sum() / S1, (a1**4).sum() / (12 * S1)


fmt = lambda v: "  ".join(f"{x:8.4f}" for x in v)
print("   cuboctahedral shell            ", fmt(koide_sums(cubo)), "  exact 5/6, 113/192, 19/480, 41/2880")
n = rng.normal(size=(2_000_000, 3)); n /= np.linalg.norm(n, axis=1)[:, None]
p = n @ b; q = n @ xi; a1 = p * q; a2 = 0.5 * p**2 * (1 - q**2); a3 = -0.5 * p**3 * q * (1 - q**2)
S1 = 12 * (a1**2).mean()
iso = (S1, 12 * (2 * a1 * a3 + a2**2).mean() / S1, 12 * (a1**2 * a2).mean() / S1, 12 * (a1**4).mean() / (12 * S1))
print("   twelve neighbours, no angular order", fmt(iso), "  exact 4/5, 2/7, 1/7, 1/84")
for name, P in (("icosahedral, orientation-averaged ", ico), ("cuboctahedral, randomly rotated   ", cubo)):
    vals = np.array([koide_sums(Rotation.random(random_state=i).apply(P)) for i in range(2000)])
    print(f"   {name}", fmt(vals.mean(0)), " +-", fmt(vals.std(0)))
g3, g4, need = -7.0, 23.5, 0.6272
for name, (S1, W, C3, C4) in (("cuboctahedral", koide_sums(cubo)), ("no angular order", iso)):
    roll = -0.019 if name == "cuboctahedral" else 0.0        # Cosserat rolling tilt, zero on an isotropic shell
    eta = W + g3 * C3 + g4 * C4 + roll
    print(f"   eta4 at (gamma3, gamma4) = (-7, 23.5) on the {name} shell: {eta:+.3f} (requirement {need}); gamma4 needed: {(need - W - g3*C3 - roll)/C4:.1f}")
