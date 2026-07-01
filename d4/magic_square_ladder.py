"""
Why the matter symmetry is capped at E6: the ladder of the magic square.

The vacuum is a supersolid, so it carries two order parameters. The crystal is
octonionic; that slot is fixed, and its derivation algebra is always
G2 = Der(O), dimension 14. The condensate is the second slot, and the
Freudenthal-Tits magic square makes the whole exceptional series a function of
that second algebra alone:

    condensate = R  (no phase)        ->  F4   (dim 52)
    condensate = C  (a complex scalar)->  E6   (dim 78)   <- our vacuum
    condensate = H  (an SU(2) triple) ->  E7   (dim 133)
    condensate = O  (octonionic)      ->  E8   (dim 248)

The vacuum's condensate is a single complex scalar: one amplitude and one
phase, the field psi = sqrt(rho_s) exp(i theta). Its division algebra is C, so
the group is E6, and it cannot be E7 or E8. To climb the ladder the vacuum
would need a three-component (quaternionic, SU(2)) condensate for E7, or an
octonionic one as rich as the crystal for E8. It has neither. So E6 is a
ceiling on the matter side, a second and independent reason to the no-order-5
argument that already blocks E8 on the force side.

This script makes that ladder concrete. For each division algebra K' it builds
the Jordan algebra J3(K') from a Cayley-Dickson multiplication, computes the
derivation algebra Der(J3(K')) from scratch as the span of the commutators of
the multiplication operators, and confirms the dimensions

    Der(J3(R)) = so(3)  = 3
    Der(J3(C)) = su(3)  = 8
    Der(J3(H)) = sp(3)  = 21
    Der(J3(O)) = f4     = 52.

It then feeds these into the Tits construction

    dim T(O, K') = dim Der(O) + dim Der(J3(K')) + dim Im(O) * dim J3(K')_0
                 = 14 + dim Der(J3(K')) + 7 * (dim J3(K') - 1)

and recovers 52, 78, 133, 248 = dim F4, E6, E7, E8. No exceptional group is
imported; only the two-factor arithmetic of the division algebras is used.
Companion to condensate_f4_to_e6.py, which handles the F4 -> E6 step in full.
"""

import numpy as np

P1 = 2147483647   # 2**31 - 1, prime, for exact rank over F_p
P2 = 2147483629   # prime

# ---------------------------------------------------------------------------
# Cayley-Dickson division algebras R, C, H, O as R^(2^n).
#   product:  (a,b)(c,d) = (a c - conj(d) b, d a + b conj(c))
#   conjugate: conj(a,b) = (conj(a), -b)
# The construction is verified below (norm multiplicativity, and the associator
# pattern: R,C,H associative; O alternative but non-associative).
# ---------------------------------------------------------------------------
def cd_conj(u):
    n = len(u)
    if n == 1:
        return u.copy()
    h = n // 2
    return np.concatenate([cd_conj(u[:h]), -u[h:]])

def cd_mul(u, v):
    n = len(u)
    if n == 1:
        return np.array([u[0] * v[0]], dtype=float)
    h = n // 2
    a, b, c, d = u[:h], u[h:], v[:h], v[h:]
    left = cd_mul(a, c) - cd_mul(cd_conj(d), b)
    right = cd_mul(d, a) + cd_mul(b, cd_conj(c))
    return np.concatenate([left, right])

def cd_check(dim, name):
    rng = np.random.default_rng(dim)
    e = np.eye(dim)
    # identity element is e0
    for i in range(dim):
        assert np.allclose(cd_mul(e[0], e[i]), e[i]) and np.allclose(cd_mul(e[i], e[0]), e[i])
    # imaginary units square to -1
    for i in range(1, dim):
        assert np.allclose(cd_mul(e[i], e[i]), -e[0])
    # norm multiplicativity |uv| = |u||v|
    for _ in range(200):
        u, v = rng.standard_normal(dim), rng.standard_normal(dim)
        assert np.isclose(np.linalg.norm(cd_mul(u, v)),
                          np.linalg.norm(u) * np.linalg.norm(v)), name
    # associativity holds for R,C,H; fails for O
    assoc = True
    for _ in range(200):
        u, v, w = rng.standard_normal(dim), rng.standard_normal(dim), rng.standard_normal(dim)
        if not np.allclose(cd_mul(cd_mul(u, v), w), cd_mul(u, cd_mul(v, w))):
            assoc = False
            break
    return assoc

# ---------------------------------------------------------------------------
# J3(K'): 3x3 K'-Hermitian matrices, dim = 3 + 3*D (D = dim K').
# Coordinates: [t0,t1,t2 (real diag), off01(D), off12(D), off20(D)].
# Layout    [ t0    o01   conj(o20) ]
#           [ conj(o01)  t1   o12   ]
#           [ o20   conj(o12)  t2   ]
# ---------------------------------------------------------------------------
class Jordan:
    def __init__(self, D):
        self.D = D
        self.dimJ = 3 + 3 * D

    def coords_to_matrix(self, v):
        D = self.D
        M = np.zeros((3, 3, D))
        M[0, 0, 0], M[1, 1, 0], M[2, 2, 0] = v[0], v[1], v[2]
        o01 = v[3:3 + D]; o12 = v[3 + D:3 + 2 * D]; o20 = v[3 + 2 * D:3 + 3 * D]
        M[0, 1] = o01;        M[1, 0] = cd_conj(o01)
        M[1, 2] = o12;        M[2, 1] = cd_conj(o12)
        M[2, 0] = o20;        M[0, 2] = cd_conj(o20)
        return M

    def matrix_to_coords(self, M):
        D = self.D
        v = np.zeros(self.dimJ)
        v[0], v[1], v[2] = M[0, 0, 0], M[1, 1, 0], M[2, 2, 0]
        v[3:3 + D] = M[0, 1]
        v[3 + D:3 + 2 * D] = M[1, 2]
        v[3 + 2 * D:3 + 3 * D] = M[2, 0]
        return v

    def matmul(self, A, B):
        D = self.D
        C = np.zeros((3, 3, D))
        for i in range(3):
            for k in range(3):
                acc = np.zeros(D)
                for j in range(3):
                    acc = acc + cd_mul(A[i, j], B[j, k])
                C[i, k] = acc
        return C

    def product(self, u, w):
        X, Y = self.coords_to_matrix(u), self.coords_to_matrix(w)
        S = 0.5 * (self.matmul(X, Y) + self.matmul(Y, X))
        return self.matrix_to_coords(S)

    def basis(self):
        E = []
        for i in range(self.dimJ):
            v = np.zeros(self.dimJ); v[i] = 1.0
            E.append(v)
        return E

    def L_int(self, u):
        """Multiplication operator L_u : Y -> u o Y, scaled by 2 (integer)."""
        cols = [self.product(u, e) for e in self.basis()]
        M = 2.0 * np.array(cols).T
        Mi = np.rint(M).astype(np.int64)
        assert abs(M - Mi).max() < 1e-9
        return Mi

# ---------------------------------------------------------------------------
# Exact rank over F_p via incremental Gaussian elimination (as in the companion).
# ---------------------------------------------------------------------------
class ModEchelon:
    def __init__(self, p):
        self.p = p; self.rows = []; self.cols = []

    def add(self, vec):
        p = self.p
        v = np.mod(vec.astype(np.int64), p)
        for r, c in zip(self.rows, self.cols):
            if v[c]:
                v = np.mod(v - v[c] * r, p)
        nz = np.nonzero(v)[0]
        if nz.size == 0:
            return False
        c = int(nz[0])
        inv = pow(int(v[c]), p - 2, p)
        self.rows.append(np.mod(v * inv, p)); self.cols.append(c)
        return True

    def rank(self):
        return len(self.rows)

def span_rank(mats, p):
    E = ModEchelon(p)
    for A in mats:
        E.add(A.flatten().astype(np.int64))
    return E.rank()

def svd_rank(mats, tol=1e-8):
    M = np.array([A.flatten().astype(float) for A in mats])
    s = np.linalg.svd(M, compute_uv=False)
    kept = s[s > tol * s[0]]
    return len(kept), (kept[-1] if kept.size else float("nan")), \
        (s[s <= tol * s[0]][0] if (s <= tol * s[0]).any() else 0.0)

def derivation_dim(D):
    """dim Der(J3(K')) from the span of all [L_X, L_Y]."""
    J = Jordan(D)
    Ls = [J.L_int(e) for e in J.basis()]
    comms = []
    for i in range(len(Ls)):
        for j in range(i + 1, len(Ls)):
            comms.append(Ls[i] @ Ls[j] - Ls[j] @ Ls[i])
    d1 = span_rank(comms, P1)
    d2 = span_rank(comms, P2)
    dsvd, keep, drop = svd_rank(comms)
    assert d1 == d2 == dsvd, f"rank disagreement at D={D}: {d1},{d2},{dsvd}"
    return d1, keep, drop

# ---------------------------------------------------------------------------
# Run the ladder.
# ---------------------------------------------------------------------------
ALGEBRAS = [(1, "R", "so(3)", 3, "F4", 52),
            (2, "C", "su(3)", 8, "E6", 78),
            (4, "H", "sp(3)", 21, "E7", 133),
            (8, "O", "f4", 52, "E8", 248)]

CONDENSATE = {1: "no phase (real)",
              2: "complex scalar, one amplitude + one phase  <- our vacuum",
              4: "quaternionic (SU(2) triple)",
              8: "octonionic"}

DIM_G2 = 14      # dim Der(O), the crystal's fixed contribution

def line():
    print("-" * 74)

print(__doc__.strip().splitlines()[0])
line()
print("Cayley-Dickson algebras (built and self-checked):")
for (D, name, *_ ) in ALGEBRAS:
    assoc = cd_check(D, name)
    tag = "associative" if assoc else "non-associative (alternative)"
    print(f"  {name:2s}  dim {D:2d}   {tag}")
line()
print("Der(J3(K')) from the span of [L_X, L_Y], and the Tits total:")
print(f"  crystal slot fixed: dim Der(O) = dim G2 = {DIM_G2}")
print()
header = f"  {'K':2s} {'condensate':<46s} {'DerJ3':>6s} {'expect':>7s} {'Tits':>5s} {'grp':>4s}"
print(header)
all_ok = True
for (D, name, algname, expect_der, grp, expect_dim) in ALGEBRAS:
    der, keep, drop = derivation_dim(D)
    J = Jordan(D)
    tits = DIM_G2 + der + 7 * (J.dimJ - 1)     # 7 = dim Im(O)
    ok = (der == expect_der and tits == expect_dim)
    all_ok &= ok
    mark = "ok" if ok else "MISMATCH"
    print(f"  {name:2s} {CONDENSATE[D]:<46s} {der:>6d} {algname:>7s} {tits:>5d} {grp:>4s}  {mark}")
line()
print("Reading the ladder:")
print("  The crystal fixes the octonionic slot (G2 = 14, always).")
print("  The condensate fixes the second slot, and with it the whole group.")
print("  A complex scalar condensate -> C -> E6. That is the vacuum, and the")
print("  ceiling: a quaternionic condensate would give E7, an octonionic one")
print("  E8, and the vacuum's single complex order parameter is neither.")
print(f"\n  ladder verified end to end: {all_ok}")
