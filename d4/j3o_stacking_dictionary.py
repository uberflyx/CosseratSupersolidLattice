"""
Is the three-layer stacking cell the physical seat of J3(O)?

The vacuum-symmetry chapter reaches the continuous E6 as the reduced structure
algebra of the exceptional Jordan algebra J3(O), with F4 the crystal-alone
symmetry and the condensate's compression supplying the enlargement (companion
script condensate_f4_to_e6.py). The one premise that script left open was
physical: does anything in the vacuum actually carry the J3(O) structure, or
is the algebra only laid beside the lattice?

This script tests the identification that closes that premise. The candidate
seat is the crystal's own three-layer stacking cell:

    - the three DIAGONAL entries of the Hermitian 3x3 octonionic matrix are
      the breathing amplitudes of the three stacking layers A, B, C;
    - the three OFF-DIAGONAL octonions are the three inter-layer slip
      registries, one per bond sector (A-B, B-C, C-A), exactly the three
      Shockley-partial types the framework identifies with the colour charges;
    - Hermitian conjugation is slip reversal (the B-to-A registry is the
      conjugate of A-to-B);
    - the stacking rotation A -> B -> C is the cyclic permutation of the
      matrix, which must therefore be a Jordan automorphism;
    - the crystal's Spin(8), the symmetry that fixes the layer frame, must
      see the three registries as its three inequivalent octets 8_v, 8_s,
      8_c, cycled by that same rotation (triality).

Counting alone matches (3 + 3x8 = 27), but counting is cheap. The checks
below test the structure:

  CHECK 1  The layer cycle is a Jordan automorphism. It permutes the three
           diagonal idempotents cyclically and carries the three registry
           slots around with them, so the stacking rotation IS the triality
           automorphism realised on the algebra.
  CHECK 2  The Peirce decomposition relative to one layer's idempotent
           lands exactly on coordinate slots: eigenvalue 1 = that layer's
           own amplitude (dim 1), eigenvalue 1/2 = the two registries that
           touch the layer (dim 16, the SO(10) spinor of matter), and
           eigenvalue 0 = the far registry plus the other two amplitudes
           (dim 10). Choosing the layer a defect sits on IS choosing the
           idempotent, and the matter 16 is the pair of registries the
           defect disturbs.
  CHECK 3  The stabiliser of the full layer frame inside F4 (derivations
           annihilating all three idempotents) has dimension 28 = so(8);
           it leaves each registry slot invariant; and the three slot
           representations are pairwise INEQUIVALENT irreducible octets
           (intertwiner test). The layer cycle conjugates this so(8) to
           itself while permuting the slots: triality, exhibited on the
           physical registries.
  CHECK 4  The cell's volume form (the Jordan determinant) is the UNIQUE
           cubic invariant of the E6 action on the 27, up to scale. The
           nullspace of the E6 action on Sym^3(27*), dimension 3654, is
           computed head on and comes out exactly one-dimensional, with
           the determinant spanning it.
  CHECK 5  The Jordan product is RECONSTRUCTED from that volume form. The
           sharp (adjugate) product defined by T(x # y, z) = 3 N(x, y, z),
           with N the polarised determinant and T the trace form, returns
           the Jordan product through the classical Springer identity
               x # y = x o y - (1/2) tr(x) y - (1/2) tr(y) x
                        + (1/2)( tr(x) tr(y) - T(x, y) ) 1,
           and the trace form itself is the Hessian of the determinant at
           the identity. So once the 27 coordinates, the volume form, and
           the undeformed cell (the identity element) are given, the
           multiplication is not a further postulate: it is recovered by
           polarisation.

Together: the linear seat, the symmetry action, the Peirce split, the
uniqueness of the invariant, and the reconstruction of the product. What the
identification then still owes the physics is only the explicit weight map,
which defect fills which weight of 8_s + 8_c; the frame that map lands in is
fixed here.

Octonion and Jordan machinery is shared with condensate_f4_to_e6.py (same
conventions: coordinates 0-2 diagonal; slot x = entry (2,3), coords 3-10;
slot y = entry (3,1), coords 11-18; slot z = entry (1,2), coords 19-26).
Everything is built from scratch; no E6, F4 or so(8) is imported anywhere.
The heavy step (CHECK 4) is a dense 3654 x 3654 eigenproblem assembled from
sparse generator actions; it runs in well under a minute in plain
numpy/scipy, so no numba JIT is needed here.
"""

import numpy as np
from scipy import sparse

rng = np.random.default_rng(7)
DIM_J = 27

# ---------------------------------------------------------------------------
# 1. Octonions (same table as condensate_f4_to_e6.py, self-checked).
# ---------------------------------------------------------------------------
def _build_octonion_table():
    """Structure constants from the Fano triples (e_a e_b = e_c cyclically)."""
    triples = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6),
               (2, 5, 7), (3, 4, 7), (3, 6, 5)]
    C = np.zeros((8, 8, 8))
    C[0, :, :] = np.eye(8)          # 1 * e_b = e_b
    C[:, 0, :] = np.eye(8)          # e_a * 1 = e_a
    for a in range(1, 8):
        C[a, a, 0] = -1.0           # e_a^2 = -1
    for (a, b, c) in triples:
        for (i, j, k) in ((a, b, c), (b, c, a), (c, a, b)):
            C[i, j, k] = 1.0
            C[j, i, k] = -1.0
    return C

OCT = _build_octonion_table()

def omul(u, v):
    return np.einsum('a,b,abk->k', u, v, OCT)

def oconj(u):
    w = u.copy()
    w[1:] = -w[1:]
    return w

# self-check: composition norm |uv| = |u||v| (the defining octonion property)
for _ in range(20):
    u, v = rng.standard_normal(8), rng.standard_normal(8)
    assert abs(np.linalg.norm(omul(u, v))
               - np.linalg.norm(u) * np.linalg.norm(v)) < 1e-10

# ---------------------------------------------------------------------------
# 2. J3(O): coordinates <-> Hermitian octonionic 3x3 matrices, Jordan product.
# ---------------------------------------------------------------------------
def coords_to_matrix(v):
    a, b, c = v[0], v[1], v[2]
    x, y, z = v[3:11], v[11:19], v[19:27]
    M = np.zeros((3, 3, 8))
    M[0, 0, 0] = a; M[1, 1, 0] = b; M[2, 2, 0] = c
    M[0, 1] = z; M[1, 0] = oconj(z)
    M[1, 2] = x; M[2, 1] = oconj(x)
    M[2, 0] = y; M[0, 2] = oconj(y)
    return M

def matrix_to_coords(M):
    v = np.zeros(DIM_J)
    v[0], v[1], v[2] = M[0, 0, 0], M[1, 1, 0], M[2, 2, 0]
    v[3:11] = M[1, 2]
    v[11:19] = M[2, 0]
    v[19:27] = M[0, 1]
    return v

def matmul_oct(A, B):
    C = np.zeros((3, 3, 8))
    for i in range(3):
        for k in range(3):
            acc = np.zeros(8)
            for j in range(3):
                acc += omul(A[i, j], B[j, k])
            C[i, k] = acc
    return C

def jp(u, v):
    """Jordan product X o Y = (XY + YX)/2, in 27-coordinate form."""
    X, Y = coords_to_matrix(u), coords_to_matrix(v)
    return matrix_to_coords(0.5 * (matmul_oct(X, Y) + matmul_oct(Y, X)))

E = np.eye(DIM_J)                       # basis e_0 .. e_26 of J as columns
ONE = np.zeros(DIM_J); ONE[:3] = 1.0    # identity element diag(1,1,1)
SLOT = {'diag': [0, 1, 2],
        'x': list(range(3, 11)),        # entry (2,3): the B-C registry
        'y': list(range(11, 19)),       # entry (3,1): the C-A registry
        'z': list(range(19, 27))}       # entry (1,2): the A-B registry

def tr(u):
    return u[0] + u[1] + u[2]

def T(u, v):
    """Trace bilinear form T(u,v) = tr(u o v)."""
    return tr(jp(u, v))

def L_op(u):
    """Multiplication operator L_u : v -> u o v as a 27x27 matrix."""
    return np.array([jp(u, e) for e in E]).T

def line(ch='-'):
    print(ch * 74)

print("The three-layer stacking cell as J3(O): five structural checks")
line('=')

# ---------------------------------------------------------------------------
# CHECK 1: the layer cycle A -> B -> C is a Jordan automorphism that cycles
#          the idempotents and carries the registry slots with them.
# ---------------------------------------------------------------------------
P = np.zeros((3, 3, 8))                 # cyclic permutation as octonion matrix
P[0, 2, 0] = P[1, 0, 0] = P[2, 1, 0] = 1.0   # row i <- old row i-1
PT = np.transpose(P, (1, 0, 2))

def sigma(v):
    """The stacking rotation on the cell: M -> P M P^T."""
    return matrix_to_coords(matmul_oct(matmul_oct(P, coords_to_matrix(v)), PT))

SIG = np.array([sigma(e) for e in E]).T    # sigma as a 27x27 matrix

auto_ok = all(np.allclose(sigma(jp(u, v)), jp(sigma(u), sigma(v)))
              for u, v in (rng.standard_normal((2, DIM_J)) for _ in range(30)))
idem_cycle = (np.allclose(SIG @ E[:, 0], E[:, 1])
              and np.allclose(SIG @ E[:, 1], E[:, 2])
              and np.allclose(SIG @ E[:, 2], E[:, 0]))

def slot_image(name):
    img = SIG[:, SLOT[name]]
    for tgt in ('x', 'y', 'z'):
        sub = img[SLOT[tgt], :]
        if np.linalg.matrix_rank(sub, tol=1e-9) == 8 and \
           np.allclose(np.delete(img, SLOT[tgt], axis=0), 0):
            return tgt
    return '?'

cycle = {s: slot_image(s) for s in ('x', 'y', 'z')}
print("CHECK 1  the stacking rotation on the cell")
print(f"  Jordan automorphism (30 random pairs):     {auto_ok}")
print(f"  cycles the layer idempotents E1->E2->E3:   {idem_cycle}")
print(f"  carries the registry slots around:         "
      f"x->{cycle['x']}, y->{cycle['y']}, z->{cycle['z']}"
      f"   (a 3-cycle: {sorted(cycle.values()) == ['x', 'y', 'z'] and len(set(cycle.values())) == 3})")
assert auto_ok and idem_cycle
line()

# ---------------------------------------------------------------------------
# CHECK 2: Peirce decomposition w.r.t. one layer = coordinate slots.
# ---------------------------------------------------------------------------
L1 = L_op(E[:, 0])                       # multiplication by the layer-A idempotent
evals, evecs = np.linalg.eigh(0.5 * (L1 + L1.T))
def peirce_space(lam):
    return evecs[:, np.abs(evals - lam) < 1e-9]

pe1, pe12, pe0 = peirce_space(1.0), peirce_space(0.5), peirce_space(0.0)

def support(cols):
    return sorted(set(np.nonzero(np.abs(cols) > 1e-9)[0]))

print("CHECK 2  Peirce split relative to layer A's idempotent")
print(f"  dims (1, 1/2, 0) = ({pe1.shape[1]}, {pe12.shape[1]}, {pe0.shape[1]})"
      f"   (want 1, 16, 10: {(pe1.shape[1], pe12.shape[1], pe0.shape[1]) == (1, 16, 10)})")
s1, s12, s0 = support(pe1), support(pe12), support(pe0)
print(f"  eigenvalue-1 support   = layer A amplitude:        {s1 == [0]}")
print(f"  eigenvalue-1/2 support = registries y and z only:  "
      f"{s12 == sorted(SLOT['y'] + SLOT['z'])}")
print(f"  eigenvalue-0 support   = far registry x + B, C:    "
      f"{s0 == sorted([1, 2] + SLOT['x'])}")
assert (pe1.shape[1], pe12.shape[1], pe0.shape[1]) == (1, 16, 10)
assert s12 == sorted(SLOT['y'] + SLOT['z'])
line()

# ---------------------------------------------------------------------------
# F4 as the span of commutators of multiplication operators (needed below).
# ---------------------------------------------------------------------------
def independent(mats, tol=1e-8):
    """Greedy SVD-based extraction of a linearly independent subset."""
    kept, rows = [], np.zeros((0, DIM_J * DIM_J))
    for M in mats:
        r = M.ravel()[None, :]
        cand = np.vstack([rows, r])
        if np.linalg.matrix_rank(cand, tol=tol) > rows.shape[0]:
            kept.append(M)
            rows = cand
    return kept

L_all = [L_op(e) for e in E]
comms = [L_all[i] @ L_all[j] - L_all[j] @ L_all[i]
         for i in range(DIM_J) for j in range(i + 1, DIM_J)]
F4 = independent(comms)
assert len(F4) == 52, f"dim F4 = {len(F4)}"

# E6 = F4  (+)  multiplications by TRACELESS elements (26 of them)
traceless = [E[:, 0] - E[:, 1], E[:, 1] - E[:, 2]] + \
            [E[:, k] for k in SLOT['x'] + SLOT['y'] + SLOT['z']]
E6 = F4 + [L_op(u) for u in traceless]
assert len(independent(E6)) == 78, "dim E6 != 78"

# ---------------------------------------------------------------------------
# CHECK 3: the layer-frame stabiliser in F4 is so(8), the slots are its three
#          inequivalent octets, and the layer cycle permutes them: triality.
# ---------------------------------------------------------------------------
# stabiliser = derivations killing all three idempotents
A = np.zeros((3 * DIM_J, 52))
for a, D in enumerate(F4):
    A[:, a] = np.concatenate([D @ E[:, 0], D @ E[:, 1], D @ E[:, 2]])
_, s, Vt = np.linalg.svd(A)
null = Vt[np.sum(s > 1e-9):, :]
STAB = [sum(c * D for c, D in zip(row, F4)) for row in null]
print("CHECK 3  the crystal's Spin(8) on the three registries")
print(f"  dim of layer-frame stabiliser in F4 = {len(STAB)}   "
      f"(want 28 = dim so(8): {len(STAB) == 28})")

def block(D, rows, cols):
    return D[np.ix_(rows, cols)]

inv_ok = all(np.allclose(block(D, SLOT[a], SLOT[b]), 0)
             for D in STAB for a in ('x', 'y', 'z') for b in ('x', 'y', 'z')
             if a != b)
print(f"  each registry slot invariant under it:      {inv_ok}")

def hom_dim(a, b):
    """dim of intertwiners T with T rho_a(D) = rho_b(D) T for all D."""
    rows = []
    for D in STAB:
        Ra, Rb = block(D, SLOT[a], SLOT[a]), block(D, SLOT[b], SLOT[b])
        # vec(T Ra - Rb T) = (Ra^T (x) I - I (x) Rb) vec(T)
        rows.append(np.kron(Ra.T, np.eye(8)) - np.kron(np.eye(8), Rb))
    _, sv, _ = np.linalg.svd(np.vstack(rows))
    return int(np.sum(sv < 1e-8)) + 64 - len(sv) if len(sv) < 64 else int(np.sum(sv < 1e-8))

homs = {(a, b): hom_dim(a, b) for a in ('x', 'y', 'z') for b in ('x', 'y', 'z')}
diag_ok = all(homs[(a, a)] == 1 for a in ('x', 'y', 'z'))
off_ok = all(homs[(a, b)] == 0 for a in ('x', 'y', 'z') for b in ('x', 'y', 'z') if a != b)
print(f"  each slot an irreducible octet (Hom = 1):   {diag_ok}")
print(f"  the three octets pairwise inequivalent:     {off_ok}"
      f"   (the 8_v, 8_s, 8_c of triality)")
norm_ok = all(np.linalg.matrix_rank(
    np.vstack([np.array([Dp.ravel()]),
               np.array([D.ravel() for D in STAB])]), tol=1e-7) == 28
    for Dp in (SIG @ D @ np.linalg.inv(SIG) for D in STAB))
print(f"  layer cycle conjugates so(8) to itself:     {norm_ok}"
      f"   (so it acts as the OUTER triality on the slots)")
assert len(STAB) == 28 and inv_ok and diag_ok and off_ok and norm_ok
line()

# ---------------------------------------------------------------------------
# The determinant (volume form) via the Jordan Newton identity.
# ---------------------------------------------------------------------------
def det(u):
    u2 = jp(u, u)
    u3 = jp(u, u2)
    s1, s2, s3 = tr(u), tr(u2), tr(u3)
    return (s1**3 - 3 * s1 * s2 + 2 * s3) / 6.0

assert abs(det(ONE) - 1.0) < 1e-12
v = np.zeros(DIM_J); v[0], v[1], v[2] = 2.0, 3.0, 5.0
assert abs(det(v) - 30.0) < 1e-9          # diagonal: det = abc

def N3(x, y, z):
    """Symmetric trilinear polarisation of det."""
    return (det(x + y + z) - det(x + y) - det(x + z) - det(y + z)
            + det(x) + det(y) + det(z)) / 6.0

# ---------------------------------------------------------------------------
# CHECK 4: the volume form is the UNIQUE cubic invariant of E6 on the 27.
# ---------------------------------------------------------------------------
print("CHECK 4  uniqueness of the invariant cubic (the volume form)")
mono = [(i, j, k) for i in range(DIM_J) for j in range(i, DIM_J)
        for k in range(j, DIM_J)]
NM = len(mono)
mono_id = {m: n for n, m in enumerate(mono)}
midx = np.array(mono)                      # NM x 3
print(f"  dim Sym^3(27) = {NM}   (want 3654: {NM == 3654})")

M = np.zeros((NM, NM))
for g in E6:
    rows_all, cols_all, vals_all = [], [], []
    for p in range(3):
        keep = np.delete(midx, p, axis=1)              # NM x 2
        for b in range(DIM_J):
            trip = np.sort(np.column_stack(
                [keep, np.full(NM, b)]), axis=1)
            rows = np.fromiter((mono_id[tuple(t)] for t in map(tuple, trip)),
                               dtype=np.int64, count=NM)
            rows_all.append(rows)
            cols_all.append(np.arange(NM))
            vals_all.append(g[midx[:, p], b])
    Ag = sparse.coo_matrix(
        (np.concatenate(vals_all),
         (np.concatenate(rows_all), np.concatenate(cols_all))),
        shape=(NM, NM)).tocsr()
    M += (Ag.T @ Ag).toarray()

evals4 = np.linalg.eigvalsh(M)
nnull = int(np.sum(evals4 < 1e-8 * evals4[-1]))
print(f"  nullspace of the E6 action on Sym^3:  dim = {nnull}   "
      f"(want 1: {nnull == 1})")
print(f"  spectral gap: lambda_2 / lambda_1_nonzero-scale = "
      f"{evals4[1] / evals4[-1]:.3e} vs {evals4[0] / evals4[-1]:.3e}")

# the determinant's coefficient vector spans that nullspace
cdet = np.zeros(NM)
for n, (i, j, k) in enumerate(mono):
    mult = 6 if i < j < k else (1 if i == j == k else 3)
    cdet[n] = mult * N3(E[:, i], E[:, j], E[:, k])
resid = float(np.sqrt(cdet @ M @ cdet) / np.linalg.norm(cdet))
print(f"  determinant lies in the nullspace (residual): {resid:.2e}   "
      f"(want ~0: {resid < 1e-7})")
assert nnull == 1 and resid < 1e-7
line()

# ---------------------------------------------------------------------------
# CHECK 5: the Jordan product is reconstructed from the volume form.
# ---------------------------------------------------------------------------
print("CHECK 5  the product recovered from the volume form")
G = np.array([[T(E[:, i], E[:, j]) for j in range(DIM_J)]
              for i in range(DIM_J)])
Ginv = np.linalg.inv(G)

def sharp(x, y):
    """x # y defined by T(x#y, z) = 3 N3(x, y, z) for all z."""
    rhs = np.array([3.0 * N3(x, y, E[:, k]) for k in range(DIM_J)])
    return Ginv @ rhs

springer_ok = True
for _ in range(25):
    x, y = rng.standard_normal(DIM_J), rng.standard_normal(DIM_J)
    lhs = sharp(x, y)
    rhs = jp(x, y) - 0.5 * tr(x) * y - 0.5 * tr(y) * x \
        + 0.5 * (tr(x) * tr(y) - T(x, y)) * ONE
    springer_ok &= np.allclose(lhs, rhs, atol=1e-8)
print(f"  Springer identity  x#y = x o y - (1/2) tr(x) y - (1/2) tr(y) x")
print(f"                       + (1/2)(tr x tr y - T(x,y)) 1 :   {springer_ok}")

# the trace form itself is the Hessian of det at the identity:
# 6 N3(1, x, y) = tr(x) tr(y) - T(x, y)   =>   T is encoded in det too.
hess_ok = True
for _ in range(25):
    x, y = rng.standard_normal(DIM_J), rng.standard_normal(DIM_J)
    hess_ok &= abs(6 * N3(ONE, x, y) - (tr(x) * tr(y) - T(x, y))) < 1e-8
print(f"  trace form = Hessian of det at the identity:      {hess_ok}")
print("  => given the 27 coordinates, the volume form, and the undeformed")
print("     cell (the identity), the Jordan multiplication follows by")
print("     polarisation. The product is reconstructed, not postulated.")
assert springer_ok and hess_ok
line('=')

print("VERDICT")
print("  The three-layer stacking cell carries J3(O) structurally: three")
print("  layer amplitudes on the diagonal, three octonionic registries on")
print("  the off-diagonal, slip reversal as conjugation, and the stacking")
print("  rotation as the triality automorphism that cycles both the layers")
print("  and the three inequivalent Spin(8) octets living on the registries.")
print("  The cell's volume form is the unique cubic invariant of E6 and")
print("  already contains the Jordan product by polarisation. The premise")
print("  condensate_f4_to_e6.py flagged as open is closed at the level of")
print("  structure; what remains for the physics is the explicit weight map,")
print("  which defect fills which weight of 8_s + 8_c.")
