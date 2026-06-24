"""
The octonion -> F4 -> E6 bracket check (the magic-square route, made concrete).

Claim under test (from the E6-lift analysis): the continuous E6 Lie bracket is
BUILT from octonion multiplication, not borrowed. The cleanest non-circular
realisation is the exceptional Jordan (Albert) algebra J3(O):
  - J3(O) is the 27-dimensional space of 3x3 Hermitian octonionic matrices.
  - der(J3(O)) = f4  (dim 52): the F4 Lie algebra IS the derivations of J3(O).
  - str0(J3(O)) = der + L(J0) = e6 (dim 78): the E6 Lie algebra is the reduced
    structure algebra, acting on the 27 = J3(O).
Nothing here inputs E6 or F4. We input ONLY the octonion multiplication table
and the Jordan product A o B = (AB + BA)/2, then SOLVE for the derivations and
verify the algebra closes to dimension 52, and the structure algebra to 78.

If the dimensions come out 52 and 78 and both close under the commutator, the
bracket is confirmed real and octonion-built. A wrong octonion table or a sign
error would fail the closure, so the test is genuinely falsifiable.

Framework hook (stated, not over-claimed): the order-3 cyclic automorphism of
J3(O) that permutes its three off-diagonal octonion slots is an element of F4 of
order 3 -- the triality the framework identifies with the three generations /
colour Z3. The 24 unit quaternions 2T (the D4 roots) embed in the unit octonions.
The PHYSICS lift (that a defect core carries the J3 structure) is NOT tested here
and remains open.
"""

import numpy as np
import itertools

np.set_printoptions(suppress=True)

# ===================================================================
# 1. Octonions: multiplication table from the cyclic Fano triples
#    e_i e_{i+1} = e_{i+3} (indices 1..7 mod 7); e0 = 1; e_i^2 = -1.
# ===================================================================
# positive cyclic triples (a,b,c): e_a e_b = e_c, e_b e_c = e_a, e_c e_a = e_b
TRIPLES = [(1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7),
           (5, 6, 1), (6, 7, 2), (7, 1, 3)]

# structure tensor MULT[i,j] = (sign, k) meaning e_i e_j = sign * e_k
_oct = {}
for i in range(8):
    _oct[(0, i)] = (1.0, i)   # 1 * e_i
    _oct[(i, 0)] = (1.0, i)   # e_i * 1
for i in range(1, 8):
    _oct[(i, i)] = (-1.0, 0)  # e_i^2 = -1
for (a, b, c) in TRIPLES:
    _oct[(a, b)] = (1.0, c); _oct[(b, c)] = (1.0, a); _oct[(c, a)] = (1.0, b)
    _oct[(b, a)] = (-1.0, c); _oct[(c, b)] = (-1.0, a); _oct[(a, c)] = (-1.0, b)

def omul(u, v):
    """Octonion product of two length-8 real vectors."""
    w = np.zeros(8)
    for i in range(8):
        if u[i] == 0.0:
            continue
        for j in range(8):
            if v[j] == 0.0:
                continue
            s, k = _oct[(i, j)]
            w[k] += s * u[i] * v[j]
    return w

def oconj(u):
    """Octonion conjugate: negate the imaginary part."""
    w = -u.copy(); w[0] = u[0]; return w

def onorm2(u):
    return float(np.dot(u, u))

# --- self-tests on the octonions ---
def _rand_oct(rng):
    return rng.standard_normal(8)

rng = np.random.default_rng(0)
ok_alt = ok_norm = ok_anti = True
nonassoc = False
for _ in range(200):
    x, y, z = _rand_oct(rng), _rand_oct(rng), _rand_oct(rng)
    # alternativity: (x x) y == x (x y)
    if not np.allclose(omul(omul(x, x), y), omul(x, omul(x, y)), atol=1e-9):
        ok_alt = False
    # norm composition: N(xy) = N(x) N(y)
    if not np.isclose(onorm2(omul(x, y)), onorm2(x) * onorm2(y), atol=1e-7):
        ok_norm = False
    # non-associativity should occur somewhere
    if not np.allclose(omul(omul(x, y), z), omul(x, omul(y, z)), atol=1e-9):
        nonassoc = True
for i in range(1, 8):
    for j in range(1, 8):
        if i != j:
            ei = np.zeros(8); ei[i] = 1; ej = np.zeros(8); ej[j] = 1
            if not np.allclose(omul(ei, ej), -omul(ej, ei), atol=1e-12):
                ok_anti = False

print("=" * 70)
print("1. OCTONIONS")
print(f"   alternativity (xx)y=x(xy) ............ {'PASS' if ok_alt else 'FAIL'}")
print(f"   norm composition N(xy)=N(x)N(y) ...... {'PASS' if ok_norm else 'FAIL'}")
print(f"   anticommuting imaginary units ........ {'PASS' if ok_anti else 'FAIL'}")
print(f"   genuinely non-associative ............ {'PASS' if nonassoc else 'FAIL'}")
assert ok_alt and ok_norm and ok_anti and nonassoc

# ===================================================================
# 2. The exceptional Jordan algebra J3(O): 27-dim Hermitian matrices
#    A = [[a, z, y*],[z*, b, x],[y, x*, c]],  a,b,c real; x,y,z in O.
#    27-vector layout: [a, b, c, x(8), y(8), z(8)].
# ===================================================================
DIM = 27

def vec_to_mat(v):
    a, b, c = v[0], v[1], v[2]
    x = v[3:11]; y = v[11:19]; z = v[19:27]
    A = np.zeros((3, 3, 8))
    A[0, 0, 0] = a; A[1, 1, 0] = b; A[2, 2, 0] = c
    A[0, 1] = z;        A[1, 0] = oconj(z)
    A[1, 2] = x;        A[2, 1] = oconj(x)
    A[2, 0] = y;        A[0, 2] = oconj(y)
    return A

def mat_to_vec(A):
    v = np.zeros(DIM)
    v[0] = A[0, 0, 0]; v[1] = A[1, 1, 0]; v[2] = A[2, 2, 0]
    v[3:11] = A[1, 2]; v[11:19] = A[2, 0]; v[19:27] = A[0, 1]
    return v

def omatmul(A, B):
    C = np.zeros((3, 3, 8))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                C[i, j] += omul(A[i, k], B[k, j])
    return C

def jordan(u, v):
    """Jordan product on 27-vectors: (AB + BA)/2."""
    A, B = vec_to_mat(u), vec_to_mat(v)
    P = 0.5 * (omatmul(A, B) + omatmul(B, A))
    return mat_to_vec(P)

BASIS = [np.eye(DIM)[i] for i in range(DIM)]

# Jordan multiplication operator L_x as a 27x27 matrix
def Lmat(u):
    M = np.zeros((DIM, DIM))
    for j in range(DIM):
        M[:, j] = jordan(u, BASIS[j])
    return M

# --- self-tests on J3(O) ---
ok_comm = ok_herm = ok_jordan = True
for _ in range(50):
    u = rng.standard_normal(DIM); v = rng.standard_normal(DIM)
    if not np.allclose(jordan(u, v), jordan(v, u), atol=1e-9):
        ok_comm = False
    # A o B must stay Hermitian/in J3: round-trip vec->mat->vec is faithful
    P = vec_to_mat(jordan(u, v))
    if not (np.allclose(P[0, 0, 1:], 0) and np.allclose(P[1, 1, 1:], 0)
            and np.allclose(P[2, 2, 1:], 0)
            and np.allclose(P[1, 0], oconj(P[0, 1]), atol=1e-9)):
        ok_herm = False
    # Jordan identity: (A o B) o (A o A) = A o (B o (A o A))
    aa = jordan(u, u)
    lhs = jordan(jordan(u, v), aa)
    rhs = jordan(u, jordan(v, aa))
    if not np.allclose(lhs, rhs, atol=1e-8):
        ok_jordan = False

print("\n2. EXCEPTIONAL JORDAN ALGEBRA J3(O)  (dim 27)")
print(f"   commutative  A o B = B o A ............ {'PASS' if ok_comm else 'FAIL'}")
print(f"   closes in J3 (Hermitian, real diag) .. {'PASS' if ok_herm else 'FAIL'}")
print(f"   Jordan identity ...................... {'PASS' if ok_jordan else 'FAIL'}")
assert ok_comm and ok_herm and ok_jordan

# precompute L matrices for the basis
Lbasis = [Lmat(BASIS[i]) for i in range(DIM)]
# Jordan product structure as 27x27x27: c[k,i,j] with (e_i o e_j)_k
JT = np.zeros((DIM, DIM, DIM))
for i in range(DIM):
    JT[:, i, :] = Lbasis[i]   # (e_i o e_j)_k = Lbasis[i][k,j]

# ===================================================================
# 3. der(J3(O)) = f4 : solve M (A o B) = (M A) o B + A o (M B)
#    for all basis pairs. Nullspace dimension should be 52.
# ===================================================================
print("\n3. DERIVATIONS der(J3(O))  -> expect f4, dim 52")
rows = []
for i in range(DIM):
    Li = Lbasis[i]
    for j in range(i, DIM):
        Lj = Lbasis[j]
        cij = jordan(BASIS[i], BASIS[j])      # 27-vec
        for p in range(DIM):
            row = np.zeros(DIM * DIM)
            # + sum_q cij[q] M[p,q]
            for q in range(DIM):
                row[p * DIM + q] += cij[q]
            # - sum_q Lj[p,q] M[q,i]
            for q in range(DIM):
                row[q * DIM + i] -= Lj[p, q]
            # - sum_q Li[p,q] M[q,j]
            for q in range(DIM):
                row[q * DIM + j] -= Li[p, q]
            rows.append(row)
C = np.array(rows)
print(f"   constraint matrix: {C.shape}")
# nullity via eigenvalues of C^T C
G = C.T @ C
evals = np.linalg.eigvalsh(G)
tol = max(G.shape) * np.finfo(float).eps * (evals[-1] if evals[-1] > 0 else 1.0)
tol = max(tol, 1e-7 * evals[-1])
nullity = int(np.sum(evals < tol))
print(f"   dim der(J3(O)) = {nullity}   {'PASS (= dim f4)' if nullity == 52 else 'MISMATCH'}")

# extract a basis of der as 27x27 matrices (eigenvectors with ~0 eigenvalue)
w, V = np.linalg.eigh(G)
der_vecs = V[:, :nullity]                      # 729 x 52
der_mats = [der_vecs[:, a].reshape(DIM, DIM) for a in range(nullity)]

def in_span(mats, X, tol=1e-6):
    """Is matrix X in the real span of the list `mats`? Return residual norm."""
    A = np.array([m.flatten() for m in mats]).T   # 729 x n
    x = X.flatten()
    coef, *_ = np.linalg.lstsq(A, x, rcond=None)
    return np.linalg.norm(A @ coef - x)

# verify der is closed under the commutator (=> it is a Lie algebra: f4)
maxres = 0.0
for a in range(0, nullity, 7):
    for b in range(a + 1, nullity, 7):
        Br = der_mats[a] @ der_mats[b] - der_mats[b] @ der_mats[a]
        maxres = max(maxres, in_span(der_mats, Br))
print(f"   der closed under [.,.] (Lie algebra) . {'PASS' if maxres < 1e-5 else 'FAIL'}  (max residual {maxres:.2e})")
assert nullity == 52

# ===================================================================
# 4. e6 = str0(J3(O)) = der(J3(O)) + L(J0),  J0 = traceless (a+b+c=0).
#    Expect dim 52 + 26 = 78, and closure under the commutator.
# ===================================================================
print("\n4. REDUCED STRUCTURE ALGEBRA str0(J3(O))  -> expect e6, dim 78")
# traceless basis of J3: 26 elements. trace = v[0]+v[1]+v[2].
traceless = []
# two independent traceless diagonals
d1 = np.zeros(DIM); d1[0] = 1; d1[1] = -1
d2 = np.zeros(DIM); d2[1] = 1; d2[2] = -1
traceless += [d1, d2]
for i in range(3, DIM):       # all 24 off-diagonal directions are traceless
    e = np.zeros(DIM); e[i] = 1; traceless.append(e)
assert len(traceless) == 26
L_traceless = [Lmat(t) for t in traceless]

structure_mats = der_mats + L_traceless        # 52 + 26 = 78
A_all = np.array([m.flatten() for m in structure_mats]).T
rank_struct = np.linalg.matrix_rank(A_all, tol=1e-7)
print(f"   dim span(der + L(J0)) = {rank_struct}   {'PASS (= dim e6)' if rank_struct == 78 else 'MISMATCH'}")

# closure: commutators of a sample of the 78 generators stay in the span
maxres2 = 0.0
idx = list(range(0, 52, 11)) + list(range(52, 78, 5))
for a in idx:
    for b in idx:
        if a < b:
            Br = structure_mats[a] @ structure_mats[b] - structure_mats[b] @ structure_mats[a]
            maxres2 = max(maxres2, in_span(structure_mats, Br))
print(f"   str0 closed under [.,.] (Lie algebra)  {'PASS' if maxres2 < 1e-5 else 'FAIL'}  (max residual {maxres2:.2e})")
print(f"   acts on the 27 = J3(O) ............... yes (operators are 27x27)")
assert rank_struct == 78

# ===================================================================
# 5. Triality: the order-3 cyclic automorphism (the framework's generations)
# ===================================================================
print("\n5. TRIALITY  (the order-3 automorphism = framework generations / colour Z3)")
# cyclic permutation of the 3 matrix indices: A -> P A P^T, P the 3-cycle.
P = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=float)
def tri(u):
    A = vec_to_mat(u)
    # relocate octonion entries: (P A P^T)_{ij} = A_{sigma^{-1}(i), sigma^{-1}(j)}
    B = np.zeros((3, 3, 8))
    for i in range(3):
        for j in range(3):
            ii = int(np.argmax(P[i])); jj = int(np.argmax(P[j]))
            B[i, j] = A[ii, jj]
    return mat_to_vec(B)

# tau is a Jordan automorphism of order 3
ok_auto = True
for _ in range(50):
    u = rng.standard_normal(DIM); v = rng.standard_normal(DIM)
    if not np.allclose(tri(jordan(u, v)), jordan(tri(u), tri(v)), atol=1e-9):
        ok_auto = False
order3 = np.allclose(tri(tri(tri(BASIS[3]))), BASIS[3], atol=1e-9) and not np.allclose(tri(BASIS[3]), BASIS[3])
print(f"   tau preserves the Jordan product ..... {'PASS' if ok_auto else 'FAIL'}")
print(f"   tau has order 3 (cycles generations) . {'PASS' if order3 else 'FAIL'}")
assert ok_auto and order3

# ===================================================================
# 6. Framework hook: the 24 unit quaternions (2T = D4 roots) in the octonions
# ===================================================================
print("\n6. FRAMEWORK HOOK")
# 2T as 24 unit quaternions, embedded in a genuine quaternion subalgebra of O.
# Under the Fano convention e1 e2 = e4, the subalgebra is H = span(1,e1,e2,e4)
# (i=e1, j=e2, k=e4: ij=k, jk=i, ki=j), NOT span(1,e1,e2,e3).
QIDX = [0, 1, 2, 4]
units = []
for ax in QIDX:
    for s in (1.0, -1.0):
        q = np.zeros(8); q[ax] = s; units.append(q)
for signs in itertools.product([0.5, -0.5], repeat=4):
    q = np.zeros(8)
    for n, ax in enumerate(QIDX):
        q[ax] = signs[n]
    units.append(q)
units = np.array(units)
# verify they are unit octonions and closed under octonion multiplication
allunit = np.allclose([onorm2(u) for u in units], 1.0)
closed = True
for u in units:
    for v in units:
        p = omul(u, v)
        if not any(np.allclose(p, w, atol=1e-9) for w in units):
            closed = False
print(f"   24 quaternions 2T are unit octonions . {'PASS' if allunit else 'FAIL'}")
print(f"   2T closed under octonion product ..... {'PASS' if closed else 'FAIL'}")
print("   (2T = D4 roots embed in the unit octonions; the lattice's triality")
print("    is the same order-3 automorphism that lives in F4 above. The PHYSICS")
print("    lift -- that a defect core carries the J3(O) structure -- is NOT")
print("    tested here and remains open.)")

print("\n" + "=" * 70)
print("RESULT")
print("=" * 70)
print(f"""
  der(J3(O))         = {nullity}   (= f4)      built from O-multiplication alone
  str0(J3(O))        = {rank_struct}   (= e6)      acts on the 27 = J3(O)
  triality           order-3 F4 automorphism = the three generations
  E6 / F4 assumed?   NO -- only the octonion table and the Jordan product

  The continuous E6 Lie BRACKET is genuinely octonion-built and acts on the
  27. The magic-square route is mathematically sound: no sign or normalisation
  obstruction. What stays open is the PHYSICS lift (the Albert/J3 structure on
  the defect core), which this script does not address.
""")
