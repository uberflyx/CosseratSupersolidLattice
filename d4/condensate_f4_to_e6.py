"""
Does the condensate supply the step from F4 to E6?

The vacuum-symmetry chapter reaches the matter symmetry E6 as the reduced
structure algebra of the exceptional (Albert) Jordan algebra J3(O), with
F4 = Der(J3(O)) as the crystal-alone automorphism symmetry. The open question
flagged in the closing section (Route 3) is the PHYSICS: which degree of freedom
promotes the crystal's F4 to the full E6?

This script tests the claim that the supersolid's SECOND order parameter, the
condensate, is exactly that degree of freedom. The condensate is a complex
scalar (amplitude + phase) that lives in the compression channel: it changes the
local density, i.e. it scales volume, and it carries a U(1) phase.

Algebra of the claim
--------------------
  - F4 = Der(J3(O)), dim 52. The derivations fix the identity element and the
    cubic norm (a volume-like invariant on the algebra). These are the
    crystal-alone symmetries.
  - E6 = str_0(J3(O)) = F4 (+) L(J0), dim 78, where L(J0) is multiplication by a
    TRACELESS Jordan element. The extra 26 directions are 'scalings': the coset
    E6/F4 is the 26 of F4. To grow F4 into E6 you must adjoin a multiplication
    (scaling) operator.
  - The condensate compresses the three ABC stacking layers. Resolved on those
    three layers a compression is a DIAGONAL Jordan element; its traceless part
    is a layer-DIFFERENTIAL compression, an element of L(J0) = E6/F4. Its trace
    part is the uniform breathing L(1) = identity, the central U(1) that
    accompanies E6 (the phase). So the two real dimensions of the complex
    condensate map onto (differential compression -> coset) + (phase -> centre).

The sharp, falsifiable test
---------------------------
Start with F4 (52 generators) plus ONE condensate generator
    C = L_{diag(1, 1, -2)}        (two layers squeezed, one relaxed)
and close under the commutator. Because the 26 of F4 is irreducible, one
generator is enough: a single round of [F4, C] must sweep out the entire
26-dimensional coset, and the algebra must saturate at exactly dim E6 = 78.
The script verifies this. It thereby shows that a single breathing mode of the
condensate is sufficient to grow the crystal's F4 into E6.

Method
------
Everything is built from scratch. No E6 or F4 is imported at any step:
    octonion table  ->  J3(O) and its Jordan product
                    ->  the 27 multiplication operators L_X (27x27 real matrices)
                    ->  F4 as the span of all [L_X, L_Y]
                    ->  E6 as F4 (+) L(J0)
                    ->  the single-generator closure test.
Dimensions are counted by exact rank over two prime fields (no floating-point
ambiguity), and cross-checked by a floating-point SVD whose singular-value gap
is printed so the reader can see the count is unambiguous. The calculation is
small (well under a second), so no numba/JIT is used; the bottleneck is human
verification, not compute.

Companion to the vacuum-symmetry chapter (Route 3 of the open-questions section)
and to affine_e6_gravity_node.py, which handles the finite McKay side.
"""

import numpy as np

# Two primes for exact rank over F_p. Agreement across both (and with the SVD
# cross-check) makes an accidental rank drop through a bad prime vanishingly
# unlikely. Both fit in int64 arithmetic with room for the products below.
P1 = 2147483647   # 2**31 - 1, prime
P2 = 2147483629   # prime

# ---------------------------------------------------------------------------
# 1. The octonions O, built from the Fano-plane multiplication rule.
#    Basis e0 = 1 (real unit) and e1..e7 (imaginary units).
#    The seven Fano lines are the base triple (1,2,4) shifted cyclically mod 7.
# ---------------------------------------------------------------------------
FANO_LINES = [tuple(((k + s - 1) % 7) + 1 for k in (0, 1, 3)) for s in range(7)]
# -> (1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3)

def _build_octonion_table():
    """Return T with T[i][j] = (k, sign) meaning e_i e_j = sign * e_k."""
    T = [[(0, 0) for _ in range(8)] for _ in range(8)]
    for i in range(8):
        T[0][i] = (i, 1)      # 1 * e_i = e_i
        T[i][0] = (i, 1)      # e_i * 1 = e_i
    for i in range(1, 8):
        T[i][i] = (0, -1)     # e_i^2 = -1
    for (a, b, c) in FANO_LINES:
        # cyclic: ab=c, bc=a, ca=b ; anticommuting partners get the minus sign
        T[a][b] = (c, 1);  T[b][a] = (c, -1)
        T[b][c] = (a, 1);  T[c][b] = (a, -1)
        T[c][a] = (b, 1);  T[a][c] = (b, -1)
    return T

OCT = _build_octonion_table()

def omul(u, v):
    """Octonion product of two length-8 real vectors."""
    r = np.zeros(8)
    for i in range(8):
        ui = u[i]
        if ui == 0.0:
            continue
        for j in range(8):
            vj = v[j]
            if vj == 0.0:
                continue
            k, s = OCT[i][j]
            r[k] += s * ui * vj
    return r

def oconj(u):
    """Octonion conjugate: negate the imaginary part."""
    c = -u.copy()
    c[0] = u[0]
    return c

def _octonion_self_check():
    rng = np.random.default_rng(0)
    e = np.eye(8)
    # e_i^2 = -1 for imaginary units
    for i in range(1, 8):
        assert np.allclose(omul(e[i], e[i]), -e[0])
    # norm multiplicativity |uv| = |u||v|
    for _ in range(200):
        u, v = rng.standard_normal(8), rng.standard_normal(8)
        assert np.isclose(np.linalg.norm(omul(u, v)),
                          np.linalg.norm(u) * np.linalg.norm(v))
    # alternative: (uu)v = u(uv)
    for _ in range(200):
        u, v = rng.standard_normal(8), rng.standard_normal(8)
        assert np.allclose(omul(omul(u, u), v), omul(u, omul(u, v)))
    # genuinely non-associative: exhibit one failure
    a, b, c = e[1], e[2], e[3]
    assert not np.allclose(omul(omul(a, b), c), omul(a, omul(b, c)))

_octonion_self_check()

# ---------------------------------------------------------------------------
# 2. The Albert algebra J3(O): 3x3 Hermitian octonionic matrices.
#    Coordinates (27): [a, b, c, x(8), y(8), z(8)] laid out as
#        [ a   z   y~ ]
#        [ z~  b   x  ]
#        [ y   x~  c  ]
#    with ~ the octonion conjugate (Hermitian: M_ji = conj(M_ij)).
# ---------------------------------------------------------------------------
DIM_J = 27

def coords_to_matrix(v):
    """27 real coordinates -> 3x3 array of octonions (each length 8)."""
    a, b, c = v[0], v[1], v[2]
    x, y, z = v[3:11], v[11:19], v[19:27]
    M = np.zeros((3, 3, 8))
    M[0, 0, 0] = a; M[1, 1, 0] = b; M[2, 2, 0] = c
    M[0, 1] = z;        M[1, 0] = oconj(z)
    M[1, 2] = x;        M[2, 1] = oconj(x)
    M[2, 0] = y;        M[0, 2] = oconj(y)
    return M

def matrix_to_coords(M):
    """Hermitian 3x3 octonion matrix -> 27 real coordinates (with checks)."""
    assert abs(M[0, 0, 1:]).max() < 1e-9    # diagonal entries real
    assert abs(M[1, 1, 1:]).max() < 1e-9
    assert abs(M[2, 2, 1:]).max() < 1e-9
    assert np.allclose(M[1, 0], oconj(M[0, 1]))   # Hermitian
    assert np.allclose(M[2, 1], oconj(M[1, 2]))
    assert np.allclose(M[0, 2], oconj(M[2, 0]))
    v = np.zeros(DIM_J)
    v[0], v[1], v[2] = M[0, 0, 0], M[1, 1, 0], M[2, 2, 0]
    v[3:11] = M[1, 2]     # x
    v[11:19] = M[2, 0]    # y
    v[19:27] = M[0, 1]    # z
    return v

def matmul_oct(A, B):
    """Product of two 3x3 octonion matrices (octonion entry arithmetic)."""
    C = np.zeros((3, 3, 8))
    for i in range(3):
        for k in range(3):
            acc = np.zeros(8)
            for j in range(3):
                acc += omul(A[i, j], B[j, k])
            C[i, k] = acc
    return C

def jordan_product(u, v):
    """Jordan product X o Y = (XY + YX)/2 in 27-coordinate form."""
    X, Y = coords_to_matrix(u), coords_to_matrix(v)
    S = 0.5 * (matmul_oct(X, Y) + matmul_oct(Y, X))
    return matrix_to_coords(S)

# Basis of J: 3 diagonal units, then 24 off-diagonal (slot x/y/z) x (octonion 0..7)
def _basis_of_J():
    B = []
    for d in range(3):                       # diag(1,0,0), (0,1,0), (0,0,1)
        v = np.zeros(DIM_J); v[d] = 1.0; B.append(v)
    for slot, base in (("x", 3), ("y", 11), ("z", 19)):
        for k in range(8):
            v = np.zeros(DIM_J); v[base + k] = 1.0; B.append(v)
    return B

JBASIS = _basis_of_J()
assert len(JBASIS) == DIM_J

ONE = np.zeros(DIM_J); ONE[0] = ONE[1] = ONE[2] = 1.0   # Jordan identity diag(1,1,1)

def _jordan_self_check():
    rng = np.random.default_rng(1)
    for _ in range(50):
        u, v = rng.standard_normal(DIM_J), rng.standard_normal(DIM_J)
        assert np.allclose(jordan_product(u, v), jordan_product(v, u))   # commutative
        assert np.allclose(jordan_product(ONE, u), u)                    # 1 o u = u

_jordan_self_check()

# ---------------------------------------------------------------------------
# 3. Multiplication operators L_X : Y -> X o Y, as 27x27 real matrices.
#    Scaled by 2 so every entry is an integer (the 1/2 in the Jordan product
#    is cleared), which lets the rank counts below be done in exact integer
#    arithmetic. Scaling by a nonzero constant does not change any rank.
# ---------------------------------------------------------------------------
def L_operator(u):
    cols = [jordan_product(u, e) for e in JBASIS]     # columns X o E_mu
    return np.array(cols).T                            # 27x27

def L_int(u):
    M = 2.0 * L_operator(u)
    Mi = np.rint(M).astype(np.int64)
    assert abs(M - Mi).max() < 1e-9                    # genuinely half-integer
    return Mi

L_basis_int = [L_int(e) for e in JBASIS]               # 27 operators, integer

# L is injective (L_X = 0  <=>  X o 1 = X = 0), so dim L(J) = 27.
assert sum(np.any(A) for A in L_basis_int) == DIM_J

def bracket(A, B):
    """Operator commutator [A,B] = AB - BA (ordinary 27x27 matrix product)."""
    return A @ B - B @ A

# ---------------------------------------------------------------------------
# Exact-rank bookkeeping over F_p, via incremental Gaussian elimination.
# Each 27x27 integer operator is flattened to a length-729 integer vector.
# ---------------------------------------------------------------------------
class ModEchelon:
    """Maintain a row-echelon basis over F_p; add(v) returns True iff v is new."""
    def __init__(self, p):
        self.p = p
        self.rows = []      # reduced pivot vectors (mod p)
        self.cols = []      # pivot column of each row

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
        inv = pow(int(v[c]), p - 2, p)      # Fermat inverse (p prime)
        v = np.mod(v * inv, p)
        self.rows.append(v)
        self.cols.append(c)
        return True

    def rank(self):
        return len(self.rows)

def flat(A):
    return A.flatten().astype(np.int64)

def span_rank(mats, p):
    E = ModEchelon(p)
    for A in mats:
        E.add(flat(A))
    return E.rank()

def independent_mats(mats, p):
    """Greedily return a subset of `mats` forming a basis of their span (F_p)."""
    E = ModEchelon(p)
    keep = []
    for A in mats:
        if E.add(flat(A)):
            keep.append(A)
    return keep

def svd_rank(mats, tol=1e-8):
    """Floating-point rank plus the singular-value gap, for cross-checking."""
    M = np.array([A.flatten().astype(float) for A in mats])
    s = np.linalg.svd(M, compute_uv=False)
    kept = s[s > tol * s[0]]
    dropped = s[s <= tol * s[0]]
    smallest_kept = kept[-1] if kept.size else float("nan")
    largest_dropped = dropped[0] if dropped.size else 0.0
    return len(kept), smallest_kept, largest_dropped

# ---------------------------------------------------------------------------
# 4. F4 = Der(J3(O)) as the span of all commutators [L_X, L_Y].
# ---------------------------------------------------------------------------
commutators = []
for i in range(DIM_J):
    for j in range(i + 1, DIM_J):
        commutators.append(bracket(L_basis_int[i], L_basis_int[j]))

dim_f4_p1 = span_rank(commutators, P1)
dim_f4_p2 = span_rank(commutators, P2)
dim_f4_svd, f4_keep, f4_drop = svd_rank(commutators)

f4_basis = independent_mats(commutators, P1)          # 52 operator matrices

# ---------------------------------------------------------------------------
# 5. E6 = str_0(J3(O)) = F4 (+) L(J0), J0 = traceless Jordan elements.
#    Traceless basis: the 24 off-diagonal units (already traceless) plus two
#    traceless diagonals h1 = diag(1,-1,0), h2 = diag(0,1,-1).
# ---------------------------------------------------------------------------
h1 = np.zeros(DIM_J); h1[0] = 1.0; h1[1] = -1.0
h2 = np.zeros(DIM_J); h2[1] = 1.0; h2[2] = -1.0
J0_elements = [e for e in JBASIS[3:]] + [h1, h2]       # 24 + 2 = 26 traceless
assert len(J0_elements) == 26
L_J0 = [L_int(e) for e in J0_elements]

e6_spanning = f4_basis + L_J0
dim_e6_p1 = span_rank(e6_spanning, P1)
dim_e6_p2 = span_rank(e6_spanning, P2)
dim_e6_svd, e6_keep, e6_drop = svd_rank(e6_spanning)

e6_basis = independent_mats(e6_spanning, P1)           # 78 operator matrices

# Closure: every [A,B] of E6 basis elements stays inside the 78-dim span.
E_ref = ModEchelon(P1)
for A in e6_basis:
    E_ref.add(flat(A))
closure_leaks = 0
for i in range(len(e6_basis)):
    for j in range(i + 1, len(e6_basis)):
        test = ModEchelon(P1)
        for r, c in zip(E_ref.rows, E_ref.cols):       # seed with the reduced basis
            test.rows.append(r.copy()); test.cols.append(c)
        if test.add(flat(bracket(e6_basis[i], e6_basis[j]))):
            closure_leaks += 1

# ---------------------------------------------------------------------------
# 6. The sharp test: F4 + one condensate generator, closed under one bracket.
#    C = L_{diag(1,1,-2)} is a single differential-compression (breathing) mode.
# ---------------------------------------------------------------------------
Xc = np.zeros(DIM_J); Xc[0] = 1.0; Xc[1] = 1.0; Xc[2] = -2.0   # diag(1,1,-2)
C = L_int(Xc)

# (a) C is a genuine scaling: it lies in the E6/F4 coset, not in F4 itself.
in_f4          = (span_rank(f4_basis + [C], P1) == len(f4_basis))
in_e6          = (span_rank(e6_basis + [C], P1) == len(e6_basis))

# (b) The 26 is F4-irreducible: the F4-MODULE generated by C (repeated action
#     of F4 on C) fills the whole coset. E6 = F4 (+) 26 is a symmetric-pair
#     grading -- [F4, coset] stays in the coset -- so the module lives entirely
#     in the 26 dimensions, and reaching all 26 shows the coset is one
#     irreducible piece: a single breathing mode already knows about all of it.
module = ModEchelon(P1)
module.add(flat(C))
frontier = [C]
while frontier:
    grown = []
    for s in frontier:
        for g in f4_basis:
            b = bracket(g, s)
            if module.add(flat(b)):
                grown.append(b)
    frontier = grown
module_dim = module.rank()

# (c) Full Lie closure of <F4, C> lands on exactly E6 and goes no further,
#     because once the coset is full, [coset, coset] closes back into F4.
gen = list(f4_basis) + [C]
prev = -1
cur = span_rank(gen, P1)
rounds = 0
while cur != prev and rounds < 6:
    basis_now = independent_mats(gen, P1)
    new = [bracket(A, B) for i, A in enumerate(basis_now)
           for B in basis_now[i + 1:]]
    gen = basis_now + new
    prev, cur = cur, span_rank(gen, P1)
    rounds += 1
dim_generated = cur

# (d) The phase: uniform breathing L(1) is the central U(1) that sits OUTSIDE
#     E6. It is the identity operator, commutes with everything, and adjoining
#     it lifts str_0 = E6 (78) to str = E6 x U(1) (79).
L_one = L_int(ONE)
is_identity_map = np.array_equal(L_one, 2 * np.eye(DIM_J, dtype=np.int64))
commutes_with_all = all(not np.any(bracket(L_one, A)) for A in e6_basis)
dim_str = span_rank(e6_basis + [L_one], P1)

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def line():
    print("-" * 72)

print(__doc__.strip().splitlines()[0])
line()
print("Octonions O and Albert algebra J3(O) built and self-checked.")
print(f"  dim J3(O) = {DIM_J}   (3 real diagonal + 3 octonion off-diagonal)")
line()
print("F4 = Der(J3(O)) = span of all [L_X, L_Y]:")
print(f"  dim over F_{P1} : {dim_f4_p1}")
print(f"  dim over F_{P2} : {dim_f4_p2}")
print(f"  dim by SVD       : {dim_f4_svd}   "
      f"(smallest kept {f4_keep:.3g}, largest dropped {f4_drop:.2g})")
print(f"  => F4 recovered, dim {dim_f4_p1} = dim F4 (52).")
line()
print("E6 = str_0(J3(O)) = F4 (+) L(J0), J0 = traceless Jordan elements:")
print(f"  dim over F_{P1} : {dim_e6_p1}")
print(f"  dim over F_{P2} : {dim_e6_p2}")
print(f"  dim by SVD       : {dim_e6_svd}   "
      f"(smallest kept {e6_keep:.3g}, largest dropped {e6_drop:.2g})")
print(f"  closure: {closure_leaks} bracket(s) left the 78-dim span "
      f"(0 => E6 is a Lie algebra).")
print(f"  => E6 recovered, dim {dim_e6_p1} = dim E6 (78); coset E6/F4 has "
      f"dim {dim_e6_p1 - dim_f4_p1} (the 26 of F4).")
line()
print("SHARP TEST -- one condensate generator C = L_diag(1,1,-2):")
print(f"  C lies inside F4 ?           {in_f4}   (want False: C is a scaling)")
print(f"  C lies inside E6 ?           {in_e6}   (want True:  C is in E6/F4)")
print(f"  F4-module generated by C   : {module_dim}"
      f"   (= 26 => the coset is one irreducible piece)")
print(f"  full closure <F4, C>       : {dim_generated}"
      f"   (converged in {rounds} round(s))")
verdict = (in_f4 is False and in_e6 is True and
           module_dim == 26 and dim_generated == 78)
print(f"  => one breathing mode + F4 generates all of E6: {verdict}")
line()
print("THE PHASE -- uniform breathing L(1), the central U(1):")
print(f"  L(1) is the identity map ?   {is_identity_map}")
print(f"  L(1) commutes with all E6 ?  {commutes_with_all}   (central)")
print(f"  dim str = E6 + L(1)        : {dim_str}   (= 79 = dim E6 + 1)")
print(f"  => the phase is the U(1) accompanying E6; str = E6 x U(1).")
line()
print("CONCLUSION")
print("  The crystal alone realises F4 = Der(J3(O)), the symmetries that fix")
print("  the octonionic volume form. Adjoining the condensate's compression")
print("  (a single scaling generator) grows F4 into E6 = str_0(J3(O)), and its")
print("  phase supplies the accompanying central U(1). One complex order")
print("  parameter, two real dimensions, is exactly the (O, C) ingredient the")
print("  Freudenthal-Tits square needs for E6. What stays open is unchanged:")
print("  whether a defect core physically carries the J3(O) structure at all.")
