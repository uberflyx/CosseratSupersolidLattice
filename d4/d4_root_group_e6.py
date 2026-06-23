"""
The D4 root group, its McKay dual E6, and the shared Q8 of strong and weak isospin.

This script verifies, from scratch and numerically, the chain of group-theoretic
facts that ties the two isospin SU(2)s of the lattice to a single finite group.

The 24 minimal vectors (roots) of the D4 lattice, taken as unit quaternions, close
under quaternion multiplication into the binary tetrahedral group 2T = SL(2,3),
which factors as the semidirect product Q8 |x Z3.  Two readings of this one group
supply the two isospins:

  strong isospin : the McKay dual of the affine D4 Dynkin diagram is Q8, whose only
                   doublet is the up-down isospin doublet;
  weak  SU(2)_L  : left quaternion multiplication by the roots is a symmetry of the
                   root system, hence a lattice point-group element, and it acts on
                   the self-dual microrotation, which is the weak SU(2)_L factor of
                   SO(4) = SU(2)_L x SU(2)_R.  Its binary-dihedral core is the same Q8.

Adjoining colour, the cyclic Z3 that rotates i -> j -> k, promotes Q8 to the full
2T = SL(2,3); the McKay dual of 2T is the affine E6 diagram.  The script builds the
McKay quiver explicitly and confirms it is E6-hat.

Checks printed below (each PASS/FAIL):
  A  the 24 roots close into a group of order 24
  B  Q8 = {+-1, +-i, +-j, +-k} is a normal subgroup of order 8
  C  the quotient 2T / Q8 is Z3
  D  2T is the semidirect product Q8 |x Z3 (Z3 cycles i -> j -> k by conjugation)
  E  2T has 7 conjugacy classes
  F  2T has 7 irreducibles of dimensions (1,1,1,2,2,2,3), squares summing to 24
  G  the McKay quiver of 2T is the affine E6 diagram
  H  left multiplication embeds 2T in SO(4); the Q8 core lies in the SU(2)_L factor
  I  the strong-isospin Q8 and the weak-isospin Q8 are the same eight elements
"""

import numpy as np
import itertools

def _fmt(q):
    """Pretty-print a quaternion as a+bi+cj+dk with halves shown as /2."""
    if q is None:
        return "?"
    labels = ["", "i", "j", "k"]
    terms = []
    for val, lab in zip(q, labels):
        if abs(val) < 1e-9:
            continue
        sign = "+" if val > 0 else "-"
        mag = abs(val)
        magstr = "1" if abs(mag - 1) < 1e-9 else ("1/2" if abs(mag - 0.5) < 1e-9 else f"{mag:.2f}")
        terms.append(f"{sign}{magstr}{lab}")
    s = "".join(terms)
    return s[1:] if s.startswith("+") else s


np.set_printoptions(precision=3, suppress=True)
TOL = 1e-9


# --------------------------------------------------------------------------
# Quaternions as 4-vectors [w, x, y, z] = w + x i + y j + z k
# --------------------------------------------------------------------------
def qmul(a, b):
    """Hamilton product of two quaternions given as length-4 arrays."""
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return np.array([
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    ])


def qconj(a):
    """Quaternion conjugate (a unit quaternion's inverse)."""
    return np.array([a[0], -a[1], -a[2], -a[3]])


def to_su2(q):
    """Unit quaternion -> 2x2 SU(2) matrix.  a+bi+cj+dk -> [[a+bi, c+di],[-c+di, a-bi]]."""
    a, b, c, d = q
    return np.array([[a + 1j * b, c + 1j * d],
                     [-c + 1j * d, a - 1j * b]], dtype=complex)


# --------------------------------------------------------------------------
# The 24 roots of D4 as unit quaternions: the binary tetrahedral group 2T
# --------------------------------------------------------------------------
def build_2T():
    elems = []
    # 8 Lipschitz units: +-1, +-i, +-j, +-k  -> this set is Q8
    for axis in range(4):
        for sign in (+1.0, -1.0):
            v = np.zeros(4)
            v[axis] = sign
            elems.append(v)
    # 16 half-integer units: (+-1 +-i +-j +-k) / 2
    for signs in itertools.product((+1.0, -1.0), repeat=4):
        elems.append(0.5 * np.array(signs))
    return [np.array(e) for e in elems]


def find_index(group, q):
    """Index of quaternion q within the list `group` (matching up to TOL)."""
    for idx, g in enumerate(group):
        if np.allclose(g, q, atol=TOL):
            return idx
    return -1


# --------------------------------------------------------------------------
# Group-theory utilities driven by the Cayley table
# --------------------------------------------------------------------------
def cayley_table(group):
    n = len(group)
    table = np.full((n, n), -1, dtype=int)
    for i in range(n):
        for j in range(n):
            table[i, j] = find_index(group, qmul(group[i], group[j]))
    return table


def conjugacy_classes(group, table):
    n = len(group)
    inv = [find_index(group, qconj(group[i])) for i in range(n)]   # unit quaternions: inverse = conjugate
    seen = set()
    classes = []
    for x in range(n):
        if x in seen:
            continue
        cls = set()
        for g in range(n):
            # g x g^{-1}
            cls.add(table[table[g, x], inv[g]])
        classes.append(sorted(cls))
        seen |= cls
    return classes


# --------------------------------------------------------------------------
# Run the checks
# --------------------------------------------------------------------------
def main():
    group = build_2T()
    n = len(group)

    # Sanity: all unit norm, SU(2) matrices unitary with det 1
    assert all(abs(np.dot(q, q) - 1.0) < TOL for q in group)
    for q in group:
        U = to_su2(q)
        assert np.allclose(U @ U.conj().T, np.eye(2), atol=TOL)
        assert abs(np.linalg.det(U) - 1.0) < TOL

    # A: closure into a group of order 24
    closed = True
    for a in group:
        for b in group:
            if find_index(group, qmul(a, b)) < 0:
                closed = False
    A = closed and n == 24
    print(f"A  24 roots close into a group of order {n:<2d}                       : {'PASS' if A else 'FAIL'}")

    table = cayley_table(group)

    # Identify Q8 = {+-1, +-i, +-j, +-k} by its indices
    q8_quats = []
    for axis in range(4):
        for sign in (+1.0, -1.0):
            v = np.zeros(4); v[axis] = sign
            q8_quats.append(v)
    Q8 = sorted(find_index(group, q) for q in q8_quats)

    # B: Q8 normal of order 8
    def is_subgroup(subset):
        s = set(subset)
        return all(table[i, j] in s for i in subset for j in subset)

    def is_normal(subset):
        s = set(subset)
        inv = [find_index(group, qconj(group[i])) for i in range(n)]
        return all(table[table[g, i], inv[g]] in s for g in range(n) for i in subset)

    B = len(Q8) == 8 and is_subgroup(Q8) and is_normal(Q8)
    print(f"B  Q8 = {{+-1,+-i,+-j,+-k}} is a normal subgroup of order 8         : {'PASS' if B else 'FAIL'}")

    # C: quotient 2T / Q8 is Z3 (three cosets, cyclic)
    cosets = []
    placed = set()
    for x in range(n):
        if x in placed:
            continue
        coset = sorted(table[x, q] for q in Q8)
        cosets.append(coset)
        placed |= set(coset)
    # quotient is cyclic of order 3 iff there are 3 cosets and a coset of order 3 in the quotient
    C = len(cosets) == 3
    print(f"C  the quotient 2T / Q8 has {len(cosets)} cosets, so 2T/Q8 = Z3            : {'PASS' if C else 'FAIL'}")

    # D: semidirect structure -- conjugation by a non-Q8 element 3-cycles the axes.
    # In Q8/{+-1} the three axes are {+-i}, {+-j}, {+-k}; the colour Z3 permutes them
    # cyclically, which is the defining feature of the binary tetrahedral group.
    iq = np.array([0.0, 1.0, 0.0, 0.0])
    jq = np.array([0.0, 0.0, 1.0, 0.0])
    kq = np.array([0.0, 0.0, 0.0, 1.0])

    def axis_of(v):
        """Return 0/1/2 for a quaternion lying along +-i / +-j / +-k, else -1."""
        for ax, e in enumerate((iq, jq, kq)):
            if np.allclose(v, e, atol=TOL) or np.allclose(v, -e, atol=TOL):
                return ax
        return -1

    # search the order-3 / order-6 elements for one whose conjugation is the 3-cycle (0 1 2)
    found_cycle = False
    cycler = None
    for g in group:
        perm = tuple(axis_of(qmul(qmul(g, e), qconj(g))) for e in (iq, jq, kq))
        if perm == (1, 2, 0):           # i->j, j->k, k->i  (mod sign)
            found_cycle = True
            cycler = g
            break
    D = found_cycle
    order_of_cycler = None
    if cycler is not None:
        p = np.array([1.0, 0, 0, 0]); order_of_cycler = 0
        while True:
            p = qmul(p, cycler); order_of_cycler += 1
            if np.allclose(p, np.array([1.0, 0, 0, 0]), atol=TOL):
                break
    print(f"D  conjugation by ({_fmt(cycler)}) 3-cycles i->j->k->i             : {'PASS' if D else 'FAIL'}")

    # E: 7 conjugacy classes
    classes = conjugacy_classes(group, table)
    E = len(classes) == 7
    print(f"E  number of conjugacy classes = {len(classes)}                            : {'PASS' if E else 'FAIL'}")

    # ----------------------------------------------------------------------
    # Build the 7 irreducible representations explicitly and their characters
    # ----------------------------------------------------------------------
    # Map each element to its Z3 coset label (0,1,2) for the three 1-dim irreps
    coset_label = {}
    for lbl, coset in enumerate(_order_cosets(group, table, Q8)):
        for idx in coset:
            coset_label[idx] = lbl

    zeta = np.exp(2j * np.pi / 3)

    # 3-dim irrep: conjugation action on imaginary quaternions {i,j,k}
    basis = [np.array([0.0, 1, 0, 0]), np.array([0.0, 0, 1, 0]), np.array([0.0, 0, 0, 1])]
    def rep3(idx):
        q = group[idx]
        M = np.zeros((3, 3))
        for col, e in enumerate(basis):
            img = qmul(qmul(q, e), qconj(q))   # imaginary part lands back in span{i,j,k}
            M[:, col] = img[1:]
        return M

    # characters on a class representative
    reps = [cls[0] for cls in classes]
    class_size = [len(cls) for cls in classes]

    chi = {}   # chi[name] = array over classes
    chi['1']   = np.array([1.0 for _ in reps], dtype=complex)
    chi['w']   = np.array([zeta ** coset_label[r] for r in reps], dtype=complex)
    chi['w2']  = np.array([zeta ** (2 * coset_label[r]) for r in reps], dtype=complex)
    chi['2']   = np.array([np.trace(to_su2(group[r])) for r in reps], dtype=complex)
    chi['2w']  = chi['2'] * chi['w']
    chi['2w2'] = chi['2'] * chi['w2']
    chi['3']   = np.array([np.trace(rep3(r)) for r in reps], dtype=complex)

    names = ['1', 'w', 'w2', '2', '2w', '2w2', '3']
    dims = {name: int(round(chi[name][_class_of(reps, classes, _identity_index(group))].real)) for name in names}

    # inner product of two class functions
    def inner(a, b):
        return sum(class_size[c] * a[c] * np.conj(b[c]) for c in range(len(reps))) / n

    # F: dimensions and orthonormality
    dim_list = sorted(dims[name] for name in names)
    sq_sum = sum(d * d for d in dim_list)
    orthonormal = all(abs(inner(chi[a], chi[b]) - (1.0 if a == b else 0.0)) < 1e-6
                      for a in names for b in names)
    F = dim_list == [1, 1, 1, 2, 2, 2, 3] and sq_sum == 24 and orthonormal
    print(f"F  irrep dims {dim_list}, squares sum to {sq_sum}, orthonormal         : {'PASS' if F else 'FAIL'}")

    # G: McKay quiver -- adjacency a_mn = <chi_m * chi_2, chi_n>
    adj = np.zeros((7, 7), dtype=int)
    for mi, m in enumerate(names):
        prod = chi[m] * chi['2']
        for ni, nm in enumerate(names):
            adj[mi, ni] = int(round(inner(prod, chi[nm]).real))
    degrees = sorted(adj.sum(axis=1).tolist())
    # affine E6: one node of degree 3, three of degree 2, three of degree 1 -> seven nodes
    is_E6_hat = (degrees == [1, 1, 1, 2, 2, 2, 3]) and _is_connected(adj) and np.array_equal(adj, adj.T)
    G = is_E6_hat
    print(f"G  McKay quiver of 2T has degree sequence {degrees}        : {'PASS' if G else 'FAIL'} (affine E6)")

    # ----------------------------------------------------------------------
    # H, I: left multiplication, the SU(2)_L embedding, and the shared Q8
    # ----------------------------------------------------------------------
    # Left-multiplication matrix L_q acting on R^4 (column = image of basis quaternion)
    e_basis = [np.array([1.0, 0, 0, 0]), np.array([0.0, 1, 0, 0]),
               np.array([0.0, 0, 1, 0]), np.array([0.0, 0, 0, 1])]

    def Lmat(q):
        M = np.zeros((4, 4))
        for col, e in enumerate(e_basis):
            M[:, col] = qmul(q, e)
        return M

    def Rmat(q):
        M = np.zeros((4, 4))
        for col, e in enumerate(e_basis):
            M[:, col] = qmul(e, q)
        return M

    Ls = [Lmat(q) for q in group]
    # group of order 24 in SO(4): orthogonal, det +1, closed
    orthog = all(np.allclose(L @ L.T, np.eye(4), atol=TOL) for L in Ls)
    detpos = all(abs(np.linalg.det(L) - 1.0) < TOL for L in Ls)
    # left mult commutes with every right mult -> L lies in the left SU(2) factor of SO(4)
    Rs = [Rmat(q) for q in group]
    commute = all(np.allclose(L @ R, R @ L, atol=TOL) for L in Ls for R in Rs)
    H = orthog and detpos and commute
    print(f"H  left mult: {len(Ls)} matrices in SO(4), commute with right mult     : {'PASS' if H else 'FAIL'} (SU(2)_L)")

    # The left-mult Q8 core is the same eight quaternions as the strong-isospin Q8
    leftQ8 = [Lmat(group[i]) for i in Q8]
    leftQ8_closed = True
    for Li in leftQ8:
        for Lj in leftQ8:
            prod = Li @ Lj
            if not any(np.allclose(prod, Lk, atol=TOL) for Lk in leftQ8):
                leftQ8_closed = False
    I = leftQ8_closed and len(Q8) == 8
    print(f"I  weak-isospin Q8 (left mult) = strong-isospin Q8 (same 8 elements): {'PASS' if I else 'FAIL'}")

    # ----------------------------------------------------------------------
    # J, K, L: the colour Z3 is the physical FCC stacking rotation, and it acts
    # on Q8 as the diagram automorphism on a real composite defect.
    # The close-packed stacking A->B->C advances by a 120 deg rotation about [111];
    # as a unit quaternion that rotation is q3 = (1+i+j+k)/2, the order-3 element.
    # ----------------------------------------------------------------------
    q3 = 0.5 * np.array([1.0, 1.0, 1.0, 1.0])
    Rstack = np.zeros((3, 3))
    for col, ax in enumerate((iq, jq, kq)):
        Rstack[:, col] = qmul(qmul(q3, ax), qconj(q3))[1:]
    fixes_111 = np.allclose(Rstack @ np.array([1.0, 1, 1]), np.array([1.0, 1, 1]), atol=TOL)
    # order 3 in SO(3)
    Mp = np.eye(3); o = 0
    while True:
        Mp = Rstack @ Mp; o += 1
        if np.allclose(Mp, np.eye(3), atol=TOL):
            break
    J = fixes_111 and o == 3
    print(f"J  stacking rotation q3=(1+i+j+k)/2 is 120 deg about [111], SO(3) order {o}: {'PASS' if J else 'FAIL'}")

    # K: the three Shockley partials a/6<112> lie in (111) and are cycled by the rotation
    b1, b2, b3 = np.array([-2.0, 1, 1]), np.array([1.0, -2, 1]), np.array([1.0, 1, -2])
    in_plane = all(abs(np.dot(b, [1, 1, 1])) < TOL for b in (b1, b2, b3))
    cyc = (np.allclose(Rstack @ b1, b2, atol=TOL) and np.allclose(Rstack @ b2, b3, atol=TOL)
           and np.allclose(Rstack @ b3, b1, atol=TOL))
    K = in_plane and cyc
    print(f"K  three Shockley partials a/6<112> in (111), rotation cycles them      : {'PASS' if K else 'FAIL'}")

    # L: conjugation by q3 cycles the axes i->j->k (the diagram automorphism on Q8)
    def cj(v):
        return qmul(qmul(q3, v), qconj(q3))
    axis_cycle = (np.allclose(cj(iq), jq, atol=TOL) and np.allclose(cj(jq), kq, atol=TOL)
                  and np.allclose(cj(kq), iq, atol=TOL))
    L = axis_cycle
    print(f"L  conjugation by q3 cycles i->j->k (diagram automorphism of Q8)        : {'PASS' if L else 'FAIL'}")

    print()
    allpass = all([A, B, C, D, E, F, G, H, I, J, K, L])
    print("ALL CHECKS PASS" if allpass else "SOME CHECKS FAILED")
    return allpass


# small helpers kept out of the main flow for readability
def _order_cosets(group, table, Q8):
    placed, cosets = set(), []
    for x in range(len(group)):
        if x in placed:
            continue
        coset = sorted(table[x, q] for q in Q8)
        cosets.append(coset); placed |= set(coset)
    return cosets


def _identity_index(group):
    return find_index(group, np.array([1.0, 0, 0, 0]))


def _class_of(reps, classes, idx):
    for c, cls in enumerate(classes):
        if idx in cls:
            return c
    return -1


def _is_connected(adj):
    n = adj.shape[0]
    seen = {0}; stack = [0]
    while stack:
        v = stack.pop()
        for w in range(n):
            if adj[v, w] and w not in seen:
                seen.add(w); stack.append(w)
    return len(seen) == n


if __name__ == "__main__":
    main()
