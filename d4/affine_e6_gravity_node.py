"""
Is gravity the affine node of E6-hat?

The vacuum-symmetry chapter reaches the affine diagram E6-hat by McKay duality on
the root group 2T = SL(2,3), then DELETES the affine node to get the finite Lie
group E6, and concludes "gravity stands outside E6." This script tests a sharper
reading: that the deleted node is gravity itself.

Two independent symmetries label one lattice.
  - O_h (point group) sorts the FORCES: each irrep is one deformation channel.
  - 2T  (root group)  sorts the MATTER: its McKay dual is the affine E6 diagram.

They are linked through the rotational tetrahedral group T = A4, which is at once
    a subgroup of O  (so O_h irreps restrict to it) and
    the quotient 2T/{+-1} (so its irreps ARE the non-spinorial McKay nodes).
So every force channel, restricted to A4 and lifted, lands on one McKay node.

The script:
  1. builds 2T as 24 unit quaternions, constructs all 7 irreps EXPLICITLY,
     verifies the character table by orthogonality;
  2. builds the McKay graph (tensor with the natural 2-dim rep), confirms it is
     affine E6, and identifies the affine node (trivial rep);
  3. branches each O_h force channel  O_h -> O=S4 -> A4 -> McKay node;
  4. prints the force/gravity -> node correspondence.

No physics input beyond "A1g = compression = gravity" (established in the
monograph). Everything else is character theory.
"""

import numpy as np
import itertools

w = np.exp(2j * np.pi / 3)  # primitive cube root of unity

# ---------------------------------------------------------------------------
# 1. The root group 2T = binary tetrahedral group, as 24 unit quaternions.
#    8 Lipschitz units  {+-1,+-i,+-j,+-k}  and  16 Hurwitz units (+-1+-i+-j+-k)/2.
# ---------------------------------------------------------------------------
def quaternions_2T():
    units = []
    for axis in range(4):                      # the 8 axis units +-1,+-i,+-j,+-k
        for sign in (1.0, -1.0):
            q = np.zeros(4)
            q[axis] = sign
            units.append(q)
    for s in itertools.product([0.5, -0.5], repeat=4):   # the 16 half-integer units
        units.append(np.array(s, dtype=float))
    return units

def qmul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return np.array([
        a0*b0 - a1*b1 - a2*b2 - a3*b3,
        a0*b1 + a1*b0 + a2*b3 - a3*b2,
        a0*b2 - a1*b3 + a2*b0 + a3*b1,
        a0*b3 + a1*b2 - a2*b1 + a3*b0,
    ])

def q_to_su2(q):
    """Natural 2-dim (spin-1/2) representation: quaternion -> SU(2) matrix."""
    a, b, c, d = q
    return np.array([[a + 1j*b, c + 1j*d],
                     [-c + 1j*d, a - 1j*b]])

G = quaternions_2T()
assert len(G) == 24

# group multiplication table (indices) and inverses, by nearest-quaternion match
def find(q):
    for idx, g in enumerate(G):
        if np.allclose(g, q, atol=1e-9):
            return idx
    raise ValueError("not a group element")

mult = [[find(qmul(G[i], G[j])) for j in range(24)] for i in range(24)]
inv = [next(j for j in range(24) if mult[i][j] == find(np.array([1.0, 0, 0, 0]))) for i in range(24)]

# conjugacy classes
seen = [False]*24
classes = []
for i in range(24):
    if seen[i]:
        continue
    orb = set()
    for g in range(24):
        orb.add(mult[mult[g][i]][inv[g]])
    for x in orb:
        seen[x] = True
    classes.append(sorted(orb))
classes.sort(key=lambda c: (len(c), c[0]))
print(f"2T has {len(classes)} conjugacy classes, sizes {[len(c) for c in classes]}")
reps = [c[0] for c in classes]          # one representative element per class

# ---------------------------------------------------------------------------
# Build all 7 irreps explicitly, then read characters on the class reps.
# ---------------------------------------------------------------------------
def char_of_matrix_rep(matrices):
    return np.array([np.trace(matrices[r]) for r in reps])

# trivial and the two other 1-dim reps, via 2T -> 2T/Q8 = Z3.
# Q8 = the 8 axis-unit quaternions (indices); the Z3 grade is which coset.
Q8 = set(range(24))
axis_idx = set()
for idx, g in enumerate(G):
    if sum(abs(x) for x in g) == 1:
        axis_idx.add(idx)
# coset grade in Z3: pick generator t = (1+i+j+k)/2, grade(g) = which power of t-coset
t = find(np.array([0.5, 0.5, 0.5, 0.5]))
def z3_grade(g):
    # g * Q8 coset: test membership of g, g*t^{-1}, g*t^{-2} in Q8
    for n in range(3):
        gg = g
        for _ in range(n):
            gg = mult[gg][inv[t]]
        if gg in axis_idx:
            return n
    raise ValueError
grades = [z3_grade(g) for g in range(24)]

one  = [np.array([[1.0]]) for _ in range(24)]
onew  = [np.array([[w**grades[g]]]) for g in range(24)]
onew2 = [np.array([[w**(2*grades[g])]]) for g in range(24)]

# natural 2-dim (faithful), and its two twists by the 1-dim characters
two  = [q_to_su2(G[g]) for g in range(24)]
two_w  = [two[g] * w**grades[g] for g in range(24)]
two_w2 = [two[g] * w**(2*grades[g]) for g in range(24)]

# 3-dim: symmetric square of the natural 2-dim is the spin-1 rep (dim 3)
def sym2(M):
    a, b = M[0, 0], M[0, 1]
    c, d = M[1, 0], M[1, 1]
    # basis e1^2, e1e2(sym), e2^2
    return np.array([
        [a*a, 2*a*b, b*b],
        [a*c, a*d + b*c, b*d],
        [c*c, 2*c*d, d*d],
    ])
three = [sym2(two[g]) for g in range(24)]

irreps = {
    "1   (rho0, trivial)": one,
    "1'  (rho1, colour w)": onew,
    "1'' (rho2, colour w2)": onew2,
    "2   (rho3, gen 1)": two,
    "2'  (rho4, gen 2)": two_w,
    "2'' (rho5, gen 3)": two_w2,
    "3   (rho6, vector)": three,
}
names = list(irreps)
chars = {n: char_of_matrix_rep(m) for n, m in irreps.items()}

# verify it is a complete, correct character table by orthogonality
sizes = np.array([len(c) for c in classes])
def inner(a, b):
    return np.sum(sizes * a * np.conj(b)) / 24
ortho_ok = True
for i, ni in enumerate(names):
    for j, nj in enumerate(names):
        val = inner(chars[ni], chars[nj])
        expected = 1.0 if i == j else 0.0
        if not np.isclose(val, expected, atol=1e-9):
            ortho_ok = False
print(f"character table orthonormal (7 distinct irreps, dims {[int(round(chars[n][0].real)) for n in names]}): {ortho_ok}")
print(f"sum of squares of dims = {sum(int(round(chars[n][0].real))**2 for n in names)} (should be 24)")

# ---------------------------------------------------------------------------
# 2. McKay graph: connect R--S with multiplicity <R (x) 2_natural, S>.
# ---------------------------------------------------------------------------
nat = chars["2   (rho3, gen 1)"]
print("\nMcKay graph (tensor with natural 2-dim); entry = edge multiplicity:")
print(f"{'':24}" + "".join(f"{n.split('(')[0].strip():>6}" for n in names))
mckay = {}
for ni in names:
    row = []
    for nj in names:
        m = inner(chars[ni] * nat, chars[nj])
        row.append(int(round(m.real)))
    mckay[ni] = row
    print(f"{ni:24}" + "".join(f"{x:6d}" for x in row))

# affine node = trivial; its neighbours; node "marks" = dims
dims = {n: int(round(chars[n][0].real)) for n in names}
print("\nNode dims (= Coxeter/Kac marks of affine E6): "
      + ", ".join(f"{n.split('(')[0].strip()}:{dims[n]}" for n in names))
print("Three tips (mark 1) =", [n.split('(')[0].strip() for n in names if dims[n] == 1])
print("Centre (mark 3)     =", [n.split('(')[0].strip() for n in names if dims[n] == 3])
print("Mids (mark 2)        =", [n.split('(')[0].strip() for n in names if dims[n] == 2])
print("=> affine node deleted to reach finite E6 is the TRIVIAL rep rho0.")

# ---------------------------------------------------------------------------
# 3. Branch O_h force channels  O_h -> O=S4 -> A4 -> McKay node.
#    A4 irreps {triv, 1_w, 1_w2, 3} are exactly the non-spinorial 2T nodes
#    {rho0, rho1, rho2, rho6}.  S4 -> A4 branching by class fusion.
# ---------------------------------------------------------------------------
# S4 character table (classes e,(12),(12)(34),(123),(1234); sizes 1,6,3,8,6)
S4 = {
    "A1":  np.array([1,  1,  1,  1,  1]),
    "A2":  np.array([1, -1,  1,  1, -1]),
    "E":   np.array([2,  0,  2, -1,  0]),
    "T1":  np.array([3,  1, -1,  0, -1]),
    "T2":  np.array([3, -1, -1,  0,  1]),
}
# A4 classes: e,(12)(34),(123),(132); both 3-cycle classes fuse from S4 (123).
# restrict S4 char to A4 classes -> (e, (12)(34), (123), (132)):
def s4_to_a4(chi):
    return np.array([chi[0], chi[2], chi[3], chi[3]])
A4 = {
    "triv": np.array([1, 1, 1, 1]),
    "1_w":  np.array([1, 1, w, w**2]),
    "1_w2": np.array([1, 1, w**2, w]),
    "3":    np.array([3, -1, 0, 0]),
}
a4_sizes = np.array([1, 3, 4, 4])
def a4_decomp(chi):
    out = []
    for nm, ch in A4.items():
        m = np.sum(a4_sizes * chi * np.conj(ch)) / 12
        if round(m.real) != 0:
            out.append((nm, int(round(m.real))))
    return out

node_of_a4 = {"triv": "rho0 (affine node)", "1_w": "rho1 (tip)",
              "1_w2": "rho2 (tip)", "3": "rho6 (centre)"}
force_of_ohirrep = {
    "A1g": "GRAVITY (compression / breathing)",
    "A2u": "pseudoscalar A2u (eta' anomaly)",
    "Eg":  "nuclear tensor (f2 1270)",
    "T1u": "electromagnetism (photon)",
    "T2g": "strong",
    "T1g": "weak",
}
oh_to_o = {"A1g": "A1", "A2u": "A2", "Eg": "E", "T1u": "T1", "T2g": "T2", "T1g": "T1"}

print("\nForce channel  ->  O_h irrep  ->  A4 content  ->  McKay node(s):")
for oh, force in force_of_ohirrep.items():
    a4c = a4_decomp(s4_to_a4(S4[oh_to_o[oh]]))
    nodes = " + ".join(f"{m}x {node_of_a4[nm]}" for nm, m in a4c)
    print(f"  {oh:4} {force:34} -> {nodes}")

print("\nThe three spinorial nodes rho3,rho4,rho5 (the generations) are NOT reached")
print("by any deformation channel: integer-spin O_h modes restrict to integer-spin")
print("A4 irreps only. Matter lives in the double cover; force/gravity do not.")

# ---------------------------------------------------------------------------
# 4. The colour Z3 cyclically relates the three spinor nodes (and the 3 tips).
#    Tensoring by the colour character chi_w should send rho3->rho4->rho5->rho3.
# ---------------------------------------------------------------------------
def decomp(chi):
    return {n: int(round(inner(chi, chars[n]).real)) for n in names}
def single(chi):
    d = decomp(chi); return next(n for n, m in d.items() if m == 1)
chiw = chars["1'  (rho1, colour w)"]
cycle2 = [single(chars[n] * chiw) for n in
          ["2   (rho3, gen 1)", "2'  (rho4, gen 2)", "2'' (rho5, gen 3)"]]
cycle1 = [single(chars[n] * chiw) for n in
          ["1   (rho0, trivial)", "1'  (rho1, colour w)", "1'' (rho2, colour w2)"]]
print("\nColour character (x) spinor nodes cycles them:",
      " -> ".join(c.split('(')[0].strip() for c in cycle2))
print("Colour character (x) tip nodes cycles them:    ",
      " -> ".join(c.split('(')[0].strip() for c in cycle1))
print("Centre rho6 fixed by colour:",
      single(chars["3   (rho6, vector)"] * chiw).split('(')[0].strip() == "3")
