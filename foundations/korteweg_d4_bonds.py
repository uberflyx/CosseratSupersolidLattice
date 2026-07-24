"""Korteweg gradient coefficient on the D4 lattice: 12-bond slice vs all 24 bonds.

Checks, exactly (Fraction arithmetic) and numerically:
  1. The 12 compact-crossing bonds of D4 are diagonals: each carries a spatial
     displacement, and together they carry 1/3 of the total spatial second moment.
  2. The FCC 12-bond fourth moment is anisotropic (so the slice-only beta_latt
     needs an angular average and is direction-dependent at the +/-14% level).
  3. The full 24-bond (24-cell) fourth moment is exactly isotropic (spherical
     design property), so beta_latt is a single number with no averaging:
     beta_latt = c^2 l^2 / 24, replacing c^2 l^2 / 20.
  4. Downstream: f_s = 5/6 (was 4/5), v_2 = 4c exactly (was sqrt(40/3) c),
     Debye-Waller exponent 2W = 12.9/(1-f_s) strengthens the Lorentz
     suppression by ~5.6 decades.
  5. Independence of the compact-axis choice: any 3D hyperplane projection of
     the 24-cell gives the same second and fourth moments.
"""

from fractions import Fraction as Fr
from itertools import permutations, product
import numpy as np

# ── D4 minimal vectors: all permutations of (+-1, +-1, 0, 0), squared length 2 ──
def d4_vectors():
    vecs = set()
    for pos in permutations(range(4), 2):
        for s1, s2 in product((1, -1), repeat=2):
            v = [0, 0, 0, 0]
            v[pos[0]], v[pos[1]] = s1, s2
            vecs.add(tuple(v))
    return sorted(vecs)

V = d4_vectors()
assert len(V) == 24, "D4 kissing number is 24"

# Integer units: |delta|^2 = 2  <->  physical bond length l  =>  1 unit^2 = l^2/2.
UNIT2 = Fr(1, 2)   # l^2 per integer unit^2

# Split by the compact axis e4
in_slice = [v for v in V if v[3] == 0]
crossing = [v for v in V if v[3] != 0]
assert (len(in_slice), len(crossing)) == (12, 12)

# ── 1. spatial second moments ──
def spatial_M2(vecs):
    """Sum over bonds of the spatial (first three components) second-moment tensor."""
    M = [[Fr(0)] * 3 for _ in range(3)]
    for v in vecs:
        for i in range(3):
            for j in range(3):
                M[i][j] += Fr(v[i] * v[j]) * UNIT2
    return M

M_slice, M_cross = spatial_M2(in_slice), spatial_M2(crossing)
tr = lambda M: sum(M[i][i] for i in range(3))
print("spatial second moment (units of l^2):")
print(f"  in-slice  Sum|d_s|^2 = {tr(M_slice)}   tensor diag {[M_slice[i][i] for i in range(3)]}")
print(f"  crossing  Sum|d_s|^2 = {tr(M_cross)}   tensor diag {[M_cross[i][i] for i in range(3)]}")
print(f"  crossing share of total = {tr(M_cross) / (tr(M_slice) + tr(M_cross))}  (claim: 1/3)")

# ── 2 & 3. fourth moments along arbitrary spatial directions ──
def S4_along(vecs, khat):
    """Sum over bonds of (delta_spatial . khat)^4, in units of l^4."""
    total = Fr(0)
    for v in vecs:
        d = sum(Fr(v[i]) * khat[i] for i in range(3))
        total += d ** 4 * UNIT2 ** 2
    return total

def S2_along(vecs, khat):
    total = Fr(0)
    for v in vecs:
        d = sum(Fr(v[i]) * khat[i] for i in range(3))
        total += d ** 2 * UNIT2
    return total

dirs = {
    "[100]": (Fr(1), Fr(0), Fr(0)),
    "[110]/sqrt2": None,   # irrational; handled numerically below
    "[111]/sqrt3": None,
}

print("\nfourth moment Sum (d.khat)^4 (units of l^4):")
print(f"  {'direction':<14} {'12-bond slice':>14} {'24-bond D4':>12}")
k100 = (Fr(1), Fr(0), Fr(0))
print(f"  {'[100]':<14} {str(S4_along(in_slice, k100)):>14} {str(S4_along(V, k100)):>12}")

# numerical for irrational directions
def S4_num(vecs, khat):
    k = np.asarray(khat, float); k /= np.linalg.norm(k)
    A = np.array([v[:3] for v in vecs], float) * np.sqrt(0.5)  # physical units of l
    return float(np.sum((A @ k) ** 4))

for name, k in (("[110]", (1, 1, 0)), ("[111]", (1, 1, 1)), ("[210]", (2, 1, 0))):
    print(f"  {name:<14} {S4_num(in_slice, k):>14.6f} {S4_num(V, k):>12.6f}")

# random directions: verify exact isotropy of the 24-bond set
rng = np.random.default_rng(1)
vals24 = [S4_num(V, rng.normal(size=3)) for _ in range(2000)]
vals12 = [S4_num(in_slice, rng.normal(size=3)) for _ in range(2000)]
print(f"\n  12-bond S4 over random khat: min {min(vals12):.4f}, max {max(vals12):.4f}, "
      f"mean {np.mean(vals12):.4f}  (angular avg 12/5 = 2.4; spread = anisotropy)")
print(f"  24-bond S4 over random khat: min {min(vals24):.4f}, max {max(vals24):.4f}  "
      f"(exactly 3: spherical-design isotropy)")

# ── beta_latt from the ratio of fourth to second moments ──
# omega^2 = (2K_s/m0)[ (1/4) S2 k^2 - (1/48) S4 k^4 ],  calibrated so the k^2 term is c^2 k^2
# =>  beta_latt / c^2 = (1/12) * S4 / S2   (per unit l^2)
def beta(vecs, khat_exact):
    s2, s4 = S2_along(vecs, khat_exact), S4_along(vecs, khat_exact)
    return Fr(1, 12) * s4 / s2

print("\nbeta_latt / (c^2 l^2):")
print(f"  12-bond, angular-averaged S4=12/5, S2=4 : {Fr(1,12) * Fr(12,5) / 4}   (the printed 1/20)")
print(f"  24-bond, exact any direction           : {beta(V, k100)}   (the corrected value)")

# ── downstream ──
beta_QM = Fr(1, 4)
for label, b in (("12-bond", Fr(1, 20)), ("24-bond", Fr(1, 24))):
    ratio = b / beta_QM
    fs = 1 - ratio
    fn = 1 - fs
    v2sq = Fr(8, 3) / fn          # v2^2/c^2 = C11/(mu f_n), C11 = (8/3) mu
    W = 12.9 / float(1 - fs)      # Debye-Waller exponent 2W = 12.9/(1-f_s)
    print(f"\n  {label}: beta_latt/beta_QM = {ratio},  f_s = {fs},  f_n = {fn}")
    print(f"           v2 = sqrt({v2sq}) c = {float(v2sq)**0.5:.4f} c")
    print(f"           2W = {W:.1f}  ->  Bragg suppression e^-2W = 10^-{W/np.log(10):.1f}")

d = (12.9/ (1/6.) - 12.9/(1/5.)) / np.log(10)
print(f"\n  Lorentz-suppression gain from f_s = 4/5 -> 5/6: {d:.2f} decades (stronger)")

# ── 5. projection independence: random 3D hyperplanes through the 24-cell ──
print("\nprojection independence (random compact-axis choices):")
for trial in range(3):
    # random unit 4-vector as the compact axis; project bonds onto its orthogonal 3-space
    n = rng.normal(size=4); n /= np.linalg.norm(n)
    A = np.array(V, float) * np.sqrt(0.5)          # physical units
    P = np.eye(4) - np.outer(n, n)                 # projector onto the spatial 3-space
    B = A @ P.T
    # orthonormal basis of the 3-space
    Q, _ = np.linalg.qr(np.c_[n, rng.normal(size=(4, 3))])
    E3 = Q[:, 1:]                                  # 4x3, spans the spatial hyperplane
    C = A @ E3                                     # spatial components, shape (24,3)
    S2 = np.einsum('bi,bj->ij', C, C)
    k = rng.normal(size=3); k /= np.linalg.norm(k)
    S4 = float(np.sum((C @ k) ** 4))
    print(f"  axis {trial+1}: Sum d_s d_s = {S2[0,0]:.4f} I (off-diag max "
          f"{np.max(np.abs(S2 - np.diag(np.diag(S2)))):.1e}),  S4(khat) = {S4:.6f}")
print("  -> second moment 6 l^2 I and fourth moment 3 l^4, whatever axis is compact.")
