"""
The condensate charge is the Peirce grading of the stacking cell.

Two results from opposite ends of the monograph meet here. The superfluid
chapter identifies the U(1) that completes Spin(10) to E6 as the condensate
phase. The symmetry chapter identifies the 27 of E6 as the stacking cell,
whose Peirce decomposition with respect to the defect's layer has eigenvalues
lambda = 1 (the layer's own amplitude), 1/2 (the sixteen coordinates of the
two touching registries, the matter), and 0 (the ten it does not touch).

Under E6 -> SO(10) x U(1)_psi the 27 branches as

    27 = 16(q=+1) + 10(q=-2) + 1(q=+4)

(charges in the standard integer normalisation).  This script verifies the
exact affine identity

    q_psi = 6 * lambda_Peirce - 2

on all three blocks, and checks the two structural facts that pin it: the
charge is traceless over the 27 (as a generator of E6 must be), and the
Peirce eigenvalues sum to the rank-weighted trace (sum lambda = 9 = 3 * 3).

Physical reading: the condensate charge of a cell coordinate measures its
contact with the defect's own layer.  Full contact (the layer amplitude)
carries +4; half contact (the two touching slip registries, the matter 16)
carries +1; no contact (the far registry and the other two amplitudes)
carries -2.  The U(1) the superfluid supplies is not an extra label bolted
onto the cell; it IS the cell's contact grading, the same splitting that
makes the 16 the matter.  Geometry and charge are one function.
"""

from fractions import Fraction as F

# Peirce eigenvalue and dimension of each block of the 27.
blocks = [
    ("layer amplitude (1)",        F(1),      1),
    ("touching registries (16)",   F(1, 2),  16),
    ("far registry + rest (10)",   F(0),     10),
]

# Standard E6 -> SO(10) x U(1)_psi charges of 27 = 16 + 10 + 1.
q_standard = {"layer amplitude (1)": F(4),
              "touching registries (16)": F(1),
              "far registry + rest (10)": F(-2)}

line = "-" * 74
print(__doc__.strip().splitlines()[0])
print(line)
print(f"  {'block':<28s} {'dim':>4s} {'lambda':>7s} {'6*lambda-2':>11s} {'q_psi':>6s}")
ok = True
trace_q = F(0)
trace_lam = F(0)
for name, lam, dim in blocks:
    q_pred = 6 * lam - 2
    q_std = q_standard[name]
    ok &= (q_pred == q_std)
    trace_q += dim * q_std
    trace_lam += dim * lam
    print(f"  {name:<28s} {dim:>4d} {str(lam):>7s} {str(q_pred):>11s} {str(q_std):>6s}")
assert ok, "affine identity fails"
assert trace_q == 0, "U(1) generator is not traceless"
assert trace_lam == 9, "Peirce trace is not 9"
print(line)
print(f"  identity q_psi = 6*lambda - 2 holds on all three blocks: PASS")
print(f"  tracelessness: sum(dim * q) = {trace_q} (must be 0):       PASS")
print(f"  Peirce trace:  sum(dim * lambda) = {trace_lam} = 3 x rank:    PASS")
print(line)
print("VERDICT")
print("  The U(1) that lifts Spin(10) to E6, the condensate phase, is the")
print("  Peirce grading of the stacking cell: q = 6*lambda - 2.  A cell")
print("  coordinate's condensate charge is its contact with the defect's")
print("  layer.  The superfluid's completing charge and the geometry that")
print("  selects the matter 16 are the same function, so the two chapters")
print("  are welded by one line.")
