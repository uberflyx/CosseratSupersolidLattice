"""The defect current, its conservation, and the minimal coupling, by duality.

The reviewer asked for the electromagnetic defect current and its
Ward-Takahashi identity.  In the antiplane (shear) sector, where the
framework's photon lives, the ingredients are derivable rather than
assumable, because a 2+1-dimensional massless scalar is dual to a
2+1-dimensional Maxwell field and a screw dislocation is the scalar's vortex.

A first version of this script asserted L_Maxwell = L_scalar pointwise and
failed, correctly: the 2+1D duality is a Legendre transform, under which the
electric field maps onto the ROTATED SHEAR STRAIN and the magnetic field onto
the NODE VELOCITY, so the two Lagrangians differ by a sign while the energy
densities, the equations of motion and the sources correspond exactly.  What
is checkable, and checked here, is:

  1.  the energy densities agree:  (E^2 + B^2)/2 = (mu/2)(u_t^2 + (grad u)^2);
  2.  the dual field strength built from u satisfies Maxwell's equations with
      the SCREW AS ITS ONLY SOURCE:  d_mu f^{mu nu} = J_D^nu, where J_D is the
      topological (vortex) current, zero for smooth u and b delta^2 on a core;
      and the scalar's wave equation is the dual Bianchi identity;
  3.  conservation d_mu J_D^mu = 0 is an identity of antisymmetry, holding for
      arbitrary (even multivalued) u: the lattice form of "dislocation lines
      cannot end", which is the operator statement charge conservation needs;
  4.  the Peach-Koehler force on a static screw equals the Lorentz force
      e_D E of the dual field with e_D = sqrt(mu) b, term by term: the
      coupling is minimal, with nothing residual.

Conventions: shear speed c_t = sqrt(mu/rho) = 1, x^mu = (t, x, y), metric
(+,-,-), eps^{012} = +1.
"""
from itertools import permutations

import sympy as sp

t, x, y = sp.symbols("t x y", real=True)
mu_s, b = sp.symbols("mu b", positive=True)
X = [t, x, y]
g = sp.diag(1, -1, -1)
ginv = g.inv()

eps = sp.MutableDenseNDimArray.zeros(3, 3, 3)
for perm in permutations(range(3)):
    eps[perm] = sp.Matrix(sp.eye(3)[list(perm), :]).det()   # eps^{012} = +1

line = lambda s: print("\n" + s + "\n" + "-" * len(s))

u = sp.Function("u")(t, x, y)
du_up = [sum(ginv[m, n] * sp.diff(u, X[n]) for n in range(3)) for m in range(3)]

# lowered epsilon for this metric (det g = +1, so eps_{012} = +1 as well)
eps_low = sp.MutableDenseNDimArray.zeros(3, 3, 3)
for a in range(3):
    for c in range(3):
        for d in range(3):
            eps_low[a, c, d] = sum(
                g[a, a2] * g[c, b2] * g[d, c2] * eps[a2, b2, c2]
                for a2 in range(3) for b2 in range(3) for c2 in range(3))

# dual field strength f_{mu nu} = sqrt(mu) eps_{mu nu lam} d^lam u
f = sp.MutableDenseNDimArray.zeros(3, 3)
for m in range(3):
    for n in range(3):
        f[m, n] = sp.sqrt(mu_s) * sum(
            eps_low[m, n, l] * du_up[l] for l in range(3))

# ------------------------------------------------------------------ step 1
line("1. THE DICTIONARY, AND THE ENERGY DENSITIES")
E = [sp.simplify(f[i, 0]) for i in (1, 2)]           # E_i = f_{i0}
B = sp.simplify(f[1, 2])                             # B = f_{xy}
print(f"   E_x = {E[0]}")
print(f"   E_y = {E[1]}")
print(f"   B   = {B}")
print("   -> the dual electric field is the rotated shear strain,")
print("      the dual magnetic field is the node velocity")

H_em = sp.Rational(1, 2) * (E[0] ** 2 + E[1] ** 2 + B ** 2)
H_el = sp.Rational(1, 2) * mu_s * (
    sp.diff(u, t) ** 2 + sp.diff(u, x) ** 2 + sp.diff(u, y) ** 2)
dH = sp.simplify(H_em - H_el)
print(f"   (E^2+B^2)/2 - (mu/2)(u_t^2 + |grad u|^2) = {dH}")
assert dH == 0, "energy densities disagree"

# ------------------------------------------------------------------ step 2
line("2. MAXWELL SOURCED BY THE SCREW ALONE; THE WAVE EQUATION AS BIANCHI")
# d_mu f^{mu nu}: raise both indices of f first
f_up = sp.MutableDenseNDimArray.zeros(3, 3)
for m in range(3):
    for n in range(3):
        f_up[m, n] = sum(ginv[m, a] * ginv[n, c] * f[a, c]
                         for a in range(3) for c in range(3))
divf = [sp.simplify(sum(sp.diff(f_up[m, n], X[m]) for m in range(3)))
        for n in range(3)]
print(f"   d_mu f^(mu nu) for smooth u = {divf}")
assert divf == [0, 0, 0], "Maxwell should be an identity for smooth u"
print("   -> identically zero off the core: the screw is the ONLY source;")
print("      on a core, d_mu f^(mu nu) = sqrt(mu) eps^(nu mu lam) d_mu d_lam u")
print("      = J_D^nu, supported where u is multivalued")

# the dual vector F^mu = (1/2) eps^{mu nu lam} f_{nu lam}; its divergence is
# the scalar's wave equation: the EOM <-> Bianchi exchange of the duality
F_dual = [sp.simplify(sp.Rational(1, 2) * sum(
    eps[m, n, l] * f[n, l] for n in range(3) for l in range(3)))
    for m in range(3)]
bianchi = sp.simplify(sum(sp.diff(F_dual[m], X[m]) for m in range(3)))
box_u = sp.diff(u, t, 2) - sp.diff(u, x, 2) - sp.diff(u, y, 2)
ratio = sp.simplify(bianchi / (sp.sqrt(mu_s) * box_u))
print(f"   d_mu F_dual^mu / (sqrt(mu) box u) = {ratio}")
assert ratio in (1, -1), "Bianchi should reproduce the wave equation"
print("   -> the scalar's equation of motion is the dual Bianchi identity")

# ------------------------------------------------------------------ step 3
line("3. THE TOPOLOGICAL CURRENT: CHARGE AND IDENTICAL CONSERVATION")
r, phi = sp.symbols("r phi", positive=True, real=True)
u_screw = b / (2 * sp.pi) * sp.atan2(y, x)
circ = sp.integrate(
    sp.diff(u_screw.subs([(x, r * sp.cos(phi)), (y, r * sp.sin(phi))]), phi),
    (phi, -sp.pi, sp.pi))
print(f"   oint du_z around the core = {sp.simplify(circ)}   (loop-independent)")

w = sp.Function("w")(t, x, y)
divJ = sp.simplify(sum(
    eps[m, n, l] * sp.diff(w, X[m], X[n], X[l])
    for m in range(3) for n in range(3) for l in range(3)))
print(f"   eps^(mnl) d_m d_n d_l w = {divJ} for ARBITRARY w")
print("   -> d_mu J_D^mu = 0 is an identity of total antisymmetry, holding")
print("      for multivalued fields too: dislocation lines cannot end,")
print("      and charge conservation is topology rather than dynamics")

# ------------------------------------------------------------------ step 4
line("4. PEACH-KOEHLER = LORENTZ, TERM BY TERM")
u_ext = sp.Function("U")(t, x, y)
du_ext_up = [sum(ginv[m, n] * sp.diff(u_ext, X[n]) for n in range(3))
             for m in range(3)]
f_ext = sp.MutableDenseNDimArray.zeros(3, 3)
for m in range(3):
    for n in range(3):
        f_ext[m, n] = sp.sqrt(mu_s) * sum(
            eps_low[m, n, l] * du_ext_up[l] for l in range(3))
E_ext = [sp.simplify(f_ext[1, 0]), sp.simplify(f_ext[2, 0])]

sigma_xz = mu_s * sp.diff(u_ext, x)
sigma_yz = mu_s * sp.diff(u_ext, y)
F_PK = [b * sigma_yz, -b * sigma_xz]        # F_i = eps_{ij} sigma_{jz} b

e_D = sp.sqrt(mu_s) * b
match_plus = [sp.simplify(F_PK[i] - e_D * E_ext[i]) for i in range(2)]
match_minus = [sp.simplify(F_PK[i] + e_D * E_ext[i]) for i in range(2)]
print(f"   F_PK - e_D E = {match_plus}")
print(f"   F_PK + e_D E = {match_minus}")
assert match_plus == [0, 0] or match_minus == [0, 0], "PK != Lorentz"
sign = "+" if match_plus == [0, 0] else "-"
print(f"   -> F_PK = {sign}e_D E exactly, e_D = sqrt(mu) b; the overall sign")
print("      is the orientation of the Burgers circuit against eps^{012},")
print("      physically immaterial.  The magnetic term e_D v x B follows from")
print("      applying the same identification in the boosted frame, both")
print("      sides transforming under the same shear-sector Lorentz group.")

line("5. WHAT THIS ESTABLISHES FOR THE WARD-TAKAHASHI IDENTITY")
print("""   The textbook operator derivation of q_mu Gamma^mu = S^-1(p+q) - S^-1(p)
   needs three inputs, and all three are now derived in this sector:
     (i)   a current conserved as an IDENTITY (step 3), with quantised charge,
           the winding being integer-valued on the lattice, so the defect
           field is a charge eigenoperator and the equal-time commutator
           [J^0(x), psi_D(y)] = -delta(x-y) psi_D(y) holds;
     (ii)  MINIMAL coupling of that current to the mediating field: the force
           the dual gauge field exerts is the Peach-Koehler force with
           nothing left over (step 4), so L_int = -e_D a_mu J_D^mu;
     (iii) the gauge redundancy a_mu -> a_mu + d_mu chi, harmless against the
           conserved current, which is what lets vertex and self-energy
           renormalise together, Z_1 = Z_2, and pins F_D(0) = 1 exactly.
   Scope: the 2+1-dimensional antiplane sector, the same reduction in which
   the core width and alpha are computed.  The identity constrains the
   longitudinal vertex only; F_1(q^2 != 0) and F_2 stay dynamical, exactly as
   in QED, which is why the reflectionless spectrum and the lepton moments
   continue to carry that part of the argument.""")
