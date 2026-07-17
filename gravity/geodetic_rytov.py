"""Geodetic precession as Rytov rotation: symbolic holonomy and the GP-B number.

The spin axis of a defect gyroscope parallel-transports along its ray through
the graded-index lattice (the Rytov-Vladimirskii law). For the conformal
spatial metric dl^2 = n^2 (dr^2 + r^2 dphi^2), transport once around a circular
orbit of radius r rotates the axis by

    Delta_alpha_space = -2 pi r d(ln n)/dr = 2 pi eps      (eps = GM/c^2 r)

and the Thomas rotation of the accumulated boosts adds pi*eps, giving the
de Sitter rate Omega_geo = (3/2) eps Omega. This script verifies the holonomy
symbolically from the conformal connection and evaluates Gravity Probe B.
"""

import numpy as np
import sympy as sp

# --- symbolic holonomy --------------------------------------------------------
# Conformal 2D metric g_ij = n(r)^2 delta_ij in Cartesian coordinates.
# Transport equation: dV^i/dphi + Gamma^i_jk t^j V^k = 0 along x(phi) = r(cos, sin).
r_s, phi = sp.symbols('r phi', positive=True)
x, y = r_s*sp.cos(phi), r_s*sp.sin(phi)
n = sp.Function('n')(sp.sqrt(x**2 + y**2))

X = [x, y]
g = [[n**2 if i == j else 0 for j in range(2)] for i in range(2)]
ginv = [[1/n**2 if i == j else 0 for j in range(2)] for i in range(2)]

# Christoffels of the conformal metric, evaluated on the circle
coords = sp.symbols('X0 X1')
nn = sp.Function('n')(sp.sqrt(coords[0]**2 + coords[1]**2))
Gamma = [[[sp.simplify(
    sp.Rational(1,2)*sum(1/nn**2 * (1 if i == l else 0) * (
        sp.diff(nn**2*(1 if l == j else 0), coords[k]) +
        sp.diff(nn**2*(1 if l == k else 0), coords[j]) -
        sp.diff(nn**2*(1 if j == k else 0), coords[l]))
        for l in range(2)))
    for k in range(2)] for j in range(2)] for i in range(2)]

t = [sp.diff(x, phi), sp.diff(y, phi)]        # tangent dx/dphi
V1, V2 = sp.Function('V1')(phi), sp.Function('V2')(phi)
V = [V1, V2]
sub = {coords[0]: x, coords[1]: y}
eqs = [sp.Eq(sp.diff(V[i], phi),
             sp.simplify(-sum(Gamma[i][j][k].subs(sub) * t[j] * V[k]
                              for j in range(2) for k in range(2))))
       for i in range(2)]

# Specialise to n = 1 + a/rho (a = GM/c^2), keep leading order in a
a = sp.symbols('a', positive=True)
rho = sp.symbols('rho', positive=True)
nfun = 1 + a/rho
eqs_lin = [sp.Eq(e.lhs, sp.series(e.rhs.subs(sp.Function('n')(rho), nfun)
           .subs(rho, r_s).doit(), a, 0, 2).removeO()) for e in eqs]

# Solve the linearised transport over one revolution with V(0) = (1, 0)
sol = sp.dsolve(eqs_lin, [V1, V2], ics={V1.subs(phi, 0): 1, V2.subs(phi, 0): 0})
V1f = sp.series(sol[0].rhs.subs(phi, 2*sp.pi), a, 0, 2).removeO()
V2f = sp.series(sol[1].rhs.subs(phi, 2*sp.pi), a, 0, 2).removeO()
alpha = sp.simplify(sp.atan2(V2f, V1f).rewrite(sp.atan))
print(f"rotation after one orbit (leading order in eps = a/r): {sp.simplify(V2f)}")
print("expected: 2*pi*a/r  (i.e. 2 pi eps, prograde)")

# --- Gravity Probe B ----------------------------------------------------------
c    = 2.99792458e8
GM_E = 3.986004418e14
r_gp = 6.371e6 + 6.42e5          # 642 km altitude
eps  = GM_E/(c**2 * r_gp)
Omega = np.sqrt(GM_E/r_gp**3)
Om_geo = 1.5 * eps * Omega        # rad/s
mas_yr = Om_geo * 3.156e7 * np.degrees(1) * 3600e3
print(f"\nGP-B: eps = {eps:.3e}, Omega_geo = {Om_geo:.3e} rad/s "
      f"= {mas_yr:.0f} mas/yr (measured 6602 +/- 18)")

# Earth-Moon around the Sun: the de Sitter precession of the lunar orbit
GM_S = 1.32712440018e20
AU   = 1.495978707e11
eps_s = GM_S/(c**2*AU)
Om_dS = 1.5*eps_s*np.sqrt(GM_S/AU**3) * 3.156e7 * np.degrees(1)*3600e3
print(f"Earth-Moon de Sitter: {Om_dS:.1f} mas/yr (LLR-confirmed ~19.2)")
