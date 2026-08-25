"""Angular harmonic content of the screw and edge dislocation fields.

The paper labels these fields by multipole order in three places and the labels
disagree.  This resolves the disagreement by decomposing each field explicitly,
in both the Cartesian and the polar component basis, since the harmonic content
is basis-dependent and that is the likely origin of the confusion.

Line along z throughout.  Screw: b along z.  Edge: b along x, glide plane normal
along y.  Isotropic linear elasticity, Hirth & Lothe conventions.
"""
import sympy as sp

r, th, b, mu, nu = sp.symbols('r theta b mu nu', positive=True)

def harmonics(expr, mmax=4):
    """Fourier content of expr(theta) over [0, 2pi], as {(kind, m): coeff}."""
    out = {}
    a0 = sp.simplify(sp.integrate(expr, (th, 0, 2*sp.pi)) / (2*sp.pi))
    if sp.simplify(a0) != 0:
        out[('const', 0)] = sp.simplify(a0)
    for m in range(1, mmax + 1):
        c = sp.simplify(sp.integrate(expr*sp.cos(m*th), (th, 0, 2*sp.pi))/sp.pi)
        s = sp.simplify(sp.integrate(expr*sp.sin(m*th), (th, 0, 2*sp.pi))/sp.pi)
        if sp.simplify(c) != 0:
            out[('cos', m)] = sp.simplify(c)
        if sp.simplify(s) != 0:
            out[('sin', m)] = sp.simplify(s)
    return out

def show(title, fields):
    print(f"\n{title}")
    for name, expr in fields.items():
        h = harmonics(sp.simplify(expr))
        terms = ", ".join(f"{k[0]}({k[1]}theta)" if k[0] != 'const' else "const"
                          for k in h) or "zero"
        print(f"  {name:22s} -> {terms}")

# ---------------- SCREW ----------------------------------------------------
# u_z = b*theta/(2 pi); the only strains are the anti-plane pair.
show("SCREW, Cartesian strain components (x r, in units of b/(4 pi))", {
    "eps_xz": -sp.sin(th),
    "eps_yz":  sp.cos(th),
})
show("SCREW, polar strain components (x r, in units of b/(4 pi))", {
    "eps_rz":   sp.cos(th)*sp.cos(th) + sp.sin(th)*sp.sin(th) - 1,   # = 0
    "eps_thz":  sp.Integer(1),
})

# ---------------- EDGE -----------------------------------------------------
# Polar stresses: sigma_rr = sigma_thth = -D sin(th)/r,  sigma_rth = D cos(th)/r
show("EDGE, polar stress components (x r/D)", {
    "sigma_rr":   -sp.sin(th),
    "sigma_thth": -sp.sin(th),
    "deviator (sigma_rr - sigma_thth)/2": sp.Integer(0),
    "sigma_rth":   sp.cos(th),
})

# Cartesian stresses, from the standard plane-strain solution.
D = 1  # overall scale irrelevant to the harmonic content
sxx = -sp.sin(th)*((1 - 2*nu) + 2*sp.cos(th)**2)
syy =  sp.sin(th)*((3 - 2*nu) - 2*sp.cos(th)**2)
sxy =  sp.cos(th)*((1 - 2*nu) + 2*sp.sin(th)**2)
show("EDGE, Cartesian stress components (x r, general nu)", {
    "sigma_xx": sxx, "sigma_yy": syy, "sigma_xy": sxy,
})
show("EDGE, Cartesian stress components at the incompressible point nu=1/2", {
    "sigma_xx": sxx.subs(nu, sp.Rational(1,2)),
    "sigma_yy": syy.subs(nu, sp.Rational(1,2)),
    "sigma_xy": sxy.subs(nu, sp.Rational(1,2)),
})

# Displacement field of the edge: the ramp plus genuine second harmonics.
show("EDGE, displacement (x 2 pi / b), harmonic part only", {
    "u_x  (drop theta ramp)": sp.sin(2*th)/(4*(1-nu)),
    "u_y  (drop ln r)":       -sp.cos(2*th)/(4*(1-nu)),
})
