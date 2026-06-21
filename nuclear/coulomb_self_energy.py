"""
The 3/5 in the nuclear Coulomb self-energy, from first principles.

Cosserat Supersolid Lattice framework. The Coulomb term of the mass formula reads
the screw-charge field as Maxwell's, so the self-energy of a uniformly charged drop
is E = (3/5) Q^2 / (4 pi eps0 R). This builds that drop shell by shell and checks
the 3/5 numerically; it is a geometric integral, not a fit.
"""
import numpy as np
# numerical check of E = (3/5) Q^2/(4 pi eps0 R) by direct shell integration
# work in units 1/(4 pi eps0) = 1, Q=1, R=1; uniform density rho = Q/((4/3) pi R^3)
R, Q = 1.0, 1.0
rho = Q/((4/3)*np.pi*R**3)
r = np.linspace(1e-9, R, 200000)
qr = (4/3)*np.pi*r**3*rho            # charge enclosed at radius r
dq = 4*np.pi*r**2*rho               # shell charge per dr
phi = np.where(r>0, qr/r, 0.0)      # potential at surface of inner sphere (1/4pieps0=1)
E_num = np.trapezoid(phi*dq, r)         # dE = phi dq, integrated
print(f"numerical self-energy  E = {E_num:.6f}  (units Q^2/4pi eps0 R)")
print(f"analytic 3/5           = {3/5:.6f}")
print(f"ratio                  = {E_num/(3/5):.6f}\n")

# the 2/5 the user may be thinking of: moment of inertia of a solid sphere
# I = integral r_perp^2 dm ; for uniform sphere = (2/5) M R^2 -- a DIFFERENT integral
M=1.0
# I = (2/3) integral_0^R r^2 dm_shell with dm = (3M/R^3) r^2 dr ... do it:
dm = (M/((4/3)*np.pi*R**3))*4*np.pi*r**2
I_num = np.trapezoid((2/3)*r**2*dm, r)   # <r_perp^2>=(2/3)r^2 averaged over a shell
print(f"moment of inertia I    = {I_num:.6f} M R^2   (this is the 2/5, a different integral)")
print(f"analytic 2/5           = {2/5:.6f}")
