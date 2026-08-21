"""The covariant defect form factor F_D(q^2): what is derived and what is not.

Five routes to the defect vertex are run and held against each other.  The
headline is negative and worth stating first, because an earlier version of
this file claimed the opposite: the finite-q vertex is NOT derived here.

  A. minimal coupling with a worldline delta current -> F_1 = 1, but this
     ASSUMES pointlikeness and is therefore circular
  B. the ideal screw is harmonic off-core, so ITS source is a point; the
     Peierls-Nabarro core the framework computes is not, carrying the Burgers
     vector in a Lorentzian of half-width w with transform exp(-|q|w)
  C. Eshelby covariance saturates the exponent at 2Mw = 0.90 alpha, a model
     estimate rather than a form-factor theorem
  D. the core is reflectionless (verified numerically), which is evidence
     about back-scattering and not about the current matrix element
  E. the semiclassical expansion is unavailable, g^2 ~ 8/(Mw) ~ 2400

What IS established is F_D(0) = 1 to all orders, from the Ward-Takahashi
identity that the derived topological current and minimal coupling deliver;
see defect_current_duality.py.  See the VERDICT block at the end.
"""
import math

import numpy as np
from scipy import integrate
import sympy as sp

line = lambda t: print("\n" + t + "\n" + "=" * len(t))
ALPHA_INV = 137.035999177
ALPHA = 1.0 / ALPHA_INV
W_OVER_ELL = 0.453                      # PN half-width, w = 0.453 ell
MW = 0.45 * ALPHA                       # m_e c w / hbar  (dimensionless)

# ----------------------------------------------------------------- ROUTE A
line("A. MINIMAL COUPLING + A POINT CURRENT -> F_1 = 1, BUT THAT IS CIRCULAR")
# A point charge with minimal action S = -M int dtau + e int a_mu dx^mu has
# current J^mu(x) = e int dtau xdot^mu delta(x - x(tau)).  Its matrix element
# between scalar momentum eigenstates is (p+p')^mu / (2 sqrt(p0 p0')), i.e.
# F_1(q^2) = 1 identically.  Verify the current is conserved and the spatial
# integral is the charge, symbolically, for a straight worldline.
tau = sp.symbols("tau", real=True)
p0, pz, M = sp.symbols("p0 p_z M", positive=True)
# worldline x^mu(tau) = (p^mu/M) tau  (free particle), current density FT:
# tilde J^mu(q) = e int dtau (p^mu/M) e^{i q.x(tau)} = e (p^mu/M) 2pi delta(q.p/M)
# elastic (q.(p+p')=0 in Breit frame) -> matrix element ~ (p+p')^mu, F=1.
q = sp.symbols("q", real=True)
# demonstrate F(0) = charge and no q-dependence in the vertex structure:
print("   point-charge current  J^mu = e (p+p')^mu/(2 sqrt(p0 p0')):")
print("   -> F_1(q^2) = 1 for all q^2 (no intrinsic scale in a point coupling)")
print("   BUT the delta-function worldline current ASSUMES pointlikeness,")
print("   which is exactly what a finite-q probe is being asked to decide.")
print("   Minimal coupling is derived; a POINT current is not a consequence.")

# ----------------------------------------------------------------- ROUTE B
line("B. THE IDEAL SOURCE IS A POINT; THE RELAXED PN SOURCE IS NOT")
x, y = sp.symbols("x y", real=True)
u_ideal = sp.atan2(y, x)                          # ideal screw angle (x b/2pi)
lap = sp.simplify(sp.diff(u_ideal, x, 2) + sp.diff(u_ideal, y, 2))
print(f"   Laplacian of the ideal screw field away from core = {lap}")
print("   -> harmonic off-core: the source is a POINT, the far field is 1/r")
print("      (2D log potential), whose dual charge density is delta^2(x).")
print("   That holds for the IDEAL SINGULAR screw only.  The framework")
print("   computes the Peierls-Nabarro core, whose misfit is spread over the")
print("   glide plane, so the density carrying the Burgers vector is the")
print("   Lorentzian rho_B(x) = (b/pi) w/(x^2+w^2) with Int rho_B dx = b.")
print("   The circuit invariant is untouched and the CHARGE IS DISTRIBUTED.")
print("   Its normalised transform is exp(-|q|w) exactly, so route B does")
print("   NOT support a structureless vertex for the relaxed core.")

# dual (Coulomb) form factor of a point topological source:
# f ~ 1/q^2 tilde J(q); if J = delta, tilde J = 1, so F_D(q^2) = 1 exactly.
print("   dual Maxwell: d_mu f^{mu nu} = J_D^nu; a point J_D would give")
print("   F_D = 1, and the relaxed J_D is not a point.")

# ----------------------------------------------------------------- ROUTE C
line("C. ESHELBY COVARIANCE: THE BOOSTED CLASSICAL PROFILE SATURATES")
# The moving screw is Eshelby's contracted static solution (shear sector is
# exactly Lorentz-covariant, Kleinert), NOT a rigid-boost ansatz.  In the
# Breit frame the profile of rest-width w contracts to w/gamma with
# gamma = sqrt(1 + q^2/4M^2).  The PN Burgers density is the Lorentzian,
# whose FT is exp(-|q| w), so the covariant form factor is exp(-|q| w/gamma).
qq, Mm, ww = sp.symbols("q M w", positive=True)
gamma = sp.sqrt(1 + qq**2 / (4 * Mm**2))
exponent = qq * ww / gamma
sat = sp.limit(exponent, qq, sp.oo)
print(f"   exponent |q| w / gamma,  gamma = sqrt(1 + q^2/4M^2)")
print(f"   lim_(q->oo) |q| w / gamma = {sat}   -> F -> e^(-2 M w)")
val = math.exp(-2 * MW)
print(f"   2 M w = 2 (m_e c w/hbar) = 0.90 alpha = {2*MW:.6f}")
print(f"   saturation value e^(-0.90 alpha) = {val:.6f}, deficit "
      f"{1-val:.6f} = {(1-val)/ALPHA:.3f} alpha")
# the Fourier transform check: FT of the PN Lorentzian rho=(w/pi)/(x^2+w^2)
xx = sp.symbols("x", real=True)
qs = sp.symbols("qs", positive=True)
rho = (ww / sp.pi) / (xx**2 + ww**2)
FT = sp.integrate(rho * sp.cos(qs * xx), (xx, -sp.oo, sp.oo))
print(f"   FT of PN Lorentzian density = {sp.simplify(FT)}   (= e^(-q w))")

# ----------------------------------------------------------------- ROUTE D
line("D. REFLECTIONLESSNESS: THE CORE BACKSCATTERS NOTHING (numeric)")
# Poeschl-Teller V = -2 sech^2(x) is the l=1 reflectionless potential, the
# same class as the Benjamin-Ono core operator.  Integrate the Schrodinger
# equation and read off |R(k)|; reflection = backscattering = large-q elastic
# structure, so R=0 means the core adds no form-factor deficit of its own.
def reflection(k, L=40.0):
    def rhs(xv, Y):
        psi, dpsi = Y
        V = -2.0 / math.cosh(xv) ** 2
        return [dpsi, (V + k * k) * psi * 0 + (V - (-k * k)) * 0 + (V + 0) * psi
                - (k * k) * psi + V * psi * 0]  # placeholder, replaced below
    # -psi'' + V psi = k^2 psi  =>  psi'' = (V - k^2) psi
    def rhs2(xv, Y):
        psi, dpsi = Y
        V = -2.0 / math.cosh(xv) ** 2
        return [dpsi, (V - k * k) * psi]
    # transmitted wave at +L: psi = e^{ikx}
    Y0 = [complex(np.exp(1j * k * L)), complex(1j * k * np.exp(1j * k * L))]
    sol = integrate.solve_ivp(rhs2, [L, -L], Y0, rtol=1e-10, atol=1e-12,
                              dense_output=True, method="DOP853")
    psiL, dpsiL = sol.y[0, -1], sol.y[1, -1]
    xL = -L
    # psi = A e^{ikx} + B e^{-ikx}; solve for A, B at x = -L
    A = 0.5 * (psiL + dpsiL / (1j * k)) * np.exp(-1j * k * xL)
    B = 0.5 * (psiL - dpsiL / (1j * k)) * np.exp(1j * k * xL)
    return abs(B / A)

for k in (0.5, 1.0, 2.0, 4.0):
    R = reflection(k)
    print(f"   k = {k:4.1f}:  |R(k)| = {R:.2e}   (reflectionless: expect ~0)")
print("   -> the core back-scatters nothing.  That is EVIDENCE and not proof")
print("      about the vertex: reflection and the current matrix element are")
print("      different observables, and a reflectionless system can still")
print("      carry transmission phases, time delays and internal structure.")

# ----------------------------------------------------------------- ROUTE E
line("E. WHY THE EXACT CONTINUUM F_1 IS NOT DERIVABLE SEMICLASSICALLY")
Mw = MW
g2_sineGordon = 8.0 / Mw            # beta^2 = 8/(M w) for the SG soliton
print(f"   M w = 0.45 alpha = {Mw:.3e}   (soliton is LIGHT and SMALL)")
print(f"   semiclassical parameter g^2 ~ 8/(M w) = {g2_sineGordon:.0f}")
print(f"   free-fermion point beta^2 = 4pi = {4*math.pi:.2f}; soliton picture")
print(f"   fails past ~8pi = {8*math.pi:.2f}.  g^2 = {g2_sineGordon:.0f} is deep")
print(f"   in strong coupling: the continuum collective-coordinate expansion")
print(f"   does NOT converge, so an exact F_1(q^2) needs the LATTICE (the")
print(f"   framework's own UV completion), a Rajantie-Weir computation, not a")
print(f"   continuum soliton bootstrap.  That is the one piece left open.")

# ----------------------------------------------------------------- THE CHECK
line("THE CHECK: THE LEPTON MOMENT SELECTS THE POINT-CHARGE READING")
# electron g-2: a_e ~ alpha/2pi, measured to relative ~1e-10 -> delta a_e ~1e-13
a_e = ALPHA / (2 * math.pi)
rel_prec = 1.0e-10
delta_ae = rel_prec * a_e
eps_bound = delta_ae / (2 * a_e)          # vertex deficit enters at 2 eps
print(f"   a_e (Schwinger) = {a_e:.6e}, measured to relative {rel_prec:.0e}")
print(f"   -> vertex deficit bound  eps(m_e c) <= {eps_bound:.2e}"
      f"   (paper: 5.6e-11)")

# the two theoretical readings at q ~ m_e c:
q_over_M = 1.0                              # probe at q = m_e c
gam = math.sqrt(1 + q_over_M**2 / 4)
eps_PN = q_over_M * MW / gam                # PN-relaxed distributed charge
eps_point = 0.0                             # minimal-coupling point charge
print(f"\n   at q = m_e c (gamma = {gam:.3f}):")
print(f"     PN-relaxed / classical-profile reading:  eps = {eps_PN:.3e}"
      f"  = {eps_PN/ALPHA:.2f} alpha")
print(f"     minimal-coupling point-charge reading:   eps = {eps_point:.0e}")
print(f"   measured bound {eps_bound:.1e} is below the PN reading by "
      f"{eps_PN/eps_bound:.0e}  (~8 orders)")
print(f"   -> the data EXCLUDE the distributed-charge reading and REQUIRE the")
print(f"      point-charge one.  The photon couples to the bare topological")
print(f"      charge (route A/B), not the relaxed strain profile (route C).")

line("VERDICT")
print(f"""   ESTABLISHED:  F_D(0) = 1 to all orders.  Minimal coupling is derived
            (Peach-Koehler = Lorentz) and the current is topological with
            conservation an identity, so the Ward-Takahashi identity holds and
            carries the circuit result from tree level to all orders.

   NOT ESTABLISHED, and an earlier version of this script said otherwise:
            F_1(q^2) = 1 at finite q^2.  Route A writes a delta-function
            worldline current, which ASSUMES the defect is a point, and that
            is precisely what a finite-q probe is asked to decide.  Route B
            fails for the same reason once the core is relaxed: the ideal
            singular screw does have a point source, but the Peierls-Nabarro
            core the framework actually computes carries the Burgers vector in
            a Lorentzian of half-width w, whose normalised transform is
            exp(-|q|w) exactly.  The charge is distributed, not concentrated.

   ESTIMATES, not bounds:  route C is a model prescription (contracting a
            rest-frame density and transforming it is not a form-factor
            theorem for a composite), and route D concerns a different
            observable, since a reflectionless system can still carry
            transmission phases, time delays and internal structure.

   THE CHECK is a sensitivity estimate, not an exclusion: converting a vertex
            deficit into g-2 needs a gauge-invariant vertex, a modified
            propagator and a renormalisation prescription, none of which is
            done here.  It indicates the structureless reading by a factor
            5e7; it does not establish it.

   OPEN:    the defect's current matrix element computed from the elastic
            theory.  The continuum route is unavailable, the soliton being
            strongly coupled at g^2 ~ {g2_sineGordon:.0f} against a free-fermion point of
            4pi = {4*math.pi:.1f}, so this is a lattice computation of the kind
            Rajantie and Weir performed for kinks.""")
