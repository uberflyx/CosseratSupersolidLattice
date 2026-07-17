"""Energy bookkeeping for gravitational waves in the Cosserat lattice.

Question: which lattice channel carries the energy of a gravitational wave,
and does the h-field's effective stiffness reproduce the GR flux normalisation
c^3 omega^2 h^2 / (16 pi G)?

Three candidate channels for a wave of measured metric strain h at frequency omega:
  (a) bare microrotation (rocking) channel: energy density ~ (1/2) rho j omega^2 phi^2
  (b) bare shear strain channel:            energy density ~ (1/2) mu h^2
  (c) collective (19-node tunnelling) channel with stiffness c^4/(16 pi G),
      i.e. the Einstein-Hilbert normalisation, which the lattice G-derivation
      should supply as mu / (16 pi K alpha^19).

The script tests the identity  c^2 / (16 pi G rho l^2)  ==  1 / (16 pi K alpha^19)
which holds iff hbar/(m0 c) == l == r_e, i.e. the node Compton wavelength equals
the classical electron radius. It then evaluates saturation of channels (a),(b)
near a merger, and computes the Hulse-Taylor / double-pulsar GR decay numbers
as the target for the quadrupole derivation.

CODATA 2022 / PDG values from the project files.
"""

import numpy as np

# --- constants (CODATA 2022) -------------------------------------------------
c      = 2.99792458e8          # m/s
G_exp  = 6.67430e-11           # m^3 kg^-1 s^-2
hbar   = 1.054571817e-34       # J s
alpha  = 7.2973525643e-3       # fine-structure constant
m_e    = 9.1093837139e-31      # kg
r_e    = 2.8179403205e-15      # m, classical electron radius
Msun   = 1.98892e30            # kg

# --- lattice parameters ------------------------------------------------------
m0   = m_e / alpha             # node mass
ell  = r_e                     # lattice length scale
rho  = m0 / ell**3             # lattice density
mu   = rho * c**2              # shear modulus
j    = ell**2                  # microinertia
K    = (1 + 1/np.pi) * (1 - 17*alpha/18)   # prefactor in the G formula

G_lat = K * (hbar * c / m0**2) * alpha**19
print(f"lattice G   = {G_lat:.6e}   (exp {G_exp:.6e}, ratio {G_lat/G_exp:.6f})")

# --- identity check: node Compton wavelength vs r_e --------------------------
lam0 = hbar / (m0 * c)
print(f"hbar/(m0 c) = {lam0:.6e} m ;  r_e = {r_e:.6e} m ;  ratio = {lam0/r_e:.8f}")

# --- the stiffness identity --------------------------------------------------
lhs = c**2 / (16*np.pi * G_lat * rho * ell**2)     # required enhancement of gradient modulus
rhs = 1.0 / (16*np.pi * K * alpha**19)             # inverse tunnelling suppression
print(f"c^2/(16 pi G rho l^2) = {lhs:.6e}")
print(f"1/(16 pi K alpha^19)  = {rhs:.6e}   ratio = {lhs/rhs:.8f}")

# --- GR flux for reference wave ----------------------------------------------
def flux_GR(h, omega):
    """time-averaged flux, one polarisation, W/m^2"""
    return c**3 * omega**2 * h**2 / (32*np.pi*G_exp)

# --- channel energy densities for a wave of metric strain h ------------------
def e_rocking_max(omega, phi_max=1.0):
    return 0.5 * rho * j * omega**2 * phi_max**2   # J/m^3, saturated rocking

def e_shear(h):
    return 0.5 * mu * h**2

# GW150914-like numbers: h ~ 1e-21 at f=100 Hz at D=410 Mpc; h ~ 0.1 near source
for label, h, f, note in [("at Earth", 1e-21, 100.0, ""),
                          ("near merger", 0.1, 300.0, "(r ~ 100 km)")]:
    om = 2*np.pi*f
    F  = flux_GR(h, om)
    e  = F / c
    print(f"\n{label} {note}: h={h:.1e}, f={f:.0f} Hz")
    print(f"  GR energy density      = {e:.3e} J/m^3")
    print(f"  rocking channel MAX    = {e_rocking_max(om):.3e} J/m^3  "
          f"-> {'SATURATES, cannot carry it' if e_rocking_max(om) < e else 'ok'}")
    print(f"  bare shear at strain h = {e_shear(h):.3e} J/m^3  "
          f"(ratio to GR: {e_shear(h)/e:.3e})")

# frequency at which the BARE channels happen to cross the GR normalisation
om_star = np.sqrt(16*np.pi*G_exp*rho)
print(f"\nbare-channel crossing frequency f* = {om_star/2/np.pi:.1f} Hz "
      "(mid LIGO band: why the l-gradient dictionary looked plausible)")

# --- Hulse-Taylor and double pulsar: the target numbers ----------------------
def pbdot_GR(m1, m2, Pb, e):
    """Peters (1964) orbital period decay for an eccentric binary."""
    M  = m1 + m2
    mu_r = m1*m2 / M
    n  = 2*np.pi / Pb
    a  = (G_exp*M / n**2)**(1/3)
    fe = (1 + 73/24*e**2 + 37/96*e**4) / (1 - e**2)**3.5
    # dE/dt = -32/5 G^4 m1^2 m2^2 M / (c^5 a^5) * f(e);  Pbdot = -(3/2)(Pb/E)... standard:
    return -192*np.pi/5 * (2*np.pi*G_exp/Pb)**(5/3) * m1*m2 / M**(1/3) / c**5 * fe

# PSR B1913+16 (Weisberg & Huang 2016): m1=1.438, m2=1.390 Msun, Pb=27906.98 s, e=0.6171
pb1913 = pbdot_GR(1.438*Msun, 1.390*Msun, 27906.980895, 0.6171340)
print(f"\nPSR B1913+16: Pbdot_GR = {pb1913:.4e}  (published GR value -2.4025e-12; "
      f"observed/GR = 0.9983 +/- 0.0016)")
# PSR J0737-3039 (Kramer et al 2021): m1=1.338185, m2=1.248868, Pb=8834.534998 s, e=0.087777
pb0737 = pbdot_GR(1.338185*Msun, 1.248868*Msun, 8834.534998, 0.0877775)
print(f"PSR J0737-3039A/B: Pbdot_GR = {pb0737:.4e}  (published -1.247920e-12; "
      f"observed/GR = 0.999963 +/- 0.000063)")

# =============================================================================
# Symbolic verifications promised in the monograph text
# =============================================================================
import sympy as sp

# --- (1) Planck-force identity: c^4/G == mu l^2 / (K alpha^19) ---------------
lhs_pf = c**4 / G_lat
rhs_pf = mu * ell**2 / (K * alpha**19)
print(f"\nPlanck force c^4/G      = {lhs_pf:.6e} N")
print(f"mu l^2 / (K alpha^19)   = {rhs_pf:.6e} N   ratio = {lhs_pf/rhs_pf:.10f}")

# --- (2) TT-projector sphere average on symmetric traceless tensors ----------
# <Lambda_ij,kl A_ij A_kl> = (2/5) A_ij A_ij, verified by symbolic integration.
th, ph = sp.symbols('theta phi', real=True)
n = sp.Matrix([sp.sin(th)*sp.cos(ph), sp.sin(th)*sp.sin(ph), sp.cos(th)])
# generic symmetric traceless tensor
a11, a12, a13, a22, a23 = sp.symbols('a11 a12 a13 a22 a23', real=True)
A = sp.Matrix([[a11, a12, a13], [a12, a22, a23], [a13, a23, -a11-a22]])
P = sp.eye(3) - n*n.T
integrand = sp.S(0)
for i in range(3):
    for j_ in range(3):
        for k in range(3):
            for l in range(3):
                Lam = P[i,k]*P[j_,l] - sp.Rational(1,2)*P[i,j_]*P[k,l]
                integrand += Lam * A[i,j_] * A[k,l]
avg = sp.integrate(sp.integrate(integrand*sp.sin(th), (th, 0, sp.pi)),
                   (ph, 0, 2*sp.pi)) / (4*sp.pi)
AA = sum(A[i,j_]**2 for i in range(3) for j_ in range(3))
ratio = sp.simplify(avg / AA)
print(f"sphere average of TT projector on traceless symmetric tensors: {ratio} (expect 2/5)")

# --- (3) Peters-Mathews eccentricity factor f(e) by orbit average ------------
# Instantaneous quadrupole luminosity on a Keplerian ellipse, averaged over one
# period, over the circular-orbit result: f(e) = (1 + 73/24 e^2 + 37/96 e^4)/(1-e^2)^{7/2}
ecc, nu = sp.symbols('e nu', positive=True)
# Peters (1964) Eq. (5.4): dE/dt as a function of true anomaly nu, in units of
# the circular-orbit luminosity at semi-major axis a:
# relative to the circular-orbit luminosity at the same semi-major axis a
dEdt = (1 + ecc*sp.cos(nu))**4 / (1 - ecc**2)**5 * \
       ((1 + ecc*sp.cos(nu))**2 + sp.Rational(1,12)*ecc**2*sp.sin(nu)**2)
# time average: dt = (1-e^2)^{3/2}/(1+e cos nu)^2 dnu / (2 pi) over one period
favg = sp.integrate(dEdt * (1 - ecc**2)**sp.Rational(3,2)
                    / (1 + ecc*sp.cos(nu))**2, (nu, 0, 2*sp.pi)) / (2*sp.pi)
favg = sp.simplify(sp.expand_trig(favg))
target = (1 + sp.Rational(73,24)*ecc**2 + sp.Rational(37,96)*ecc**4) \
         / (1 - ecc**2)**sp.Rational(7,2)
print(f"orbit-averaged f(e) - target = {sp.simplify(favg - target)} (expect 0)")
