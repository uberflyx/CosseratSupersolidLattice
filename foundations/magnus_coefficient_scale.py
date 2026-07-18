"""A3: the physical Magnus coefficient gamma and the phonon-edge energy.

Context. The core's transverse fluctuation doublet (glide leg a, stacking leg b)
carries a gyroscopic Magnus term because the electron line is also a condensate
vortex. In circular combinations psi_pm = psi_a -+ i psi_b the coupled pair
decouples into an omega-dependent core well/hill:

    (H + sign * gamma * omega * beta(x)) psi_pm = omega^2 psi_pm,   sign = -+1,

with H = |k| + U(x) the established core-fluctuation operator (continuum edge
normalised to 1), beta(x) = 1/(1 + (x/w)^2) the dimensionless vortex-circulation
localisation profile (unit peak at the core), and gamma the dimensionless Magnus
coefficient this script derives.

Two deliverables:
  (1) gamma = omega_M / omega_edge, the ratio of the vortex Magnus (cyclotron)
      gyration frequency to the edge frequency; parameter-free up to the vortex
      core-mass convention (the standard log ambiguity), gamma = 2 f_s ell / w.
  (2) the edge energy in MeV under the SHEAR-WAVE reading hbar c / w. CAVEAT
      (honest, unresolved): the core operator's own dispersion is
      omega^2 ~ (mu_bar |k| + V'') / rho_eff, and until rho_eff is derived from
      the PN energy functional the MeV anchor is provisional. The dimensionless
      structure (order-one gamma, one bound circular sense) does not depend on
      the anchor; every MeV below does.

Then the Magnus scan is extended to the derived gamma range to read the depth
of the bound circular mode.

All framework relations used are already in the monograph:
  kappa = h/m0 = 2*pi*c*ell         (circulation quantum, sec:condensate_vortices)
  rho_s = (4/5) rho                 (superfluid fraction, sec:fs_bootstrap)
  rho   = m0 / ell^3                 (working mass density)
  hbar c / ell = m0 c^2 = m_e c^2 / alpha   (bootstrap; ell = r_e)
  continuum edge of the core operator at |k| = 1/w   (sec:core_fluctuations)
  transverse fluctuations are shear waves, omega = c k  (shear branch)
"""

import numpy as np
import sympy as sp

# ----------------------------------------------------------------------
# Part 1: symbolic reduction of gamma and the edge frequency
# ----------------------------------------------------------------------
c, ell, w, m0, h, hbar, fs = sp.symbols('c ell w m0 h hbar f_s', positive=True)

kappa = h / m0                 # circulation quantum
kappa = kappa.subs(h, 2*sp.pi*hbar)
kappa = kappa.subs(hbar, m0*c*ell)          # bootstrap hbar = m0 c ell -> kappa = 2 pi c ell
rho = m0 / ell**3
rho_s = fs * rho

# vortex Magnus (cyclotron) frequency omega_M = rho_s kappa / m_star,
# bare core-tube mass per unit length m_star = pi rho w^2
m_star = sp.pi * rho * w**2
omega_M = sp.simplify(rho_s * kappa / m_star)

# transverse phonon edge: the core operator's continuum edge sits at |k| = 1/w;
# shear branch omega = c k  =>  omega_edge = c / w
omega_edge = c / w

gamma = sp.simplify(omega_M / omega_edge)

print("=== Part 1: symbolic ===")
print("kappa           =", kappa, "   (= 2 pi c ell)")
print("omega_M         =", omega_M, "   (bare core mass m* = pi rho w^2)")
print("omega_edge      =", omega_edge)
print("gamma           =", gamma, "   (= 2 f_s ell / w)")
print()

# ----------------------------------------------------------------------
# Part 2: numbers
# ----------------------------------------------------------------------
mec2 = 0.51099895069        # MeV, CODATA 2022
alpha_inv = 137.035999177   # CODATA 2022
m0c2 = mec2 * alpha_inv     # MeV; also = hbar c / ell
fs_val = 4.0/5.0

w_par = 0.452               # glide leg, units of ell
w_perp = 0.639              # stacking leg, units of ell

print("=== Part 2: numbers ===")
print(f"m0 c^2 = m_e c^2 / alpha = hbar c / ell = {m0c2:.4f} MeV")
for name, wv in [("glide (a)", w_par), ("stacking (b)", w_perp)]:
    E_edge = m0c2 / wv                       # hbar c / w = (hbar c/ell)/(w/ell)
    g_bare = 2*fs_val / wv                    # gamma = 2 f_s ell / w, bare core mass
    print(f"  {name:14s} w/ell={wv:.3f}   edge hbar c/w = {E_edge:7.2f} MeV   "
          f"gamma_bare = {g_bare:5.3f}")

# log-enhanced vortex mass: m* = pi rho w^2 [1 + 2 ln(R/w)], R = L4 = sqrt(6) ell
L4 = np.sqrt(6.0)
print("\n  vortex-mass log correction, cut at R = L4 = sqrt(6) ell:")
for name, wv in [("glide (a)", w_par), ("stacking (b)", w_perp)]:
    Lam = 1.0 + 2.0*np.log(L4/wv)            # backflow enhancement factor
    g_log = (2*fs_val/wv) / Lam
    print(f"  {name:14s} ln(R/w)={np.log(L4/wv):.3f}  enhancement={Lam:.3f}  "
          f"gamma_log = {g_log:5.3f}")

# ----------------------------------------------------------------------
# Part 3: Magnus scan extended to the derived gamma range; read E_minus
# (reuses the established core operator H = |k| + (x^2-1)/(x^2+1), edge at 1)
# ----------------------------------------------------------------------
N, Lbox = 1536, 250.0
x = (np.arange(N) - N/2) * (2*Lbox/N)
dx = x[1]-x[0]
k = np.fft.fftfreq(N, d=dx)*2*np.pi
Fm = np.fft.fft(np.eye(N), axis=0)
Lam_op = (np.fft.ifft(np.abs(k)[:,None]*Fm, axis=0)).real
Lam_op = 0.5*(Lam_op+Lam_op.T)
Hc = Lam_op + np.diag((x**2-1.0)/(x**2+1.0))
evc, vecc = np.linalg.eigh(Hc)
psi_th = vecc[:,1]                     # odd threshold state
beta = 1.0/(x**2+1.0)                  # dimensionless localisation, unit peak

def circular_branch(gamma, sign, E0=1.0, iters=60):
    E = E0
    for _ in range(iters):
        womega = np.sqrt(max(E, 1e-12))
        Heff = Hc + sign*gamma*womega*np.diag(beta)
        ev, vec = np.linalg.eigh(Heff)
        ov = (vec.T @ psi_th)**2
        j = np.argmax(ov)
        Enew = ev[j]
        if abs(Enew-E) < 1e-11:
            E = Enew; break
        E = 0.5*(E+Enew)
    return E, ov[j]

print("\n=== Part 3: bound circular mode across the derived gamma range ===")
print(" gamma   E_minus   sqrt(E_-)   ov      -> bound-mode energy (MeV), glide / stacking")
for gamma in [0.4, 0.57, 0.66, 0.81, 1.6, 2.5, 3.54]:
    Em, ov = circular_branch(gamma, -1)
    s = np.sqrt(max(Em, 0.0))
    Eg = s * (m0c2/w_par)
    Es = s * (m0c2/w_perp)
    print(f" {gamma:4.2f}   {Em:+.4f}   {s:.4f}    {ov:.3f}    {Eg:6.1f} / {Es:6.1f}")
