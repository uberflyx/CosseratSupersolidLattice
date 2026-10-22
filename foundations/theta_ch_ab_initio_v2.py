#!/usr/bin/env python3
"""
theta_ch_ab_initio_v2.py: the checks behind the lattice chirality theta_ch = alpha^2/(2 pi)
=========================================================================================
Companion to "The Cosserat Supersolid" monograph (M. A. Cox, 2026).
Repository: https://github.com/uberflyx/CosseratSupersolidLattice

What it computes
----------------
Route A. The static microrotation of a screw dislocation.  The microrotation
  amplitude G(r) solves the Cosserat equilibrium

      G'' + G'/r - (1/r^2 + q^2) G = -q^2 f(r),     f(r) = r/(r^2 + w^2),

  with q^2 = 2 kappa_c/gamma.  The strain profile f is the macrorotation,
  curl u = (b/2 pi) f(r) rhat, and the microrotation is phi = (b/4 pi) G(r) rhat,
  so G = f is full tracking (phi = curl(u)/2).  The equation is solved by its
  Green's function,

      G(r) = q^2 [ K1(qr) int_0^r I1(qs) f(s) s ds + I1(qr) int_r^inf K1(qs) f(s) s ds ],

  which is regular at the origin and decays to f at infinity, and cross-checked
  with a boundary-value solver.  It is evaluated at the physical core half-width
  w = W d, with W = ln(1/alpha)/(2 pi) and d = l/sqrt(3) the misfit period, so
  w = 0.45211 l.  Reported: the tracking ratio G/f at the origin, at 1/q and at
  3w; the radius where G/f reaches 0.9; where the mismatch norm int (f-G)^2 r dr
  sits; the virial identity; the overlap angle with two core weights; and the
  static energy partition.

  Energy partition.  The antiplane energy density is
      (mubar/2)|grad u3|^2 + (kappa_c/4)|curl u - 2 phi|^2 + (gamma/2)|grad phi|^2,
  with mubar = mu + kappa_c/2 the true shear modulus.  At the stationary phi the
  rotational part (the last two terms) equals (kappa_c/4) int curl u.(curl u - 2 phi),
  so per unit length it is E_rot = (mubar b^2/4 pi) lambda_rot with
      lambda_rot = (kappa_c/(2 mubar)) int_0^inf f (f - G) r dr,
  a finite number, while the classical part within radius R is
      E_cl(R) = (mubar b^2/4 pi) (1/2)[ln(1 + R^2/w^2) - R^2/(R^2 + w^2)],
  which grows as ln(R/w).

Route B. The cross-coupling term -alpha^2/pi of the alpha equation, tested by
  removing it from the complete equation (Scheme A on-site resummation plus the
  inter-valley term summed over all shells) and solving self-consistently.

Route C. The index-structure theorem (the static hemitropic overlap
  eps_(ij) kappa_(ij) vanishes identically for the anti-plane screw), checked by
  building both tensors on a grid, and the constitutive normalisation
  theta_ch = alpha^2 N^2/2 = alpha^2/(2 pi).

Units: lattice spacing l = 1, Eringen shear modulus mu = 1, gamma = mu l^2,
kappa_c = 2 mu/(pi - 2) (the rolling coupling N^2 = 1/pi).

Requires numpy and scipy.
"""

import numpy as np
from scipy.integrate import solve_bvp, quad
from scipy.special import ive, kve, kv
from scipy.optimize import brentq

# -- Constants ---------------------------------------------------------------
ALPHA_INV_EXP = 137.035999177          # CODATA 2022
alpha_exp = 1.0 / ALPHA_INV_EXP
N2 = 1.0 / np.pi
mu = 1.0                               # Eringen shear modulus
kappa_c = 2 * N2 * mu / (1 - 2 * N2)   # = 2 mu/(pi - 2)
gamma_c = mu * 1.0                     # gamma = mu l^2
mu_bar = mu + kappa_c / 2              # true (long-range) shear modulus
q = np.sqrt(2 * kappa_c / gamma_c)     # screening wavevector
d = 1.0 / np.sqrt(3.0)                 # misfit period, partial hop
W = np.log(ALPHA_INV_EXP) / (2 * np.pi)
w = W * d                              # physical core half-width
m_e_MeV = 0.51099895069                # CODATA 2022

sep = "=" * 72
sep2 = "-" * 72

print(sep)
print("  CHECKS BEHIND theta_ch = alpha^2/(2 pi)")
print("  M. A. Cox, The Cosserat Supersolid (2026)")
print(sep)
print(f"\n  alpha^-1 = {ALPHA_INV_EXP}   N^2 = 1/pi = {N2:.10f}")
print(f"  kappa_c/mu = 2/(pi-2) = {kappa_c:.6f}   mubar/mu = {mu_bar:.6f}   gamma = {gamma_c:.1f} mu l^2")
print(f"  q = sqrt(2 kappa_c/gamma) = {q:.6f}/l     1/q = {1/q:.6f} l")
print(f"  w = W d = {W:.6f} x {d:.6f} l = {w:.6f} l     qw = {q*w:.6f}")


# ==========================================================================
# ROUTE A: the static microrotation profile
# ==========================================================================
print(f"\n{sep}\nROUTE A: Cosserat microrotation profile G(r)\n{sep}")


def f_strain(r):
    return r / (r**2 + w**2)


def green(r):
    """G(r) and G'(r) from the Green's function (exponentially scaled Bessels)."""
    a = quad(lambda s: ive(1, q*s) * np.exp(q*(s - r)) * f_strain(s) * s,
             0.0, r, limit=400, epsabs=1e-15, epsrel=1e-12)[0]
    b = quad(lambda s: kve(1, q*s) * np.exp(-q*(s - r)) * f_strain(s) * s,
             r, np.inf, limit=400, epsabs=1e-15, epsrel=1e-12)[0]
    x = q * r
    G = q*q * (kve(1, x)*a + ive(1, x)*b)
    dK1 = q * (-kve(0, x) - kve(1, x)/x)     # d/dr K1(qr), scaled
    dI1 = q * (ive(0, x) - ive(1, x)/x)      # d/dr I1(qr), scaled
    Gp = q*q * (dK1*a + dI1*b)
    return G, Gp


ratio = lambda r: green(r)[0] / f_strain(r)

# Origin: G/f -> (a^2/2) int_0^inf K1(x) x^2/(x^2 + a^2) dx with a = qw
aw = q * w
ratio0_closed = aw**2 / 2 * quad(lambda x: kv(1, x) * x*x / (x*x + aw*aw), 0, np.inf, limit=400)[0]
print(f"\n  G/f at the origin: closed integral {ratio0_closed:.4f}, Green's function at r = 1e-4 l: {ratio(1e-4):.4f}")
print(f"  G/f at r = 1/q = {1/q:.3f} l: {ratio(1/q):.4f}")
print(f"  G/f at r = 3w  = {3*w:.3f} l: {ratio(3*w):.4f}")
r90 = brentq(lambda r: ratio(r) - 0.9, 0.2, 10.0, xtol=1e-12)
print(f"  G/f = 0.9 at r = {r90:.4f} l = {r90/w:.3f} w = {r90*q:.3f}/q")

# Grid for the integrals
R_out = 60.0
r = np.concatenate([np.geomspace(1e-5, 1e-2, 200)[:-1], np.linspace(1e-2, 4.0, 1600),
                    np.linspace(4.0, R_out, 1200)[1:]])
GG = np.array([green(ri) for ri in r])
G, Gp = GG[:, 0], GG[:, 1]
f = f_strain(r)
tr = np.trapezoid

# Boundary-value cross-check: regular at the origin (G ~ r), G -> f at R_out
ra = 1e-4
def ode(rr, y):
    return np.vstack([y[1], -y[1]/rr + (1/rr**2 + q*q)*y[0] - q*q*f_strain(rr)])
def bc(ya, yb):
    return np.array([ya[1] - ya[0]/ra, yb[0] - f_strain(R_out)])
r_init = np.concatenate([np.linspace(ra, 4.0, 400), np.linspace(4.05, R_out, 200)])
sol = solve_bvp(ode, bc, r_init, np.vstack([f_strain(r_init), np.gradient(f_strain(r_init), r_init)]),
                tol=1e-8, max_nodes=200000)
dev = max(abs(sol.sol(rr)[0] / green(rr)[0] - 1) for rr in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0])
print(f"\n  Boundary-value solver: converged {sol.success}; largest relative deviation from the")
print(f"  Green's function at r = 0.01 to 5 l: {dev:.1e}")

# Virial identity: gamma int (G'^2 + G^2/r^2) r dr = 2 kappa_c int (f G - G^2) r dr
I_curv = tr((Gp**2 + (G/r)**2) * r, r) + 1.0 / R_out**2      # analytic tail, G -> 1/r
I_GfmG = tr(G * (f - G) * r, r)
vir = gamma_c * I_curv / (2 * kappa_c * I_GfmG) - 1
print(f"  Virial identity: relative error {vir:.1e}")

# Mismatch localisation
mis = (f - G)**2 * r
I_mis = tr(mis, r)
for k in (1, 2, 3):
    m = r <= k / q
    print(f"  Mismatch norm int (f-G)^2 r dr inside {k}/q = {k/q:.3f} l: {tr(mis[m], r[m]) / I_mis * 100:.1f}%")

# Core-weighted overlap angle (the global cosine diverges with the box)
def sin2(weight):
    c = tr(weight*f*G*r, r) / np.sqrt(tr(weight*f*f*r, r) * tr(weight*G*G*r, r))
    return 1 - c*c
print(f"  Core-weighted sin^2(theta): misfit weight w^2/(r^2+w^2): {sin2(w*w/(r*r + w*w)):.3f};"
      f" unrelaxed-rotation weight 1 - G/f: {sin2(1 - G/f):.3f}")

# Energy partition
lam_identity = kappa_c / (2 * mu_bar) * quad(lambda s: f_strain(s) * (f_strain(s) - green(s)[0]) * s,
                                             0, np.inf, limit=400, epsabs=1e-13)[0]
lam_parts = (kappa_c / (2 * mu_bar)) * I_mis + (gamma_c / (4 * mu_bar)) * I_curv
print(f"\n  Rotational energy E_rot = (mubar b^2/4pi) lambda_rot:")
print(f"    lambda_rot from (kappa_c/2mubar) int f(f-G) r dr       = {lam_identity:.5f}")
print(f"    lambda_rot from lag + curvature terms separately       = {lam_parts:.5f}")
print(f"    kappa_c/(2 mubar) = 1/(pi-1) = {kappa_c/(2*mu_bar):.6f}")
for R in (20.0, 100.0, 1000.0):
    bracket = 0.5 * (np.log(1 + R*R/(w*w)) - R*R/(R*R + w*w))
    print(f"    R = {R:6.0f} l: E_cl bracket {bracket:.4f}, E_rot/E_cl = {lam_identity/bracket:.4f}")
print(f"  The static share falls as 1/ln R and is not the modulus ratio N^2 = {N2:.4f}.")


# ==========================================================================
# ROUTE B: the cross-coupling term in the alpha equation
# ==========================================================================
print(f"\n{sep}\nROUTE B: Cross-coupling in the alpha equation\n{sep}")

bare = np.exp(np.pi**2 / 2)
iv_all = lambda a: 6*a/np.pi**2 * (1/np.sqrt(1 - 4*a*a) - 1/(1 - a*a))   # all shells


def series_full(a):
    return bare - 2 - a - a/(np.pi*(1 - a)) - iv_all(a)


def series_no_cross(a):
    return bare - 2 - a - a/np.pi - iv_all(a)


a_full = brentq(lambda a: 1/a - series_full(a), 1e-4, 0.1, xtol=1e-20, rtol=4*np.finfo(float).eps)
a_noc = brentq(lambda a: 1/a - series_no_cross(a), 1e-4, 0.1, xtol=1e-20, rtol=4*np.finfo(float).eps)
res_full = (1/a_full - ALPHA_INV_EXP) / ALPHA_INV_EXP * 1e9
res_noc = (1/a_noc - ALPHA_INV_EXP) / ALPHA_INV_EXP * 1e9
shift = 1/a_noc - 1/a_full
print(f"\n  Complete equation: alpha^-1 = {1/a_full:.12f}  ({res_full:+.4f} ppb)")
print(f"  Without -alpha^2/pi and its higher powers: alpha^-1 = {1/a_noc:.12f}  ({res_noc:+.3f} ppb)")
print(f"  Shift = {shift:.4e};  alpha^2/pi = {a_full**2/np.pi:.4e};  difference {abs(shift/(a_full**2/np.pi) - 1)*100:.2f}%")
print(f"  Residual worsens by a factor {res_noc/res_full:.3g}")


# ==========================================================================
# ROUTE C: index structure and constitutive normalisation
# ==========================================================================
print(f"\n{sep}\nROUTE C: index structure and normalisation\n{sep}")

# Build u = (0, 0, u3(x, y)) and phi = (phi1, phi2, 0) for the screw on a grid and
# contract the Cosserat strain eps_ji = u_j,i - e_jik phi_k with kappa_ji = phi_j,i.
eps_lc = np.zeros((3, 3, 3))
eps_lc[0, 1, 2] = eps_lc[1, 2, 0] = eps_lc[2, 0, 1] = 1
eps_lc[0, 2, 1] = eps_lc[2, 1, 0] = eps_lc[1, 0, 2] = -1
h = 1e-4
rng = np.random.default_rng(1)
pts = rng.uniform(-3, 3, size=(200, 2))
pts = pts[np.hypot(pts[:, 0], pts[:, 1]) > 0.05]


def fields(x, y):
    rr = np.hypot(x, y)
    th = np.arctan2(y, x)
    u3 = (1/(2*np.pi)) * th * rr*rr / (rr*rr + w*w)          # a regularised screw, b = 1
    g = np.interp(rr, r, G) / (4*np.pi)                     # phi = (b/4pi) G rhat
    return np.array([0.0, 0.0, u3]), np.array([g*x/rr, g*y/rr, 0.0])


worst = 0.0
for x0, y0 in pts:
    grad_u = np.zeros((3, 3)); grad_phi = np.zeros((3, 3))   # [j, i] = d_i of component j
    for i, (dx, dy) in enumerate([(h, 0), (0, h)]):
        up, pp = fields(x0 + dx, y0 + dy); um, pm = fields(x0 - dx, y0 - dy)
        grad_u[:, i] = (up - um) / (2*h); grad_phi[:, i] = (pp - pm) / (2*h)
    _, phi0 = fields(x0, y0)
    eps = grad_u - np.einsum('jik,k->ji', eps_lc, phi0)       # z-derivatives vanish
    kap = grad_phi
    sym = lambda A: 0.5*(A + A.T)
    val = np.sum(sym(eps) * sym(kap))
    scale = np.sqrt(np.sum(eps**2) * np.sum(kap**2))
    worst = max(worst, abs(val) / scale)
print(f"\n  eps_(ij) kappa_(ij) / (|eps| |kappa|), largest over {len(pts)} points: {worst:.1e}")
print("  The strain has a 3 on one index, the curvature has none, so every product vanishes.")

theta_ch = a_full**2 / (2*np.pi)
print(f"\n  Series cross term alpha^2/pi = alpha^2 N^2; constitutive matching 2 mu theta_ch b^2 = (alpha^2/pi) mu b^2")
print(f"  theta_ch = alpha^2/(2 pi) = {theta_ch:.4e}   check 2 pi theta_ch/alpha^2 = {2*np.pi*theta_ch/a_full**2:.12f}")

m1_meV = a_full**3 * m_e_MeV / (4*np.pi**2) * 1e9
print(f"\n  Neutrino scale m1 = theta_ch^2 m0 = alpha^3 m_e/(4 pi^2) = {m1_meV:.2f} meV")
print("  (the Z3 eigenvalues and the mass sum built on this scale are not computed here)")
