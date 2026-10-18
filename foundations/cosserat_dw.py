#!/usr/bin/env python3
"""
cosserat_dw.py - zero-point displacement and Debye-Waller chain of the Cosserat vacuum
=====================================================================================
Companion script to "The Cosserat Supersolid" monograph (M. A. Cox).

What it computes, in the framework's own parametrisation:

1. The six-branch zero-point mean-square displacement <u^2> of a node.
   Moduli in units of the relaxed shear modulus mubar = rho c^2:
       N^2 = 1/pi,  mu_tot = mubar/(1 - N^2),  kappa_c = 2 N^2 mu_tot,
       mu (Eringen) = mu_tot - kappa_c,  gamma = mu ell^2,  J = gamma/c^2,
       c_L/c = sqrt((3 - 4N^2)/(1 - N^2))   (isotropic D4 constants).
   The transverse 2x2 problem (displacement u, microrotation phi) is diagonalised
   at each k, and each mode adds its displacement share times hbar/(2 m0 omega).
   Two closed forms are checked against the diagonalisation:
       <u^2>_LA = (k_D l)^2/(8 sqrt2 pi^2) (c/c_L) l^2
       <u^2>_TA = ln(1 + j (k_D l)^2/4) / (sqrt2 pi^2 j) l^2,   j = J/(rho l^2)
   (the acoustic transverse branch is exactly omega = c k at J = gamma/c^2, with
   displacement share 1/(1 + j k^2 l^2/4)).
2. The nearest-neighbour correlation C(l), the bond relative width sigma_rel, and
   the split longitudinal/transverse bond widths.
3. The Debye-Waller chain at the smallest in-plane D4* vector |G_min| = 2 pi sqrt2/l:
   2W_ac, the product with the core form factor, the 10 GeV velocity anomaly, the
   counterflow-bounded third factor exp(-2W_ac w f_s/(1-f_s)) with w = (k_c/k_D)^2,
   the GRB 090510 time delay, and the f_s floors; plus the Leggett ceiling.

Units: lengths in l (the lattice spacing), speeds in c, hbar/m0 = c l.
Run: python3 cosserat_dw.py
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

pi = np.pi
# ---------------------------------------------------------------- medium
N2 = 1/pi
mu_tot = 1/(1 - N2)
kap = 2*N2*mu_tot
mu = mu_tot - kap
gam = mu                      # gamma = mu l^2
J = gam                       # J = gamma / c^2   (j = J/(rho l^2))
cL = np.sqrt((3 - 4*N2)/(1 - N2))
Vp = 1/np.sqrt(2)             # primitive cell of the close-packed slice
kD = (6*pi**2/Vp)**(1/3)
pf = Vp/(2*pi**2)             # per-node mode weight V_p k^2 dk / (2 pi^2)


def transverse(k):
    """Return [(omega, displacement share)] for the acoustic and optical roots."""
    K = np.array([[(mu + kap)*k**2, -1j*kap*k], [1j*kap*k, gam*k**2 + 2*kap]])
    Mh = np.diag([1.0, 1/np.sqrt(J)])           # M^(-1/2), M = diag(rho, J)
    w2, V = np.linalg.eigh(Mh @ K @ Mh)
    out = []
    for i in range(2):
        v = Mh @ V[:, i]                         # M-normalised eigenvector
        out.append((np.sqrt(max(w2[i], 0.0)), abs(v[0])**2))
    return out


def u2_branches():
    LA = quad(lambda k: pf*k**2/(2*cL*k), 0, kD)[0]
    TA = quad(lambda k: 2*pf*k**2*transverse(k)[0][1]/(2*max(transverse(k)[0][0], 1e-12)), 1e-6, kD, limit=200)[0]
    TO = quad(lambda k: 2*pf*k**2*transverse(k)[1][1]/(2*transverse(k)[1][0]), 1e-6, kD, limit=200)[0]
    return LA, TA, TO


def weights(k):
    (wa, sa), (wo, so) = transverse(k)
    gL = 1/(2*cL*k)
    gT = sa/(2*max(wa, 1e-12)) + so/(2*wo)
    return gL, gT


def j0(z): return np.sin(z)/z
def j1(z): return np.sin(z)/z**2 - np.cos(z)/z


def correlation_and_widths(u2):
    den = quad(lambda k: k**2*(weights(k)[0] + 2*weights(k)[1]), 1e-6, kD, limit=200)[0]
    num = quad(lambda k: k**2*(weights(k)[0] + 2*weights(k)[1])*j0(k), 1e-6, kD, limit=200)[0]
    C = num/den
    srel = np.sqrt(2/3*u2*(1 - C))
    F0 = quad(lambda k: k**2*(weights(k)[0] + 2*weights(k)[1])/3, 1e-6, kD, limit=200)[0]
    Fpar = quad(lambda k: k**2*(weights(k)[0]*(j0(k) - 2*j1(k)/k) + weights(k)[1]*2*j1(k)/k), 1e-6, kD, limit=200)[0]
    Fperp = quad(lambda k: k**2*(weights(k)[0]*j1(k)/k + weights(k)[1]*(j0(k) - j1(k)/k)), 1e-6, kD, limit=200)[0]
    spar = np.sqrt(2*pf*(F0 - Fpar)); sperp = np.sqrt(2*pf*(F0 - Fperp))
    return C, srel, spar, sperp


# ---------------------------------------------------------------- DW chain
m_e, alpha = 0.51099895000, 7.2973525693e-3  # MeV and fine-structure constant, CODATA 2018
m0 = m_e/alpha                               # node mass, MeV
E_BZ = pi*np.sqrt(2)*m0
w_core = (pi/4)/np.sqrt(3)                   # PN core half-width, units of l
Gmin = 2*pi*np.sqrt(2)
F2 = np.exp(-2*Gmin*w_core)
nshell, fermi, dt_bound = 6, 7e-19, 0.86
z_grb, E1, E2 = 0.903, 100.0, 31000.0        # MeV, observed photon pair
H0, Om = 67.4, 0.315
H0_s = H0/3.0856775814913673e19


def dv(E, S): return nshell*S*(E*E_BZ)**2/(E**4 + E_BZ**4)


def log_S(fs, twoW, w, S12):
    """Natural log of the total suppression, kept in log space so f_s -> 1 cannot underflow."""
    return np.log(S12) - twoW*w*fs/(1 - fs)


def delay(S):
    f = lambda z: (dv(E1*(1+z), S) - dv(E2*(1+z), S))/((1+z)*H0_s*np.sqrt(Om*(1+z)**3 + 1 - Om))
    return abs(quad(f, 0, z_grb, limit=200)[0])


if __name__ == "__main__":
    print(f"medium: mu={mu:.5f} kappa_c={kap:.5f} mubar={mu + kap/2:.5f} j=J/(rho l^2)={J:.5f} "
          f"omega0={np.sqrt(2*kap/J):.4f} c/l c_L={cL:.4f} c k_D l={kD:.4f}")
    LA, TA, TO = u2_branches(); u2 = LA + TA + TO
    print(f"\n<u^2> six-branch: LA {LA:.4f} + TA {TA:.4f} + TO {TO:.4f} = {u2:.4f} l^2  (Lindemann {np.sqrt(u2):.3f})")
    print(f"  closed forms: LA {(kD**2/(8*np.sqrt(2)*pi**2))/cL:.4f}, TA {np.log(1 + J*kD**2/4)/(np.sqrt(2)*pi**2*J):.4f}")
    lin = kD**2/(8*np.sqrt(2)*pi**2)
    print(f"  three-branch Debye value (all at c): {3*lin:.4f}; linear branches with c_L: {lin*(1/cL + 2):.4f}")
    dev = max(abs(transverse(k)[0][0] - k) for k in np.linspace(0.01, kD, 400))
    print(f"  check: acoustic transverse branch max |omega - c k| over the Debye sphere = {dev:.1e} c/l")
    C, srel, spar, sperp = correlation_and_widths(u2)
    print(f"  C(l) = {C:.4f}, sigma_rel = {srel:.4f} l, sigma_par = {spar:.4f} l, sigma_perp = {sperp:.4f} l, ratio {sperp/spar:.4f}")

    u2a = 0.4045                                   # adopted value
    twoW = Gmin**2*u2a/3
    S12 = F2*np.exp(-twoW)
    w = 1/kD**2
    print(f"\nDW chain at <u^2> = {u2a}: 2W_ac = {twoW:.3f}, e^-2W = {np.exp(-twoW):.2e}, |F|^2 = {F2:.2e}, product {S12:.2e}")
    print(f"  dv/c at 10 GeV (first two factors) = {dv(1e4, S12):.2e} against {fermi:.0e}")
    print(f"  {'f_s':>7} {'S_total':>10} {'dv/c 10GeV':>11} {'GRB dt (s)':>11}")
    for fs in (5/6, 0.90, 0.95, 0.96, 0.97, 0.98, 0.99):
        S = S12*np.exp(-twoW*w*fs/(1 - fs))
        print(f"  {fs:7.4f} {S:10.2e} {dv(1e4, S):11.2e} {delay(S):11.2e}")
    for ww, lab in ((1/kD**2, "k_c = 1/l"), ((2/kD)**2, "k_c = 2/l")):
        fv = brentq(lambda fs: log_S(fs, twoW, ww, S12) + np.log(dv(1e4, 1.0)/fermi), 0.5, 0.999)
        ft = brentq(lambda fs: log_S(fs, twoW, ww, S12) + np.log(delay(1.0)/dt_bound), 0.5, 0.999)
        print(f"  floor ({lab}): velocity f_s > {fv:.4f}, time delay f_s > {ft:.4f}")
    G1sq = 6*pi**2                                 # |G|^2 l^2 of the first [111] harmonic
    delta = np.exp(-G1sq*u2a/6)
    print(f"\nLeggett ceiling: first [111] harmonic delta = {delta:.4e}, f_s <= 1 - 2 delta^2 = {1 - 2*delta**2:.5f}")
