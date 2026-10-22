"""Core half-widths of the screw and edge dislocation from the Cosserat
Peierls-Nabarro balance.

Both characters of one glide system solve the same variational balance,

    integral_0^inf (K(k)/k) exp(-2 k w) dk  =  mu_bar / (2 d_111),

with the Frenkel misfit, whose pinning on the right belongs to the glide plane
and is common to both.  Only the elastic line kernel K(k) differs.

Screw (anti-plane).  The exact projection of the Cosserat bulk stiffness onto
the glide line,

    K_s(k) = (mu_bar/2)|k| + (kappa_c/4) k^2 / sqrt(k^2 + q^2),
    q^2 = 2 kappa_c / gamma.

Edge (plane strain).  The same energy, with the slip written as a plastic
distortion and the displacement split along and across each bulk wavevector K.
The longitudinal channel carries the dilatation and gives the classical Poisson
factor; the transverse channel carries the slip's rotation, where the shear
stiffness 2 mu_bar and the relaxed rotational stiffness
g(K) = kappa_c K^2/(K^2 + q^2) act in series.  Per bulk mode,

    E(K)/B^2 = mu_bar kx^2 ky^2 / ((1-nu) K^4)
             + [2 mu_bar g/(2 mu_bar + g)] kx^4 / K^4,

and integrating over ky gives

    K_e(k) = mu_bar |k| / (2(1-nu)) + (gamma/2) |k|^3 (1 - |k|/sqrt(k^2 + p^2)),
    p^2 = (mu_bar/mu_tot) q^2.

At long wavelength K_e is 1/(1-nu) times the screw's long-range kernel, the
classical Poisson factor; at the core's wavelengths it is NOT 1/(1-nu) times
K_s, because the rotational piece carries no Poisson factor and switches on as
|k|^3.  The Poisson ratio is the crystalline skeleton's, nu = lambda/(2(lambda
+ mu_bar)) = (1-2N^2)/(4-6N^2) = 0.17385, with lambda = C_12.

Checks, in order:
  (1) both closed-form kernels against a direct minimisation of the Cosserat
      energy mode by mode, followed by the ky quadrature;
  (2) the Cauchy limits (kappa_c = 0, nu = 1/4) analytically:
      screw w = d_111/2 = d/sqrt(2), edge w = d_111/(2(1-nu));
  (3) the screw width at the vacuum point, w/d = 0.83664;
then the edge width, the edge-to-screw ratio, the result of scaling the whole
screw kernel by 1/(1-nu) for comparison, and the edge width carried to the
adopted screw width w = d ln(1/alpha)/(2 pi) by two different ways of closing
the screw's residual.

Units: lengths in ell (nearest-neighbour spacing), moduli in the Eringen shear
modulus mu = 1.  d = ell/sqrt(3) is the partial (misfit) period.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import struve, y1

# Vacuum-point Cosserat constants (Eringen units, mu = 1)
N2 = 1.0 / np.pi                       # Cosserat coupling number
KAPPA_C = 2.0 * N2 / (1.0 - 2.0 * N2)  # = 2/(pi-2) = 1.752
MU = 1.0
MU_TOT = MU + KAPPA_C                  # lattice-scale shear modulus, 2.752
MU_BAR = MU + KAPPA_C / 2.0            # true (long-range) shear modulus, 1.876
GAMMA = MU * 1.0                       # curvature modulus, gamma = mu ell^2
LAMBDA = (1.0 - 2.0 * N2) * MU_TOT     # C_12 = mu at the vacuum point
NU = LAMBDA / (2.0 * (LAMBDA + MU_BAR))  # 0.17385
Q2 = 2.0 * KAPPA_C / GAMMA             # 3.504 / ell^2
P2 = (MU_BAR / MU_TOT) * Q2            # 2.389 / ell^2

D111 = np.sqrt(2.0 / 3.0)              # {111} interplanar spacing
D_PARTIAL = 1.0 / np.sqrt(3.0)         # Shockley partial period
RHS = MU_BAR / (2.0 * D111)            # Frenkel misfit pinning, 1.149 mu/ell

ALPHA = 1.0 / 137.035999177            # CODATA 2022
W_ADOPTED = np.log(1.0 / ALPHA) / (2.0 * np.pi)   # adopted screw w/d = 0.78308


# --------------------------------------------------------------------------
# Closed-form kernels
# --------------------------------------------------------------------------
def kernel_screw(k):
    return 0.5 * MU_BAR * k + 0.25 * KAPPA_C * k * k / np.sqrt(k * k + Q2)


def kernel_edge(k):
    return (MU_BAR * k / (2.0 * (1.0 - NU))
            + 0.5 * GAMMA * k**3 * (1.0 - k / np.sqrt(k * k + P2)))


# --------------------------------------------------------------------------
# Check (1): direct minimisation of the Cosserat energy, mode by mode
# --------------------------------------------------------------------------
def _mode_energy(kx, ky, disp, rot, B, character, lam=LAMBDA, mubar=MU_BAR,
                 kc=KAPPA_C, gam=GAMMA):
    """Energy of one bulk mode for the centrosymmetric energy
    W = (lam/2) tr(e)^2 + mubar e_(ij) e_(ij) + (kc/2) e_[ij] e_[ij]
        + (gam/2) |grad phi|^2,
    e_ij = d_j u_i - eps_ijk phi_k - beta*_ij, with the Fourier convention
    d_j u_i -> k_j v_i (real amplitudes).  The slip B sits in beta*_xy (edge)
    or beta*_zy (screw)."""
    kv = np.array([kx, ky, 0.0])
    eps = np.zeros((3, 3, 3))
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1.0
    eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1.0
    e = np.outer(disp, kv) - np.einsum("ijk,k->ij", eps, rot)
    if character == "edge":
        e[0, 1] -= B
    else:
        e[2, 1] -= B
    sym = 0.5 * (e + e.T)
    asym = 0.5 * (e - e.T)
    return (0.5 * lam * np.trace(e)**2 + mubar * np.sum(sym * sym)
            + 0.5 * kc * np.sum(asym * asym)
            + 0.5 * gam * (kx * kx + ky * ky) * np.sum(rot * rot))


def mode_energy_min(kx, ky, character, **moduli):
    """Minimum over displacement and microrotation at unit slip, as 2E/B^2."""
    if character == "edge":
        def pack(x, B):
            return np.array([x[0], x[1], 0.0]), np.array([0.0, 0.0, x[2]]), B
    else:
        def pack(x, B):
            return np.array([0.0, 0.0, x[0]]), np.array([x[1], x[2], 0.0]), B
    # quadratic form in (x1, x2, x3, B) by polarisation
    n = 4
    Qm = np.zeros((n, n))
    def f(z):
        d_, r_, B = pack(z[:3], z[3])
        return _mode_energy(kx, ky, d_, r_, B, character, **moduli)
    I = np.eye(n)
    for i in range(n):
        for j in range(n):
            Qm[i, j] = 0.5 * (f(I[i] + I[j]) - f(I[i]) - f(I[j]))
    Axx, AxB, ABB = Qm[:3, :3], Qm[:3, 3], Qm[3, 3]
    Emin = ABB - AxB @ np.linalg.solve(Axx, AxB)
    return 2.0 * Emin


def kernel_direct(k, character, **moduli):
    """K(k) = int dky/(2 pi) 2E/B^2, with ky = k tan(t) mapping [0, inf) to [0, pi/2)."""
    def integrand(t):
        ky = k * np.tan(t)
        return mode_energy_min(k, ky, character, **moduli) * k / np.cos(t)**2
    val, _ = quad(integrand, 0.0, 0.5 * np.pi, limit=200, epsabs=1e-11, epsrel=1e-9)
    return 2.0 * val / (2.0 * np.pi)       # even in ky


# --------------------------------------------------------------------------
# The balance
# --------------------------------------------------------------------------
def lhs_screw(w):
    """Closed form: mu_bar/(4w) + (kappa_c/4) q [(pi/2)(H_1(2qw) - Y_1(2qw)) - 1]."""
    q = np.sqrt(Q2)
    a = 2.0 * q * w
    return MU_BAR / (4.0 * w) + 0.25 * KAPPA_C * q * ((np.pi / 2.0) * (struve(1, a) - y1(a)) - 1.0)


def lhs(kernel, w, scale=1.0):
    val, _ = quad(lambda k: scale * kernel(k) / k * np.exp(-2.0 * k * w),
                  0.0, np.inf, limit=400, epsabs=1e-14, epsrel=1e-12)
    return val


def solve(f, rhs=RHS):
    return brentq(lambda w: f(w) - rhs, 0.05, 5.0, xtol=1e-14)


if __name__ == "__main__":
    print(f"kappa_c = {KAPPA_C:.5f}, mu_bar = {MU_BAR:.5f}, mu_tot = {MU_TOT:.5f}, "
          f"lambda = {LAMBDA:.5f}")
    print(f"nu = {NU:.5f}, 1/(1-nu) = {1/(1-NU):.5f}, q = {np.sqrt(Q2):.5f}, "
          f"p = {np.sqrt(P2):.5f}  [1/ell]")

    print("\n(1) closed-form kernels against direct mode-by-mode minimisation")
    for k in (0.05, 0.5, 1.0, 3.0, 10.0):
        ks, ke = kernel_direct(k, "screw"), kernel_direct(k, "edge")
        print(f"  k = {k:5.2f}: screw {kernel_screw(k)/ks - 1:+.1e}   "
              f"edge {kernel_edge(k)/ke - 1:+.1e}   "
              f"edge/[screw/(1-nu)] = {kernel_edge(k)/(kernel_screw(k)/(1-NU)):.4f}")

    print("\n(2) Cauchy limit (kappa_c = 0, lambda = mu_bar = 1, so nu = 1/4)")
    cauchy = dict(lam=1.0, mubar=1.0, kc=0.0, gam=GAMMA)
    for k in (0.3, 2.0):
        print(f"  k = {k}: direct screw / (|k|/2) = {kernel_direct(k, 'screw', **cauchy)/(k/2):.8f}, "
              f"direct edge / (|k|/(2(1-1/4))) = {kernel_direct(k, 'edge', **cauchy)/(k/1.5):.8f}")
    # with K = |k|/2 or |k|/(2(1-nu)), the balance reads 1/(4w) [/(1-nu)] = 1/(2 d_111)
    print(f"  so screw w = d_111/2 = {D111/2:.5f} ell = {D111/2/D_PARTIAL:.5f} d "
          f"(1/sqrt2 = {1/np.sqrt(2):.5f}), "
          f"edge w = d_111/(2(1-nu)) = {D111/1.5:.5f} ell, the classical Peierls widths")

    w_s = solve(lhs_screw)
    print(f"\n(3) Cosserat screw:  w = {w_s:.5f} ell = {w_s / D_PARTIAL:.5f} d "
          f"(alpha chapter: 0.83664 d)")
    print(f"    quadrature check: {solve(lambda w: lhs(kernel_screw, w)):.5f} ell")

    w_e = solve(lambda w: lhs(kernel_edge, w))
    print(f"\nCosserat edge:   w = {w_e:.5f} ell = {w_e / D_PARTIAL:.5f} d")
    print(f"edge/screw ratio = {w_e / w_s:.5f} (long-wave factor 1/(1-nu) = {1/(1-NU):.5f})")
    w_u = solve(lambda w: lhs(kernel_screw, w, 1.0 / (1.0 - NU)))
    print(f"for comparison, screw kernel x 1/(1-nu) at every k: w = {w_u:.5f} ell, "
          f"ratio {w_u / w_s:.5f}")

    # Carry the ratio to the adopted screw width.
    ws_adopt = W_ADOPTED * D_PARTIAL
    print(f"\nadopted screw width w = {W_ADOPTED:.5f} d = {ws_adopt:.5f} ell")
    print(f"ratio x adopted: w_edge = {w_e / w_s * ws_adopt:.5f} ell "
          f"= {w_e / w_s * W_ADOPTED:.5f} d")
    # closure A: a stiffer misfit S (equivalently a uniformly softer kernel)
    S = brentq(lambda S: solve(lhs_screw, S * RHS) - ws_adopt, 0.9, 1.3, xtol=1e-12)
    wA = solve(lambda w: lhs(kernel_edge, w), S * RHS)
    print(f"closure A, misfit stiffened by S = {S:.5f}: edge w = {wA:.5f} ell")
    # closure B: soften only the rotational stiffening of both kernels by s
    def ks_s(k, s): return 0.5 * MU_BAR * k + s * 0.25 * KAPPA_C * k * k / np.sqrt(k * k + Q2)
    def ke_s(k, s): return (MU_BAR * k / (2.0 * (1.0 - NU))
                            + s * 0.5 * GAMMA * k**3 * (1.0 - k / np.sqrt(k * k + P2)))
    s = brentq(lambda s: solve(lambda w: lhs(lambda k: ks_s(k, s), w)) - ws_adopt,
               0.0, 1.0, xtol=1e-10)
    wB = solve(lambda w: lhs(lambda k: ke_s(k, s), w))
    print(f"closure B, rotational stiffening scaled by s = {s:.4f}: edge w = {wB:.5f} ell")

    ELL_FM = 2.8179403205                      # lattice spacing = r_e [fm]
    print(f"\nedge transverse half-width L_b = {w_e / w_s * ws_adopt * ELL_FM:.3f} fm")
