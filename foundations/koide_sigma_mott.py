#!/usr/bin/env python3
"""
koide_sigma_mott.py — the Koide amplitude derived at the Mott point
=====================================================================

Closes the last open parameter of the charged-lepton sector: the SCP
dressing amplitude sigma, previously selected by the Koide target, is
derived from the quantum mechanics of the D4 compact ring at the Mott
critical point. Zero free parameters.

Chain (no Koide input anywhere):
  1. The spatial pair-potential cage cannot localise a node at the Mott
     point (the Einstein SCP map has no crystalline fixed point): the
     pair potential does not set sigma.
  2. The compact ring confines by topology regardless of well depth.
     The exact ring Schroedinger solve, in the deep-tunnelling regime
     the Mott condition forces, gives a Wannier width pinned by the
     ring geometry: sigma_W ~ d4/2 = l/(2 sqrt 2) = 0.354 l, insensitive
     to hbar (x2.6) and to the well depth (x4) at the few-percent level.
  3. D4 local isotropy (all 24 roots project with component 1/sqrt2 on
     every coordinate axis; Zener A = 1) extends the same amplitude to
     all four directions.
  4. The pair potential dressed at sigma = sigma_W gives (exact
     Hermite-convolution derivatives) gamma3_eff, gamma4_eff; the D4
     eta4 lattice sums then close the 0.46% Koide gap:
     eta4 = 0.61-0.62 against the required 0.624, i.e. the
     FCC-normalised gamma4 = 17.5 the monograph predicted.

Inputs: c, hbar, m_e through alpha and the Morse constraint a*l = 7/3
(gamma3 = -7 from the gravitational sector). Natural units l = 1,
V''(l) = 1, m0 = 1; the Mott condition hbar = m0 c l with
c^2 = S2 V'' l^2/(2 m0) gives hbar = sqrt(S2/2) (D4 sums: S2 = 1).

References:
  Cox, "The Cosserat Supersolid" (Koide chapter, SCP and D4 sections)
  Nosanow, Phys. Rev. 146, 120 (1966)   (self-consistent phonon method)
  Koide, Phys. Lett. B 120, 161 (1983)  (the charged-lepton relation)

Mitchell A. Cox, 2026
University of the Witwatersrand, Johannesburg
"""

import math
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from numpy.polynomial.laguerre import laggauss
from scipy.linalg import eigh

# ════════════════════════════════════════════════════════════════════
#  PHYSICAL CONSTANTS OF THE CHAIN (natural units l = 1, V'' = 1, m0 = 1)
# ════════════════════════════════════════════════════════════════════

A_MORSE = 7.0 / 3.0                 # Morse range: a*l = 7/3 (gamma3 = -7)
D_MORSE = 1.0 / (2 * A_MORSE**2)    # well depth 9/98 from V''(l) = 1
GAMMA3_BARE = -3 * A_MORSE          # = -7 (gravitational sector)
GAMMA4_BARE = 7 * A_MORSE**2        # = 343/9 = 38.1

D4_STEP = 1 / math.sqrt(2)          # compact component of a temporal root
L_RING = 3 * D4_STEP                # compact circumference (root metric)

ALPHA_INV = 137.035999              # for the bare PN overlap w/d
LOG_A = math.log(ALPHA_INV)         # ln(1/alpha) = 4.92024
W_OVER_D = LOG_A / (2 * math.pi)    # = 0.78308

# C_eff: strain-weighted overlap difference, exact closed form in alpha.
# The core's peak strain is eps_p = b/(2 pi w) = 1/ln(1/alpha) since the
# hop transfers one Shockley partial (b = d).  The self element samples
# the core's own strain intensity; the transfer element samples the
# symmetrised branch mean (eps_A^2 + eps_B^2)/2, forced by the localised
# valley configurations of the tight-binding ansatz.  All integrals are
# elementary Lorentzian moments (residues); the difference collapses to
#   C_eff = pi^2 (7 pi^2 + 11 L^2) / (16 L^2 (pi^2 + L^2)^2),  L = ln(1/alpha).
C_EFF = (math.pi**2 * (7 * math.pi**2 + 11 * LOG_A**2)
         / (16 * LOG_A**2 * (math.pi**2 + LOG_A**2)**2))   # = 0.0073585

def morse(r):
    """Morse pair potential, V''(1) = 1, minimum at r = 1."""
    r = np.asarray(r, dtype=float)
    return D_MORSE * (np.exp(np.clip(-2 * A_MORSE * (r - 1), -200, 200))
                      - 2 * np.exp(np.clip(-A_MORSE * (r - 1), -200, 200)))

# ════════════════════════════════════════════════════════════════════
#  D4 GEOMETRY AND THE eta4 LATTICE SUMS
# ════════════════════════════════════════════════════════════════════

def d4_roots():
    """The 24 D4 roots (+-e_i +- e_j)/sqrt2, i < j <= 4 (bond length 1)."""
    roots = []
    for i in range(4):
        for j in range(i + 1, 4):
            for si in (1, -1):
                for sj in (1, -1):
                    v = np.zeros(4)
                    v[i], v[j] = si, sj
                    roots.append(v / math.sqrt(2))
    return np.array(roots)

def eta4_sums(roots):
    """Bond-stretch lattice sums for anti-plane shear {111}<110>."""
    b_hat = np.array([-1, 1, 0, 0]) / math.sqrt(2)
    xi_hat = np.array([1, 1, -2, 0]) / math.sqrt(6)
    a1 = (roots @ b_hat) * (roots @ xi_hat)
    a2 = 0.5 * (roots @ b_hat)**2 * (1 - (roots @ xi_hat)**2)
    S2 = float(np.sum(a1**2))
    return {"S2": S2,
            "W4": float(np.sum(a2**2) / S2),
            "C3": float(np.sum(a1**2 * a2) / S2),
            "C4": float(np.sum(a1**4) / (12 * S2))}

# ════════════════════════════════════════════════════════════════════
#  1. THE SPATIAL PAIR CAGE HAS NO CRYSTALLINE FIXED POINT
# ════════════════════════════════════════════════════════════════════

def V_dressed(R, sigma, ndim=3, nq=64):
    """<V(|R e_x + u|)>, u Gaussian, per-component width sigma, ndim dims."""
    if sigma < 1e-10:
        return morse(np.atleast_1d(R))
    xi, wh = hermegauss(nq)
    wz = wh / math.sqrt(2 * math.pi)
    k = ndim - 1
    tl, wl = laggauss(nq)
    qw = wl * tl**(k / 2 - 1) / math.gamma(k / 2)
    q = 2 * sigma**2 * tl
    R = np.atleast_1d(np.asarray(R, dtype=float))
    Rz = R[:, None, None] + sigma * xi[None, :, None]
    r = np.maximum(np.sqrt(Rz**2 + q[None, None, :]), 0.02)
    return np.sum(morse(r) * wz[None, :, None] * qw[None, None, :], axis=(1, 2))

def spatial_cage_map(sigma, hbar, roots3, ndim=3):
    """One step of the Einstein SCP map sigma -> sigma' for the 12-bond
    spatial cage at physical spacing 1: Phi'' = sum V''_eff c^2 + V'_eff (1-c^2)."""
    s_rel = math.sqrt(2) * sigma          # uncorrelated relative smearing
    h = 5e-3
    xs = np.array([1 - 2*h, 1 - h, 1.0, 1 + h, 1 + 2*h])
    ys = V_dressed(xs, s_rel, ndim)
    V1 = (ys[3] - ys[1]) / (2 * h)
    V2 = (ys[3] - 2 * ys[2] + ys[1]) / h**2
    e = np.zeros(roots3.shape[1]); e[0] = 1
    c2 = (roots3 @ e)**2
    Phi = float(np.sum(V2 * c2 + V1 * (1 - c2)))
    if Phi <= 0:
        return None
    return math.sqrt(hbar / (2 * math.sqrt(Phi)))

def demonstrate_no_spatial_fixed_point():
    roots = d4_roots()
    spatial = roots[np.abs(roots[:, 3]) < 1e-12][:, :3]
    hbar = math.sqrt(5.0 / 6.0 / 2.0)     # FCC sums, Mott condition
    print("  Einstein SCP map for the 12-bond spatial pair cage "
          "(Mott hbar = %.4f):" % hbar)
    s = 0.05
    for it in range(12):
        s_new = spatial_cage_map(s, hbar, spatial)
        if s_new is None:
            print(f"    s = {s:.3f} -> cage curvature <= 0 (no restoring force)")
            break
        print(f"    s = {s:.3f} -> {s_new:.3f}")
        if s_new > 0.9:
            print("    -> runaway delocalisation: no crystalline fixed point")
            break
        s = s_new

# ════════════════════════════════════════════════════════════════════
#  2. THE COMPACT RING: EXACT SOLVE AND GEOMETRIC SATURATION
# ════════════════════════════════════════════════════════════════════

def cage_potential(x4, depth_scale=1.0):
    """Periodic ring cage from the 12 temporal Morse bonds (6 per layer),
    each with spatial offset 1/sqrt2 (root metric)."""
    x4 = np.asarray(x4, dtype=float)
    V = np.zeros_like(x4)
    for x_layer in (D4_STEP, 2 * D4_STEP):
        dx = x4 - x_layer
        dx -= L_RING * np.round(dx / L_RING)
        V += 6 * depth_scale * morse(np.sqrt(0.5 + dx**2))
    return V

def ring_solve(hbar, depth_scale=1.0, N=1200):
    """Exact Schroedinger solve on the periodic ring; returns the
    maximally localised Wannier width and tight-binding parameters."""
    x = np.linspace(0, L_RING, N, endpoint=False)
    dx = x[1] - x[0]
    V = cage_potential(x, depth_scale)
    V0 = V.min(); i0 = int(np.argmin(V))
    K = hbar**2 / (2 * dx**2)
    idx = np.arange(N)
    H = np.zeros((N, N))
    H[idx, idx] = 2 * K + (V - V0)
    H[idx, (idx + 1) % N] = -K
    H[idx, (idx - 1) % N] = -K
    ev, evec = eigh(H, subset_by_index=[0, 3])
    psi = evec[:, :3].copy()
    for j in range(3):
        if psi[i0, j] < 0:
            psi[:, j] *= -1
    w = psi.sum(axis=1) / math.sqrt(3)
    w /= math.sqrt(float((w**2).sum()) * dx)
    xs = x - x[i0]
    xs -= L_RING * np.round(xs / L_RING)
    p = w**2
    m1 = float((xs * p).sum()) * dx
    m2 = float((xs * xs * p).sum()) * dx
    sigma_W = math.sqrt(m2 - m1 * m1)
    # tight-binding parameters and adjacent-site Wannier overlap
    E0 = ev[0]; E1 = 0.5 * (ev[1] + ev[2])
    t_tb = (E1 - E0) / 3
    eps_tb = (E0 + 2 * E1) / 3
    w_shift = np.interp((x - D4_STEP) % L_RING, x, w)
    C12 = float((w * w_shift).sum()) * dx
    barrier = float(cage_potential(np.array([D4_STEP / 2]))[0]) - V0
    zpe = 0.5 * hbar * math.sqrt(6.0)     # V''_cage = 6 V''
    frac = float(p[np.abs(xs) < D4_STEP / 2].sum()) * dx
    return dict(sigma_W=sigma_W, t=t_tb, eps=eps_tb, C12=C12,
                barrier=barrier, zpe=zpe, frac_local=frac)

def free_ring_wannier():
    """Free-particle floor: Wannier width of the 3 lowest plane waves."""
    N = 6000
    x = np.linspace(0, L_RING, N, endpoint=False); dx = x[1] - x[0]
    w = (1 + 2 * np.cos(2 * np.pi * x / L_RING))
    w /= math.sqrt(float((w**2).sum()) * dx)
    xs = x - L_RING * np.round(x / L_RING)
    return math.sqrt(float(((xs**2) * (w**2)).sum()) * dx)

# ════════════════════════════════════════════════════════════════════
#  3. EXACT DERIVATIVES OF THE DRESSED PAIR POTENTIAL
# ════════════════════════════════════════════════════════════════════

def _hermite_prob(n, x):
    if n == 0:
        return np.ones_like(x)
    Hm2, Hm1 = np.ones_like(x), x
    for k in range(2, n + 1):
        Hm2, Hm1 = Hm1, x * Hm1 - (k - 1) * Hm2
    return Hm1

def dressed_props(sigma, ndim=3, nq=80, Vfun=morse):
    """Dressed minimum d_eff and gamma3_eff, gamma4_eff by exact
    differentiation under the Gaussian: V^(n)(R) = sigma^-n
    <He_n(u/sigma) F(R+u)>, spectrally accurate (no polynomial fitting)."""
    k = ndim - 1
    tl, wl = laggauss(64)
    qw = wl * tl**(k / 2 - 1) / math.gamma(k / 2)
    q = 2 * sigma**2 * tl
    xi, wh = hermegauss(nq)
    wz = wh / math.sqrt(2 * math.pi)

    def F_perp(z):
        z = np.atleast_1d(z)
        r = np.maximum(np.sqrt(z[:, None]**2 + q[None, :]), 1e-3)
        return Vfun(r) @ qw

    def dn(R, n):
        Hn = _hermite_prob(n, xi)
        return float(F_perp(R + sigma * xi) @ (wz * Hn)) / sigma**n

    lo, hi = 1.0, 2.4
    flo = dn(lo, 1)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if flo * dn(mid, 1) <= 0:
            hi = mid
        else:
            lo = mid
            flo = dn(lo, 1)
    d = 0.5 * (lo + hi)
    V2, V3, V4 = dn(d, 2), dn(d, 3), dn(d, 4)
    return d, V2, d * V3 / V2, d * d * V4 / V2

# ════════════════════════════════════════════════════════════════════
#  4. THE KOIDE CHAIN
# ════════════════════════════════════════════════════════════════════

def koide_chain(sigma, sums_d4, ndim=3):
    """sigma -> dressed gammas -> eta4 (D4 sums) -> t/eps -> Q."""
    d, V2, g3, g4 = dressed_props(sigma, ndim)
    eta4 = sums_d4["W4"] + g3 * sums_d4["C3"] + g4 * sums_d4["C4"]
    te_bare = 4 * W_OVER_D**2 / (1 + 4 * W_OVER_D**2)
    te = te_bare * (1 - eta4 * C_EFF)
    Q = 1.0 / 3.0 + 2.0 * te**2 / 3.0
    # FCC-normalised equivalent quartic (the monograph's gamma4 = 17.5)
    g4_fcc = (eta4 - (0.668 - 7 * 0.03967 - 0.019)) / 0.0142
    return dict(d_eff=d, g3=g3, g4=g4, eta4=eta4, te=te, Q=Q, g4_fcc=g4_fcc)

# ════════════════════════════════════════════════════════════════════
#  MAIN
# ════════════════════════════════════════════════════════════════════

def ceff_verify(n=200000):
    """Quadrature check of the closed-form C_eff and the rejected
    alternative constructions (coherent-average and saddle-field)."""
    delta = 2 * math.pi / LOG_A
    eps_p2 = 1 / LOG_A**2
    x = np.linspace(-40, 40, n)
    Lo = 1 / (1 + x**2)
    Ls = 1 / (1 + (x - delta)**2)
    Lm = 1 / (1 + (x - delta / 2)**2)
    den = np.trapezoid(Lo * Ls, x)
    mean_branch = np.trapezoid(Lo * Ls * (Lo**2 + Ls**2) / 2, x) / den
    coh_avg = np.trapezoid(Lo * Ls * ((Lo + Ls) / 2)**2, x) / den
    saddle = np.trapezoid(Lo * Ls * Lm**2, x) / den
    return {"C_quad": (0.625 - mean_branch) * eps_p2,
            "C_coherent": (0.625 - coh_avg) * eps_p2,
            "C_saddle": (0.625 - saddle) * eps_p2}

if __name__ == "__main__":
    print("=" * 72)
    print("  THE KOIDE AMPLITUDE FROM THE MOTT POINT (zero free parameters)")
    print("=" * 72)

    v = ceff_verify()
    print(f"\n  C_eff closed form = {C_EFF:.8f}  (quadrature: {v['C_quad']:.8f})")
    print(f"  rejected constructions: coherent-average {v['C_coherent']:.6f}, "
          f"saddle-field {v['C_saddle']:.6f}")

    roots = d4_roots()
    sums = eta4_sums(roots)
    print(f"\n  D4 lattice sums: S2 = {sums['S2']:.4f}, W4 = {sums['W4']:.4f}, "
          f"C3 = {sums['C3']:.4f}, C4 = {sums['C4']:.4f}")
    hbar_d4 = math.sqrt(sums["S2"] / 2)
    print(f"  Mott condition hbar = m0 c l -> hbar = sqrt(S2/2) = {hbar_d4:.4f}")

    print(f"\n{'-' * 72}\n  STEP 1: THE SPATIAL PAIR CAGE DOES NOT SET sigma\n{'-' * 72}")
    demonstrate_no_spatial_fixed_point()

    print(f"\n{'-' * 72}\n  STEP 2: THE COMPACT RING (topology confines)\n{'-' * 72}")
    print(f"  Ring: step d4 = {D4_STEP:.4f} l, circumference L = {L_RING:.4f} l")
    print(f"  Geometric anchor d4/2 = {D4_STEP/2:.4f} l; "
          f"free-ring floor = {free_ring_wannier():.4f} l")
    central = ring_solve(hbar_d4)
    print(f"\n  Exact solve at Mott hbar = {hbar_d4:.4f}:")
    print(f"    sigma_W = {central['sigma_W']:.4f} l   "
          f"(= {central['sigma_W']/(D4_STEP/2):.3f} x d4/2)")
    print(f"    ZPE/barrier = {central['zpe']/central['barrier']:.2f}  "
          f"(deep tunnelling); localisation = {central['frac_local']:.0%}")
    print(f"    band t/eps = {central['t']/central['eps']:.3f}; "
          f"Wannier overlap C12 = {central['C12']:+.3f}")

    print("\n  Insensitivity scan (hbar x2.6, well depth x4):")
    band = []
    for hb in (math.sqrt(5.0/12), hbar_d4, 1.0):
        for ds in (0.5, 1.0, 2.0):
            r = ring_solve(hb, ds)
            band.append(r["sigma_W"])
            print(f"    hbar = {hb:.3f}, depth x{ds:3.1f}: "
                  f"sigma_W = {r['sigma_W']:.4f} l")
    sig_lo = min(ring_solve(h)["sigma_W"] for h in (math.sqrt(5/12), hbar_d4, 1.0))
    sig_hi = max(ring_solve(h)["sigma_W"] for h in (math.sqrt(5/12), hbar_d4, 1.0))

    print(f"\n{'-' * 72}\n  STEP 3: D4 ISOTROPY -> sigma equal in all four directions\n{'-' * 72}")
    e4 = np.abs(roots[:, 3])
    print(f"  Every root projects on every axis with component 0 or "
          f"{np.max(e4):.4f} = 1/sqrt2: the projected site spacing is d4 in "
          f"all four directions;\n  Zener A = 1 makes the cage isotropic. "
          f"The dressing amplitude is sigma_W per component.")

    print(f"\n{'-' * 72}\n  STEP 4: THE KOIDE CHAIN (exact derivatives, 3D dressing)\n{'-' * 72}")
    te_bare = 4 * W_OVER_D**2 / (1 + 4 * W_OVER_D**2)
    eta_req = (te_bare - 1 / math.sqrt(2)) / te_bare / C_EFF
    print(f"  bare t/eps = {te_bare:.5f}; gap to 1/sqrt2 = "
          f"{(te_bare - 1/math.sqrt(2))/te_bare*100:.3f}%; "
          f"required eta4 = {eta_req:.3f}")
    print(f"\n  {'sigma':>8s}{'gamma3':>9s}{'gamma4':>9s}{'eta4':>8s}"
          f"{'g4(FCC)':>9s}{'t/eps':>10s}{'Q':>10s}")
    for s, tag in ((sig_lo, "band low "), (D4_STEP/2, "d4/2     "),
                   (central["sigma_W"], "Mott hbar"), (sig_hi, "band high")):
        r = koide_chain(s, sums)
        print(f"  {s:8.4f}{r['g3']:9.2f}{r['g4']:9.1f}{r['eta4']:8.3f}"
              f"{r['g4_fcc']:9.1f}{r['te']:10.5f}{r['Q']:10.5f}   {tag}")
    print(f"\n  Koide requirement: t/eps = {1/math.sqrt(2):.5f}, "
          f"Q_obs = 0.666661(7)")
    r = koide_chain(central["sigma_W"], sums)
    print(f"  At the Mott-hbar solve: the quartic correction covers "
          f"{r['eta4']/eta_req*100:.0f}% of the bare gap;")
    print(f"  the FCC-normalised quartic is gamma4 = {r['g4_fcc']:.1f} "
          f"(FCC-bookkeeping requirement: 18.0).")

    print(f"\n{'-' * 72}\n  APPENDIX: exact-derivative SCP table (replaces fitted table)\n{'-' * 72}")
    print(f"  {'s':>5s}{'d_eff':>8s}{'V2eff':>8s}{'g3eff':>8s}{'g4eff':>8s}{'g4/g3^2':>9s}")
    print(f"  {0.0:5.2f}{1.0:8.3f}{1.0:8.3f}{GAMMA3_BARE:8.2f}{GAMMA4_BARE:8.1f}"
          f"{GAMMA4_BARE/GAMMA3_BARE**2:9.3f}")
    for s in (0.10, 0.20, 0.25, 0.30, 0.35):
        d, V2, g3, g4 = dressed_props(s)
        print(f"  {s:5.2f}{d:8.3f}{V2:8.3f}{g3:8.2f}{g4:8.1f}{g4/g3**2:9.3f}")
