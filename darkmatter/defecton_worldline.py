"""Defecton dispersion from the D4 worldline propagator, against the tight-binding band.

The vacancy is a point defect in the four-dimensional crystal, so its worldline is a
directed path on the 24 bonds of D4. Summing paths with fugacity z per step gives the
lattice propagator G(k) = 1/(1 - z S(k)), S(k) = sum_delta exp(i k.delta), and the
particle sits at the pole. Splitting k = (k_s, k_4) and continuing k_4 -> i q turns the
pole condition into a dispersion relation for the spatial wavevector k_s.

Checks, numerically and against the analytic expansion:
  1. The pole condition at k_s = 0 fixes the rest mass m from the fugacity z.
  2. Expanding to O(k_s^2) gives the inertial mass m*, and the ratio is
         m*/m = 3 sinh(x) / [ x (2 + cosh(x)) ],   x = q_0 l / sqrt(2) = m/(sqrt(2) m_0),
     which is 1 + O(x^4): m*/m = 1 - x^4/180 + ...
     The fugacity cancels from the ratio, so the equivalence principle holds to
     (m/m_0)^4/720 whatever the hop amplitude.
  3. The full dispersion tracks E = sqrt(m^2 c^4 + hbar^2 c^2 k^2), i.e. the worldline
     register is relativistic, not a non-relativistic band.
  4. Contrast: the equal-time tight-binding band on the same shell gives
     m* = m_0/(6 alpha) with rest energy m_0(1 - 24 alpha), hence m*/m_g ~ 28.
"""

import numpy as np
from itertools import permutations, product
from scipy.optimize import brentq

ALPHA = 1 / 137.035999177          # CODATA 2022 fine structure constant
M0_MEV = 0.51099895069 / ALPHA     # node mass m_0 = m_e/alpha [MeV]

# ── the 24 minimal vectors of D4, in units where the bond length is l = 1 ──
def d4_bonds():
    v = set()
    for p in permutations(range(4), 2):
        for s1, s2 in product((1, -1), repeat=2):
            w = [0, 0, 0, 0]
            w[p[0]], w[p[1]] = s1, s2
            v.add(tuple(w))
    return np.array(sorted(v), float) * np.sqrt(0.5)   # |delta| = 1

D = d4_bonds()
IN_SLICE = D[D[:, 3] == 0]          # 12 equal-time bonds
CROSSING = D[D[:, 3] != 0]          # 12 diagonals crossing the compact axis


def S12(ks):
    """Structure factor of the 12 in-slice bonds at spatial wavevector ks."""
    return float(np.sum(np.cos(IN_SLICE[:, :3] @ ks)))


def A(ks):
    """Spatial part of the crossing-bond structure factor (the k_4 = 0 envelope)."""
    return float(np.sum(np.cos(CROSSING[:, :3] @ ks))) / np.cosh(0.0) / 1.0


def cosh_q(ks, z):
    """Solve 1 = z[S12 + cosh(q l/sqrt2) * A] for cosh(q l/sqrt2).

    The crossing sum factorises as 2 cos(k4 a) * 2 sum_i cos(k_i a) with a = l/sqrt2,
    so continuing k4 -> i q replaces the cosine by a cosh and leaves A untouched.
    """
    a_env = 4.0 * float(np.sum(np.cos(np.sqrt(0.5) * ks)))   # = 4 sum_i cos(a k_i)
    return (1.0 / z - S12(ks)) / a_env


def q_of_k(ks, z):
    """Euclidean pole location q (in units 1/l) at spatial wavevector ks."""
    c = cosh_q(np.asarray(ks, float), z)
    if c < 1.0:
        return np.nan                       # below the mass shell: no real pole
    return np.sqrt(2.0) * np.arccosh(c)


def mass_from_z(z):
    """Rest mass in units of m_0: m/m_0 = q_0 l."""
    return q_of_k(np.zeros(3), z)


def z_from_mass(m_over_m0):
    """Invert: the fugacity that puts the rest mass at m_over_m0."""
    x = m_over_m0 / np.sqrt(2.0)
    return 1.0 / (12.0 * (1.0 + np.cosh(x)))


def m_star_numeric(z, h=1e-4):
    """Inertial mass from the curvature of E(k) at k = 0, in units of m_0.

    E = hbar c q, and m* = hbar^2 k^2 / (2 (E - E_0)) in the small-k limit.
    In units where l = 1 and hbar/(c l) = m_0: E/m_0 c^2 = q l.
    """
    q0 = q_of_k(np.zeros(3), z)
    qh = q_of_k(np.array([h, 0.0, 0.0]), z)
    return h ** 2 / (2.0 * (qh - q0))


def m_star_over_m_analytic(x):
    """3 sinh x / [x (2 + cosh x)] — the fugacity has cancelled."""
    return 3.0 * np.sinh(x) / (x * (2.0 + np.cosh(x)))


if __name__ == "__main__":
    print("Defecton dispersion: worldline register on D4\n" + "=" * 62)

    # ── 1. the band register, for contrast ──────────────────────────────────
    m_band_star = 1.0 / (6.0 * ALPHA)          # m*/m_0 from the 24-bond band curvature
    m_band_rest = 1.0 - 24.0 * ALPHA           # band-bottom rest mass / m_0
    print("\n[1] equal-time tight-binding band (all 24 bonds as spatial hops)")
    print(f"    m*   = m_0/(6 alpha) = {m_band_star:8.2f} m_0 = {m_band_star*M0_MEV/1e3:.3f} GeV")
    print(f"    m_g  = m_0(1-24 alpha) = {m_band_rest:6.4f} m_0 = {m_band_rest*M0_MEV:.1f} MeV")
    print(f"    m*/m_g = {m_band_star/m_band_rest:.1f}   (equivalence principle violated)")

    # ── 2. the worldline register ───────────────────────────────────────────
    print("\n[2] worldline register: pole of G(k) = 1/(1 - z S(k))")
    print(f"    critical fugacity z_c = 1/24 = {1/24:.6f} (massless)")
    print(f"    {'m/m_0':>8} {'z':>10} {'z/z_c':>8} {'m*/m (num)':>12} {'m*/m (ana)':>12} {'1-x^4/180':>12}")
    for m_over in (0.05, 0.2, 0.5, m_band_rest, 1.0, 2.0):
        z = z_from_mass(m_over)
        x = m_over / np.sqrt(2.0)
        num = m_star_numeric(z) / m_over
        print(f"    {m_over:8.4f} {z:10.6f} {z*24:8.4f} {num:12.8f} "
              f"{m_star_over_m_analytic(x):12.8f} {1 - x**4/180:12.8f}")

    # ── 3. is the dispersion relativistic? ──────────────────────────────────
    m_over = m_band_rest
    z = z_from_mass(m_over)
    print(f"\n[3] full dispersion at m = {m_over:.4f} m_0 against sqrt(m^2 + k^2)")
    print(f"    {'k l':>8} {'E/m_0c^2 (lattice)':>20} {'relativistic':>14} {'rel. error':>12}")
    for kl in (0.05, 0.1, 0.2, 0.4, 0.8):
        E = q_of_k(np.array([kl, 0.0, 0.0]), z)
        Erel = np.hypot(m_over, kl)
        print(f"    {kl:8.3f} {E:20.8f} {Erel:14.8f} {abs(E/Erel-1):12.2e}")

    # isotropy of the worldline dispersion
    kl = 0.3
    dirs = {"[100]": (1, 0, 0), "[110]": (1, 1, 0), "[111]": (1, 1, 1), "[210]": (2, 1, 0)}
    print(f"\n    isotropy at |k|l = {kl}:")
    for name, d in dirs.items():
        u = np.array(d, float); u /= np.linalg.norm(u)
        print(f"      {name:>6}: E/m_0c^2 = {q_of_k(kl*u, z):.10f}")

    # ── 4. the defecton number ──────────────────────────────────────────────
    x = m_band_rest / np.sqrt(2.0)
    dev = 1.0 - m_star_over_m_analytic(x)
    print(f"\n[4] defecton at m_g = {m_band_rest*M0_MEV:.1f} MeV:")
    print(f"    fugacity z = {z_from_mass(m_band_rest):.6f} = {z_from_mass(m_band_rest)*24:.4f} z_c")
    print(f"    m*/m_g = {m_star_over_m_analytic(x):.6f}")
    print(f"    equivalence-principle violation = {dev:.2e}  (band register: "
          f"{m_band_star/m_band_rest - 1:.0f})")
    print(f"    leading form (m/m_0)^4/720 = {(m_band_rest)**4/720:.2e}")
