"""The hopping amplitude samples the lattice at finite wavevector, and that
changes the rotational share of h/A.

The ring Hamiltonian carries an on-site energy A and a hop h, and the spectrum's
shape depends only on their ratio.  The leading-order form h/A = (1 + N^2)/2
uses the bare Cosserat coupling number N^2 = 1/pi for the rotational channel.
N^2 is a rolling-CONTACT quantity, defined at a point contact, so it is the
short-wavelength limit of the rotational share and not its value at every scale.

In a Cosserat medium the microrotation is screened over 1/q, so the medium's
response is wavevector dependent.  Eliminating the microrotation at fixed
displacement (see derive_mu_eff below) gives

    mu_eff(k) = mu + kappa k^2 / (k^2 + q^2),

running from mu at long wavelength, where the rotation is screened out, to
mu + kappa at short wavelength, where it participates fully.  The rotational
share of the response at wavevector k is therefore

    N^2_eff(k) = kappa x / (mu + kappa x),      x = k^2/(k^2 + q^2),

which recovers N^2 only as k -> infinity.

A and h sit at opposite ends of this.  A is the far-field self-energy, where the
microrotation is screened and contributes nothing, so A takes the k -> 0 limit
and is uncorrected.  The hop proceeds by a Shockley partial step through the
stacking channel, of length d = ell/sqrt(3), so h samples the lattice at
k = 2 pi / d.  That wavevector is not fitted: it is the step the hop takes.
"""

import numpy as np
import sympy as sp

# --------------------------------------------------------------------------
# Inputs, all previously fixed
# --------------------------------------------------------------------------
ALPHA = 1 / 137.035999177
M_E = 0.51099895069          # MeV
N2 = 1 / np.pi               # Cosserat coupling number, rolling-contact value
Q = 1.87                     # Cosserat screening wavevector [1/ell]
D_PARTIAL = 1 / np.sqrt(3)   # Shockley partial period [ell]

# NuFIT 6.1, normal ordering, IC24 with SK atmospheric data
DM21, S21 = 7.537, 0.10      # 1e-5 eV^2
DM31, S31 = 2.511, 0.021     # 1e-3 eV^2


def derive_mu_eff():
    """Symbolic check of mu_eff(k), by relaxing the microrotation."""
    k, q, mu, kap, u = sp.symbols('k q mu kappa u', positive=True)
    omega = -sp.I * k * u / 2                      # macrorotation
    phi = omega / (1 + k**2 / q**2)                # relaxed microrotation
    gam = 4 * kap / q**2                           # curvature modulus from q
    W = (sp.Rational(1, 2) * mu * k**2 * u**2
         + 2 * kap * sp.Abs(phi - omega)**2
         + sp.Rational(1, 2) * gam * k**2 * sp.Abs(phi)**2)
    return sp.simplify(sp.expand(2 * W / (k**2 * u**2)))


def rotational_share(k):
    """N^2_eff(k): the fraction of the response carried by the rotation."""
    kap_over_mu = N2 / (1 - N2)                    # = 1/(pi - 1)
    x = k**2 / (k**2 + Q**2)
    return kap_over_mu * x / (1 + kap_over_mu * x)


def h_over_A(k=None, dress=True):
    """h/A = (1/2)(1 + rotational share), rotational leg dressed.

    The translational leg stays bare: it is already inside the self-consistent
    alpha, which is solved against the full periodic potential V_PN(u), so
    dressing it would double-book.  The rotational leg is not, because V_PN
    knows nothing of the microrotation, so it competes and carries (1 - alpha).
    """
    n2 = N2 if k is None else rotational_share(k)
    return 0.5 * (1 + n2 * ((1 - ALPHA) if dress else 1.0))


def spectrum(hA, scale=1.0):
    """Neutrino masses [meV] from the Z_3 ring, and the splitting ratio."""
    delta = N2**2 * (1 - ALPHA)**2                 # chiral phase, per-leg dressed
    phase = 2 * np.pi / 3 + delta
    lv = np.sort([1 + 2 * hA * np.cos(phase + 2 * np.pi * j / 3)
                  for j in range(3)])
    m1 = ALPHA**3 * M_E / (4 * np.pi**2) * 1e9     # meV
    m = m1 * lv / lv[0] * scale
    R = (m[1]**2 - m[0]**2) / (m[2]**2 - m[0]**2)
    return m, R


def main():
    print("mu_eff(k) =", derive_mu_eff(), "\n")

    R_meas = DM21 / (DM31 * 100)
    sig = R_meas * np.hypot(S21 / DM21, S31 / DM31)
    print(f"measured R = {R_meas:.6f} +- {sig:.6f}\n")

    k_hop = 2 * np.pi / D_PARTIAL                  # the Shockley step
    print(f"{'case':38s} {'k/q':>6s} {'N2eff/N2':>9s} {'h/A':>9s} "
          f"{'R':>10s} {'pull':>7s}")
    for lab, k in [("k -> oo   (leading order)", None),
                   ("2 pi / ell  (lattice period)", 2 * np.pi),
                   ("2 pi / d    (Shockley step)", k_hop)]:
        hA = h_over_A(k)
        _, R = spectrum(hA)
        share = 1.0 if k is None else rotational_share(k) / N2
        print(f"{lab:38s} {np.inf if k is None else k/Q:6.2f} {share:9.4f} "
              f"{hA:9.6f} {R:10.6f} {(R - R_meas)/sig:+7.2f}")

    # With the shape settled, whatever is left must be a common scale factor.
    hA = h_over_A(k_hop)
    m, R = spectrum(hA)
    w21 = np.sqrt(DM21 / ((m[1]**2 - m[0]**2) * 1e-1)) - 1
    w31 = np.sqrt(DM31 / ((m[2]**2 - m[0]**2) * 1e-3)) - 1
    print(f"\nresidual after the shape is corrected:")
    print(f"  scale wanted by Dm21 = {w21:+.4%}")
    print(f"  scale wanted by Dm31 = {w31:+.4%}")
    print(f"  agreement between them: {abs(w21 - w31):.2e}"
          f"  -> a pure scale residual, as it must be if the shape is right")
    print(f"  the common shift is {w21/ALPHA:.2f} alpha, against a budget of"
          f" 2 x 2 alpha on theta_ch")

    m_s, _ = spectrum(hA, 1 + w21)
    print(f"\nwith that scale applied (NOT adopted; theta_ch NLO is uncomputed):")
    print(f"  m = {np.round(m_s, 4)} meV,  sum = {m_s.sum():.2f} meV")
    print(f"  Dm21 = {(m_s[1]**2 - m_s[0]**2)*1e-1:.4f}  (meas {DM21})")
    print(f"  Dm31 = {(m_s[2]**2 - m_s[0]**2)*1e-3:.4f}  (meas {DM31})")

    # How sharply does k pick itself out?  This is a test, not a robust fit.
    print("\nsensitivity: R against the assumed hop wavevector")
    for f in (0.8, 0.9, 1.0, 1.1, 1.2):
        _, R = spectrum(h_over_A(f * k_hop))
        print(f"  k = {f:.1f} x 2pi/d : R = {R:.6f}, pull {(R-R_meas)/sig:+.2f}")


if __name__ == "__main__":
    main()
