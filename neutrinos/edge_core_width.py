"""Core half-width of the edge dislocation from the Cosserat Peierls-Nabarro balance.

The screw (electron) core half-width follows from the canonical PN equilibrium

    integral_0^inf  (K(k)/k) exp(-2 k w) dk  =  mu_bar / (2 d_111),

with K(k) the anti-plane Cosserat surface kernel.  The neutrino is built from
edge-dislocation modes, and in plane strain the elastic kernel of the slip
plane carries the classical Poisson enhancement 1/(1 - nu).  The vacuum is
near-incompressible (nu -> 1/2), so the edge kernel is 2 K(k) to all relevant
digits.  The Frenkel misfit stiffness on the right is a property of the glide
plane and is common to both characters.

The script verifies the two known limits before producing the new number:
  (i)  Cauchy screw limit (N^2 = 0): w/d = 1/sqrt(2) analytically;
  (ii) Cosserat screw at the vacuum point N^2 = 1/pi: w = 0.454 ell
       (w/d = 0.78616 with d = ell/sqrt(3) the partial period),
       the value quoted in the alpha chapter.
It then solves the same equilibrium with the doubled kernel for the edge.

Units: lengths in ell (nearest-neighbour spacing), moduli in the Eringen
shear modulus mu.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

# Vacuum-point Cosserat constants (Eringen units, mu = 1)
N2 = 1.0 / np.pi                       # Cosserat coupling number
KAPPA_C = 2.0 * N2 / (1.0 - 2.0 * N2)  # = 2/(pi-2) = 1.752
MU = 1.0
MU_TOT = MU + KAPPA_C                  # short-wavelength shear modulus, 2.752
MU_BAR = MU + KAPPA_C / 2.0            # long-wavelength shear modulus, 1.876
GAMMA = MU * 1.0                       # curvature modulus, gamma = mu ell^2
P2 = KAPPA_C * (2 * MU + KAPPA_C) / ((MU + KAPPA_C) * GAMMA)  # 2.389 / ell^2
Q2 = (MU_TOT / MU_BAR) * P2                                   # 3.504 / ell^2

D111 = np.sqrt(2.0 / 3.0)              # {111} interplanar spacing, 0.8165 ell
D_PARTIAL = 1.0 / np.sqrt(3.0)         # Shockley partial period, 0.5774 ell
RHS = MU_BAR / (2.0 * D111)            # Frenkel misfit stiffness, 1.149 mu/ell


def kernel(k, factor=1.0):
    """Cosserat surface kernel K(k); factor=2 gives the edge (nu = 1/2)."""
    return factor * 0.5 * MU_TOT * k * (k * k + P2) / (k * k + Q2)


def lhs(w, factor=1.0):
    """Elastic restoring integral of the PN equilibrium at core half-width w."""
    val, _ = quad(lambda k: kernel(k, factor) / k * np.exp(-2.0 * k * w),
                  0.0, np.inf, limit=200)
    return val


def solve_width(factor=1.0, rhs=RHS):
    """Core half-width w (in ell) balancing elastic pull against misfit pinning."""
    return brentq(lambda w: lhs(w, factor) - rhs, 0.05, 5.0, xtol=1e-12)


if __name__ == "__main__":
    # Check (i): Cauchy screw limit, kernel k/2, analytic w = mu/(4 RHS_c)
    rhs_cauchy = 1.0 / (2.0 * D111)
    w_c = brentq(lambda w: 0.25 / w - rhs_cauchy, 0.05, 5.0, xtol=1e-14)
    print(f"Cauchy screw:  w = {w_c:.5f} ell = {w_c / D_PARTIAL:.5f} d "
          f"(analytic 1/sqrt2 = {1/np.sqrt(2):.5f})")

    # Check (ii): Cosserat screw at the vacuum point
    w_s = solve_width(factor=1.0)
    print(f"Cosserat screw: w = {w_s:.5f} ell = {w_s / D_PARTIAL:.5f} d "
          f"(alpha chapter: 0.454 ell, 0.78616 d)")

    # New number: Cosserat edge (nu = 1/2, kernel doubled)
    w_e = solve_width(factor=2.0)
    print(f"Cosserat edge:  w = {w_e:.5f} ell = {w_e / D_PARTIAL:.5f} d")
    print(f"edge/screw ratio = {w_e / w_s:.5f} (classical value 1/(1-nu) = 2)")

    ELL_FM = 2.8179403205                      # lattice spacing = r_e [fm]
    print(f"\nneutrino transverse half-width L_b = {w_e * ELL_FM:.3f} fm "
          f"(was quoted 0.78 ell = {0.78 * ELL_FM:.2f} fm)")


def make_figure(path="edge_core_balance.pdf"):
    """Plot the PN balance: elastic restoring integrals against the misfit line."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    w_grid = np.linspace(0.12, 1.6, 300)
    lhs_screw = np.array([lhs(w, 1.0) for w in w_grid])
    lhs_edge = np.array([lhs(w, 2.0) for w in w_grid])
    w_s, w_e = solve_width(1.0), solve_width(2.0)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(w_grid, lhs_screw, color="#1f77b4", lw=1.8,
            label=r"screw kernel $\hat K(k)$")
    ax.plot(w_grid, lhs_edge, color="#922b21", lw=1.8,
            label=r"edge kernel $2\hat K(k)$  ($\nu \to 1/2$)")
    ax.axhline(RHS, color="0.25", lw=1.2, ls="--",
               label=r"misfit pinning $\bar\mu/(2d_{111})$")
    for w0, c in [(w_s, "#1f77b4"), (w_e, "#922b21")]:
        ax.plot([w0], [RHS], "o", color=c, ms=6, zorder=5)
        ax.vlines(w0, 0, RHS, color=c, lw=0.8, ls=":")
    ax.annotate(r"$w_{\rm screw} = 0.454\,\ell$", (w_s, RHS),
                xytext=(w_s + 0.05, RHS + 0.55), fontsize=10, color="#1f77b4")
    ax.annotate(r"$w_{\rm edge} = 0.860\,\ell$", (w_e, RHS),
                xytext=(w_e + 0.05, RHS + 0.55), fontsize=10, color="#922b21")
    ax.set_xlabel(r"core half-width $w$  [$\ell$]")
    ax.set_ylabel(r"elastic restoring integral  [$\mu/\ell$]")
    ax.set_xlim(0.12, 1.6); ax.set_ylim(0, 5.2)
    ax.legend(frameon=False, fontsize=9.5)
    fig.tight_layout()
    fig.savefig(path)
    print(f"wrote {path}")


if __name__ == "__main__":
    make_figure()
