#!/usr/bin/env python3
"""
crystallisation_baryogenesis.py
===============================================================================
Order-of-magnitude estimate of the cosmic baryon asymmetry from the
crystallisation of the four-dimensional vacuum lattice.

Mechanism (the three Sakharov conditions, all met by the lattice):

  1. Baryon-number violation is geometric.  Above the deconfinement
     temperature T_c the compact (stacking) direction is disordered, the Z_3
     stacking symmetry is restored, and the winding number that defines baryon
     number is undefined.  Below T_c the winding is a topological invariant and
     baryon number is exact.  The crystallisation at T_c freezes the windings
     in: the net frozen winding IS the net baryon number.  No sphaleron and no
     grand-unified decay is needed, and there is no wash-out afterwards because
     a topological invariant cannot relax.

  2. Departure from equilibrium is the transition itself, which is first order
     (confirmed by the FCC 3-state Potts Monte Carlo: latent heat, two-state
     signal).  Bubble nucleation drives the freeze-in out of equilibrium.

  3. C and CP violation is the propagation chirality Im(Sigma_13) ~ alpha*k_4,
     a time-odd optical activity of the compact direction that biases the
     forward winding A->B->C (a baryon) over the reverse C->B->A (an
     antibaryon).  The bias is geometric, so its SIGN is fixed by the stacking
     handedness: the framework predicts matter over antimatter, not merely an
     asymmetry of unknown direction, and ties that sign to the chiral-magnetic
     charge correlator measured in heavy-ion collisions.

The net baryon-to-photon ratio is

        eta_B  ~  theta_eff * kappa * f_eff

with theta_eff the CP bias per frozen winding, kappa a geometric prefactor
(winding-carrying quark density over entropy), and f_eff the net efficiency of
the transition (the fraction of the biased winding that survives the bubble
wall).  theta_eff and kappa are fixed by the framework; f_eff and the exact
power of alpha in theta_eff are the two open pieces, set respectively by the
bubble-wall transport problem and the multi-period Dyson resummation of the
transfer matrix.  This script makes the estimate explicit, brackets it across
the physical range of freeze-in wavevectors (crossover to zone boundary),
attaches the early-universe epoch, and maps how the result depends on those
two open inputs.

No numba: every quantity here is closed-form arithmetic plus one light 1-D
integral and a small parameter scan, so a JIT compiler would buy nothing.

All inputs are framework quantities; none is fitted to the observed eta_B.
===============================================================================
"""

import numpy as np
from scipy.special import zeta

# ----------------------------------------------------------------------------
# Framework inputs (no fitting).  ell = r_e is the sole dimensionful scale.
# ----------------------------------------------------------------------------
ALPHA   = 1.0 / 137.035999084          # fine-structure constant
M_E     = 0.51099895069                # electron rest energy [MeV]
M0      = M_E / ALPHA                  # node rest energy = hbar c / ell  [MeV]

TC      = 156.0                        # deconfinement temperature [MeV] (Polyakov/Potts)
TGEOM   = M0 / np.sqrt(6.0)            # geometric temperature, hbar c/(sqrt6 ell) [MeV]
K4ELL_X = 1.47                         # compact wavevector at the chirality crossover
K4ELL_ZB = np.pi                       # zone-boundary wavevector (lattice scale, k_4 ell = pi)

# relativistic degrees of freedom just above T_c
G_QUARK = 12.0                         # u,d quark dof carrying the winding (2 flavour x 2 spin x 3 colour)
G_STAR  = 17.25                        # total relativistic dof (quarks, gluons, photons, leptons)

# observations
ETA_OBS  = 6.1e-10                     # baryon-to-photon ratio
NGAMMA_OVER_S = 1.0 / 7.04             # photon-to-entropy ratio today
M_PL     = 1.22e22                     # Planck mass [MeV] (sets the radiation-era expansion rate)

# Standard Model benchmark for context
DELTA_CKM_EFF = 1.0e-20                # effective Jarlskog CP measure


# ----------------------------------------------------------------------------
# Building blocks
# ----------------------------------------------------------------------------
def cp_phase(power=2.0, k4ell=K4ELL_X):
    """CP bias per frozen winding.

    theta_eff is the square of the chiral mixing amplitude v_mix ~ alpha*(k4 ell),
    so theta_eff ~ (alpha k4 ell)^2 = alpha^2 * (k4 ell)^2 at leading order.
    'power' carries the open question of the multi-period resummation: if the
    per-period chirality partly cancels, the effective power of alpha rises
    above 2 and theta_eff shrinks.
    """
    return (ALPHA ** power) * (k4ell ** 2)


def geometric_prefactor(g_quark=G_QUARK, g_star=G_STAR):
    """kappa = n_quark / s.

    The windings are supplied by the quark population of equilibrium number
    density n_q = (3 zeta(3)/4 pi^2) g_q T^3 (the 3/4 is the fermion factor),
    and the comoving asymmetry is carried as the ratio to the entropy density
    s = (2 pi^2/45) g_* T^3.  The common T^3 cancels, so kappa is a pure number.
    """
    n_over_s = (zeta(3) / np.pi**2) * 0.75 * g_quark / ((2 * np.pi**2 / 45.0) * g_star)
    return n_over_s


def eta_B(theta_eff, kappa, f_eff=1.0):
    """Baryon-to-photon ratio. Y_B = n_B/s is conserved; convert to eta_B today."""
    Y_B = theta_eff * kappa * f_eff          # baryon-to-entropy ratio at freeze-in
    return Y_B / NGAMMA_OVER_S


def window_integral(t_geom, t_c=TC):
    """Cooling-window weight for the chirality, theta_eff(T) ~ T^2, normalised to T_c^3.

    Used only to show that the L_4 correction (which moves the lower edge T_geom)
    is immaterial to the asymmetry: the integral barely changes.
    """
    return (t_c**3 - t_geom**3) / (3.0 * t_c**3)


def cp_phase_bracket(power=2.0):
    """CP bias at the two physical reference wavevectors.

    The single-point estimate uses the chirality crossover, k_4 ell = 1.47.
    But the windings freeze in at T_c on the scale of the partial-dislocation
    core, of size ell, so the freeze-in samples wavevectors up to the zone
    boundary k_4 ell = pi.  theta_eff is therefore a bracket, not a point, and
    so are eta_B and the efficiency the data demand.
    """
    return cp_phase(power, K4ELL_X), cp_phase(power, K4ELL_ZB)


def hubble_time(T, g_star=G_STAR):
    """Age of the radiation-dominated universe at temperature T [MeV], in seconds.

    From H = 1.66 sqrt(g_*) T^2 / M_Pl and t = 1/(2H), with the numerical
    prefactor folded in: t[s] = 2.42 g_*^{-1/2} (T/MeV)^{-2}.  Fixes the epoch
    at which crystallisation baryogenesis occurs.
    """
    return 2.42 / (np.sqrt(g_star) * T**2)


# ----------------------------------------------------------------------------
# Main estimate
# ----------------------------------------------------------------------------
def main():
    kappa = geometric_prefactor()
    theta = cp_phase(power=2.0)
    eta_naive = eta_B(theta, kappa, f_eff=1.0)
    overshoot = eta_naive / ETA_OBS
    f_required = ETA_OBS / eta_naive

    print("=" * 70)
    print("Crystallisation baryogenesis: order-of-magnitude estimate")
    print("=" * 70)
    print(f"  node energy m0 c^2            = {M0:8.3f} MeV")
    print(f"  geometric temperature T_geom = {TGEOM:8.3f} MeV   (= m0 c^2 / sqrt6)")
    print(f"  deconfinement      T_c       = {TC:8.3f} MeV")
    print()
    print(f"  CP bias per winding  theta_eff = (alpha*{K4ELL_X})^2 = {theta:.3e}")
    print(f"  geometric prefactor  kappa     = n_q/s             = {kappa:.4f}")
    print()
    print(f"  eta_B (f_eff = 1, naive)       = {eta_naive:.3e}")
    print(f"  observed eta_B                 = {ETA_OBS:.3e}")
    print(f"  naive / observed (overshoot)   = {overshoot:.2e}")
    print(f"  -> required net efficiency f_eff = {f_required:.2e}")
    print()
    print(f"  framework CP supply  theta_eff ~ {theta:.1e}")
    print(f"  Standard Model CP    delta_CKM ~ {DELTA_CKM_EFF:.0e}")
    print(f"  framework exceeds SM by        ~ {theta/DELTA_CKM_EFF:.0e}")
    print("  (The SM undershoots the observed value; the framework overshoots,")
    print("   which is the recoverable direction: efficiency suppresses, it")
    print("   cannot enhance beyond the available CP supply.)")
    print()

    # Sensitivity to the two open inputs ------------------------------------
    # The exact power of alpha in theta_eff is set by the (unfinished) Dyson
    # resummation; the efficiency f_eff by the (unsolved) bubble-wall transport.
    # For each candidate power, report the f_eff that reproduces the observed value.
    print("  Required efficiency vs the open resummation power of alpha:")
    print("    p   theta_eff      eta_B(f=1)    f_eff to match observed")
    for p in (2.0, 2.5, 3.0, 3.5, 4.0):
        th = cp_phase(power=p)
        en = eta_B(th, kappa, f_eff=1.0)
        fr = ETA_OBS / en
        print(f"    {p:>3}  {th:.3e}   {en:.3e}   {fr:.2e}")
    print("  (A softening toward alpha^3 lifts the required efficiency into the")
    print("   10^-3 range typical of a first-order transition, which is more")
    print("   natural than the 10^-6 demanded at alpha^2.)")
    print()

    # L_4 dependence: the corrected circumference moves only the window edge.
    w_old = window_integral(33.0)          # old L_4 = 3 ell / sqrt2
    w_new = window_integral(TGEOM)         # corrected L_4 = sqrt6 ell
    print(f"  window weight, T_geom = 33.0 MeV : {w_old:.4f}")
    print(f"  window weight, T_geom = {TGEOM:4.1f} MeV : {w_new:.4f}")
    print(f"  -> L_4 correction shifts eta_B by {100*(w_new-w_old)/w_old:+.1f}%  (immaterial)")
    print()

    # CP-bias bracket: freeze-in samples wavevectors up to the lattice scale --
    th_x, th_zb = cp_phase_bracket(power=2.0)
    en_x, en_zb = eta_B(th_x, kappa), eta_B(th_zb, kappa)
    f_x, f_zb = ETA_OBS / en_x, ETA_OBS / en_zb
    print("  CP-bias bracket (crossover k4l=1.47  ->  zone boundary k4l=pi):")
    print(f"    theta_eff      : {th_x:.2e}  ->  {th_zb:.2e}")
    print(f"    eta_B (f=1)    : {en_x:.2e}  ->  {en_zb:.2e}")
    print(f"    required f_eff : {f_x:.2e}  ->  {f_zb:.2e}")
    print()

    # Early-universe timing: when does the freeze-in happen? -----------------
    t_c    = hubble_time(TC)
    t_geom = hubble_time(TGEOM)
    print("  Early-universe timing (radiation era):")
    print(f"    freeze-in at T_c      = {TC:5.1f} MeV  ->  t = {t_c*1e6:6.1f} us")
    print(f"    chirality off at Tgeom= {TGEOM:5.1f} MeV  ->  t = {t_geom*1e3:6.2f} ms")
    print(f"    window crossed over   ~ {(t_geom - t_c)*1e3:.2f} ms; asymmetry frozen at t ~ {t_c*1e6:.0f} us")
    print("    (topologically protected thereafter: no wash-out to the present)")
    print("=" * 70)
    print("  Prediction with teeth: a definite sign (matter over antimatter),")
    print("  fixed by the stacking handedness A->B->C and correlated with the")
    print("  sign of the heavy-ion chiral-magnetic charge correlator.")
    print("  Open pieces: the alpha-power (Dyson resummation) and f_eff")
    print("  (bubble-wall transport).")
    print("=" * 70)


if __name__ == "__main__":
    main()
