#!/usr/bin/env python3
"""
front_velocity_ceiling.py
=========================
The crystallisation-front velocity from lattice mechanics: v_front = c/sqrt(3).

This script carries the three numerical checks behind the derivation in the
cosmology chapter ("The growing crystal: expansion without stretching"):

1. Kinetic capacity. The Wilson-Frenkel attachment capacity is the Debye
   attempt rate times the hop distance, omega_D * ell = pi c (exact, from
   Theta_D = pi m0 and the bootstrap m0 c^2 = hbar c / ell). It exceeds every
   relativistic ceiling, so the front is transport-limited, not
   kinetics-limited (Turnbull's collision-limited regime at maximal
   undercooling).

2. Channel census and the Chapman-Jouguet collapse. Of the exhaust's three
   longitudinal channels, only the massless transverse-phonon (radiation) gas
   carries the bulk pressure P = u/3; first sound is thermally unpopulated
   and second sound is a counterflow (temperature, not pressure). A steady,
   unsupported detonation with the products at rest behind it is exactly
   sonic with respect to them (Steinhardt 1982), and the rigid lattice
   downstream (no flow) collapses the standard Jouguet speed v_J(alpha) to
   its alpha -> 0 value c/sqrt(3) at every transition strength alpha. The
   flowing-downstream formula (Espinosa et al. 2010) is evaluated for
   comparison.

3. Bath temperature. Depositing the latent heat (one node mass m0 per node)
   into the gapless bath (two transverse branches at c, plus second sound at
   3.65 c) gives a newborn temperature of ~131 MeV, below the 156-220 MeV
   melting band: the front is never heat-choked, and second sound clears the
   latent heat backward at 3.65 c, ~6x faster than the front lays it down.

Units: c = 1; energies in m0 c^2 = 70.025 MeV; kB = hbar = 1; ell = 1.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

m0 = 70.025            # MeV, node rest energy m0 = m_e / alpha
v2 = 3.65              # second sound speed [c]
CS = 1.0 / np.sqrt(3)  # radiation-gas sound speed [c]

# ---------------------------------------------------------------------------
# 1. kinetic capacity
# ---------------------------------------------------------------------------
print("=== 1. Wilson-Frenkel kinetic capacity ===")
v_kin = np.pi   # omega_D * ell = (pi c / ell) * ell
print(f"v_kin = omega_D * ell = pi c = {v_kin:.4f} c  >>  c/sqrt(3) = {CS:.4f} c")
print("kinetics cannot limit the front; it saturates the transport ceiling\n")

# ---------------------------------------------------------------------------
# 2. channel census and the Chapman-Jouguet collapse
# ---------------------------------------------------------------------------
print("=== 2. exhaust channels and the CJ ceiling ===")
print(f"first sound  ~1e20 c : carries pressure, thermal weight (c/v)^3 ~ 1e-60")
print(f"second sound {v2} c : counterflow, carries temperature, not pressure")
print(f"radiation gas        : P = u/3  ->  c_s = 1/sqrt(3) = {CS:.6f} c")
# isotropy check: <cos^2 theta> = 1/3
cos2 = quad(lambda th: np.cos(th)**2 * np.sin(th), 0, np.pi)[0] / 2.0
print(f"isotropy check: <cos^2 theta> = {cos2:.6f} (the 1/3 of P = u/3)")

def v_jouguet_flowing(alpha):
    """Standard relativistic Jouguet speed for a FLOWING downstream
    (Steinhardt 1982; Espinosa, Konstandin, No & Servant 2010), cs^2 = 1/3."""
    return (np.sqrt(alpha * (2.0 + 3.0 * alpha)) + 1.0) / (np.sqrt(3.0) * (1.0 + alpha))

print("\nflowing-downstream Jouguet speed vs the rigid-downstream value:")
for a in (1e-4, 0.03, 0.1, 0.3, 1.0):
    print(f"  alpha = {a:6g}:  v_J(flowing) = {v_jouguet_flowing(a):.4f} c")
print(f"  rigid downstream (lattice, no flow): v_front = {CS:.4f} c at every alpha\n")

# ---------------------------------------------------------------------------
# 3. newborn bath temperature
# ---------------------------------------------------------------------------
print("=== 3. bath temperature from latent heat = m0 per node ===")

def u_branch(T, v, npol=1):
    """Debye energy per node volume for a gapless branch of speed v [c],
    Debye cutoff omega_D = pi (c/ell); units m0 = ell = 1."""
    xD = np.pi / T
    integ = quad(lambda x: x**3 / np.expm1(x), 0, xD)[0]
    return npol / (2 * np.pi**2 * v**3) * T**4 * integ

def u_total(T):
    return u_branch(T, 1.0, npol=2) + u_branch(T, v2, npol=1)

T_bath = brentq(lambda T: u_total(T) - 1.0, 0.05, np.pi)
share_2nd = u_branch(T_bath, v2, 1) / u_total(T_bath)
print(f"T_bath = {T_bath:.3f} m0 = {T_bath * m0:.1f} MeV  (melting band 156-220 MeV)")
print(f"second-sound share of the bath energy: {share_2nd * 100:.2f}%")
print(f"heat clearing: v2 / v_front = {v2 / CS:.2f} (second sound outruns the front)")
print("\nconclusion: v_front = c/sqrt(3) exactly at leading order;")
print("massive (optical) branches can only lower it, so it is an upper bound.")
