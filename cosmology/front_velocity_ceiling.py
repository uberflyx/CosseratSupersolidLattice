#!/usr/bin/env python3
"""
front_velocity_ceiling.py
=========================
The crystallisation-front velocity from lattice mechanics.

This script carries the numerical checks behind the derivation in the cosmology
chapter ("The growing crystal and the cosmological redshift"):

1. Kinetic capacity. The Wilson-Frenkel attachment capacity is the Debye
   attempt rate times the hop distance, omega_D * ell = pi c (exact, from
   Theta_D = pi m0 and the bootstrap m0 c^2 = hbar c / ell). It exceeds every
   relativistic ceiling, so the front is transport-limited, not
   kinetics-limited (Turnbull's collision-limited regime).

2. Newborn bath temperature. Depositing the latent heat (one node mass m0 per
   node) into the Debye-cut gapless bath gives the newborn temperature, with
   and without second sound (1.592 c) in the bath.

3. The Chapman-Jouguet ceiling. With the crystal and condensate held still
   (rigid downstream), the pressure rarefaction that fixes the Jouguet point
   travels through the gas of thermal quanta. Its sound speed is
       u^2 = T S^2 / (rho_n C_v),  rho_n,i = 4 eps_i / (3 v_i^2),
   which for branches with eps_i ~ T^4 reduces to
       u^2 = (1/3) sum(eps_i) / sum(eps_i / v_i^2).
   Transverse quanta alone give c/sqrt(3); second sound in the bath lifts it.
   First sound is thermally unpopulated (weight (c/v1)^3 ~ 1e-60). A steady,
   unsupported detonation with the products at rest is sonic with respect to
   them (Steinhardt 1982), so v_front = u. The flowing-downstream Jouguet
   speed (Espinosa et al. 2010) is evaluated for comparison.

Units: c = 1; energies in m0 c^2 = 70.025 MeV; kB = hbar = 1; ell = 1.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

m0 = 70.025            # MeV, node rest energy m0 = m_e / alpha
v2 = 1.591558          # second sound speed [c], sqrt((3-4N2)/(1-N2)) at N2=1/pi
CS = 1.0 / np.sqrt(3)  # sound speed of a gas of quanta at c [c]

# ---------------------------------------------------------------------------
# 1. kinetic capacity
# ---------------------------------------------------------------------------
print("=== 1. Wilson-Frenkel kinetic capacity ===")
v_kin = np.pi   # omega_D * ell = (pi c / ell) * ell
print(f"v_kin = omega_D * ell = pi c = {v_kin:.4f} c  >>  c/sqrt(3) = {CS:.4f} c")
print("kinetics cannot limit the front; it saturates the transport ceiling\n")

# ---------------------------------------------------------------------------
# 2. newborn bath temperature
# ---------------------------------------------------------------------------
print("=== 2. bath temperature from latent heat = m0 per node ===")


def u_branch(T, v, npol=1):
    """Debye energy per node volume for a gapless branch of speed v [c],
    frequency cutoff omega_D = pi (c/ell); units m0 = ell = 1."""
    xD = np.pi / T
    integ = quad(lambda x: x**3 / np.expm1(x), 0, xD)[0]
    return npol / (2 * np.pi**2 * v**3) * T**4 * integ


def bath(with_second_sound):
    tot = lambda T: u_branch(T, 1.0, 2) + (u_branch(T, v2, 1) if with_second_sound else 0.0)
    T = brentq(lambda T: tot(T) - 1.0, 0.05, np.pi)
    share = u_branch(T, v2, 1) / tot(T) if with_second_sound else 0.0
    return T, share


T_ss, share_ss = bath(True)
T_tr, _ = bath(False)
print(f"two transverse + second sound: T = {T_ss:.3f} m0 = {T_ss * m0:.1f} MeV, "
      f"second-sound share {share_ss * 100:.2f}%")
print(f"two transverse only:           T = {T_tr:.3f} m0 = {T_tr * m0:.1f} MeV")
print("(melting band: stacking transition 156.2 MeV, crystallisation 220 MeV)\n")

# ---------------------------------------------------------------------------
# 3. the Chapman-Jouguet ceiling
# ---------------------------------------------------------------------------
print("=== 3. exhaust sound speed and the CJ ceiling ===")


def u_gas(shares_speeds):
    """Sound speed of a quasiparticle gas: u^2 = (1/3) sum eps_i / sum(eps_i/v_i^2)."""
    num = sum(s for s, v in shares_speeds)
    den = sum(s / v**2 for s, v in shares_speeds)
    return np.sqrt(num / (3.0 * den))


cos2 = quad(lambda th: np.cos(th)**2 * np.sin(th), 0, np.pi)[0] / 2.0
print(f"isotropy check: <cos^2 theta> = {cos2:.6f} (the 1/3 of P = u/3)")
print(f"first sound ~1e20 c: thermal weight (c/v)^3 ~ 1e-60, no rarefaction")
print(f"transverse quanta alone:              u = {u_gas([(1.0, 1.0)]):.4f} c")
for s in (share_ss, 0.25):
    print(f"second sound at {s * 100:4.1f}% of the bath:   u = {u_gas([(1 - s, 1.0), (s, v2)]):.4f} c")


def v_jouguet_flowing(alpha):
    """Standard relativistic Jouguet speed for a FLOWING downstream
    (Steinhardt 1982; Espinosa, Konstandin, No & Servant 2010), cs^2 = 1/3."""
    return (np.sqrt(alpha * (2.0 + 3.0 * alpha)) + 1.0) / (np.sqrt(3.0) * (1.0 + alpha))


print("\nflowing-downstream Jouguet speed vs the rigid-downstream value:")
for a in (1e-4, 0.03, 0.1, 0.3, 1.0):
    print(f"  alpha = {a:6g}:  v_J(flowing) = {v_jouguet_flowing(a):.4f} c")
print(f"  rigid downstream (lattice, no flow): v_front = u at every alpha")
print(f"heat clearing: v2 / (c/sqrt3) = {v2 / CS:.2f} (second sound, if present, outruns the front)")
print("\nconclusion: v_front = c/sqrt(3) with transverse quanta alone, 0.60-0.63 c if second")
print("sound shares the bath; massive (optical) branches can only lower it.")
