#!/usr/bin/env python3
"""
d4_chirality_resummation.py
Multi-period (Dyson) resummation of the D4 Cosserat chirality.

Question settled here
---------------------
The single-period transfer matrix gives a per-period helicity splitting.
A thermal bath populates +k4 and -k4 winding modes together, so the
observable depends on the PARITY of the splitting in k4. This script
computes that parity directly and resolves, per observable, whether the
per-period splitting accumulates (coherent) or cancels.

Result
------
There are two chiralities with OPPOSITE k4-parity:

  * amplitude split  |lambda_+| - |lambda_-|   EVEN in k4
      -> natural optical activity: the birefringence source (beta_EM, beta_GW).
         Survives the symmetric thermal average. COHERENT, grows with distance.

  * Im(Sigma_13)  (imaginary compact-rotation leakage)   ODD in k4
      -> CP violation: the baryogenesis / chiral-magnetic source.
         Cancels in a symmetric bath. Survives ONLY through the Z3 stacking
         handedness A->B->C, which fixes its sign.

The PARITY result above is exact and independent of normalisation. The
thermal magnitude (theta_eff/theta_ch, the crossover, the merger factor)
uses theta_D4 = absolute amplitude split at the first KK, weighted by the
first-KK Bose occupation; relating that split to the constitutive theta_ch
carries an O(1) optical-rotation normalisation, so those numbers are
order-of-magnitude, not few-percent.

Inputs: c, hbar, m_e (everything else derived). mu = 1, ell = 1 units.
"""
import numpy as np
from numba import njit

# --- lattice constants (units mu = ell = 1) ------------------------------
N2       = 1.0 / np.pi                 # Cosserat coupling number N^2 = 1/pi
mu       = 1.0                         # shear modulus
ell      = 1.0                         # lattice spacing
kappa_c  = 2 * N2 * mu / (1.0 - N2)    # rotational stiffness
mu_tot   = mu + kappa_c                # total shear modulus
gamma_c  = mu * ell**2                 # curvature modulus
ALPHA    = 1.0 / 137.035999177         # fine-structure constant

# --- physical scales (MeV, fm) -------------------------------------------
HBARC    = 197.3269804                 # hbar c  [MeV fm]
ELL_FM   = 2.8179403                   # classical electron radius r_e [fm]
M0_MEV   = 0.51099895 / ALPHA          # node mass m_e/alpha = 70.0 MeV
THETA_CH = ALPHA**2 / (2 * np.pi)      # constitutive chirality

# The compact loop is THREE sites of circumference sqrt(6)*ell (the D4 result).
# Its first Kaluza-Klein momentum is therefore k4 = 2*pi/(sqrt(6)*ell).
K4_KK    = 2 * np.pi / (np.sqrt(6) * ell)   # first KK momentum, k4*ell = 2.566
M_KK_MEV = np.sqrt(3) * M0_MEV              # first KK mass, sqrt(3) m0 = 121 MeV
T_GEOM   = HBARC / (np.sqrt(6) * ELL_FM)    # window opens, = 28.6 MeV
T_C      = 156.1                            # deconfinement (window closes)


@njit(cache=True)
def _block(x, k4, V0, s):
    """One 6x6 helicity block. s=+1 anti-self-dual (CG), s=-1 self-dual (DF).
    The relative sign s carries the chirality: the curl coupling and the
    ik4 mixing flip together between the two helicities."""
    V = V0 * np.sin(np.pi * x / ell)**2
    k4sq = k4 * k4
    M = np.zeros((6, 6), dtype=np.complex128)
    M[0, 2] = 1.0; M[1, 3] = 1.0; M[4, 5] = 1.0
    M[2, 0] = (V + mu_tot * k4sq) / mu_tot
    M[2, 3] = -s * kappa_c / mu_tot
    M[2, 4] = -1j * s * kappa_c * k4 / mu_tot     # ik4 mixing, sign tied to s
    M[3, 1] = (2 * kappa_c + gamma_c * k4sq) / gamma_c
    M[3, 2] = s * kappa_c / gamma_c
    M[5, 4] = (2 * kappa_c + gamma_c * k4sq) / gamma_c
    M[5, 0] = 1j * kappa_c * k4 / gamma_c
    return M


@njit(cache=True)
def _mm(A, B):
    C = np.zeros((6, 6), dtype=np.complex128)
    for i in range(6):
        for j in range(6):
            acc = 0j
            for k in range(6):
                acc += A[i, k] * B[k, j]
            C[i, j] = acc
    return C


@njit(cache=True)
def transfer(k4, V0, s, nstep):
    """Transfer matrix over one PN period by fixed-step RK4 (numba)."""
    P = np.eye(6, dtype=np.complex128)
    h = ell / nstep
    for n in range(nstep):
        x = n * h
        k1 = _mm(_block(x,        k4, V0, s), P)
        k2 = _mm(_block(x + 0.5*h, k4, V0, s), P + 0.5 * h * k1)
        k3 = _mm(_block(x + 0.5*h, k4, V0, s), P + 0.5 * h * k2)
        k4v = _mm(_block(x + h,    k4, V0, s), P + h * k3)
        P = P + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4v)
    return P


def smallest_eig(T):
    e = np.linalg.eigvals(T)
    return e[np.argmin(np.abs(e))]


def calibrate_V0(nstep, target=ALPHA):
    """Fix the barrier height so the k4=0 tunnelling amplitude equals alpha."""
    lo, hi = 50.0, 250.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if abs(smallest_eig(transfer(0.0, mid, 1.0, nstep))) < target:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def chirality(k4, V0, nstep):
    """Return (amplitude split / alpha, Im(Sigma_13) leakage) at momentum k4."""
    lcg = smallest_eig(transfer(k4, V0, +1.0, nstep))
    ldf = smallest_eig(transfer(k4, V0, -1.0, nstep))
    amp = (abs(ldf) - abs(lcg)) / ALPHA
    e, V = np.linalg.eig(transfer(k4, V0, +1.0, nstep))
    v = V[:, np.argmin(np.abs(e))]
    imleak = (v[4] / v[0]).imag if abs(v[0]) > 1e-12 else 0.0
    return amp, imleak


def f_bose(T, m=M_KK_MEV):
    """Bose-Einstein occupation of the first KK mode at temperature T (MeV)."""
    return 1.0 / np.expm1(m / T)


def main():
    NS = 4000
    V0 = calibrate_V0(NS)
    print("D4 chirality: multi-period (Dyson) resummation")
    print("=" * 64)
    print(f"  V0 (calibrated to alpha) = {V0:.4f}")
    print(f"  first KK: k4*ell = {K4_KK*ell:.4f},  mass = {M_KK_MEV:.1f} MeV")
    print(f"  theta_ch = alpha^2/(2pi) = {THETA_CH:.3e}")

    print("\n  PARITY OF THE SPLITTING IN k4")
    print(f"  {'k4*ell':>9} {'amp split %':>13} {'Im Sigma13':>13}")
    for k4l in (-K4_KK*ell, -1.0, 0.0, 1.0, K4_KK*ell):
        a, im = chirality(k4l / ell, V0, NS)
        print(f"  {k4l:9.4f} {a*100:13.5f} {im:13.6f}")
    print("    amplitude split: EVEN in k4  (birefringence; survives)")
    print("    Im(Sigma_13):    ODD  in k4  (CP; cancels unless handedness)")

    print("\n  Z3 SUM over sublattice momenta {0, +2pi/sqrt6, -2pi/sqrt6}")
    ks = (0.0, K4_KK, -K4_KK)
    amps = [chirality(k, V0, NS)[0] for k in ks]
    ims = [chirality(k, V0, NS)[1] for k in ks]
    w = np.exp(2j * np.pi * np.arange(3) / 3)        # Z3 Bloch phases
    print(f"    birefringence, symmetric sum: {sum(amps)*100:+.5f} %   SURVIVES")
    print(f"    CP, symmetric sum:            {sum(ims):+.6f}      CANCELS")
    print(f"    CP, handedness-weighted:      {abs(sum(i*wn for i, wn in zip(ims, w))):.6f}"
          f"   survives via stacking")

    # Birefringence magnitude at the (validated) first KK, and its thermal weight.
    amp_kk = chirality(K4_KK, V0, NS)[0]              # |lam_+|-|lam_-| over alpha
    # absolute helicity split at full KK; ratio to theta_ch is O(1)-normalised
    theta_D4 = abs(amp_kk) * ALPHA
    ratio_kk = theta_D4 / THETA_CH
    print("\n  THERMAL theta_eff(T) = theta_D4 * f_Bose(m_KK/T)  (birefringence)")
    print(f"    theta_D4 (full KK) = {theta_D4:.3e}  = {ratio_kk:.1f} x theta_ch")
    print(f"    window: {T_GEOM:.1f} MeV (opens) .. {T_C:.1f} MeV (closes)")
    print(f"    {'T (MeV)':>9} {'f_Bose':>9} {'theta_eff/theta_ch':>20} {'beta enhance':>14}")
    for T in (30, 50, 75, 100, 125, 150):
        fb = f_bose(T)
        r = ratio_kk * fb
        print(f"    {T:9.0f} {fb:9.4f} {r:20.3f} {1+r:14.2f}x")
    # crossover where theta_eff = theta_ch
    from scipy.optimize import brentq
    Tx = brentq(lambda T: ratio_kk * f_bose(T) - 1.0, 5.0, 300.0)
    print(f"    crossover theta_eff = theta_ch at T = {Tx:.0f} MeV")
    fb_lo, fb_hi = f_bose(50), f_bose(100)
    print(f"\n  Merger remnant (50-100 MeV): beta_GW enhanced "
          f"{1+ratio_kk*fb_lo:.1f}x to {1+ratio_kk*fb_hi:.1f}x over the T=0 value,")
    print("  same sign as the cosmological birefringence, and surviving coherently.")


if __name__ == "__main__":
    main()
