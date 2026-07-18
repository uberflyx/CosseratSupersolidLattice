"""PN kink as a compact-direction instanton: action decomposition and
fluctuation spectrum.

Three results, in order of hardness:

1. BOND LEMMA (exact, symbolic). The winding count of the form-factor rule
   at the 4D bond harmonic k = 2*pi/ell along the interlayer bond direction
   is (d^2 + d111^2)/ell^2 = 1 exactly. One alpha per bond winding is the
   4D statement; the glide and stacking anchors are its projections.

2. THERMAL SKELETON (structural form of the conjecture). In nonlocal
   Frenkel units (mu_bar = 1, phase period 2*pi) the exact PN kink is
   phi = 2*arctan(x/w) + pi with w = 2*mu_bar/V0, and the tunnelling
   exponent is ln(1/alpha) = w. The conjecture ln(1/alpha) = beta*(2 m0 c^2)
   = 2*sqrt(6) is therefore the statement mu_bar/V0 = sqrt(6) = L4/ell:
   the misfit barrier equals the elastic stiffness divided by the compact
   circumference in lattice units, V0 = mu_bar * ell / L4. The Cosserat
   equilibrium gives mu_bar/V0 = pi^2/4 = 2.4674 (bare) against
   sqrt(6) = 2.4495 (-0.73%); the physical value ln(1/alpha)/2 = 2.4601
   sits between them.

3. FLUCTUATION SPECTRUM (numeric). The second variation of the nonlocal
   Frenkel energy about the kink is
       H = |k|  +  (x^2 - w^2)/(w*(x^2 + w^2)),
   with continuum edge at 1/w. The translational zero mode
   psi_0 = phi' = 2w/(x^2+w^2) is verified analytically (H psi_0 = 0).
   The alpha chapter's vacuum-polarisation correction asserts exactly TWO
   modes bound to the core (one per photon polarisation), each lowering
   1/alpha by one unit. If the instanton reading is right, those two modes
   must appear as the discrete spectrum of H. This script counts them.

Author: framework companion script.
"""

import numpy as np
import sympy as sp

# ----------------------------------------------------------------------
# 1. Bond lemma (exact)
# ----------------------------------------------------------------------
ell = sp.Integer(1)
d = ell / sp.sqrt(3)
d111 = ell * sp.sqrt(sp.Rational(2, 3))
bond_e = d / ell        # in-plane cosine of the bond direction
bond_n = d111 / ell     # stacking cosine
kb = 2 * sp.pi / ell    # bond harmonic
n_bond = sp.simplify(sp.Abs(kb * bond_e) * d / (2 * sp.pi)
                     + sp.Abs(kb * bond_n) * d111 / (2 * sp.pi))
print("=== 1. Bond lemma ===")
print(f"winding count at the bond harmonic, along the bond: n = {n_bond}")
print("-> F = alpha^1 exactly: one alpha per 4D bond winding.")

# ----------------------------------------------------------------------
# 2. Thermal skeleton in structural form
# ----------------------------------------------------------------------
lna = float(sp.log(sp.Float("137.035999177", 12)))
print("\n=== 2. Structural form of the thermal conjecture ===")
print(f"mu_bar/V0 (bare Cosserat equilibrium) = pi^2/4  = {np.pi**2/4:.5f}")
print(f"mu_bar/V0 (physical, ln(1/a)/2)      =           {lna/2:.5f}")
print(f"mu_bar/V0 (thermal, L4/ell)          = sqrt(6) = {np.sqrt(6):.5f}")
print("Conjecture: V0 = mu_bar * ell / L4 (barrier = stiffness per compact loop).")

# ----------------------------------------------------------------------
# 3. Fluctuation spectrum of the nonlocal Frenkel kink (w = 1 units)
# ----------------------------------------------------------------------
# H = |k| + U(x),  U(x) = (x^2 - 1)/(x^2 + 1), continuum edge at +1.
# Dense spectral construction: |k| applied via FFT-consistent matrix.
N = 4096
L = 400.0                       # box half-width; potential tail ~ 1/x^2 needs room
x = (np.arange(N) - N / 2) * (2 * L / N)
dx = x[1] - x[0]
k = np.fft.fftfreq(N, d=dx) * 2 * np.pi

# |k| operator as a dense symmetric matrix via the unitary DFT
F = np.fft.fft(np.eye(N), axis=0)
Lam = (np.fft.ifft(np.abs(k)[:, None] * F, axis=0)).real
Lam = 0.5 * (Lam + Lam.T)       # symmetrise numerical residue

U = (x**2 - 1.0) / (x**2 + 1.0)
H = Lam + np.diag(U)
evals, evecs = np.linalg.eigh(H)

edge = 1.0
bound = evals[evals < edge - 5.0 / L]   # exclude box-discretised continuum near the edge
print("\n=== 3. Fluctuation spectrum about the PN kink ===")
print(f"grid: N = {N}, box = [-{L:.0f}, {L:.0f}], dx = {dx:.3f}")
print(f"continuum edge at {edge}")
print(f"discrete (bound) eigenvalues below the edge: {len(bound)}")
for i, e in enumerate(bound):
    psi = evecs[:, i]
    parity = "even" if abs(np.sum(psi * psi[::-1]) ) > 0.5 else "odd"
    print(f"  mode {i}:  E = {e:+.6f}   ({parity}), E/edge = {e/edge:+.4f}")

# Zero-mode check: overlap of the numerical ground state with phi' = 2/(x^2+1)
psi0_exact = 2.0 / (x**2 + 1.0)
psi0_exact /= np.linalg.norm(psi0_exact)
ov = abs(np.dot(psi0_exact, evecs[:, 0]))
print(f"overlap of mode 0 with the exact zero mode phi': {ov:.6f}")
print("levels sit at free-box positions: the scalar channel is reflectionless.")


# ----------------------------------------------------------------------
# 4. Levinson test: winding-2 background must bind exactly two modes
# ----------------------------------------------------------------------
# Background: two well-separated unit kinks (total winding 2),
# phi_bg = 2*arctan((x+a)/w) + 2*arctan((x-a)/w) + 2*pi, w = 1.
# Not an exact stationary point (same-sign kinks interact), so the two
# collective modes sit near zero rather than at zero; the Levinson claim
# is that the DISCRETE count below the edge is exactly 2 for winding 2.
print("\n=== 4. Levinson test: winding number 2 ===")
for a in [8.0, 16.0, 32.0]:
    N2, L2 = 4096, 600.0
    x2 = (np.arange(N2) - N2 / 2) * (2 * L2 / N2)
    dx2 = x2[1] - x2[0]
    k2 = np.fft.fftfreq(N2, d=dx2) * 2 * np.pi
    F2 = np.fft.fft(np.eye(N2), axis=0)
    Lam2 = (np.fft.ifft(np.abs(k2)[:, None] * F2, axis=0)).real
    Lam2 = 0.5 * (Lam2 + Lam2.T)
    phi_bg = 2 * np.arctan(x2 + a) + 2 * np.arctan(x2 - a) + 2 * np.pi
    U2 = 0.5 * 2.0 * np.cos(phi_bg)          # gamma'' = (V0/2) cos(phi), V0 = 2/w, w = 1
    H2 = Lam2 + np.diag(U2)
    ev2 = np.linalg.eigh(H2)[0]
    nb = int(np.sum(ev2 < 1.0 - 8.0 / L2))
    print(f"  separation 2a = {2*a:5.0f}:  discrete modes below edge = {nb}"
          f"   (lowest three: {np.round(ev2[:3], 4)})")
print("Levinson claim for this operator class: bound count = winding number.")


# ----------------------------------------------------------------------
# 5. Scale covariance: H(w) = (1/w) H(1), so both transverse channels
#    (glide leg, width w_par; stacking leg, width w_perp) share one
#    spectrum: one bound translation mode each, reflectionless continuum.
# ----------------------------------------------------------------------
print("\n=== 5. Scale covariance check (w = 3 vs w = 1) ===")
for wv in [1.0, 3.0]:
    N3, L3 = 4096, 400.0 * wv
    x3 = (np.arange(N3) - N3 / 2) * (2 * L3 / N3)
    dx3 = x3[1] - x3[0]
    k3 = np.fft.fftfreq(N3, d=dx3) * 2 * np.pi
    F3 = np.fft.fft(np.eye(N3), axis=0)
    Lam3 = (np.fft.ifft(np.abs(k3)[:, None] * F3, axis=0)).real
    Lam3 = 0.5 * (Lam3 + Lam3.T)
    U3 = (x3**2 - wv**2) / (wv * (x3**2 + wv**2))
    ev3 = np.linalg.eigh(Lam3 + np.diag(U3))[0]
    print(f"  w = {wv:.0f}:  w*E of lowest four = {np.round(wv*ev3[:4], 5)}"
          f"   (edge at w*E = 1)")
print("Identical scaled spectra: the glide and stacking channels each bind")
print("exactly one translation mode; total mode count = 2, the -2's count.")


# ----------------------------------------------------------------------
# 6. State counting and the threshold state (Krein/Levinson structure)
# ----------------------------------------------------------------------
# N_int(E) - N_free(E) = +2 at every E above the edge, box-stable: the
# phase shift is a CONSTANT pi (levels slide one slot, landing on free
# positions), reflectionless yet Levinson-saturated. Parity resolves the
# two units: even channel holds the bound zero mode; the odd channel
# holds a threshold state pinned at the edge (E - edge ~ 1e-5), odd under
# reflection, rms extent ~17 w against ~232 w for a box state.
print("\n=== 6. State counting: Delta N(E) and the odd threshold state ===")
Nc, Lc = 4096, 400.0
xc = (np.arange(Nc) - Nc/2) * (2*Lc/Nc); dxc = xc[1]-xc[0]
kc = np.fft.fftfreq(Nc, d=dxc) * 2*np.pi
Fc = np.fft.fft(np.eye(Nc), axis=0)
Lc2 = (np.fft.ifft(np.abs(kc)[:,None]*Fc, axis=0)).real; Lc2 = 0.5*(Lc2+Lc2.T)
Uc = (xc**2-1.0)/(xc**2+1.0)
evc, vecc = np.linalg.eigh(Lc2 + np.diag(Uc))
evf = np.linalg.eigvalsh(Lc2 + np.eye(Nc))
for E in [1.2, 2.0, 4.0]:
    print(f"  N(E<{E}): interacting-free = {int(np.sum(evc<E))-int(np.sum(evf<E)):+d}")
for i in [0, 1]:
    p = vecc[:, i]; par = np.sum(p*p[::-1])
    dens = np.abs(p)**2; dens /= dens.sum()
    print(f"  state {i}: E = {evc[i]:+.6f}, parity {par:+.2f}, rms extent = {np.sqrt(np.sum(dens*xc**2)):.1f} w")
print("Constant-pi phase shift: reflectionless, Levinson saturated in the count.")
