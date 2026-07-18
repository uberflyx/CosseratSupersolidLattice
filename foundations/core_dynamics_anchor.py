"""The dynamical anchor of the core spectrum: rho_eff, the KK gap, and the
flavour-gated binding ratio.

Closes the structural half of A3(ii) and sharpens A4. Three results.

(1) KINETIC KERNEL (derived, quasistatic bulk slaving). The glide field
    Delta(x,t) has no mass of its own; a mode of wavevector k drags bulk
    material to depth 1/(2|k|), giving the nonlocal inertia
        M_eff(k) = rho / (4|k|).
    With the elastic kernel (mu/2)|k| and the PN valley curvature
    C = mu/(2w), the k4 = 0 dispersion is
        omega^2(k) = 2 c^2 k^2 + 2 c^2 |k| / w,
    GAPLESS at small k. The retarded check closes the loophole: a guided
    interface mode below the bulk cone would need
    (mu/2) sqrt(k^2 - omega^2/c^2) = -C < 0, impossible. So the static
    Hessian's "edge" at 1/w is a stiffness eigenvalue, not a frequency;
    the frequency continuum at k4 = 0 is the gapless bulk shear cone.

(2) NO SHARP e* (numeric). Solving the generalized problem
    H psi = omega^2 M psi at k4 = 0 (box-regularised IR), the odd
    threshold state fragments across the gapless spectrum, and the
    fragmentation deepens with box size (max single-mode overlap 0.69 at
    L = 200 falling to 0.03 at L = 250; a true bound state would hold at
    1.0 independent of box): a leaky resonance, not a particle. The electron's own
    flavour sector holds no bound excitations. Conversely the k4 = +/-1
    sector is gapped at the lattice acoustic rung
        E_1 = 3 m0 c^2 / sqrt(2) = 148.546 MeV
    (the monograph's KK rung: chain dispersion of the three-site compact
    ring at k4 d111 = 2pi/3, at the shear speed), so a core mode bound
    below E_1 in that sector is exactly stable: it cannot radiate spatial
    (k4 = 0) shear waves without violating compact quasi-momentum.
    Sharp internal excitations exist if and only if they carry flavour.
    Parameter-free inequality: the second-generation excitation must lie
    below 148.5 MeV. The muon does: m_mu / E_1 = 0.71128.

(3) BINDING RATIO (numeric, branch-continued, convention-gated). The
    Magnus-bound circular mode with gamma = 2 f_s = 8/5 (vortex core mass
    at the cell scale) and circulation healing width xi = w (the rolling
    constraint locks the condensate winding to the crystal core) gives
    sqrt(E) = 0.71112 in edge units: 0.02% from m_mu / E_1. HONESTY
    GATE: the ratio swings 0.59 (xi = ell/w_par = 2.21) to 0.85
    (xi = 0.5) across the plausible healing range, so the agreement is
    meaningful only if xi = w and the cell-scale core mass are derived,
    which is the D_4-reading construction task. Until then this is a
    target and a suggestive hit, not a derivation.
"""

import numpy as np
import sympy as sp
from scipy.linalg import sqrtm

# ----------------------------------------------------------------------
# (1) symbolic: kinetic kernel and the gapless k4 = 0 dispersion
# ----------------------------------------------------------------------
k, z, rho, Dd = sp.symbols('k z rho Delta_dot', positive=True)
T = 2 * sp.Rational(1, 2) * rho * (Dd/2)**2 * sp.integrate(sp.exp(-2*k*z),
                                                           (z, 0, sp.oo))
M_eff = sp.simplify(2*T/Dd**2)
mu, w, c, ell = sp.symbols('mu w c ell', positive=True)
om2 = sp.expand(((mu/2)*k + mu/(2*w)) / M_eff).subs(mu, rho*c**2)
print("(1) M_eff(k) =", M_eff, ";  omega^2(k) =", om2, " [gapless]")

d111 = ell*sp.sqrt(sp.Rational(2, 3))
E1_units = sp.simplify(2*(c/d111)*sp.sin(sp.pi/3) * ell/c)
m0 = 0.51099895069 * 137.035999177          # MeV (CODATA 2022)
E1 = float(E1_units)*m0
mmu = 105.6583755                            # MeV (PDG, project file)
print(f"    KK rung E_1 = {E1_units} * m0 c^2 = {E1:.4f} MeV ; "
      f"m_mu/E_1 = {mmu/E1:.5f}")

# ----------------------------------------------------------------------
# shared machinery: the established core operator
# ----------------------------------------------------------------------
N, L = 1536, 250.0
x = (np.arange(N)-N/2)*(2*L/N); dx = x[1]-x[0]
kg = np.fft.fftfreq(N, d=dx)*2*np.pi
F = np.fft.fft(np.eye(N), axis=0)
Lam = (np.fft.ifft(np.abs(kg)[:, None]*F, axis=0)).real
Lam = 0.5*(Lam+Lam.T)
Hc = Lam + np.diag((x**2-1.0)/(x**2+1.0))
evc, vecc = np.linalg.eigh(Hc)
psi0 = vecc[:, 1]                            # odd threshold (width) mode

# ----------------------------------------------------------------------
# (2) k4 = 0 generalized problem: fragmentation of the threshold state
# ----------------------------------------------------------------------
kreg = np.maximum(np.abs(kg), 2*np.pi/(2*L))
Minv = (np.fft.ifft((4*kreg)[:, None]*F, axis=0)).real
Minv = 0.5*(Minv+Minv.T)
Mih = np.real(sqrtm(Minv))
A = Mih @ Hc @ Mih; A = 0.5*(A+A.T)
om2n, phi = np.linalg.eigh(A)
psid = Mih @ phi; psid /= np.linalg.norm(psid, axis=0)
ovf = (psid.T @ psi0)**2
print(f"(2) k4=0 fragmentation: max single-mode overlap with the threshold "
      f"state = {ovf.max():.3f} (resonance, not a bound state)")

# ----------------------------------------------------------------------
# (3) branch-continued binding at gamma = 8/5, sensitivity in xi
# ----------------------------------------------------------------------
M = 360
V = vecc[:, :M]; D = np.diag(evc[:M]); t0 = V.T @ psi0

def follow(xi, gmax=1.6, steps=24):
    Bp = V.T @ np.diag(1.0/((x/xi)**2+1.0)) @ V
    target = t0.copy(); E = 1.0
    for g in np.linspace(gmax/steps, gmax, steps):
        for _ in range(50):
            om = np.sqrt(max(E, 1e-12))
            ev, U = np.linalg.eigh(D - g*om*Bp)
            j = np.argmax((U.T @ target)**2)
            En = ev[j]
            if abs(En-E) < 1e-11:
                break
            E = 0.5*(E+En)
        target = U[:, j]
    return E

print("(3) binding ratio sqrt(E) at gamma = 8/5 vs circulation width xi:")
for xi in [0.5, 1.0, 2.21]:
    E = follow(xi)
    print(f"    xi = {xi:4.2f}: sqrt(E) = {np.sqrt(max(E, 0)):.5f}"
          + ("   <- locked-profile convention; full-space ref 0.71112,"
             " target 0.71128" if xi == 1.0 else ""))
