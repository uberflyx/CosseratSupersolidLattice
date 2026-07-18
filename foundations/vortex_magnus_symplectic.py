"""Magnus/Berry term acting on the core's fluctuation spectrum.

Part (a) result (analytic, from the supersolid chapter): the two-fluid
cross-coupling is derived, g = f_n f_s K_c = (16/75) mu, while the
condensate stiffness is K_sf ~ 1e40 mu. The static condensate response
to any crystal fluctuation is therefore suppressed as g/K_sf ~ 5e-42,
and the perturbation of the threshold mode by the condensate channel is
of order (g/mu)^2 (mu/K_sf) ~ 4e-42 of the edge energy: HIERARCHY
PROTECTION, parameter-free. Yesterday's fragility was an artifact of an
artificially soft condensate channel.

Part (b): the one condensate structure that survives the freeze is the
vortex circulation, which enters the crystal channel's dynamics as a
first-order (gyroscopic/Berry) term localised where the background phase
gradient lives:

    psi_ddot + gamma (B d/dx + d/dx B) psi_dot + H_c psi = 0,
    B(x) = xi/(x^2 + xi^2)   (vortex phase gradient on the glide line),

with H_c = |k| + (x^2-1)/(x^2+1) the established crystal operator and
gamma the Magnus coefficient (scale f_s = 4/5, scanned). The gyroscopic
matrix A = gamma (B D + D B)/2 is real antisymmetric, so the linearised
system z = (psi, psi_dot), z_dot = [[0, I], [-H_c, -A]] z is Hamiltonian
in structure; frequencies come in +/- pairs, and instability would show
as complex omega.

Question: does circulation (i) leave the odd threshold mode pinned at
the edge, (ii) push it below the edge into the empty window (a genuine
in-gap state, an e*), or (iii) destabilise it?
"""

import numpy as np

N, L = 1536, 250.0
x = (np.arange(N) - N / 2) * (2 * L / N)
dx = x[1] - x[0]
k = np.fft.fftfreq(N, d=dx) * 2 * np.pi
F = np.fft.fft(np.eye(N), axis=0)
Lam = (np.fft.ifft(np.abs(k)[:, None] * F, axis=0)).real
Lam = 0.5 * (Lam + Lam.T)
Hc = Lam + np.diag((x**2 - 1.0) / (x**2 + 1.0))

evc, vecc = np.linalg.eigh(Hc)
psi_th = vecc[:, 1]                     # the odd threshold state
psi_zm = vecc[:, 0]                     # the translational zero mode

D = np.zeros((N, N))
i = np.arange(N - 1)
D[i, i + 1] = 1.0 / (2 * dx)
D[i + 1, i] = -1.0 / (2 * dx)
B = np.diag(1.0 / (x**2 + 1.0))         # xi = 1

print("gyroscopic scan: omega^2-frequencies of the modes tracked by overlap")
print(" gamma   E_th(=w^2)  ov_th    E_zm      ov_zm    max|Im w|")
for gamma in [0.0, 0.1, 0.2, 0.4, 0.8, 1.6]:
    A = 0.5 * gamma * (B @ D + D @ B)   # real antisymmetric
    M = np.block([[np.zeros((N, N)), np.eye(N)], [-Hc, -A]])
    ev, vec = np.linalg.eig(M)          # eigenvalues lambda = i*omega
    om = ev / 1j                        # omega
    # keep the omega >= 0 half; energies E = omega^2 for comparison with Hc
    sel = np.where(om.real >= -1e-9)[0]
    psis = vec[:N, sel]
    # normalise position part for overlap bookkeeping
    nrm = np.linalg.norm(psis, axis=0)
    nrm[nrm == 0] = 1.0
    psis = psis / nrm
    ov_th = np.abs(psis.conj().T @ psi_th)**2
    ov_zm = np.abs(psis.conj().T @ psi_zm)**2
    jt, jz = np.argmax(ov_th), np.argmax(ov_zm)
    Eth = (om[sel][jt]**2).real
    Ezm = (om[sel][jz]**2).real
    print(f" {gamma:4.2f}   {Eth:+.6f}  {ov_th[jt]:.3f}   {Ezm:+.6f}  {ov_zm[jz]:.3f}   {np.max(np.abs(om.imag)):.2e}")


# ----------------------------------------------------------------------
# The correct Magnus structure: antisymmetric pairing of the TWO
# transverse channels (glide leg a, stacking leg b), the vortex acting
# as an effective magnetic field on the doublet:
#   psi_a_dd + gamma B(x) psi_b_d + H psi_a = 0
#   psi_b_dd - gamma B(x) psi_a_d + H psi_b = 0
# Circular combinations psi_pm = psi_a -+ i psi_b decouple into
#   (-w^2 -+ gamma w B(x) + H) psi_pm = 0:
# an omega-dependent core well (one sense) and hill (the other).
# Solve each branch as a nonlinear-in-omega problem by iteration.
# ----------------------------------------------------------------------
Bx = 1.0 / (x**2 + 1.0)

def circular_branch(gamma, sign, E0=1.0, iters=40):
    """Solve (H + sign*gamma*w*B) psi = w^2 psi self-consistently near E0."""
    E = E0
    for _ in range(iters):
        w = np.sqrt(max(E, 1e-12))
        Heff = Hc + sign * gamma * w * np.diag(Bx)
        ev, vec = np.linalg.eigh(Heff)
        ov = (vec.T @ psi_th)**2
        j = np.argmax(ov)
        Enew = ev[j]
        if abs(Enew - E) < 1e-10:
            E = Enew; break
        E = 0.5 * (E + Enew)
    return E, ov[j]

print("\ncorrect Magnus pairing: circular branches of the threshold doublet")
print(" gamma   E_minus(bound?)  ov     E_plus(anti)   ov")
for gamma in [0.0, 0.1, 0.2, 0.4, 0.8, 1.6]:
    Em, om_ = circular_branch(gamma, -1)
    Ep, op_ = circular_branch(gamma, +1)
    print(f" {gamma:4.2f}   {Em:+.6f}      {om_:.3f}  {Ep:+.6f}     {op_:.3f}")
print("\nZero modes: at omega = 0 the Magnus term vanishes identically, so both")
print("translations stay exact zero modes; their conjugate pairing is dynamical.")


# ----------------------------------------------------------------------
# Results (2026-07-18, standalone run at N=768, L=150):
#
# Correct Magnus pairing (circular branches of the threshold doublet):
#  gamma   E_minus(bound)   ov      E_plus(anti)   ov
#  0.05    +0.98655        0.988    +1.01157      0.937
#  0.10    +0.97185        0.968    +1.02637      0.477
#  0.20    +0.94021        0.929    +1.03621      0.374
#  0.40    +0.87230        0.859    +1.07637      0.183
#  0.80    +0.73509        0.760    +1.10102      0.076
#  1.60    +0.50569        0.651    +1.02079      0.043
#
# Verdict: the circulation acts on the core's internal doublet exactly
# as an effective magnetic field. One circular sense is BOUND below the
# phonon edge, sharp and coherent (in-gap, high overlap); the opposite
# sense is anti-bound and dissolves into the continuum. The handedness
# of the bound mode is set by the circulation sense, i.e. by the sign of
# the charge. Zero modes exact at all gamma (Magnus term vanishes at
# omega = 0).
#
# Status: the physical gamma and the edge energy scale are NOT yet
# derived; observability and identification of the bound mode hang on
# both. Do not interpret before deriving them and checking the
# monograph's existing generation/internal-excitation assignments.
# ----------------------------------------------------------------------
