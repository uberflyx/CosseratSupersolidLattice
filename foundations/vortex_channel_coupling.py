"""Two-channel fluctuation problem: the electron line's crystal face
coupled to its condensate face.

Question (vortex-face conjecture): the crystal-channel fluctuation
operator holds an odd threshold state pinned at the continuum edge.
The same line is a condensate vortex. Does coupling the disregistry
channel to the condensate-phase channel (i) leave the mode pinned,
(ii) bind it below the edge (an e* prediction), or (iii) dissolve it
into a resonance decaying to electron + second sound?

Model (w = 1 units, crystal edge = 1):
  H_c = |k| + (x^2-1)/(x^2+1)          crystal channel (established)
  H_s = rho |k| + V_s(x)               condensate channel; gapless
                                        Goldstone continuum from E = 0
  C   = g * c(x),  c(x) = 2/(x^2+1)    core-localised bilinear coupling
                                        (even, so it links equal parities)
  rho = v_2/c = 3.65                   second-sound stiffness (monograph)
  V_s = 0 (pure Goldstone)  or a healing-core well -g_v*2/(x^2+xi^2)

Diagnostic: strength function of the uncoupled odd threshold state
psi_th over the coupled spectrum,
  S(E_n) = |<psi_th | n>|^2 :
one dominant eigenstate at E ~ 1  -> pinned;
dominant weight below 1           -> bound (e*);
weight spread over many states    -> dissolved, width Gamma from the
                                     strength second moment.
"""

import numpy as np

def kmat(N, L):
    x = (np.arange(N) - N / 2) * (2 * L / N)
    dx = x[1] - x[0]
    k = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    F = np.fft.fft(np.eye(N), axis=0)
    Lam = (np.fft.ifft(np.abs(k)[:, None] * F, axis=0)).real
    return x, 0.5 * (Lam + Lam.T)

N, L = 2048, 300.0
x, Lam = kmat(N, L)
Uc = (x**2 - 1.0) / (x**2 + 1.0)
Hc = Lam + np.diag(Uc)

# Uncoupled crystal spectrum: identify the odd threshold state
evc, vecc = np.linalg.eigh(Hc)
i_th = 1                                    # bound zero mode is 0; threshold is 1
psi_th = vecc[:, i_th]
par = np.sum(psi_th * psi_th[::-1])
print(f"uncoupled threshold state: E = {evc[i_th]:+.6f}, parity {par:+.2f}")

rho = 3.65
cx = 2.0 / (x**2 + 1.0)

def coupled_run(g, gv=0.0, xi=1.0):
    Vs = -gv * 2.0 / (x**2 + xi**2)
    Hs = rho * Lam + np.diag(Vs)
    C = g * np.diag(cx)
    H = np.block([[Hc, C], [C, Hs]])
    ev, vec = np.linalg.eigh(H)
    # strength function of psi_th (embedded in the crystal block)
    ov = (vec[:N, :].T @ psi_th)**2
    Ebar = np.sum(ov * ev) / np.sum(ov)
    E2 = np.sum(ov * ev**2) / np.sum(ov)
    width = np.sqrt(max(E2 - Ebar**2, 0.0))
    imax = np.argmax(ov)
    return ev[imax], ov[imax], Ebar, width, np.sum(ov[ev < 0.999])

print(f"\n rho = {rho} (second sound), V_s = 0 (pure Goldstone)")
print(" g      E(peak)    peak wt   <E>       width     wt below edge")
for g in [0.0, 0.05, 0.1, 0.2, 0.4, 0.8]:
    Ep, wp, Eb, W, below = coupled_run(g)
    print(f" {g:4.2f}   {Ep:+.5f}   {wp:.3f}    {Eb:+.5f}  {W:.5f}   {below:.4f}")

print("\n with healing-core well in the condensate channel (g_v = 0.5, xi = 1)")
print(" g      E(peak)    peak wt   <E>       width     wt below edge")
for g in [0.1, 0.4, 0.8]:
    Ep, wp, Eb, W, below = coupled_run(g, gv=0.5)
    print(f" {g:4.2f}   {Ep:+.5f}   {wp:.3f}    {Eb:+.5f}  {W:.5f}   {below:.4f}")


# ----------------------------------------------------------------------
# Physical refinement: current coupling. Crystal strain couples to the
# condensate velocity dS/dx (mass current), not to the phase itself:
# C = (g/2) [ c(x) D + D^T c(x) ], D the antisymmetric derivative.
# The symmetrised form keeps H self-adjoint; one derivative flips parity,
# so the odd crystal mode now talks to EVEN condensate states.
# ----------------------------------------------------------------------
dx = x[1] - x[0]
D = np.zeros((N, N))
idx = np.arange(N - 1)
D[idx, idx + 1] = 1.0 / (2 * dx)
D[idx + 1, idx] = -1.0 / (2 * dx)

def coupled_run_current(g, rho_v):
    Hs = rho_v * Lam
    Cc = 0.5 * g * (np.diag(cx) @ D + D.T @ np.diag(cx))
    H = np.block([[Hc, Cc], [Cc.T, Hs]])
    ev, vec = np.linalg.eigh(H)
    ov = (vec[:N, :].T @ psi_th)**2
    imax = np.argmax(ov)
    return ev[imax], ov[imax], np.sum(ov[ev < 0.999])

print("\n current (velocity) coupling, parity-flipping:")
print(" rho    g      E(peak)    peak wt   wt below edge")
for rho_v in [1.0, 3.65]:
    for g in [0.1, 0.4, 0.8]:
        Ep, wp, below = coupled_run_current(g, rho_v)
        print(f" {rho_v:4.2f}  {g:4.2f}   {Ep:+.5f}   {wp:.3f}     {below:.4f}")


# ----------------------------------------------------------------------
# Results (2026-07-18):
# - No variant binds the mode below the edge: no sharp e* in any tested
#   regime (scalar or current coupling, with or without healing core,
#   rho = 1 or 3.65). The dangerous outcome is absent.
# - Generic static coupling turns the pinned state into a narrow
#   edge-straddling resonance leaking into the condensate channel:
#   ~94% intact at g = 0.05 (scalar), ~83-87% at g = 0.1, fragmenting
#   as g grows. Centroid fixed at the edge exactly (first-moment sum
#   rule); moment width linear in g (second-moment identity).
# - The two physical coupling characters push in OPPOSITE directions:
#   density (scalar) coupling drags the peak slightly below the edge,
#   current (velocity) coupling pushes it slightly above. The two-fluid
#   combination at the physical ratio could therefore restore exact
#   pinning; that ratio must be derived from the supersolid equations,
#   not tuned.
# - The topological (Magnus/Berry) protection is untestable in a static
#   symmetric Hamiltonian by construction; the first-order-in-time
#   symplectic problem is the required next machinery.
# ----------------------------------------------------------------------
