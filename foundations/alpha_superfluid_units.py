"""
Mitch's conjecture: the -2 is 'two units into the superfluid'.

Sharp version. Per transverse channel the core operator carries TWO spectral
passengers: the bound translational zero mode (crystal) and the odd threshold
state at the edge (conjectured condensate face: it is the breathing mode, the
Magnus term binds it, and its superpotential is the vortex velocity). The
kinetic kernel |k| of the core operator is itself condensate-supplied: the
bulk-slaved inertia M_eff(k) = rho/(4|k|) is what linearises the PN energy to
|k| rather than -d^2/dx^2 (sec:core_kinetic_kernel).

DISCRIMINATING TEST. If the second unit is the condensate's, it must vanish
when the nonlocal kernel is replaced by a local crystal kernel with the same
topological content. Compare:

  (a) NONLOCAL (condensate-slaved):  H = |k|      + 1 - 2/(x^2+1)
      exact kink of the nonlocal Frenkel model, zero mode 1/(x^2+1),
      algebraic (1/x^2) tails forced by |k|.

  (b) LOCAL (bare crystal, sine-Gordon): H = -d^2/dx^2 + 1 - 2 sech^2(x)
      exact kink of the local Frenkel (sine-Gordon) model, zero mode sech(x),
      exponential tails. Poschl-Teller lambda = 1: reflectionless too.

Both are reflectionless, both have exactly one bound zero mode. The question
is the SPECTRAL SHIFT: N_int - N_free. Local Poschl-Teller is the textbook
case: shift = 1 (bound state only, no threshold companion). If the nonlocal
case gives 2, the extra unit per channel exists because and only because the
kinetic kernel is the condensate's.

SECOND TEST (winding scaling). If each unit of topological winding brings one
crystal + one condensate passenger, a winding-2 background should carry
spectral shift 4 in the nonlocal operator and 2 in the local one.

THIRD (analytic, printed for the record). The superpotential's asymptotic
winding: W = -(ln psi0)' -> p/x where p is the tail power of the zero mode.
Nonlocal: psi0 ~ 1/x^2 (BO algebraic soliton; the 1/x^2 decay is forced by
|k|), so p = 2: the superpotential is asymptotically a DOUBLE-quantum
circulation field, crystal winding + condensate winding. Local: psi0 ~ e^-x,
W -> 1 (constant, no algebraic winding at all).

Companion to alpha_core_factorisation.py and alpha_minus2_spectral.py.
"""

import numpy as np


def spectral_shift(kind, U_int, U_free, L=600.0, N=6000, window=(1.05, 4.0)):
    """N_int(E) - N_free(E) averaged over the window, plus bound/edge census."""
    x = (np.arange(N) - N // 2) * (L / N)
    k = 2.0 * np.pi * np.fft.fftfreq(N, d=L / N)
    F = np.fft.fft(np.eye(N), axis=0)
    if kind == "nonlocal":
        T = (np.fft.ifft(np.abs(k)[:, None] * F, axis=0)).real
    else:
        T = (np.fft.ifft((k**2)[:, None] * F, axis=0)).real
    T = 0.5 * (T + T.T)
    eH = np.linalg.eigvalsh(T + np.diag(U_int(x)))
    eH0 = np.linalg.eigvalsh(T + np.diag(U_free(x)))
    Es = np.linspace(*window, 50)
    dN = np.array([np.sum(eH < E) - np.sum(eH0 < E) for E in Es])
    nb = int(np.sum(eH < 1.0 - 5.0 / L))
    n_edge = int(np.sum(np.abs(eH - 1.0) < 3.0 / L))
    return dN.mean(), dN.std(), nb, n_edge, eH


print("=" * 72)
print(" TEST 1: local crystal kernel vs condensate-slaved nonlocal kernel")
print("=" * 72)

# (a) nonlocal Frenkel kink
mu, sd, nb, ne, _ = spectral_shift(
    "nonlocal",
    lambda x: (x**2 - 1.0) / (x**2 + 1.0),
    lambda x: np.ones_like(x),
)
print(f" NONLOCAL |k| + 1 - 2/(x^2+1):")
print(f"   spectral shift = {mu:.3f} +/- {sd:.3f}   bound below edge: {nb},"
      f" states pinned at edge: {ne}")

# (b) local sine-Gordon kink (Poschl-Teller lambda=1), same edge at 1
mu2, sd2, nb2, ne2, _ = spectral_shift(
    "local",
    lambda x: 1.0 - 2.0 / np.cosh(np.clip(x, -300, 300))**2,
    lambda x: np.ones_like(x),
)
print(f" LOCAL -d2/dx2 + 1 - 2 sech^2(x):")
print(f"   spectral shift = {mu2:.3f} +/- {sd2:.3f}   bound below edge: {nb2},"
      f" states pinned at edge: {ne2}")

print("""
 Verdict: if the nonlocal shift is 2 and the local shift is 1, the second
 unit per channel is supplied by the |k| kernel, i.e. by the condensate's
 bulk-slaved inertia. The crystal alone books one unit (the zero mode).
""")

print("=" * 72)
print(" TEST 2: winding scaling (two well-separated kinks)")
print("=" * 72)
a = 40.0
mu4, sd4, nb4, ne4, _ = spectral_shift(
    "nonlocal",
    lambda x: 1.0 - 2.0 / ((x - a)**2 + 1.0) - 2.0 / ((x + a)**2 + 1.0),
    lambda x: np.ones_like(x),
    L=900.0, N=6000,
)
print(f" NONLOCAL, winding 2 (separation {2*a:.0f}):")
print(f"   spectral shift = {mu4:.3f} +/- {sd4:.3f}   bound: {nb4},"
      f" edge-pinned: {ne4}")

mu5, sd5, nb5, ne5, _ = spectral_shift(
    "local",
    lambda x: 1.0 - 2.0 / np.cosh(np.clip(x - a, -300, 300))**2
                  - 2.0 / np.cosh(np.clip(x + a, -300, 300))**2,
    lambda x: np.ones_like(x),
    L=900.0, N=6000,
)
print(f" LOCAL, winding 2:")
print(f"   spectral shift = {mu5:.3f} +/- {sd5:.3f}   bound: {nb5},"
      f" edge-pinned: {ne5}")

print("""
 Verdict: nonlocal shift = 2 per winding (crystal + condensate passenger
 each), local shift = 1 per winding (crystal only), if the reading is right.
""")

print("=" * 72)
print(" TEST 3 (analytic, for the record): superpotential winding")
print("=" * 72)
print(" nonlocal zero mode 1/(x^2+1):  W = 2x/(x^2+1) -> 2/x   (winding 2)")
print(" local zero mode sech(x):       W = tanh(x)    -> 1     (no winding)")
print(" The 1/x^2 tail of the BO soliton is forced by |k|; the tail power 2")
print(" is the double circulation: one crystal quantum + one condensate")
print(" quantum, the same doubling as W = 2 v_theta.")
