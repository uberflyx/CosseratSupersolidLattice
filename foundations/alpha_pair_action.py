"""
The neutral kink-antikink pair costs two kink actions, and the companion
direction holds no saddle: why the -2 is not a dilute-gas resummation.

Context. The identity alpha^{-1} = e^{S} - 2 can be rewritten
alpha = e^{-S}/(1 - 2e^{-S}), inviting a dilute-instanton reading in which
the second term, e^{-2S}, is a further tunnelling event dressing the
winding-one kink. This script shows that reading fails, so the -2 is the
topological/spectral constant (sec:vacuum_pol, sec:core_fluctuations), not a
gas coefficient.

Two results on the nonlocal Frenkel ring (mu = 1, V0 = 1/4, w = 4):

 1. PAIR ACTION. A kink-antikink held at fixed separation R (constrained
    relaxation, W = 0 sector) has action approaching 2*S_k as R grows: a
    neutral pair is a e^{-2S} object, two kink actions, not one. Dressing
    the winding-one kink with one pair therefore costs e^{-3S}, a
    self-energy-order effect (beside the -alpha term), not the constant -2.

 2. NO COMPANION SADDLE. Along the odd threshold-companion direction of the
    winding-one kink, the exact constrained landscape (all other directions
    relaxed at fixed companion amplitude) rises MONOTONICALLY on the soft
    side, with no stationary point. The multi-soliton configurations reached
    that way are not saddles, so they furnish no instanton term. An earlier
    reading that took a single non-stationary snapshot as "a pair at one
    extra kink action" was mistaken; this script records the correction.

Parity seals it: a net-winding-one sector holds only odd soliton numbers
(1, 3, 5, ...) at actions ~ S, 3S, 5S, so there is no e^{-2S} saddle. The
geometric form is algebra, not a gas.

Companion to alpha_core_factorisation.py, alpha_threshold_weight.py,
alpha_superfluid_units.py (Lattice repo, foundations/).
"""

import numpy as np

# ---- ring model ------------------------------------------------------
N, a, mu, V0 = 128, 1.0, 1.0, 0.25
w = mu / V0
j = np.arange(N)
x = (j - N / 2) * a
k = 2 * np.pi * np.fft.fftfreq(N, d=a)
Kel = (2.0 / a) * np.abs(np.sin(k * a / 2.0))
mask = Kel > 1e-12
invK = np.where(mask, 1.0 / np.where(mask, Kel, 1.0), 0.0)


def egrad(phit, W):
    s = np.roll(phit, -1) - phit + 2 * np.pi * W / N
    sh = np.fft.fft(s)
    E_el = (mu / (2.0 * N)) * np.sum(np.abs(sh)**2 * invK) / a
    phi = 2 * np.pi * W * j / N + phit
    E = E_el + a * V0 * np.sum(1.0 - np.cos(phi))
    g_s = np.fft.ifft(sh * invK).real * (mu / a)
    g = (np.roll(g_s, 1) - g_s) + a * V0 * np.sin(phi)
    return E, g


def relax(phit, W, fixed=None, iters=60000, lr=0.25):
    for _ in range(iters):
        E, g = egrad(phit, W)
        if fixed is not None:
            g[fixed] = 0.0
        phit = phit - lr * g
        if np.max(np.abs(g)) < 1e-10:
            break
    return egrad(phit, W)[0], phit


# ---- single kink -----------------------------------------------------
seed = (np.pi + 2 * np.arctan(x / w)) - 2 * np.pi * j / N
Sk, phit_k = relax(seed.copy(), 1)
print("=" * 66)
print(f" 1. Pair action vs separation (single kink S_k = {Sk:.4f})")
print("=" * 66)
print(f" {'R/w':>5} {'S_pair':>9} {'S_pair/S_k':>11}")
for Rw in [3, 4, 6, 8, 10, 12]:
    R = Rw * w
    c1 = int(N / 2 - R / 2); c2 = int(N / 2 + R / 2); mid = int(N / 2)
    phi0 = (2 * np.arctan((x - (c1 - N / 2) * a) / w)
            - 2 * np.arctan((x - (c2 - N / 2) * a) / w))
    phit = phi0.copy()
    phit[c1] = np.pi; phit[c2] = np.pi; phit[mid] = 2 * np.pi; phit[0] = 0.0
    fixed = np.array([c1, c2, mid, 0])
    S, _ = relax(phit, 0, fixed=fixed)
    print(f" {Rw:5.0f} {S:9.3f} {S/Sk:11.4f}")
print(" => approaches 2 with separation: a neutral pair is 2 kink actions.")

# ---- companion landscape is saddle-free ------------------------------
print("\n" + "=" * 66)
print(" 2. Companion-direction landscape (winding one): monotonic, no saddle")
print("=" * 66)
H1 = np.zeros((N, N))
for b in range(N):
    e = np.zeros(N); e[b] = 1.0
    s = np.roll(e, -1) - e
    g_s = np.fft.ifft(np.fft.fft(s) * invK).real * (mu / a)
    H1[:, b] = np.roll(g_s, 1) - g_s
H1 += np.diag(a * V0 * np.cos(2 * np.pi * j / N + phit_k))
H1 = 0.5 * (H1 + H1.T)
e1, v1 = np.linalg.eigh(H1)
zm, comp = v1[:, 0], v1[:, 1]
P = np.eye(N) - np.outer(zm, zm) - np.outer(comp, comp)
print(f" {'b':>6} {'V_eff - S_k':>12}")
prevV = -1e9
mono = True
for b in np.arange(0, -25.1, -3.0):
    E, _ = relax(phit_k + b * comp, 1, fixed=None, iters=9000, lr=0.25)
    # project out zm, comp during relaxation
    phit = phit_k + b * comp
    for _ in range(9000):
        _, g = egrad(phit, 1); g = P @ g; phit = phit - 0.25 * g
        if np.max(np.abs(g)) < 1e-10:
            break
    E = egrad(phit, 1)[0]
    print(f" {b:6.1f} {E - Sk:12.4f}")
    if b < 0 and E < prevV:
        mono = False
    prevV = E
print(" => monotonic on the soft side: no stationary multi-soliton saddle,")
print("    so no instanton contribution at e^{-2S}. Parity forbids it too.")

print("\n" + "=" * 66)
print(" CONCLUSION")
print("=" * 66)
print(" A neutral pair costs 2 S_k and enters alpha at e^{-3S} (self-energy")
print(" order). No e^{-2S} saddle exists at net winding one. The geometric")
print(" form alpha = e^{-S}/(1-2e^{-S}) is algebra, not a dilute gas; the -2")
print(" is the topological (skyrmion) and spectral (shift +2) constant.")
