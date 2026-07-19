"""
[SUPERSEDED 2026-07-19 by alpha_pair_action.py.] The pair-channel reading
below OVERREACHED: it took a non-stationary snapshot of the companion
landscape as a saddle. A clean check shows a neutral kink-antikink pair costs
2 S_k (not one), so it enters alpha at e^{-3S} (self-energy order), not the
constant -2, and the companion landscape has no stationary point. The -2 is
the topological/spectral constant, not a dilute-gas resummation. Retained for
the record; do not cite as a mechanism for the -2.
"""

"""
The threshold companion opens the neutral pair channel: ring demonstration.

Context. The alpha series alpha^{-1} = e^{S} - 2 - ... reads, in resummed
form, alpha = e^{-S}/(1 - 2 e^{-S}): the bare compact-winding (registry-hop)
amplitude dressed by a geometric series whose coefficient counts the core's
two marginal threshold companions at unit weight. Parity forbids charged
instanton pairs at net winding one, so the series' second term needs a
NEUTRAL object of one kink action. This script identifies that object and
the coordinate that reaches it, on a discretised nonlocal Frenkel ring.

Model (single transverse channel, so the expected channel count is ONE):
  phi_j = 2*pi*W*j/N + phit_j on an N-site ring,
  E = (mu/2N) sum_{k != 0} |s_k|^2 / K(k) / a + a*V0*sum_j (1 - cos phi_j),
  s_j the periodic strain, K(k) the lattice |k| (the medium kernel).
W = 1 minimum: the relaxed compact kink, action S1. Its Hessian carries the
translation zero mode and the odd threshold companion (the breathing mode),
as in the continuum operator of the alpha chapter.

Results established here:
 1. LANDSCAPE. The exact constrained potential V_eff(b) along the companion
    direction (all other directions relaxed at fixed companion amplitude b)
    is one-sided: b > 0 stiffens above the Gaussian; b < 0 softens into a
    channel. One channel per companion.
 2. IDENTIFICATION. The configurations along the channel show the kink
    splitting into two kinks and one antikink, net winding one preserved:
    the companion is the nucleation coordinate of a neutral kink-antikink
    pair. Neutral insertions evade the winding-parity block.
 3. COST. The pair region opens at Delta E = 8.84 against S1 = 9.02 on this
    ring (two per cent): one extra kink action per insertion, the e^{-S}
    the resummation requires.
 4. HONEST LIMIT. The channel's excess measure over the Gaussian, computed
    by quadrature on V_eff, decays with rate ~3.3, NOT the kink rate 9.0:
    the one-dimensional reduction mixes the smooth anharmonic softening
    near b = 0 with the pair region, so the insertion WEIGHT cannot be read
    off this landscape. The weight requires the pair's own moduli
    (separation and phase) with the Dashen-Hasslacher-Neveu measure; that
    derivation remains open and is the residual of item A1.

Companion to alpha_core_factorisation.py, alpha_threshold_weight.py,
alpha_superfluid_units.py (Lattice repo, foundations/).
"""

import numpy as np

# ----------------------------------------------------------------------
# model
# ----------------------------------------------------------------------
N, a, mu, V0 = 64, 1.0, 1.0, 0.25          # w = mu/V0 = 4, L/w = 16
j = np.arange(N)
k = 2*np.pi*np.fft.fftfreq(N, d=a)
Kel = (2.0/a)*np.abs(np.sin(k*a/2.0))
mask = Kel > 1e-12
invK = np.where(mask, 1.0/np.where(mask, Kel, 1.0), 0.0)


def energy_grad(phit, W):
    """Energy and gradient in the periodic part phit; full field carries the
    winding ramp 2*pi*W*j/N."""
    s = np.roll(phit, -1) - phit + 2*np.pi*W/N
    sh = np.fft.fft(s)
    E_el = (mu/(2.0*N))*np.sum((np.abs(sh)**2)*invK)/a
    phi = 2*np.pi*W*j/N + phit
    E = E_el + a*V0*np.sum(1.0 - np.cos(phi))
    g_s = np.fft.ifft(sh*invK).real*(mu/a)
    g = (np.roll(g_s, 1) - g_s) + a*V0*np.sin(phi)
    return E, g


def relax(phit, W, iters=20000, lr=0.35, proj=None):
    for _ in range(iters):
        E, g = energy_grad(phit, W)
        if proj is not None:
            g = proj @ g
        phit = phit - lr*g
        if np.max(np.abs(g)) < 1e-11:
            break
    return energy_grad(phit, W)[0], phit


# ----------------------------------------------------------------------
# 1. relaxed compact kink and its Hessian structure
# ----------------------------------------------------------------------
x = (j - N/2)*a
w = mu/V0
seed = (np.pi + 2*np.arctan(x/w)) - 2*np.pi*j/N
S1, phit_k = relax(seed, 1)


def hessian(phit, W):
    phi = 2*np.pi*W*j/N + phit
    H = np.zeros((N, N))
    for b in range(N):
        e = np.zeros(N); e[b] = 1.0
        s = np.roll(e, -1) - e
        g_s = np.fft.ifft(np.fft.fft(s)*invK).real*(mu/a)
        H[:, b] = np.roll(g_s, 1) - g_s
    H += np.diag(a*V0*np.cos(phi))
    return 0.5*(H + H.T)


H1 = hessian(phit_k, 1)
e1, v1 = np.linalg.eigh(H1)
zm, comp = v1[:, 0], v1[:, 1]
breath = -2*x/(x**2 + w**2)
breath /= np.linalg.norm(breath)
print(f"S1 = {S1:.4f};  zero mode E = {e1[0]:+.1e};  companion E = {e1[1]:.4f}"
      f" (edge V0 = {V0}), |<companion|breathing>| = {abs(comp @ breath):.4f}")

# ----------------------------------------------------------------------
# 2-3. the exact constrained landscape along the companion, both signs
# ----------------------------------------------------------------------
P = np.eye(N) - np.outer(zm, zm) - np.outer(comp, comp)
print("\n  b     V_eff - S1     quad = lam*b^2/2")
table = []
prev = phit_k
for b in list(np.arange(0, 12.1, 1.5)) + list(np.arange(-1.5, -30.1, -1.5)):
    E, cfg = relax(prev + b*comp - (prev - phit_k)*0, 1, iters=9000, lr=0.25,
                   proj=P) if False else relax(phit_k + b*comp, 1, iters=9000,
                                               lr=0.25, proj=P)
    table.append((b, E - S1))
    print(f"{b:6.1f}  {E-S1:12.4f}  {0.5*e1[1]*b**2:12.4f}")
    if abs(b + 12.0) < 1e-9:
        phi = 2*np.pi*j/N + cfg
        s = np.roll(phi, -1) - phi; s[-1] += 2*np.pi
        print("   configuration at b = -12 (strain / max):")
        print("  ", np.round(s/np.abs(s).max(), 2)[::2])
        print(f"   net winding {s.sum()/(2*np.pi):.3f}; antikink sites "
              f"(strain < -0.05): {np.sum(s < -0.05)}; min strain {s.min():.2f}")
        print("   => two kinks and one antikink: the neutral pair, nucleated.")

# ----------------------------------------------------------------------
# 4. excess measure and its rate (the honest limit)
# ----------------------------------------------------------------------
from scipy.interpolate import interp1d
tb = np.array(sorted(table))
Vf = interp1d(tb[:, 0], tb[:, 1], kind='cubic')
bb = np.linspace(tb[0, 0], tb[-1, 0], 4001)
dV = Vf(bb); db = bb[1] - bb[0]
betas = np.linspace(0.35, 0.65, 13)
lndW = [np.log(np.sum(np.exp(-be*dV))*db - np.sqrt(2*np.pi/(be*e1[1])))
        for be in betas]
rate = -np.polyfit(betas, lndW, 1)[0]
print(f"\nexcess-measure rate over beta in [0.35, 0.65]: {rate:.2f}"
      f"  vs kink action S1 = {S1:.2f}")
print("=> the 1D reduction mixes smooth softening with the pair region;")
print("   the insertion weight needs the pair moduli (DHN measure). Open.")
