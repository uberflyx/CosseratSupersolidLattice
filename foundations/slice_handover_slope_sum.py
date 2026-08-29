#!/usr/bin/env python3
"""The slope sum S on the FCC slice under the minimal handover law.

THE ONE CALCULATION (named in open_items 2026-08-29, target on file first):
evaluate the misfit line gamma(f) of the {111} partial hop with
  - central channel: the pinned truncated Morse, a l = 7/3, D = 9/98,
    V''(l) = k_n = 1 (bond-level certified against the monograph's printed
    fault microscopy: dV(0.816) = 0.0262, dV(1.291) = 0.0223);
  - tangential channel: the lossless handover law of d4/d4_handover.py
    ported to the slice: rate-based accumulation of frustrated slip,
    stiffness k_t(r) = r0 * fade(r) * k_n with r0 = 1/(2 pi - 3) (the FCC
    rolling ratio), contact reset at opening, layer rolling angles relaxed
    explicitly (the eta_0 relief emerges, it is not input);
  - gap relaxed at every f (GSF protocol, in-plane frozen).
Then project the even part onto [1 - cos(2 pi m f/d)] and form the
absolute equilibrium factor
  P = 8 pi^2 sum_m m V_m / (1 + r0)   [V_m per wall cell, k_n units],
which the alpha inversion demands equal S_alpha = 1.001067
(foundations/misfit_slope_sum.py, committed before this build).
The same line prices the saddle (gamma_USF, per-area k_n units) for the
string-tension ledger; its closed-form confrontation needs that section's
own normalisation and is only PRINTED here, not judged.

STATUS AFTER THE FIRST RUN (2026-08-29): NOT YET A DEFENSIBLE P.
Committed as working machinery with its failure modes documented, per the
no-fudging rule: iterating knobs against the on-file target is the
forbidden move, so v1 stops where the artefacts start. What v1 established:
 1. The truncation is inferable and matters: the slice stacking fault is
    free only below the third shell, where the fault first changes
    distances (sqrt3 -> 2h = 1.633, textbook), so the truncated-shifted
    Morse needs sqrt2 < r_c < 1.633; r_c = 1.6 (the window already used in
    d4_handover's exponential fade) restores slice degeneracy by
    construction. The definitive r_c must be pinned on the D4 microscopy's
    printed remainder (-0.0051 k_n) and gamma_rigid = 0.1365 k_n, a D4
    build this file does not yet contain.
 2. Central channel confirmed barrier-poor on the slice: gamma_c,max ~
    0.006 k_n/col with the interior maximum at f/d ~ 0.25, an order below
    the tangential scale, consistent with the monograph's saddle
    statement. But the relaxed endpoint shows gamma_c(d) = -0.012: the
    gap relaxation re-engages the 2h bonds through dh < 0, reintroducing
    fault sensitivity through the back door. v2 must either pin r_c from
    D4 first or constrain dh to the physical branch.
 3. The tangential harmonics are not yet trustworthy: V_2/V_1 spans
    0.9-2.4 across fade laws where d4_handover's gamma_max was
    fade-independent to 1.00. Harmonic extraction is far more delicate
    than barrier heights; the rate-based XI accumulation needs a smoother
    integrator (finer f-grid, gradient-based minimiser, convergence in
    layers and box) before V_m mean anything.
 4. Bond-level certification stands (0.0262/0.0223 exact), and the
    absolute calibration bridge P = 8 pi^2 sum m V_m/(1+r0) is in place
    for v2.
The on-file target S_alpha = 1.001067 was not consulted in any modelling
choice above; the only comparisons run were printed after computation.
"""
import numpy as np
from scipy.optimize import minimize

# ---------- pinned contact law: truncated-shifted Morse, r_c = 1.6
# (between the second shell sqrt2 = 1.414 and the first fault-affected
# distance 2h = 1.633: the slice stacking fault is then free exactly, the
# monograph's 'free by nearest-neighbour degeneracy', and the same cutoff
# appears in d4/d4_handover.py's exponential fade window)
AM = 7/3; Dm = 1/(2*AM*AM); R = 0.5; RC = 1.6
r0 = 1/(2*np.pi - 3)
def _Vm(r):    return Dm*(1 - np.exp(-AM*(r-1)))**2
def Vmorse(r): return np.where(r < RC, _Vm(r) - _Vm(RC), 0.0)
def kn_eff(r):
    e = np.exp(-AM*(r-1)); return 2*Dm*AM*AM*(2*e*e - e)
def fade_morse(r):  return np.clip(kn_eff(r)/kn_eff(1.0), 0, None)
def fade_lin(r):    return np.clip(1-(r-1)/0.30, 0, 1)
def fade_exp(r):    return np.clip(np.exp(-2*AM*(r-1)), 0, 1)*(r < 1.6)

# ---------- {111} slice geometry
a1 = np.array([1,0,0]); a2 = np.array([0.5, np.sqrt(3)/2, 0])
h  = np.sqrt(2/3); sv = (a1+a2)/3; d = np.linalg.norm(sv)
t_hat = sv/d; w_hat = np.array([0,0,1.0])
Gt = np.outer(t_hat, w_hat) - np.outer(w_hat, t_hat)

def build(box=5, layers=3, rcpair=2.3):
    lows, ups = [], {}
    for m in range(-layers+1, 1):
        for i in range(-box, box+1):
            for j in range(-box, box+1):
                lows.append((m, i*a1 + j*a2 + m*sv + np.array([0,0,m*h])))
    for m in range(1, layers+1):
        ups[m] = m*sv + np.array([0,0,m*h])
    MU, ML, V0 = [], [], []
    for mu, pu in ups.items():
        for ml, pl in lows:
            v0 = pu - pl
            if np.linalg.norm(v0) < rcpair:
                MU.append(mu); ML.append(ml); V0.append(v0)
    return np.array(MU), np.array(ML), np.array(V0)

MU, ML, V0 = build()
ROT = [-2,-1,0,1,2]; li = {m:k for k,m in enumerate(ROT)}
IU = np.array([li.get(m,-1) for m in MU]); IL = np.array([li.get(m,-1) for m in ML])

def run(fade, r_ratio=r0, nf=61, layers_note=""):
    fs = np.linspace(0, d, nf)
    XI = np.zeros((len(V0),3)); th_p = np.zeros(len(ROT)); dh_p = 0.0; f_p = 0.0
    g = np.zeros(nf); E0 = None; x = np.zeros(1+len(ROT))
    for k,f in enumerate(fs):
        XIc = XI.copy()
        def energy(x):
            dh = x[0]; th = x[1:]
            v = V0 + f*t_hat + dh*w_hat
            r = np.linalg.norm(v, axis=1); n = v/r[:,None]
            E = Vmorse(r).sum()
            if r_ratio > 0:
                kt = r_ratio*fade(r)*kn_eff(1.0)
                dU  = (f-f_p)*t_hat + (dh-dh_p)*w_hat
                dUp = dU - (n@dU)[:,None]*n
                tu = np.where(IU>=0, th[np.maximum(IU,0)]-th_p[np.maximum(IU,0)], 0.0)
                tl = np.where(IL>=0, th[np.maximum(IL,0)]-th_p[np.maximum(IL,0)], 0.0)
                xin = XIc + dUp - R*(tu+tl)[:,None]*(n@Gt.T)
                E += 0.5*np.sum(kt*np.einsum('ij,ij->i', xin, xin))
            return E
        res = minimize(energy, x, method='Powell',
                       options={'xtol':1e-10,'ftol':1e-12,'maxiter':4000})
        x = res.x
        if E0 is None: E0 = res.fun
        g[k] = res.fun - E0
        dh = x[0]; th = x[1:]
        v = V0 + f*t_hat + dh*w_hat
        r = np.linalg.norm(v, axis=1); n = v/r[:,None]
        dU = (f-f_p)*t_hat + (dh-dh_p)*w_hat; dUp = dU - (n@dU)[:,None]*n
        tu = np.where(IU>=0, th[np.maximum(IU,0)]-th_p[np.maximum(IU,0)], 0.0)
        tl = np.where(IL>=0, th[np.maximum(IL,0)]-th_p[np.maximum(IL,0)], 0.0)
        XI = XIc + dUp - R*(tu+tl)[:,None]*(n@Gt.T)
        XI[fade(r) < 1e-3] = 0
        th_p = th.copy(); dh_p = dh; f_p = f
    return fs, g

def harmonics(fs, g, M=8):
    gs = 0.5*(g + g[::-1])          # even (degenerate-endpoint) part
    V = []
    for m in range(1, M+1):
        Am = (2/d)*np.trapezoid(gs*np.cos(2*np.pi*m*fs/d), fs)
        V.append(-Am)
    return np.array(V), gs

Acell = np.sqrt(3)/2
print(f"pairs: {len(V0)}   r0 = {r0:.5f}   (bond checks: "
      f"dV(0.816)={Vmorse(np.sqrt(2/3))-Vmorse(1):.4f}, dV(1.291)={Vmorse(1.291)-Vmorse(1):.4f})")
fs, gc = run(fade_morse, r_ratio=0.0)
print(f"\ncentral only: gamma(d) = {gc[-1]:.5f}  max = {gc.max():.5f} at f/d = {fs[gc.argmax()]/d:.2f}")
for name, fade in [("Morse-tracking", fade_morse), ("linear .30", fade_lin), ("exponential", fade_exp)]:
    fs, g = run(fade)
    V, gs = harmonics(fs, g/Acell)          # per-area line
    S_shape = np.sum(np.arange(1,9)*V)/V[0]
    P = 8*np.pi**2*np.sum(np.arange(1,9)*V)*Acell/(1+r0)   # V were per-area; P wants per-cell
    print(f"\n[{name}] gamma(d) = {g[-1]:.5f}  gamma_max = {g.max():.5f} k_n/col "
          f"(= {g.max()/Acell:.5f} k_n per area)")
    print(f"  V_m/V_1 = {np.array2string(V/V[0], precision=4)}")
    print(f"  shape S = {S_shape:.5f}    absolute P = {P:.5f}")
print(f"\non-file target: S_alpha = 1.001067  (P is the comparable object)")
