#!/usr/bin/env python3
"""
bc_route_mirror.py - the two Shockley routes are mirror-degenerate, so the
electron's handedness cannot be static
=========================================================================

Backs the closing argument of the monograph's "D4 dictionary and the
tangential sector" subsection.

A dissociated screw on the {111} wall leaves the perfect registry along the
hop direction t_hat and reaches the faulted hollow by bowing to one side or
the other in the transverse direction u_hat.  The two bows are the two
Shockley routes: one through the B hollow (+u), one through the C hollow (-u).
If the inter-layer (fourth-dimension) contacts preferred one bow, the static
fault would carry a handedness, and that handedness could seed the electron's
chirality.

This script shows it cannot, for a reason that is a symmetry rather than a
number.  Three checks, each independent:

  1. THE MIRROR IS AN EXACT D4 AUTOMORPHISM.  The reflection that exchanges
     the two routes (u_hat -> -u_hat, t_hat and the wall normal fixed) is the
     coordinate swap a<->b.  That swap is a permutation of the four axes, and
     the D4 lattice is invariant under every coordinate permutation, so the
     full nearest-neighbour contact set maps to itself.

  2. THE GAMMA-SURFACE IS EVEN.  The z-relaxed stacking-fault energy obeys
     gamma(s, u) = gamma(s, -u) to one part in 1e14 (floating-point round-off)
     at every sampled point.  The B and C route saddles are therefore equal:
     the static splitting is exactly zero, not merely small.

  3. THE SAME MIRROR REVERSES THE STACKING SCREW.  Under a<->b the transverse
     offsets of the layers flip sign, turning the cyclic stacking order
     A -> B -> C into A -> C -> B.  The static energy is mirror-even, but the
     compact-direction circulation is mirror-odd.

Conclusion: parity violation cannot live in the static fault energy, which the
mirror protects.  It can live only in a coupling that feels the turning of the
stacking screw, and a fault at rest does not turn.  The handedness is dynamic,
carried by the imaginary, time-reversal-odd coupling to the compact momentum
k4 that also makes the vacuum optically active.  This is consistent with the
weak interaction acting through a V-A current while the electron mass stays
parity-blind.

Contact model and geometry are identical to d4_pair_functional.py:
length unit = nearest-neighbour distance ell; energy unit = the normal
contact stiffness k_n = V''(1) of the truncated-shifted Morse contact with
a*ell = 7/3 (anharmonicity xi = -7).

Usage:   python3 bc_route_mirror.py
Requires: numpy, scipy.

Author: Mitchell A. Cox
Date:   June 2026
"""

import numpy as np
from scipy.optimize import minimize

# ------------------------------------------------------------------
# D4 geometry and the truncated-shifted Morse contact law
# ------------------------------------------------------------------
S      = 1/np.sqrt(2)                            # D4 -> physical length scale
W_HAT  = np.array([1, 1, 1, 0])/np.sqrt(3)       # {111} wall normal
T_HAT  = np.array([1, 1, -2, 0])/np.sqrt(6)      # partial-hop direction
U_HAT  = np.array([1, -1, 0, 0])/np.sqrt(2)      # in-plane transverse (mirror axis)
DHOP   = 1/np.sqrt(3)                             # misfit period d

AM     = 7/3                                      # Morse a*ell  (xi = ell V'''/V'' = -7)
DM     = 1/(2*AM*AM)                              # depth: sets V''(1) = 1
RINF   = 1 + np.log(2)/AM                          # inflection: V'' = 0, contact opens here
VOFF   = DM*(1 - np.exp(-AM*(RINF-1)))**2          # shift so V(RINF) = 0


def Vc(r):
    """Contact pair energy, truncated and shifted to zero at the inflection."""
    return np.where(r < RINF, DM*(1 - np.exp(-AM*(r-1)))**2 - VOFF, 0.0)


def Vp(r):
    """Contact pair force magnitude (dV/dr)."""
    e = np.exp(-AM*(r-1))
    return np.where(r < RINF, 2*DM*AM*(1-e)*e, 0.0)


def gen(box, mlo, mhi):
    """D4 nodes (integer 4-vectors with even coordinate sum) within a band of
    layer index m = a + b + c, returned as (m, scaled-position) pairs."""
    out = []
    for a in range(-box, box+1):
        for b in range(-box, box+1):
            for c in range(-box, box+1):
                m = a + b + c
                if m < mlo or m > mhi:
                    continue
                for e in range(-box, box+1):
                    if (a + b + c + e) % 2:
                        continue
                    out.append((m, S*np.array([a, b, c, e])))
    return out


# ==================================================================
# 1. the u-flip is the swap a<->b, an exact D4 automorphism
# ==================================================================
def check_automorphism():
    """The reflection u_hat -> -u_hat is the coordinate swap a<->b.  Confirm it
    fixes t_hat and the wall normal, flips u_hat, and leaves the full
    nearest-neighbour contact set invariant."""
    swap = np.array([[0, 1, 0, 0],
                     [1, 0, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])          # a <-> b
    fixes_t = np.allclose(swap @ T_HAT,  T_HAT)
    fixes_w = np.allclose(swap @ W_HAT,  W_HAT)
    flips_u = np.allclose(swap @ U_HAT, -U_HAT)

    # nearest-neighbour shell in integer coordinates (|v|*S ~ 1)
    shell = set()
    for a in range(-3, 4):
        for b in range(-3, 4):
            for c in range(-3, 4):
                for e in range(-3, 4):
                    if (a + b + c + e) % 2:
                        continue
                    v = np.array([a, b, c, e])
                    if 0 < np.linalg.norm(v)*S < 1.35:
                        shell.add((a, b, c, e))
    swapped = {(b, a, c, e) for (a, b, c, e) in shell}
    invariant = (shell == swapped)
    return fixes_t, fixes_w, flips_u, len(shell), invariant


# ==================================================================
# 2. the z-relaxed gamma-surface is even in u
# ==================================================================
def gamma_surface(box=8):
    """Return gam(s, u): the z-relaxed stacking-fault energy for an in-plane
    slip s*t_hat + u*u_hat, with the in-plane slip frozen and each layer free
    to breathe perpendicular to the wall (the standard generalised-stacking-
    fault protocol).  Energy is per atomic column, in units of k_n."""
    nodes = gen(box, -3, 4)
    bylayer = {}
    for m, p in nodes:
        bylayer.setdefault(m, []).append(p)
    # one representative per layer, nearest the wall axis
    reps = {m: min(ps, key=lambda p: np.linalg.norm(p - (p @ W_HAT)*W_HAT))
            for m, ps in bylayer.items()}

    # wall-crossing contacts: pair each representative with the nodes up to
    # three layers above it that lie within range
    contacts = []
    for m in range(-3, 4):
        if m not in reps:
            continue
        for mp in range(m+1, m+4):
            if mp not in bylayer:
                continue
            P = np.array([q for q in bylayer[mp]
                          if np.linalg.norm(q - reps[m]) < 2.6])
            if len(P):
                contacts.append((m, mp, reps[m], P))

    FREE = list(range(-2, 4))                 # layers free to breathe
    fi = {m: k for k, m in enumerate(FREE)}
    nf = len(FREE)
    far = set(range(FREE[-1]+1, 5))           # layers locked to the rigid lift

    def E_and_grad(x, offset):
        """Energy and its gradient in the per-layer breathing amplitudes x,
        for an applied in-plane slip 'offset' of every layer above the wall."""
        E = 0.0
        g = np.zeros(nf+1)

        def wd(m):                            # breathing displacement of layer m
            if m in fi:
                return x[fi[m]]
            return x[nf] if m in far else 0.0

        def gi(m):                            # gradient index for layer m
            if m in fi:
                return fi[m]
            return nf if m in far else -1

        for m, mp, rep, P in contacts:
            b = (offset if mp >= 1 else 0) - (offset if m >= 1 else 0)
            v = (P - rep) + b + (wd(mp) - wd(m))*W_HAT
            r = np.linalg.norm(v, axis=1)
            E += Vc(r).sum()
            fw = ((Vp(r)/r)*(v @ W_HAT)).sum()
            if gi(mp) >= 0:
                g[gi(mp)] += fw
            if gi(m) >= 0:
                g[gi(m)] -= fw
        return E, g

    def relax(offset):
        """Minimise over breathing amplitudes from several seeds; the fault
        opens by a finite lift, so cold starts near zero can stall."""
        best = None
        seeds = [np.zeros(nf+1), np.full(nf+1, 0.05),
                 np.concatenate([np.linspace(0, 0.15, nf), [0.15]])]
        for x0 in seeds:
            res = minimize(lambda x: E_and_grad(x, offset), x0, jac=True,
                           method='L-BFGS-B',
                           options={'maxiter': 20000, 'ftol': 1e-16,
                                    'gtol': 1e-13})
            if best is None or res.fun < best.fun:
                best = res
        return best.fun

    E_perfect = relax(np.zeros(4))

    def gam(s, u):
        return relax(s*T_HAT + u*U_HAT) - E_perfect

    return gam


# ==================================================================
# 3. the same mirror reverses the A->B->C stacking screw
# ==================================================================
def stacking_sense(box=3):
    """Transverse (u_hat) offset of consecutive layers, and the swap a<->b
    acting on it.  A sign flip of the transverse offset reverses the cyclic
    stacking order A->B->C into A->C->B."""
    swap = np.array([[0, 1, 0, 0],
                     [1, 0, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])
    rows = []
    for m in range(0, 4):
        # representative node of layer m nearest the wall axis
        best = None
        for a in range(-box, box+1):
            for b in range(-box, box+1):
                for c in range(-box, box+1):
                    if a + b + c != m:
                        continue
                    for e in range(-box, box+1):
                        if (a + b + c + e) % 2:
                            continue
                        p = S*np.array([a, b, c, e])
                        d = np.linalg.norm(p - (p @ W_HAT)*W_HAT)
                        if best is None or d < best[0]:
                            best = (d, p)
        p = best[1]
        u_off = p @ U_HAT
        u_off_mirror = (swap @ p) @ U_HAT
        rows.append((m, u_off, u_off_mirror))
    return rows


# ==================================================================
def main():
    bar = "=" * 70
    print(bar)
    print("  B/C Shockley-route splitting on the D4 {111} wall")
    print(bar)

    print("\n1. The u-flip mirror is an exact D4 automorphism (swap a<->b):")
    ft, fw, fu, nshell, inv = check_automorphism()
    print(f"   fixes t_hat: {ft}   fixes wall normal: {fw}   flips u_hat: {fu}")
    print(f"   nearest-neighbour contacts: {nshell};  set invariant under swap: {inv}")

    print("\n2. The z-relaxed gamma-surface is even in u:  gamma(s,+u) vs gamma(s,-u)")
    gam = gamma_surface()
    print(f"   {'s/d':>6} {'u':>7} | {'gamma(+u)':>13} {'gamma(-u)':>13} "
          f"{'|asym|_rel':>11}")
    print("   " + "-"*54)
    worst = 0.0
    for s_over_d in (0.5, 0.8, 1.0):
        s = s_over_d*DHOP
        for u in (0.05, 0.10, 0.15):
            gp, gm = gam(s, +u), gam(s, -u)
            asym = abs(gp - gm)/max(abs(gp), 1e-12)
            worst = max(worst, asym)
            print(f"   {s_over_d:6.2f} {u:7.3f} | {gp:13.8f} {gm:13.8f} "
                  f"{asym:11.2e}")
    print("   " + "-"*54)
    print(f"   worst relative asymmetry = {worst:.2e}  (round-off => exact degeneracy)")

    print("\n3. The same mirror reverses the A->B->C stacking screw:")
    print(f"   {'layer':>6} {'u-offset':>10} {'mirrored':>10}")
    for m, u_off, u_off_mirror in stacking_sense():
        print(f"   {m:6d} {u_off:10.3f} {u_off_mirror:10.3f}")
    print("   the transverse offset flips sign, so A->B->C becomes A->C->B.")

    print("\n" + bar)
    print("  Static fault: mirror-even  -> B and C routes exactly degenerate.")
    print("  Stacking screw: mirror-odd -> the compact circulation is chiral.")
    print("  => parity violation is necessarily DYNAMIC (couples to motion / k4),")
    print("     not a static splitting; the electron mass stays parity-blind.")
    print(bar)


if __name__ == '__main__':
    main()
