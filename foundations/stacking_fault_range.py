#!/usr/bin/env python3
"""
stacking_fault_range.py — Why the FCC slice cannot hold an electron together
============================================================================
Mitchell A. Cox, University of the Witwatersrand
Accompanies: M. A. Cox, "The Cosserat Supersolid" (2026)

The compact electron rests on the stacking-fault energy being zero on the
three-dimensional slice and positive on the D4 wall.  This script establishes
how far that contrast survives, and it is a range question rather than a
dimension question at the first step.

1. SHELL CENSUS.  An intrinsic {111} fault slides the upper half-crystal
   through one Shockley partial.  That slide leaves the twelve-fold first
   neighbour shell and the six-fold second shell of every site untouched, and
   first disturbs the twenty-four-fold third at sqrt(3) ell.  So a central
   contact reaching no further than second neighbours gives gamma_FCC = 0
   identically, which is the standard result for close-packed metals.

2. RANGE SCAN.  Opening the contact up, at the spacing where the crystal is
   unstressed and with the anharmonicity xi = -7 the gravitational sector
   fixes, the rigid slice fault energy converges to about +2e-5 k_n per site.
   Positive, so the sign the compactness argument needs survives, but some
   four thousand times below the D4 value.  Since the dissociation width runs
   as 1/gamma, that is the difference between a partial pair at 0.2 ell and
   one drawn out past a thousand bond lengths.

3. D4 INVENTORY.  Rebuilt from the root lattice: the twelve crossing bonds
   carry the whole fault energy and the in-slice bonds contribute exactly
   nothing, which is the mechanism the monograph asserts.

Units: bond length ell = 1, energies in the normal contact stiffness
k_n = V''(1).

Usage:
    python stacking_fault_range.py
"""

import itertools

import numpy as np
from scipy.optimize import brentq

# ----------------------------------------------------------------- constants
S = 1.0 / np.sqrt(2.0)
D111 = np.sqrt(2.0 / 3.0)                 # {111} interplanar spacing
E1 = np.array([1.0, 0.0])
E2 = np.array([0.5, np.sqrt(3.0) / 2.0])
SHOCKLEY = (E1 + E2) / 3.0                # partial Burgers vector, |b_p| = 1/sqrt3
W3 = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
T3 = np.array([1.0, 1.0, -2.0]) / np.sqrt(6.0)
BP3 = SHOCKLEY[0] * 0.0                    # placeholder, set below
BP = (1.0 / np.sqrt(3.0)) * T3             # the same partial, in the 3D frame


# ------------------------------------------------------------------ contacts
def morse(r, a, r0, D=1.0, rcut=None):
    """Morse pair energy, optionally truncated and shifted at rcut."""
    v = D * (1.0 - np.exp(-a * (r - r0)))**2
    if rcut is None:
        return v - D
    vs = D * (1.0 - np.exp(-a * (rcut - r0)))**2
    return np.where(r < rcut, v - vs, 0.0)


def morse_derivs(r, a, r0, D=1.0):
    """(V', V'', V''') of the untruncated Morse."""
    e = np.exp(-a * (r - r0))
    return (2 * D * a * (1 - e) * e,
            2 * D * a**2 * (2 * e**2 - e),
            2 * D * a**3 * (-4 * e**2 + e))


# =========================================================== A. the FCC fault
def fcc_shell_census(nmax=10, nlay=8):
    """Neighbour shells of a fault-plane site, before and after the slide."""
    pts = np.array([m * E1 + n * E2
                    for m in range(-nmax, nmax + 1)
                    for n in range(-nmax, nmax + 1)])
    rp, rf = [], []
    for la in range(-nlay, nlay + 1):
        dz = la * D111
        xy = pts + (la % 3) * SHOCKLEY
        rp.append(np.sqrt((xy**2).sum(1) + dz**2))
        sl = SHOCKLEY if la >= 1 else np.zeros(2)
        rf.append(np.sqrt(((xy + sl)**2).sum(1) + dz**2))
    rp, rf = np.concatenate(rp), np.concatenate(rf)
    rp, rf = rp[rp > 1e-9], rf[rf > 1e-9]
    print("  shell      R        Z(perfect)  Z(faulted)   disturbed?")
    for nm, e in [("1NN", 1.0), ("2NN", np.sqrt(2)), ("3NN", np.sqrt(3)),
                  ("4NN", 2.0), ("5NN", np.sqrt(5)), ("6NN", np.sqrt(6))]:
        cp = int((np.abs(rp - e) < 1e-6).sum())
        cf = int((np.abs(rf - e) < 1e-6).sum())
        print(f"  {nm}    {e:7.4f}     {cp:5d}       {cf:5d}        "
              f"{'no' if cp == cf else 'YES'}")


def fcc_shells(rmax=8.0):
    """FCC neighbour distances (NN spacing = 1) and multiplicities."""
    v = [np.array([1, 1, 0]) / np.sqrt(2), np.array([0, 1, 1]) / np.sqrt(2),
         np.array([1, 0, 1]) / np.sqrt(2)]
    N = int(rmax) + 3
    d = [np.linalg.norm(i * v[0] + j * v[1] + k * v[2])
         for i in range(-N, N + 1) for j in range(-N, N + 1)
         for k in range(-N, N + 1)]
    d = np.array([x for x in d if 1e-9 < x <= rmax])
    rs = np.unique(np.round(d, 9))
    return rs, np.array([int((np.abs(d - r) < 1e-7).sum()) for r in rs])


def fcc_equilibrium_contact(a_seed=2.9):
    """
    The Morse that satisfies the framework's own two conditions at once:
    zero pressure at bond length 1, and bond anharmonicity xi = -7.

    xi = -7 fixes u = exp(a(r0-1)) = (7-a)/(14-4a), leaving a single unknown a,
    which the zero-pressure condition then pins.
    """
    rs, zs = fcc_shells(8.0)

    def pressure(a):
        u = (7.0 - a) / (14.0 - 4.0 * a)
        r0 = 1.0 + np.log(u) / a
        d1, _, _ = morse_derivs(rs, a, r0)
        return float(np.sum(zs * rs * d1))

    a = brentq(pressure, 2.40, 3.45, xtol=1e-13)
    u = (7.0 - a) / (14.0 - 4.0 * a)
    r0 = 1.0 + np.log(u) / a
    _, d2, _ = morse_derivs(np.array([1.0]), a, r0)
    return a, r0, 1.0 / float(d2[0])          # D normalised so V''(1) = 1


def fcc_fault_energy(a, r0, D, rcut, nlay=18, nmax=22):
    """Rigid intrinsic-fault energy per fault-plane site."""
    pts = np.array([m * E1 + n * E2
                    for m in range(-nmax, nmax + 1)
                    for n in range(-nmax, nmax + 1)])
    dE = 0.0
    for kb in range(nlay):
        lb = -kb
        zb, rb = lb * D111, (lb % 3) * SHOCKLEY
        for la in range(1, nlay + 1):
            dz = la * D111 - zb
            if dz > rcut:
                continue
            xy = pts + (la % 3) * SHOCKLEY - rb
            rp = np.sqrt((xy**2).sum(1) + dz**2)
            rf = np.sqrt(((xy + SHOCKLEY)**2).sum(1) + dz**2)
            m = (rp < rcut) | (rf < rcut)
            dE += float((morse(rf[m], a, r0, D, rcut)
                         - morse(rp[m], a, r0, D, rcut)).sum())
    return dE


# ============================================================ B. the D_4 wall
def d4_nodes(nbox=9, ncomp=4):
    """D_4 = even-sum integer 4-vectors, scaled so the 24 minimal vectors are 1."""
    P, C = [], []
    for a in range(-nbox, nbox + 1):
        for b in range(-nbox, nbox + 1):
            for c in range(-nbox, nbox + 1):
                for e in range(-ncomp, ncomp + 1):
                    if (a + b + c + e) % 2 == 0:
                        P.append(S * np.array([a, b, c], float))
                        C.append(S * e)
    return np.array(P), np.array(C)


def d4_fault_inventory(a, r0, D, rcut=1.35, nbox=9):
    """Bond-by-bond energy change across the D_4 wall, for one column."""
    P, C = d4_nodes(nbox)
    Z = P @ W3
    tol = 1e-7
    cand = np.where((np.abs(C) < tol) & (np.abs(Z) < tol))[0]
    i0 = cand[np.argmin(np.linalg.norm(P[cand], axis=1))]
    above = Z > tol
    dsp, dc = P[above] - P[i0], C[above] - C[i0]
    rp = np.sqrt((dsp**2).sum(1) + dc**2)
    rf = np.sqrt(((dsp + BP)**2).sum(1) + dc**2)
    m = (rp < rcut) | (rf < rcut)
    dv = morse(rf[m], a, r0, D, rcut) - morse(rp[m], a, r0, D, rcut)
    return rp[m], rf[m], dv, np.abs(dc[m])

if __name__ == "__main__":
    a, r0, D = fcc_equilibrium_contact()

    print("Which shells the intrinsic {111} fault disturbs:\n")
    fcc_shell_census()
    print("\n  First and second shells preserved exactly: a contact reaching")
    print("  only that far gives gamma_FCC = 0 identically.\n")

    print("Self-consistent long-range contact (zero pressure at ell, xi = -7):")
    print(f"  a = {a:.6f}   r0 = {r0:.6f}   D = {D:.6f} k_n\n")
    print("  r_cut     gamma_FCC (k_n per fault-plane site)")
    for rc in (1.2, 1.5, 1.8, 2.2, 2.6, 3.0, 3.5, 4.0, 5.0):
        print(f"  {rc:4.2f}      {fcc_fault_energy(a, r0, D, rc):+.6e}")

    print("\nD_4 wall, bond by bond, one column:\n")
    aM, r0M, DM = 7.0 / 3.0, 1.0, 1.0 / (2 * (7.0 / 3.0) ** 2)
    rp, rf, dv, dc = d4_fault_inventory(aM, r0M, DM)
    print("  r_perfect  r_faulted    dE (k_n)   compact component")
    for x, y, z, w in sorted(zip(np.round(rp, 4), np.round(rf, 4),
                                 np.round(dv, 5), np.round(dc, 4))):
        print(f"   {x:7.4f}    {y:7.4f}   {z:+9.5f}      {w:.4f}")
    print(f"\n  crossing bonds  {dv[dc > 1e-9].sum():+.6f} k_n")
    print(f"  in-slice bonds  {dv[dc < 1e-9].sum():+.6f} k_n")
    print(f"  column total    {dv.sum():+.6f} k_n")
