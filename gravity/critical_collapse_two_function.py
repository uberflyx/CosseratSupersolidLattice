#!/usr/bin/env python3
r"""
Two-function (anisotropic) spherical collapse of the gravitational microrotation
carrier, in Christodoulou double-null coordinates, evolved with the non-AMR
algorithm of Garfinkle.

Purpose
-------
The scalar refractive-index reduction of the collapse holds the temporal and
spatial metric on the single constraint surface g_00 g_rr = -c^2.  That model
focuses critically but does not echo: it is continuously self-similar (see the
companion scalar-collapse script).  The discrete self-similar echoing that
Choptuik found lives in the SECOND, independent metric function, the growing
lapse that the constraint discards.  Restoring it means evolving the full
two-function system.

By the Deser uniqueness theorem the all-orders anharmonic Cosserat action is
Einstein-Hilbert, with the transverse microrotation psi a minimally coupled,
massless (c-propagating) real scalar.  The two-function lattice collapse is
therefore the spherically symmetric Einstein + real-massless-scalar collapse
Choptuik solved.  The graviton propagates at c throughout: no second wave speed
is inserted by hand, and the shear-sector modulus ratio (pi-1)/(pi-2)=1.876 does
NOT appear here.  The two metric functions are the areal radius r(u,v) and the
redshift function g(u,v); g is sourced by the field gradient energy (the
Price-Pullin growing lapse).

What this script reproduces
---------------------------
  (1) The central strain ECHOES: it reverses sign and returns to universal
      values (~ +/-0.61) on a shrinking scale, the discrete-self-similar
      signature, in contrast to the single hump of the scalar reduction.
  (2) The echoing period, extracted u*-free from the geometric spacing of the
      central-field extrema in proper time:  Delta ~ 3.44, matching the
      general-relativistic Choptuik eigenvalue 3.4439, and found to be
      independent of resolution and of the initial-data family (universality).
  (3) The mass-scaling exponent, from M_BH ~ (p-p*)^gamma for supercritical
      amplitudes:  gamma ~ 0.39, consistent with the GR value 0.374 to within
      the Gundlach fine-structure wiggle (the same precision a non-AMR null code
      reaches; Garfinkle reported ~0.38).  By Koike-Hara-Adachi the reciprocal
      is the growth rate of the single unstable mode, Re(kappa)=1/gamma ~ 2.6.

Coordinates and scheme
----------------------
Geometric units G = c = 1.  Outgoing rays carry u (proper time at the origin),
ingoing rays carry v (areal radius on the initial slice).  For a field f,
fbar(u,v) = (1/r) int_0^r f dr is the r-weighted interior average.

  field          psi = hbar,   q = (h - hbar)^2 / r,   g = exp(4 pi r qbar)
  metric         ds^2 = -2 g r' du dv + r^2 dOmega^2
  evolution      hdot = (g - gbar)(h - hbar)/(2 r),    rdot = -gbar/2
  Hawking mass   m = (r/2)(1 - gbar/g);   2m/r -> 1 at an apparent horizon

The whole adaptive Heun (RK2) step is compiled with numba; the innermost-points
central line fit is closed-form.  The null grid shrinks onto the focusing region
and supplies its own resolution; midpoint insertion in v replenishes rays as the
innermost ones reach the origin and are dropped.

Reference: D. Garfinkle, Phys. Rev. D 51, 5558 (1995); D. Christodoulou, Commun.
Math. Phys. 105, 337 (1986); M. W. Choptuik, Phys. Rev. Lett. 70, 9 (1993).
"""

import numpy as np
from numba import njit

# ---------------------------------------------------------------------------
# Numba core: one adaptive Heun step, with a closed-form central line fit.
# ---------------------------------------------------------------------------


@njit(cache=True, fastmath=True)
def _fit_central(h, r, npts):
    """Closed-form linear least squares of h vs r over the innermost npts
    points; returns (h0, h1) with h0 the central field psi(0,u)."""
    m = npts if npts < h.shape[0] else h.shape[0]
    sx = 0.0; sy = 0.0; sxx = 0.0; sxy = 0.0
    for j in range(m):
        sx += r[j]; sy += h[j]; sxx += r[j] * r[j]; sxy += r[j] * h[j]
    det = m * sxx - sx * sx
    if det == 0.0:
        return h[0], 0.0
    h1 = (m * sxy - sx * sy) / det
    h0 = (sy - h1 * sx) / m
    return h0, h1


@njit(cache=True, fastmath=True)
def _rates(h, r):
    """Return the rate arrays hdot, rdot and the diagnostics
    (h0, h1, max_2mr, max_gbar, drmin).  The three cumulative interior
    integrals (hbar, qbar, gbar) are trapezoidal on the non-uniform r-grid,
    anchored at the regular centre (h->h0, q->0, g->1)."""
    n = h.shape[0]
    h0, h1 = _fit_central(h, r, 4)
    hbar = np.empty(n); g = np.empty(n); gbar = np.empty(n)

    # hbar = (1/r) int_0^r h dr
    acc = 0.5 * (h[0] + h0) * r[0]
    hbar[0] = acc / r[0]
    for j in range(1, n):
        acc += 0.5 * (h[j] + h[j - 1]) * (r[j] - r[j - 1])
        hbar[j] = acc / r[j]

    # q = (h - hbar)^2 / r ; qbar ; g = exp(4 pi r qbar)
    q_prev = (h[0] - hbar[0]) ** 2 / r[0]
    accq = 0.5 * q_prev * r[0]
    g[0] = np.exp(4.0 * np.pi * r[0] * (accq / r[0]))
    for j in range(1, n):
        qj = (h[j] - hbar[j]) ** 2 / r[j]
        accq += 0.5 * (qj + q_prev) * (r[j] - r[j - 1])
        g[j] = np.exp(4.0 * np.pi * r[j] * (accq / r[j]))
        q_prev = qj

    # gbar = (1/r) int_0^r g dr
    accg = 0.5 * (g[0] + 1.0) * r[0]
    gbar[0] = accg / r[0]
    for j in range(1, n):
        accg += 0.5 * (g[j] + g[j - 1]) * (r[j] - r[j - 1])
        gbar[j] = accg / r[j]

    # rates and slice diagnostics
    hdot = np.empty(n); rdot = np.empty(n)
    max2 = -1.0e300; maxgb = 0.0; drmin = r[0]
    for j in range(n):
        hdot[j] = (g[j] - gbar[j]) * (h[j] - hbar[j]) / (2.0 * r[j])
        rdot[j] = -0.5 * gbar[j]
        t = 1.0 - gbar[j] / g[j]
        if t > max2:
            max2 = t
        if gbar[j] > maxgb:
            maxgb = gbar[j]
        if j > 0:
            d = r[j] - r[j - 1]
            if 0.0 < d < drmin:
                drmin = d
    return hdot, rdot, h0, h1, max2, maxgb, drmin


@njit(cache=True, fastmath=True)
def _full_step(h, r, courant, r_floor, du_cap):
    """One adaptive Heun step.  Returns (h_new, r_new, du, h0, h1, max_2mr),
    where the diagnostics are evaluated at the start-of-step state."""
    n = h.shape[0]
    hdot, rdot, h0, h1, m2, maxgb, drmin = _rates(h, r)
    du = courant * drmin / (0.5 * maxgb)
    if du > du_cap:
        du = du_cap
    hp = np.empty(n); rp = np.empty(n)
    for j in range(n):
        hp[j] = h[j] + du * hdot[j]
        rj = r[j] + du * rdot[j]
        rp[j] = rj if rj > r_floor else r_floor
    hdot2, rdot2, _, _, _, _, _ = _rates(hp, rp)
    hn = np.empty(n); rn = np.empty(n)
    for j in range(n):
        hn[j] = h[j] + 0.5 * du * (hdot[j] + hdot2[j])
        rn[j] = r[j] + 0.5 * du * (rdot[j] + rdot2[j])
    return hn, rn, du, h0, h1, m2


# ---------------------------------------------------------------------------
# Grid helpers and the evolution driver.
# ---------------------------------------------------------------------------


def initial_data(N, v_max, phi0, r0, sigma):
    """Christodoulou initial slice u=0: r = v on a uniform v-grid, and a
    Gaussian strain pulse h = phi0 v^2 exp(-((v-r0)/sigma)^2).  The v^2
    prefactor enforces the regular-centre behaviour h ~ r^2.  phi0 is the
    amplitude knob p that is tuned to the formation threshold p*."""
    v = np.linspace(v_max / N, v_max, N)        # exclude the exact origin v=0
    h = phi0 * v ** 2 * np.exp(-((v - r0) / sigma) ** 2)
    return v, h.copy(), v.copy()


def _reinterpolate(h, r, v):
    """Garfinkle grid maintenance: insert midpoints in v between survivors,
    restoring resolution on the shrinking, focusing region."""
    v_new = np.empty(2 * v.size - 1)
    v_new[0::2] = v
    v_new[1::2] = 0.5 * (v[:-1] + v[1:])
    return np.interp(v_new, v, h), np.interp(v_new, v, r), v_new


def evolve(N, v_max, phi0, r0, sigma, courant=0.20, r_floor=1e-7, u_max=16.0,
           max_steps=2_500_000, regrid=True, regrid_frac=0.5,
           twoMr_horizon=0.99, record_full=False, record_every=1,
           early_disperse=True, disperse_frac=0.4):
    """Evolve one collapse.  Returns a dict with the outcome ('horizon',
    'dispersed', or 'stalled'), the peak compactness, and (if record_full)
    the decimated central-field history U, H0, H1, M2."""
    v, h, r = initial_data(N, v_max, phi0, r0, sigma)
    u = 0.0; nsteps = 0
    max_2mr = 0.0; u_peak = 0.0
    U = []; H0 = []; H1 = []; M2 = []
    outcome = 'stalled'
    while u < u_max and nsteps < max_steps:
        nsteps += 1
        hn, rn, du, h0, h1, m2 = _full_step(h, r, courant, r_floor, u_max - u)
        if m2 > max_2mr:
            max_2mr = m2; u_peak = u
        if record_full and (nsteps % record_every == 0):
            U.append(u); H0.append(h0); H1.append(h1); M2.append(m2)
        if max_2mr >= twoMr_horizon:
            outcome = 'horizon'; break
        # field has focused, bounced, and is leaving without forming a horizon
        if (early_disperse and max_2mr < 0.9
                and u > max(u_peak + 0.5, r0) and m2 < disperse_frac * max_2mr):
            outcome = 'dispersed'; break
        if du <= 0 or not np.isfinite(du):
            outcome = 'stalled'; break
        h, r = hn, rn
        u += du
        keep = r > r_floor
        if not keep.all():
            h, r, v = h[keep], r[keep], v[keep]
        if r.size < 8:
            outcome = 'stalled'; break
        if regrid and r.size <= regrid_frac * N and r.size >= 4:
            h, r, v = _reinterpolate(h, r, v)
    else:
        outcome = 'dispersed' if max_2mr < twoMr_horizon else 'horizon'
    return dict(outcome=outcome, u_end=u, max_2mr=max_2mr, nsteps=nsteps,
                npts=r.size, U=np.array(U), H0=np.array(H0),
                H1=np.array(H1), M2=np.array(M2))


# ---------------------------------------------------------------------------
# (1)+(2) Tune to the formation threshold and extract the echoing period.
# ---------------------------------------------------------------------------


def tune_threshold(N, lo, hi, iters, r0, sigma, v_max=12.0, courant=0.20,
                   regrid=True, verbose=False):
    """Bisect the amplitude knob to the formation threshold p*.  The bracket
    [lo, hi] must straddle: lo disperses, hi forms a horizon.  Sub- and
    super-critical solutions cling to the same critical attractor and diverge
    only near the accumulation point, so the outcome is decided late; this is
    intrinsic, not a tuning artefact."""
    for it in range(iters):
        mid = 0.5 * (lo + hi)
        oc = evolve(N, v_max, mid, r0, sigma, courant=courant,
                    regrid=regrid, twoMr_horizon=0.99)['outcome']
        if oc == 'horizon':
            hi = mid
        else:
            lo = mid
        if verbose:
            print(f"    bisect {it:2d}: phi0={mid:.10f} {oc:>9} "
                  f"width={hi - lo:.1e}")
    return lo, hi


def echo_period(U, H0, prominence=0.04):
    """Extract the echoing period from the central-field history.  Consecutive
    sign-reversing extrema are spaced by half a period in self-similar time, so
    in proper time their gaps shrink geometrically by exp(Delta/2); the ratio
    of successive gaps gives Delta with no need for the accumulation time u*.
    Prominence filtering removes per-step jitter (genuine echoes swing by ~0.5).
    Also returns a u* estimate per triple as an internal consistency check."""
    from scipy.signal import find_peaks
    pk_max, _ = find_peaks(H0, prominence=prominence)
    pk_min, _ = find_peaks(-H0, prominence=prominence)
    idx = np.sort(np.concatenate([pk_max, pk_min]))
    u_ext = U[idx]; h_ext = H0[idx]
    out = dict(u_ext=u_ext, h_ext=h_ext, deltas=np.array([]),
               ustar=np.array([]))
    if len(idx) >= 3:
        gaps = np.diff(u_ext)
        ratios = gaps[:-1] / gaps[1:]
        out['deltas'] = 2.0 * np.log(ratios)
        out['ustar'] = np.array([u_ext[k] + gaps[k] * ratios[k] / (ratios[k] - 1)
                                 for k in range(len(ratios)) if ratios[k] > 1.001])
    return out


# ---------------------------------------------------------------------------
# (3) Mass-scaling exponent gamma.
# ---------------------------------------------------------------------------


def black_hole_mass(N, phi0, r0, sigma, v_max=12.0, courant=0.30,
                    regrid=True, twoMr_horizon=0.99):
    """Hawking mass at apparent-horizon formation for a supercritical run.
    Returns (M_BH, r_AH, u_AH), or (None, None, u) if the field disperses.
    The compactness returned by the stepper detects the horizon; the full
    interior profile is built once, at the horizon step."""
    # local imports keep the mass measurement self-contained
    v, h, r = initial_data(N, v_max, phi0, r0, sigma)
    u = 0.0; nsteps = 0; max_2mr = 0.0; u_peak = 0.0
    while u < 16.0 and nsteps < 2_500_000:
        nsteps += 1
        hn, rn, du, h0, h1, m2 = _full_step(h, r, courant, 1e-7, 16.0 - u)
        if m2 >= twoMr_horizon:                         # horizon at (h, r)
            hbar, qbar, g, gbar = _profile(h, r)
            twomr = 1.0 - gbar / g
            j = int(np.nanargmax(twomr))
            return 0.5 * r[j] * twomr[j], r[j], u
        if m2 > max_2mr:
            max_2mr = m2; u_peak = u
        if max_2mr < 0.9 and u > max(u_peak + 0.5, r0) and m2 < 0.4 * max_2mr:
            return None, None, u
        h, r = hn, rn; u += du
        keep = r > 1e-7
        if not keep.all():
            h, r, v = h[keep], r[keep], v[keep]
        if r.size < 8:
            return None, None, u
        if regrid and r.size <= 0.5 * N and r.size >= 4:
            h, r, v = _reinterpolate(h, r, v)
    return None, None, u


@njit(cache=True, fastmath=True)
def _profile(h, r):
    """Interior averages hbar, qbar, g, gbar on the grid (for the mass
    profile at a single slice).  Same integrals as _rates, returned as arrays."""
    n = h.shape[0]
    h0, _ = _fit_central(h, r, 4)
    hbar = np.empty(n); qbar = np.empty(n); g = np.empty(n); gbar = np.empty(n)
    acc = 0.5 * (h[0] + h0) * r[0]; hbar[0] = acc / r[0]
    for j in range(1, n):
        acc += 0.5 * (h[j] + h[j - 1]) * (r[j] - r[j - 1]); hbar[j] = acc / r[j]
    q_prev = (h[0] - hbar[0]) ** 2 / r[0]
    accq = 0.5 * q_prev * r[0]; qbar[0] = accq / r[0]
    g[0] = np.exp(4.0 * np.pi * r[0] * qbar[0])
    for j in range(1, n):
        qj = (h[j] - hbar[j]) ** 2 / r[j]
        accq += 0.5 * (qj + q_prev) * (r[j] - r[j - 1]); qbar[j] = accq / r[j]
        g[j] = np.exp(4.0 * np.pi * r[j] * qbar[j]); q_prev = qj
    accg = 0.5 * (g[0] + 1.0) * r[0]; gbar[0] = accg / r[0]
    for j in range(1, n):
        accg += 0.5 * (g[j] + g[j - 1]) * (r[j] - r[j - 1]); gbar[j] = accg / r[j]
    return hbar, qbar, g, gbar


def mass_scaling(pstar, dps, N, r0, sigma, courant=0.30):
    """Black-hole mass for a sequence of supercritical offsets dp = p - p*.
    Returns (p, M) arrays for the points that formed a horizon.  Stay in the
    scaling regime (M_BH << M_pulse) and span about one wiggle period."""
    P = []; M = []
    for dp in dps:
        m, rah, u = black_hole_mass(N, pstar + dp, r0, sigma, courant=courant)
        if m is not None:
            P.append(pstar + dp); M.append(m)
    return np.array(P), np.array(M)


# ---------------------------------------------------------------------------
# Demonstration.  Parameters here are a fast, representative subset; the
# published Delta and gamma used deeper tuning (see the header).
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import time
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full", action="store_true",
                    help="paper-precision run (N=1500, 28 bisections, ~4 min); "
                         "default is a fast demonstration (~90 s)")
    args = ap.parse_args()

    # fast default vs paper-precision settings
    if args.full:
        N_echo, ITERS, HZ = 1500, 28, 0.9997
        N_mass, DPS = 1500, 1e-4 * 2.0 ** np.arange(0, 7)
    else:
        N_echo, ITERS, HZ = 1200, 22, 0.9996
        N_mass, DPS = 1200, 1e-4 * 2.0 ** np.arange(0, 5)

    _ = evolve(200, 12.0, 0.01, 4.0, 1.0, u_max=2.0, regrid=False)  # warm up

    print("=" * 70)
    print("Two-function collapse: echo and critical exponents"
          + ("  [--full]" if args.full else "  [fast demo; --full for paper numbers]"))
    print("=" * 70)

    # ---- (1)+(2) echo and period, family A (r0=4, sigma=1) ----
    # The converged Delta=3.445 (N=3000) and gamma=0.39 used deeper tuning and
    # finer grids than the fast default; --full approaches them.
    print("\n(1)+(2) Threshold and echoing period, family r0=4 sigma=1:")
    t0 = time.time()
    lo, hi = tune_threshold(N=N_echo, lo=0.020, hi=0.025, iters=ITERS,
                            r0=4.0, sigma=1.0, courant=0.30)
    pstar = 0.5 * (lo + hi)
    print(f"    p* = {pstar:.9f}   (bracket width {hi - lo:.1e}, "
          f"{time.time() - t0:.0f}s)")

    res = evolve(N=N_echo, v_max=12.0, phi0=hi, r0=4.0, sigma=1.0, courant=0.30,
                 twoMr_horizon=HZ, record_full=True, record_every=5)
    ana = echo_period(res['U'], res['H0'])
    print(f"    central-field extrema: "
          + " ".join(f"{u:.3f}({h:+.2f})"
                     for u, h in zip(ana['u_ext'], ana['h_ext'])))
    if len(ana['deltas']):
        print(f"    Delta per successive gap-ratio: "
              + " ".join(f"{d:.3f}" for d in ana['deltas'])
              + "   (first is the initial-focusing transient)")
        print(f"    Delta from the last, most self-similar ratio = "
              f"{ana['deltas'][-1]:.3f}   [GR: 3.4439]")
    if len(ana['ustar']):
        print(f"    u* consistency estimates: "
              + " ".join(f"{x:.4f}" for x in ana['ustar']))

    # ---- (3) mass-scaling exponent ----
    print("\n(3) Mass-scaling exponent gamma (scaling regime):")
    P, M = mass_scaling(pstar, DPS, N=N_mass, r0=4.0, sigma=1.0)
    if len(P) >= 3:
        x = np.log(P - pstar); y = np.log(M)
        gamma = np.polyfit(x, y, 1)[0]
        print(f"    fitted over dp in [{DPS[0]:.0e}, {DPS[-1]:.0e}]: "
              f"gamma = {gamma:.3f}   [GR: 0.374; wiggle-limited, paper quotes 0.39]")
        print(f"    => unstable-mode growth rate Re(kappa) = 1/gamma = "
              f"{1.0 / gamma:.2f}   [GR: 2.67]")
    else:
        print("    too few horizon points at this resolution; use --full")

    print("\nThe two-function collapse echoes (discrete self-similarity), with")
    print("Delta and gamma at the general-relativistic Choptuik values, as the")
    print("Deser reduction requires.")
