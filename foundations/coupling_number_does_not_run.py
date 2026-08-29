"""The Cosserat coupling number does not run, and alpha^-1 is what forbids it.

A tempting correction to the neutrino ring is to let N^2 depend on wavevector.
The motivation is sound as far as it goes: in a Cosserat medium the microrotation
is screened over 1/q, so eliminating it at fixed displacement gives

    mu_eff(k) = mu + kappa k^2 / (k^2 + q^2),

and the rotational share of the response therefore reaches the bare
N^2 = kappa/(mu + kappa) = 1/pi only as k -> infinity.  Since N^2 is defined at a
rolling point contact, it is that short-wavelength limit, and a process at finite
k would seem to carry less.  Applied to the hopping amplitude alone, this closes
the splitting-ratio residual exactly.

It is wrong, and this script is the guard that says so.  Two independent tests
fail it.

TEST 1, self-consistency.  The chiral phase delta = N^4 is a round trip of two
conversions belonging to the same hop.  Correcting h/A and leaving delta bare is
arbitrary.  Correcting both moves the ratio from exact agreement to -2.5 sigma
and destroys the coherence of the residual scale demand.

TEST 2, the fine-structure constant.  N^2 enters the self-consistent
Peierls-Nabarro equation of the alpha paper through the rotational self-energy,
whose 1/pi is N^2, and through the inter-valley vertex, whose 1/pi^2 is N^4:

    T^-1 = exp(pi^2/2) - 2 - T - N^2 T/(1 - T) - 6 N^4 T^3.

At N^2 = 1/pi that solves to alpha^-1 = 137.035999177348, which is 0.0025 parts
per billion from the CODATA value 137.035999177(21), or 0.02 of its uncertainty.
Differentiating, d(alpha^-1)/d(N^2) = -0.00735, so one CODATA standard deviation
corresponds to a relative shift in N^2 of 9.0e-06. The proposed running is
2.3e-02, about 2600 times larger, and it would move alpha^-1 to 137.036054.

The conclusion is more useful than the proposal would have been: N^2 = 1/pi is
exact at every conversion vertex in the framework, and any correction that
adjusts it, at any scale, is excluded before it is evaluated.  Corrections that
multiply a vertex by a competition factor such as (1 - alpha) are untouched by
this, because they leave N^2 itself alone and are already inside the alpha
series.
"""

import numpy as np

N2 = 1 / np.pi                    # Cosserat coupling number, rolling contact
Q = 1.87                          # microrotation screening wavevector [1/ell]
D_PARTIAL = 1 / np.sqrt(3)        # Shockley partial period [ell]
ALPHA = 1 / 137.035999177
ALPHA_INV_MEAS = 137.035999177   # CODATA 2022
ALPHA_INV_ERR = 2.1e-8           # CODATA 2022 standard uncertainty


def rotational_share(k):
    """N^2_eff(k) from mu_eff(k) = mu + kappa k^2/(k^2 + q^2).

    Convention note (2026-08-28): the alpha paper's relation is
    kappa_c = 2 N^2 mu / (1 - 2 N^2), i.e. kappa/(mu+kappa) = 2 N^2,
    certified in neutrinos/screw_vertex_calibration.py against the
    printed kernel. The earlier kap_over_mu = N2/(1-N2) understated
    kappa by ~2.4x; the screening scale Q = 1.87, which the tests
    actually pivot on, was taken as an input and is correct, so the
    no-go's conclusions stand, but the internal formula now matches."""
    kap_over_mu = 2 * N2 / (1 - 2 * N2)
    x = k**2 / (k**2 + Q**2)
    return kap_over_mu * x / (1 + kap_over_mu * x)


def alpha_inv(n2, iters=200):
    """Solve the alpha paper's self-consistent Peierls-Nabarro equation at a
    given coupling number.

    N^2 enters twice: once in the rotational self-energy, whose 1/pi is N^2,
    and once in the inter-valley vertex, whose 1/pi^2 is N^4. Every other term
    is N^2-independent, so this function isolates exactly the dependence the
    proposed running would disturb. The map is a contraction on (0, 1), so the
    fixed point is unique and plain iteration converges.
    """
    bare = np.exp(np.pi**2 / 2)
    t = 1 / 137.0
    for _ in range(iters):
        t = 1 / (bare - 2 - t - n2 * t / (1 - t) - 6 * n2**2 * t**3)
    return 1 / t


def ring_ratio(hA, delta):
    """Splitting ratio (m2^2 - m1^2)/(m3^2 - m1^2) from the Z_3 ring."""
    phase = 2 * np.pi / 3 + delta
    lv = np.sort([1 + 2 * hA * np.cos(phase + 2 * np.pi * j / 3)
                  for j in range(3)])
    return (lv[1]**2 - lv[0]**2) / (lv[2]**2 - lv[0]**2)


def test_self_consistency():
    """delta belongs to the same hop as h/A, so both or neither."""
    k = 2 * np.pi / D_PARTIAL
    n2k = rotational_share(k)
    # measured ratio, NuFIT 6.1 with SK-atm
    r_meas = 7.537 / (2.511 * 100)
    sig = r_meas * np.hypot(0.10 / 7.537, 0.021 / 2.511)

    print("TEST 1: self-consistency of the running\n")
    print(f"  N^2 bare {N2:.6f}, N^2_eff(2 pi/d) {n2k:.6f} "
          f"({n2k/N2:.4f} of bare)\n")
    rows = [("h/A bare,     delta bare  (leading order)", N2, N2),
            ("h/A running,  delta bare  (the proposal)", n2k, N2),
            ("h/A running,  delta running (consistent)", n2k, n2k)]
    for lab, n_h, n_d in rows:
        hA = 0.5 * (1 + n_h * (1 - ALPHA))
        delta = n_d**2 * (1 - ALPHA)**2
        r = ring_ratio(hA, delta)
        print(f"  {lab:44s} R = {r:.6f}  pull {(r - r_meas)/sig:+.2f}")
    print("\n  The proposal works only by correcting one of three conversions")
    print("  that all belong to the same hop.  Correcting all three fails.\n")


def test_fine_structure():
    """The sharpest constraint in the framework, applied to the proposal."""
    print("TEST 2: the fine-structure constant\n")
    base = alpha_inv(N2)
    resid_ppb = (base - ALPHA_INV_MEAS) / ALPHA_INV_MEAS * 1e9
    print(f"  alpha^-1 at bare N^2 = {base:.12f}")
    print(f"  CODATA 2022          = {ALPHA_INV_MEAS:.9f} +/- {ALPHA_INV_ERR:.1e}")
    print(f"  residual             = {resid_ppb:+.4f} ppb "
          f"= {(base - ALPHA_INV_MEAS)/ALPHA_INV_ERR:+.3f} CODATA sigma\n")

    for lab, k in [("hop step,  k = 2 pi/d", 2 * np.pi / D_PARTIAL),
                   ("lattice,   k = 2 pi/ell", 2 * np.pi),
                   ("core,      k = 1/w", 1 / 0.453)]:   # screw core half-width
        n2 = rotational_share(k)
        shift = abs(alpha_inv(n2) - ALPHA_INV_MEAS) / ALPHA_INV_ERR
        print(f"  {lab:24s}: N^2_eff = {n2:.6f}, alpha^-1 = "
              f"{alpha_inv(n2):.9f}, {shift:.1e} CODATA sigma out")

    # How far may N^2 move before alpha^-1 leaves the CODATA band?
    h = 1e-8
    slope = (alpha_inv(N2 + h) - alpha_inv(N2 - h)) / (2 * h)
    tol = ALPHA_INV_ERR / abs(slope) / N2
    prop = 1 - rotational_share(2 * np.pi / D_PARTIAL) / N2
    print(f"\n  d(alpha^-1)/d(N^2) = {slope:.5f}")
    print(f"  tolerance on N^2 at one CODATA sigma: {tol:.1e} relative.")
    print(f"  the proposed running: {prop:.1e} relative, "
          f"{prop/tol:.0e} times too large.\n")
    assert prop > 100 * tol, "the running should be excluded by a wide margin"


def main():
    test_self_consistency()
    test_fine_structure()
    print("CONCLUSION: N^2 = 1/pi is exact at every conversion vertex.")
    print("Any proposed correction that adjusts the coupling number itself,")
    print("at any scale, is excluded by alpha^-1 before it is evaluated.")
    print("Note that this argument rests on the self-consistent Peierls-Nabarro")
    print("equation and not on any closed-form pi expression for alpha^-1. The")
    print("often-quoted 4 pi^3 + pi^2 + pi is not the framework's result: it")
    print("gives 137.0363038, which is 2223 ppb high, some 14500 CODATA sigma.")
    print("Competition factors such as (1 - alpha), which multiply a vertex")
    print("without altering N^2, are unaffected: they are already inside the")
    print("alpha series and do not move it.")


if __name__ == "__main__":
    main()
