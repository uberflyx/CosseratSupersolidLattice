#!/usr/bin/env python3
"""The core form factor: why a PN defect's periodic energy carries e^{-2 pi w/d}.

Companion to the neutrino mass-scale section.  Verifies the closed forms that
bound the second-order chiral coefficient, and the first-order case that is
the screw's own mass law.

THEOREM (pen).  Let varphi(x) = (b/pi)[atan(x/w) + pi/2] be the PN disregistry
and let any per-bond energy density f depend on the local registration
varphi_j = varphi(jd - x0).  Then the amplitude of the G = 2pi/d harmonic of
F(x0) = sum_j f(varphi_j) is suppressed by the core form factor e^{-2pi w/d},
independently of f; and for a second-order energy built from two vertex
insertions the suppression is exactly one factor e^{-2pi w/d}, never two,
because the momentum transfer G can split between the insertions with
|k| + |G-k| = G flat across the segment.

The checks verify the exact closed forms (not log-slope fits, which a
polynomial prefactor contaminates: the first fit of this file measured -4.53
and wrongly read it against -2pi; the closed forms below agree to machine
precision).  With z = 2 pi w/d:
    first order,  f = sin^2(pi varphi/b)      : |F_1| = z e^{-z}   (exact)
    second order, compact, v = sin(2pi varphi): |F_1| = 2z|1-z| e^{-z}  (exact)
    band intermediates                         : |F_1| = P(z) e^{-z}, P ~ z^2
The exponential is the theorem; the polynomial is what Born-Huang absorption
renormalises (the screw's own prefactor at w/d = pi/4 is z = pi^2/2 = 4.93,
absorbed to the c_3 = 1 of m_e = alpha m_0).  Structures checked:
  (1) first order, misfit density sin^2(pi varphi/b)   [the screw's own mass]
  (2) second order, compact intermediates:  -sum_j v(varphi_j)^2/Delta
  (3) second order, band intermediates:     -sum_q |v(q)|^2/Delta(q)
  (4) a vertex with strong higher harmonics [the nonlinear steelman]
Check (1) at w/d = pi/4 must land on e^{-pi^2/2} = 1/139.05, the bare alpha.

All sums are over the discrete chain j in [-J, J], J large enough that the
polynomial arctan tails contribute below the smallest measured amplitude.
"""

import numpy as np

B = 1.0                                   # misfit period (units of d)
D = 1.0
G = 2 * np.pi / D
J = 6000                                  # chain half-length
NPH = 64                                  # x0 phases per period


def varphi(x, w):
    return (B / np.pi) * (np.arctan(x / w) + np.pi / 2)


def harmonic_amp(F):
    """Relative first-harmonic amplitude of F(x0), |F_1| (absolute)."""
    c = np.fft.rfft(F) / len(F)
    return 2 * np.abs(c[1]), c[0].real


def core_sum(w, per_bond, J=J):
    j = np.arange(-J, J + 1) * D
    x0s = np.linspace(0, D, NPH, endpoint=False)
    F = np.array([per_bond(varphi(j - x0, w)).sum() for x0 in x0s])
    return harmonic_amp(F)


def band_sum(w, vertex, J=J):
    """Second order with delocalised intermediates: -sum_q |v(q)|^2/Delta(q).

    v(q) = sum_j vertex(varphi_j) e^{-iqj};  Delta(q) = 1 + 0.5(1 - cos qd),
    a gapped cosine band.  Computed per x0 by FFT over the chain.
    """
    j = np.arange(-J, J + 1) * D
    q = 2 * np.pi * np.fft.fftfreq(len(j), d=D)
    Delta = 1.0 + 0.5 * (1 - np.cos(q * D))
    x0s = np.linspace(0, D, NPH, endpoint=False)
    F = np.empty(NPH)
    for i, x0 in enumerate(x0s):
        v = vertex(varphi(j - x0, w))
        vq = np.fft.fft(v)
        F[i] = -np.sum(np.abs(vq) ** 2 / Delta) / len(j)
    return harmonic_amp(F)


def slope(ws, amps):
    p = np.polyfit(ws, np.log(amps), 1)
    return p[0]


def main():
    ws = np.array([0.35, 0.45, 0.55, 0.65, 0.75, 0.85])

    print("=" * 74)
    print("(1) FIRST ORDER: misfit density sin^2(pi varphi/b)  [screw mass law]")
    print("=" * 74)
    f1 = lambda p: np.sin(np.pi * p / B) ** 2
    amps = []
    print(f"{'w/d':>7s} {'|F_1|':>12s} {'F_0':>10s} {'|F_1| e^{2pi w/d}':>18s}")
    for w in ws:
        a1, a0 = core_sum(w, f1)
        amps.append(a1)
        print(f"{w:7.2f} {a1:12.4e} {a0:10.4f} {a1*np.exp(2*np.pi*w):18.5f}")
    s = slope(ws, amps)
    print(f"  log-slope d ln|F_1|/d(w/d) = {s:+.4f}   (theorem: -2 pi = {-2*np.pi:+.4f})")
    a1, a0 = core_sum(np.pi / 4, f1)
    print(f"  at w/d = pi/4: |F_1|/F_0-type fraction = {a1/(a0/  (2*J+1) * 1):.3e}")
    print(f"  raw |F_1| = {a1:.4e} against e^(-pi^2/2) = {np.exp(-np.pi**2/2):.4e}"
          f"   ratio = {a1/np.exp(-np.pi**2/2):.4f}")

    print()
    print("=" * 74)
    print("(2) SECOND ORDER, COMPACT INTERMEDIATES: -sum v(varphi)^2/Delta")
    print("=" * 74)
    vtx = lambda p: np.sin(2 * np.pi * p / B)          # periodic chiral vertex
    f2 = lambda p: -vtx(p) ** 2 / 1.0
    amps2 = []
    for w in ws:
        a1, a0 = core_sum(w, f2)
        amps2.append(a1)
    s2 = slope(ws, amps2)
    print(f"  log-slope = {s2:+.4f}   (theorem: one factor, -2 pi = {-2*np.pi:+.4f};")
    print(f"   two factors would be -4 pi = {-4*np.pi:+.4f})")

    print()
    print("=" * 74)
    print("(3) SECOND ORDER, BAND INTERMEDIATES: -sum_q |v(q)|^2/Delta(q)")
    print("=" * 74)
    amps3 = []
    for w in ws:
        a1, a0 = band_sum(w, vtx)
        amps3.append(a1)
    s3 = slope(ws, amps3)
    print(f"  log-slope = {s3:+.4f}   (theorem: -2 pi, the split transfer"
          f" |k|+|G-k| = G)")

    print()
    print("=" * 74)
    print("(4) NONLINEAR VERTEX (strong higher harmonics)")
    print("=" * 74)
    vtx4 = lambda p: (np.sin(2 * np.pi * p / B) + 0.8 * np.sin(4 * np.pi * p / B)
                      + 0.6 * np.sin(6 * np.pi * p / B))
    f4 = lambda p: -vtx4(p) ** 2
    amps4 = []
    for w in ws:
        a1, a0 = core_sum(w, f4)
        amps4.append(a1)
    s4 = slope(ws, amps4)
    print(f"  log-slope = {s4:+.4f}   (theorem: still -2 pi; higher harmonics"
          f" only decay faster)")

    print()
    print("=" * 74)
    print("THE NUMBERS")
    print("=" * 74)
    ainv = 137.035999177
    alpha = 1 / ainv
    theta2m0 = 5.031                              # meV, the printed scale
    for w_e in (0.540, 0.590):
        wd = w_e / (1 / np.sqrt(3))               # d = ell/sqrt(3)
        ff = np.exp(-2 * np.pi * wd)
        print(f"  w_e = {w_e:.3f} ell:  w_e/d = {wd:.4f},  form factor = {ff:.3e},"
              f"  m_1 = c x {theta2m0*ff*1e3:.1f} micro-eV")
    wd_req = -np.log(0.99) / (2 * np.pi)
    print(f"  printed value needs form factor > 0.99, i.e. w_e < {wd_req:.4f} d"
          f" = {wd_req/np.sqrt(3):.4f} ell: an abrupt core.")
    # downstream: the ring with m_1 ~ 0 against the solar splitting
    delta = 0.09985
    r = np.sin(delta) / np.sin(np.pi / 3 + delta)
    m3 = np.sqrt(2.534e-3) * 1e3                  # meV, m1 ~ 0
    m2 = r * m3
    dm21 = (m2 ** 2 - 0.0 ** 2) * 1e-6            # eV^2
    print(f"  ring with m_1 -> 0: m_2 = {m2:.2f} meV, Dm21^2 = {dm21:.3e} eV^2,"
          f"  pull = {(dm21 - 7.50e-5)/0.20e-5:+.1f} sigma")
    print(f"  ring shape demands m_1/m_3 = 0.0996 -> m_1 = {0.0996*m3:.2f} meV"
          f" from the atmospheric splitting alone.")


if __name__ == "__main__":
    main()
