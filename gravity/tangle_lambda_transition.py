"""
The lambda point of the defect gas: order parameters of the tangle's
permutation ensemble across the deconfinement transition.

Physics context (monograph, black-holes chapter; paper "Black holes as
superfluid droplets"): writing line indistinguishability into the measure
of the master action's defect path integral derives the gauged partition
function

    Z_N(x;k) = sum_{lambda |- N} (1/z_lambda)
               prod_{pair cycles c} (1 - x^{k l_c}) / (1 - x^{l_c})

from the action itself (Feynman's identical-particle construction, the
lambda-transition bookkeeping of helium, with defect lines in place of
atoms).  The class ensemble then carries two order parameters:

  rho_1     = <m_1>/N : the fixed-point fraction, the S_N Polyakov loop
                        of the Landau functional (identity condensate,
                        the black-hole/normal phase);
  f_giant   = <sum_{l > N/2} l m_l>/N : the fraction of lines living in
                        giant permutation cycles (Feynman's criterion
                        for superfluidity, the exterior record).

At x -> 0 every permutation weighs equally, E[m_l] = 1/l, and giant
cycles dominate (f_giant -> ~1/2 for the >N/2 measure, rho_1 -> 1/N):
the cycle-condensed, superfluid phase.  At x -> 1 the identity's k^d
weight crushes everything (rho_1 -> 1, f_giant -> 0): the normal phase.
The exchange sharpens with N around x_c ~ 2 ln N / ((k-1) N), the same
crossing located by the free-energy saddle analysis: the lambda point,
the Gross-Witten-Wadia condensation, and the Hawking-Page transition
are one transition under three names.

Checks and outputs:
  1. The log-space class-ensemble Z agrees with the closed form.
  2. The x -> 0 limits match the uniform-measure expectations exactly.
  3. Order-parameter curves for N = 10, 20, 40 (k = 2), with the
     numerical deconfinement points marked: fig_tangle_lambda.pdf.
"""

from math import log, lgamma, exp
import numpy as np

from tangle_matrix_model_matching import partitions, pair_cycle_lengths
from tangle_partition_exact import Z_closed
from tangle_saddle_analysis import x_deconfine


def build_class_data(N):
    """Precompute, per conjugacy class: ln(1/z_lambda), the pair-cycle
    length spectrum (as arrays of distinct lengths and counts), and the
    two order-parameter values."""
    classes = []
    for p in partitions(N):
        # ln z_lambda = sum_l [m_l ln l + ln(m_l!)]
        lnz = sum(m * log(l) + lgamma(m + 1) for l, m in p.items())
        # pair-cycle spectrum, compressed to (lengths, counts)
        spec = {}
        for l in pair_cycle_lengths(p):
            spec[l] = spec.get(l, 0) + 1
        lengths = np.array(sorted(spec), dtype=np.int64)
        counts = np.array([spec[l] for l in sorted(spec)], dtype=np.float64)
        rho1 = p.get(1, 0) / N
        fgiant = sum(l * m for l, m in p.items() if l > N / 2) / N
        classes.append((lnz, lengths, counts, rho1, fgiant))
    return classes


def ensemble(N, k, xs, classes=None):
    """Return (Z, <rho_1>, <f_giant>) arrays over the grid xs, computed
    in log space over the exact class ensemble."""
    if classes is None:
        classes = build_class_data(N)
    maxlen = max(int(c[1].max()) for c in classes)
    Zs, r1s, fgs = [], [], []
    for x in xs:
        # ln[(1 - x^{k l})/(1 - x^l)] for every length l, vectorised
        ls = np.arange(1, maxlen + 1, dtype=np.float64)
        if x <= 0:
            lnf = np.zeros(maxlen + 1)
        elif x >= 1:
            lnf = np.full(maxlen + 1, log(k))
        else:
            lnf = np.empty(maxlen + 1)
            lnf[1:] = np.log1p(-x ** (k * ls)) - np.log1p(-x ** ls)
        lnw = np.array([-lnz + float(counts @ lnf[lengths])
                        for lnz, lengths, counts, _, _ in classes])
        m = lnw.max()
        w = np.exp(lnw - m)
        Z = w.sum()
        r1 = sum(wi * c[3] for wi, c in zip(w, classes)) / Z
        fg = sum(wi * c[4] for wi, c in zip(w, classes)) / Z
        Zs.append(Z * exp(m))
        r1s.append(r1)
        fgs.append(fg)
    return np.array(Zs), np.array(r1s), np.array(fgs)


if __name__ == "__main__":
    print("1. Log-space class ensemble vs closed form, N=10, k=2:")
    cls10 = build_class_data(10)
    for x in (0.1, 0.4, 0.9):
        Z, _, _ = ensemble(10, 2, [x], cls10)
        Zc = Z_closed(10, 2, x)
        err = abs(Z[0] - Zc) / Zc
        print(f"   x={x}: rel. err. {err:.2e}  {'PASS' if err < 1e-10 else 'FAIL'}")
        assert err < 1e-10

    print("\n2. Uniform-measure (x=0) limits, N=10:")
    _, r1, fg = ensemble(10, 2, [0.0], cls10)
    print(f"   <rho_1> = {r1[0]:.6f}  (expect 1/N = {1/10:.6f})")
    print(f"   <f_giant> = {fg[0]:.6f}  (expect (N - floor(N/2))/N = 0.5)")
    assert abs(r1[0] - 0.1) < 1e-12 and abs(fg[0] - 0.5) < 1e-12

    print("\n3. Order parameters across the transition (k=2):")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.0), sharex=True)
    xs = np.linspace(0.0, 0.7, 71)
    colors = plt.cm.viridis(np.linspace(0.15, 0.75, 3))
    for N, col in zip((10, 20, 40), colors):
        cls = build_class_data(N)
        _, r1, fg = ensemble(N, 2, xs, cls)
        xc = x_deconfine(N, 2)
        ax1.plot(xs, r1, "-", color=col, lw=1.8, label=f"$N={N}$")
        ax2.plot(xs, fg, "-", color=col, lw=1.8, label=f"$N={N}$")
        for ax in (ax1, ax2):
            ax.axvline(xc, color=col, ls=":", lw=1)
        print(f"   N={N}: x_c = {xc:.4f}, rho_1(x_c) = "
              f"{ensemble(N, 2, [xc], cls)[1][0]:.3f}")
    ax1.set_xlabel(r"Boltzmann weight per winding $x$")
    ax1.set_ylabel(r"fixed-point fraction $\langle\rho_1\rangle$")
    ax1.set_title(r"(a) identity condensate (normal phase, interior)")
    ax1.legend(fontsize=8, frameon=False)
    ax2.set_xlabel(r"Boltzmann weight per winding $x$")
    ax2.set_ylabel(r"giant-cycle fraction $\langle f_{>N/2}\rangle$")
    ax2.set_title(r"(b) cycle condensate (superfluid record, exterior)")
    ax2.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig("fig_tangle_lambda.pdf")
    print("   wrote fig_tangle_lambda.pdf (dotted verticals: x_c(N))")
