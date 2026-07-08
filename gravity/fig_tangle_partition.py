"""
Figure: the exact tangle partition function.

Panel (a): low-temperature expansion coefficients dim_g of Z_N(x) for
several N, against the stabilised (N-independent) multigraph numbers:
the curves coincide exactly up to g = N/2 and peel off beyond, so the
low-temperature thermodynamics of the tangle is universal, independent
of the hole, to order x^(N/2).

Panel (b): per-pair entropy through the Hagedorn skin for a stellar
hole (2 ln N = 178, k = 2), with x(r) = exp(-2 ln N (xi/xi_*)^3): the
entropy switches on between the outer edge xi_* (first thermal winding)
and the inner edge xi_* (2 ln N)^(-1/3) (full mixing).

Output: fig_tangle_partition.pdf
"""

from math import log, exp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tangle_matrix_model_matching import molien_series
from tangle_partition_exact import per_pair_entropy

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.0))

# ---- (a) universality of the low-T coefficients -----------------------
G = 12
ref = molien_series(30, G)          # stabilised for all g <= 15
colors = plt.cm.viridis(np.linspace(0.15, 0.8, 4))
for N, col in zip((6, 8, 10, 12), colors):
    dims = molien_series(N, G)
    ax1.semilogy(range(G + 1), dims, "o-", ms=4, lw=1.2, color=col,
                 label=f"$N={N}$")
    gstar = N // 2
    ax1.axvline(gstar, color=col, lw=0.7, ls=":", alpha=0.6)
ax1.semilogy(range(G + 1), ref, "k--", lw=1.6,
             label=r"stabilised ($N\to\infty$)")
ax1.set_xlabel(r"degree $g$ (winding quanta)")
ax1.set_ylabel(r"$\dim_g$ (independent invariants)")
ax1.set_title(r"(a) low-$T$ coefficients are universal to $g = N/2$")
ax1.legend(fontsize=8, frameon=False)

# ---- (b) entropy through the skin --------------------------------------
lnN2 = 178.0                        # 2 ln N for a solar-mass hole, k = 2
xi = np.logspace(np.log10(3.0), np.log10(0.03), 400)   # xi / xi_*
s = np.array([per_pair_entropy(exp(-lnN2 * r ** 3), 2) / log(2) for r in xi])
ax2.semilogx(xi, s, "b-", lw=2)
inner = (1.0 / lnN2) ** (1.0 / 3.0)
ax2.axvline(1.0, color="0.4", ls="--", lw=1)
ax2.axvline(inner, color="0.4", ls="--", lw=1)
ax2.annotate("outer edge $\\xi_*$\n(first winding)\n$\\rho \\approx 1.4$ cm",
             xy=(1.0, 0.06), xytext=(1.6, 0.28), fontsize=8,
             arrowprops=dict(arrowstyle="->", lw=0.8))
ax2.annotate("inner edge $\\xi_*(2\\ln N)^{-1/3}$\n(full mixing)\n"
             "$\\rho \\approx 2.5$ mm",
             xy=(inner, 0.84), xytext=(0.045, 0.50), fontsize=8,
             arrowprops=dict(arrowstyle="->", lw=0.8))
ax2.set_xlabel(r"depth $\xi/\xi_*$ (redshift factor, horizon at $0$)")
ax2.set_ylabel(r"entropy per pair $s / (k_B \ln k)$")
ax2.set_title(r"(b) the skin band, $1\,M_\odot$: $x = e^{-2\ln N(\xi/\xi_*)^3}$")
ax2.set_ylim(-0.03, 1.05)
ax2.invert_xaxis()                  # falling toward the horizon, left to right

fig.tight_layout()
fig.savefig("fig_tangle_partition.pdf")
print("wrote fig_tangle_partition.pdf")
