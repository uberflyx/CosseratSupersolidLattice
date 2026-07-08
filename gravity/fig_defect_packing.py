"""
Figure: the number of ways to pack n line defects grows as exp(n^2) through
pairwise linking, while point-defect positional packing grows only as n ln n.

Panel (a): exact ln Omega versus the number of pairs n(n-1)/2 for k = 2, 3, 5
linking states per pair (Burnside orbit counts, lines unlabelled), with the
point-defect count C(N, n), N = (2n)^3, for contrast.

Panel (b): entropy per pair, [ln Omega + ln n!] / (n(n-1)/2), which the exact
finite-size formula predicts to sit at ln k.  The n! term is the relabelling
cost of indistinguishable lines and is subleading.

Output: fig_defect_packing.pdf
"""

from math import lgamma, log
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bh_defect_packing_count import count_orbits_offdiag, log_of_int

ns = list(range(4, 41, 2))
ks = (2, 3, 5)
colours = {2: "#1f77b4", 3: "#d62728", 5: "#2ca02c"}

ln_omega = {k: [log_of_int(count_orbits_offdiag(n, k)) for n in ns] for k in ks}
pairs = [n * (n - 1) / 2 for n in ns]

# Point defects: n defects on N = (2n)^3 sites
ln_point = [lgamma((2 * n) ** 3 + 1) - lgamma(n + 1) - lgamma((2 * n) ** 3 - n + 1)
            for n in ns]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.6))

for k in ks:
    ax1.plot(pairs, ln_omega[k], "o-", ms=3.5, lw=1.1, color=colours[k],
             label=rf"line defects, $k={k}$ states/pair")
ax1.plot(pairs, ln_point, "s--", ms=3.5, lw=1.1, color="0.45",
         label=r"point defects, $\binom{N}{n}$, $N=(2n)^3$")
ax1.set_xlabel(r"number of pairs $\;n(n-1)/2$")
ax1.set_ylabel(r"$\ln \Omega$")
ax1.set_title(r"(a) microstates vs pair count", fontsize=10)
ax1.legend(fontsize=7.5, frameon=False)

for k in ks:
    per_pair = [(lo + lgamma(n + 1)) / (n * (n - 1) / 2)
                for lo, n in zip(ln_omega[k], ns)]
    ax2.plot(ns, per_pair, "o-", ms=3.5, lw=1.1, color=colours[k],
             label=rf"$k={k}$")
    ax2.axhline(log(k), color=colours[k], lw=0.8, ls=":")
    ax2.text(ns[-1] + 0.6, log(k), rf"$\ln {k}$", color=colours[k],
             fontsize=8, va="center")
ax2.set_xlabel(r"number of line defects $n$")
ax2.set_ylabel(r"$[\ln\Omega + \ln n!]\,/\,[n(n-1)/2]$")
ax2.set_title(r"(b) entropy per pair $\to \ln k$ exactly", fontsize=10)
ax2.set_xlim(ns[0], ns[-1] + 5)

fig.tight_layout()
fig.savefig("fig_defect_packing.pdf")
print("wrote fig_defect_packing.pdf")
