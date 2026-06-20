"""
The nuclear asymmetry coefficient from two lattice channels.

Cosserat Supersolid Lattice framework. The asymmetry term of the semi-empirical
mass formula, a_A (N-Z)^2 / A, is the energy cost of pulling a nucleus away from
N = Z. In the lattice picture that cost has two independent parts, and this script
computes both with no parameter tuned to the answer.

  KINETIC CHANNEL.  Nucleons carry the zero-point motion of a degenerate Fermi
  gas. Surplus neutrons, barred by exclusion from the filled states, occupy higher
  momenta. For two fermion species the cost is E_F / 3, with the Fermi energy set
  by the nuclear density (the lattice inter-nucleon spacing).

  INTERACTION CHANNEL.  An unlike dock (proton-neutron) binds, like the deuteron;
  a like dock (neutron-neutron) does not, which is why there is no bound dineutron.
  When N != Z the surplus neutrons are forced into like-docks that do not hold, and
  the binding lost is the asymmetry energy. We count, on contact-maximised FCC
  clusters, the minimum number of like-docks the lattice forces (the lattice is
  geometrically frustrated, so some survive even at N = Z), and charge each forced
  like-dock the deuteron bond.

The two channels add. Result: about 18 MeV against a measured coefficient near
23 MeV, i.e. about four fifths from first principles. The shortfall is the part of
the isospin force beyond a single bond, which the like-dock count does not carry.
The interaction channel is packing-sensitive at the ~1 MeV level; the clusters here
are built by contact maximisation, the physical choice for self-bound nuclei.

Deterministic (fixed seeds) and self-contained.

Author: Mitch Cox.  Run: python asymmetry_coefficient.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from numba import njit

mpl.rcParams.update({"font.size": 10, "axes.linewidth": 0.8,
                     "figure.dpi": 140, "savefig.dpi": 140})

# Physical constants (CODATA / PDG; nucleon mass is the framework isoscalar value)
HBARC = 197.327     # MeV fm
M_N   = 938.92      # MeV, isoscalar nucleon
ELL   = 2.82        # fm, lattice constant = classical electron radius r_e
D_NUC = 1.80        # fm, inter-nucleon spacing at saturation (= 0.64 ell)
B_BOND = 2.5        # MeV, unlike-dock (deuteron-scale) binding

# The twelve FCC nearest-neighbour offsets (sites with even coordinate sum).
NN = np.array([(a, b, 0) for a in (-1, 1) for b in (-1, 1)] +
              [(a, 0, b) for a in (-1, 1) for b in (-1, 1)] +
              [(0, a, b) for a in (-1, 1) for b in (-1, 1)], dtype=np.int64)
NN_T = [tuple(d) for d in NN]


# ----------------------------------------------------------------------
# FCC packing: build a compact nucleus by maximising nearest-neighbour
# contacts (each contact is one dock, so the most-bound shape is the most
# compact). Greedy seed, then surface-swap simulated annealing.
# ----------------------------------------------------------------------
def fcc_pool(R):
    return [(x, y, z)
            for x in range(-R, R + 1)
            for y in range(-R, R + 1)
            for z in range(-R, R + 1)
            if (x + y + z) % 2 == 0]


def _contacts_in(site, chosen_set):
    return sum(1 for d in NN_T
               if (site[0] + d[0], site[1] + d[1], site[2] + d[2]) in chosen_set)


def greedy_cluster(A, pool):
    chosen = [pool[0]]
    cs = {pool[0]}
    while len(chosen) < A:
        best, best_g = None, -1
        for cand in pool:
            if cand in cs:
                continue
            g = _contacts_in(cand, cs)
            if g > best_g or (g == best_g and best is not None and
                              sum(c * c for c in cand) < sum(c * c for c in best)):
                best, best_g = cand, g
        chosen.append(best)
        cs.add(best)
    return chosen


def anneal_cluster(chosen, pool, steps=2500, T0=1.5, seed=0):
    pool_set = set(pool)
    cur = set(chosen)
    rng = np.random.default_rng(seed)

    def total_contacts(s):
        c = 0
        for p in s:
            for d in NN_T:
                if (p[0] + d[0], p[1] + d[1], p[2] + d[2]) in s:
                    c += 1
        return c // 2

    best = set(cur)
    best_c = cur_c = total_contacts(cur)
    for i in range(steps):
        T = T0 * (1 - i / steps) + 1e-3
        surf = set()
        for p in cur:
            for d in NN_T:
                q = (p[0] + d[0], p[1] + d[1], p[2] + d[2])
                if q in pool_set and q not in cur:
                    surf.add(q)
        rem = min(cur, key=lambda p: _contacts_in(p, cur))
        add = max(surf, key=lambda q: _contacts_in(q, cur)) if surf else rem
        trial = set(cur)
        trial.discard(rem)
        trial.add(add)
        if len(trial) != len(cur):
            continue
        ct = total_contacts(trial)
        if ct >= cur_c or rng.random() < np.exp((ct - cur_c) / T):
            cur, cur_c = trial, ct
            if cur_c > best_c:
                best, best_c = set(cur), cur_c
    return np.array(sorted(best), dtype=np.int64)


def packed_cluster(A):
    """Contact-maximised compact FCC cluster of A nucleons (deterministic)."""
    R = int(np.ceil(A ** (1 / 3)) + 3)
    pool = fcc_pool(R)
    pool.sort(key=lambda s: s[0] ** 2 + s[1] ** 2 + s[2] ** 2)
    return anneal_cluster(greedy_cluster(A, pool), pool)


def adjacency(coords):
    """Edge list of nearest-neighbour docks for a set of FCC sites."""
    index = {tuple(c): i for i, c in enumerate(coords)}
    edges = []
    for i, c in enumerate(coords):
        for d in NN:
            j = index.get((c[0] + d[0], c[1] + d[1], c[2] + d[2]), -1)
            if j > i:
                edges.append((i, j))
    return np.array(edges, dtype=np.int64)


# ----------------------------------------------------------------------
# Like-dock minimisation at fixed species count (fixed-magnetisation
# max-cut), by simulated annealing.
# ----------------------------------------------------------------------
@njit(cache=True)
def count_like(spin, edges):
    c = 0
    for e in range(edges.shape[0]):
        if spin[edges[e, 0]] == spin[edges[e, 1]]:
            c += 1
    return c


@njit(cache=True)
def min_like_docks(edges, A, n_up, iters, seed):
    np.random.seed(seed)
    spin = np.zeros(A, dtype=np.int64)
    order = np.arange(A)
    np.random.shuffle(order)
    for k in range(n_up):
        spin[order[k]] = 1
    best = count_like(spin, edges)
    cur = best
    for it in range(iters):
        T = 1.5 * (1.0 - it / iters) + 0.01
        ups = np.where(spin == 1)[0]
        dns = np.where(spin == 0)[0]
        i = ups[np.random.randint(len(ups))]
        j = dns[np.random.randint(len(dns))]
        spin[i] = 0
        spin[j] = 1
        c = count_like(spin, edges)
        d = c - cur
        if d <= 0 or np.random.random() < np.exp(-d / T):
            cur = c
            if cur < best:
                best = cur
        else:
            spin[i] = 1
            spin[j] = 0
    return best


def asymmetry_curve(coords, iters=8000, restarts=6):
    """Return ((N-Z), interaction asymmetry energy) and the fitted a_A for one nucleus."""
    edges = adjacency(coords)
    A = len(coords)
    base = min(min_like_docks(edges, A, A // 2, iters, s) for s in range(1, restarts + 1))
    x, y = [], []
    for dNZ in range(0, A // 2 + 1, 2):
        n_up = (A + dNZ) // 2
        lb = min(min_like_docks(edges, A, n_up, iters, s) for s in range(1, restarts + 1))
        x.append(dNZ)
        y.append((lb - base) * B_BOND)
    x = np.array(x, float)
    y = np.array(y, float)
    slope = np.sum(x ** 2 * y) / np.sum(x ** 4)      # = a_A / A
    return x, y, slope * A


def kinetic_coefficient(d=D_NUC):
    """Kinetic channel a_A = E_F/3, with E_F set by the lattice density."""
    n = d ** -3
    kF = (3 * np.pi ** 2 * n / 2) ** (1 / 3)          # two species share n
    EF = (HBARC * kF) ** 2 / (2 * M_N)
    return EF / 3, EF, kF, n


def main():
    A_list = [16, 24, 28, 40, 56]
    per_nucleus = {}
    for A in A_list:
        _, _, aA = asymmetry_curve(packed_cluster(A))
        per_nucleus[A] = aA
    a_int = np.mean(list(per_nucleus.values()))
    a_kin, EF, kF, n = kinetic_coefficient()

    print("Interaction channel (like-dock penalty, per nucleus):")
    for A, v in per_nucleus.items():
        print(f"  A = {A:>3}:  a_A = {v:5.1f} MeV")
    print(f"  mean a_A(interaction) = {a_int:.1f} MeV  "
          f"(penalty per forced like-dock = {B_BOND} MeV)")
    print(f"\nKinetic channel:  d = {D_NUC} fm -> n = {n:.3f}/fm^3, "
          f"k_F = {kF:.3f}/fm, E_F = {EF:.1f} MeV")
    print(f"  a_A(kinetic) = E_F/3 = {a_kin:.1f} MeV")
    print(f"\nTotal a_A = {a_kin:.1f} + {a_int:.1f} = {a_kin + a_int:.1f} MeV "
          f"(measured ~23 MeV)")

    # Figure: like-dock parabola for Ca-40, and the two-channel budget.
    xs, ys, a_int_fig = asymmetry_curve(packed_cluster(40))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.2))
    ax1.scatter(xs ** 2, ys, s=42, color="#1b4965", zorder=3,
                label="lattice (min like-docks)")
    xx = np.linspace(0, xs.max() ** 2, 100)
    ax1.plot(xx, (a_int_fig / 40) * xx, color="#ee6c4d", lw=1.8,
             label=(r"$E_{asym}=\frac{a_A^{int}}{A}(N-Z)^2$, "
                    f"$a_A^{{int}}$={a_int_fig:.1f} MeV"))
    ax1.set_xlabel(r"$(N-Z)^2$")
    ax1.set_ylabel("interaction asymmetry energy  (MeV)")
    ax1.set_title(r"Forced like-docks in $^{40}$Ca")
    ax1.legend(fontsize=8.5, frameon=False, loc="upper left")
    ax1.grid(alpha=0.25, lw=0.5)

    ax2.bar([0], [a_kin], width=0.55, color="#5fa8d3",
            label=f"kinetic  $E_F/3$ = {a_kin:.1f}")
    ax2.bar([0], [a_int], bottom=[a_kin], width=0.55, color="#ee6c4d",
            label=f"interaction = {a_int:.1f}")
    ax2.bar([1], [23], width=0.55, color="#cccccc",
            label=r"measured $a_A\approx23$")
    ax2.text(0, a_kin + a_int + 0.6, f"{a_kin + a_int:.1f}",
             ha="center", fontweight="bold")
    ax2.text(1, 23.6, "23", ha="center", fontweight="bold")
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(["framework", "experiment"])
    ax2.set_ylabel(r"asymmetry coefficient  $a_A$  (MeV)")
    ax2.set_title("Two derived channels vs data")
    ax2.legend(fontsize=8.5, frameon=False, loc="upper center")
    ax2.set_ylim(0, 28)
    ax2.grid(alpha=0.25, lw=0.5, axis="y")
    fig.tight_layout()
    fig.savefig("asymmetry_coefficient.pdf", bbox_inches="tight")
    print("\nsaved asymmetry_coefficient.pdf")


if __name__ == "__main__":
    main()
