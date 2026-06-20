"""
The surface-to-volume ratio of the mass formula from FCC geometry.

Cosserat Supersolid Lattice framework. The semi-empirical mass formula has a
volume term a_V A and a surface term -a_S A^(2/3). Both are built from the same
nearest-neighbour docks: an interior nucleon makes six docks (twelve FCC
neighbours, shared), a surface nucleon makes fewer. Writing the dock count of a
compact cluster as

    bonds(A) = 6 A  -  c A^(2/3),

the bulk gives a_V = 6 B_bond and the boundary deficit gives a_S = c B_bond, so

    a_S / a_V = c / 6.

The energy scale B_bond cancels: the surface-to-volume ratio is pure geometry. It
needs no knowledge of the saturation density that sets the absolute coefficients.
This script builds contact-maximised FCC clusters, counts their docks, fits c, and
compares c/6 to the measured a_S/a_V = 17.8/15.75 = 1.13.

Result: c = 7.6, so a_S/a_V = 1.27, about twelve percent above measured. A sharp
lattice surface sheds slightly more docks than the real diffuse nuclear surface;
that diffuseness is the twelve-percent correction.

Deterministic and self-contained.

Author: Mitch Cox.  Run: python surface_coefficient.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({"font.size": 10, "axes.linewidth": 0.8,
                     "figure.dpi": 140, "savefig.dpi": 140})

A_V_MEAS = 15.75    # MeV, measured volume coefficient
A_S_MEAS = 17.8     # MeV, measured surface coefficient

# Twelve FCC nearest-neighbour offsets (sites with even coordinate sum).
NN = [(a, b, 0) for a in (-1, 1) for b in (-1, 1)] + \
     [(a, 0, b) for a in (-1, 1) for b in (-1, 1)] + \
     [(0, a, b) for a in (-1, 1) for b in (-1, 1)]


def fcc_pool(R):
    return [(x, y, z)
            for x in range(-R, R + 1)
            for y in range(-R, R + 1)
            for z in range(-R, R + 1)
            if (x + y + z) % 2 == 0]


def contacts_in(site, chosen):
    return sum(1 for d in NN
               if (site[0] + d[0], site[1] + d[1], site[2] + d[2]) in chosen)


def greedy_cluster(A, pool):
    chosen = {pool[0]}
    while len(chosen) < A:
        best, best_g = None, -1
        for cand in pool:
            if cand in chosen:
                continue
            g = contacts_in(cand, chosen)
            if g > best_g or (g == best_g and best is not None and
                              sum(c * c for c in cand) < sum(c * c for c in best)):
                best, best_g = cand, g
        chosen.add(best)
    return chosen


def total_bonds(chosen):
    return sum(contacts_in(p, chosen) for p in chosen) // 2


def anneal_cluster(chosen, pool, steps=2000, T0=1.5, seed=0):
    pool_set = set(pool)
    cur = set(chosen)
    rng = np.random.default_rng(seed)
    best = set(cur)
    best_b = cur_b = total_bonds(cur)
    for i in range(steps):
        T = T0 * (1 - i / steps) + 1e-3
        surf = set()
        for p in cur:
            for d in NN:
                q = (p[0] + d[0], p[1] + d[1], p[2] + d[2])
                if q in pool_set and q not in cur:
                    surf.add(q)
        rem = min(cur, key=lambda p: contacts_in(p, cur))
        add = max(surf, key=lambda q: contacts_in(q, cur)) if surf else rem
        trial = set(cur)
        trial.discard(rem)
        trial.add(add)
        if len(trial) != len(cur):
            continue
        bt = total_bonds(trial)
        if bt >= cur_b or rng.random() < np.exp((bt - cur_b) / T):
            cur, cur_b = trial, bt
            if cur_b > best_b:
                best, best_b = set(cur), cur_b
    return best_b


def bonds_of(A):
    """Dock count of a contact-maximised compact FCC cluster of A nucleons."""
    R = int(np.ceil(A ** (1 / 3)) + 3)
    pool = fcc_pool(R)
    pool.sort(key=lambda s: s[0] ** 2 + s[1] ** 2 + s[2] ** 2)
    return anneal_cluster(greedy_cluster(A, pool), pool)


def main():
    A_list = np.array([16, 28, 40, 56, 68, 80, 100, 120, 140], float)
    bonds = np.array([bonds_of(int(A)) for A in A_list], float)
    deficit = 6 * A_list - bonds                       # docks lost at the surface
    x = A_list ** (2 / 3)
    c = np.sum(x * deficit) / np.sum(x * x)            # bonds = 6A - c A^(2/3)

    ratio_fcc = c / 6
    ratio_meas = A_S_MEAS / A_V_MEAS
    c_meas = 6 * ratio_meas
    B_bond = A_V_MEAS / 6

    print(f"{'A':>5} {'bonds':>7} {'6A':>6} {'deficit':>9} {'deficit/A^2/3':>14}")
    for A, b, d in zip(A_list, bonds, deficit):
        print(f"{A:>5.0f} {b:>7.0f} {6*A:>6.0f} {d:>9.1f} {d/A**(2/3):>14.2f}")
    print(f"\nfit  bonds(A) = 6A - c A^(2/3),  c = {c:.2f}")
    print(f"surface-to-volume ratio  a_S/a_V = c/6 = {ratio_fcc:.3f}")
    print(f"measured                 a_S/a_V       = {ratio_meas:.3f}  "
          f"(off by {100*(ratio_fcc/ratio_meas-1):+.0f}%)")
    print(f"\nwith B_bond = a_V/6 = {B_bond:.3f} MeV:  "
          f"a_S = c B_bond = {c*B_bond:.1f} MeV (measured {A_S_MEAS})")

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.scatter(x, deficit, s=46, color="#1b4965", zorder=3,
               label="FCC clusters (docks lost at surface)")
    xx = np.linspace(0, x.max() * 1.05, 100)
    ax.plot(xx, c * xx, color="#ee6c4d", lw=1.9,
            label=f"FCC geometry: $c$={c:.2f}, $a_S/a_V$={ratio_fcc:.2f}")
    ax.plot(xx, c_meas * xx, color="#777777", lw=1.6, ls="--",
            label=f"measured: $c$={c_meas:.2f}, $a_S/a_V$={ratio_meas:.2f}")
    ax.set_xlabel(r"$A^{2/3}$  (surface area)")
    ax.set_ylabel("surface dock deficit  $6A - \\mathrm{bonds}(A)$")
    ax.set_title("Surface-to-volume ratio from FCC packing")
    ax.legend(fontsize=8.8, frameon=False, loc="upper left")
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig("surface_coefficient.pdf", bbox_inches="tight")
    print("\nsaved surface_coefficient.pdf")


if __name__ == "__main__":
    main()
