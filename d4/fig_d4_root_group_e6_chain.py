"""
Figure: the D4 root group, its McKay dual E6, and the shared Q8 of the two isospins.

The schematic shows the chain that the section establishes.  The twenty-four roots of
D4, read as unit quaternions, form the binary tetrahedral group SL(2,3) = Q8 |x Z3.
Its order-eight core Q8 supplies both isospins: strong isospin as the McKay dual of the
affine D4 star (the doublet at the centre), weak SU(2)_L as the left action of the roots
on the self-dual microrotation.  The McKay dual of the full root group is the affine E6
diagram.  Numbers and group structure are produced and checked by d4_root_group_e6.py.

Output: d4_root_group_e6_chain.pdf
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

# muted, print-friendly palette
INK = "#1b1b1b"
BOX = "#f4f1ea"        # warm paper
EDGE = "#3a3a3a"
ACCENT = "#7a1f2b"     # deep red for the Q8 highlight
BLUE = "#2f4b7c"       # strong isospin
GREEN = "#2e6f4e"      # weak isospin
GOLD = "#9a6a00"       # E6


def affine_node(ax, xy, r=0.085, fill="white", lw=1.4, edge=EDGE, doubled=False):
    ax.add_patch(Circle(xy, r, facecolor=fill, edgecolor=edge, lw=lw, zorder=5))
    if doubled:
        ax.add_patch(Circle(xy, r * 0.62, facecolor="none", edgecolor=edge, lw=lw, zorder=6))


def edge(ax, a, b, lw=1.4, color=EDGE):
    ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=lw, zorder=4, solid_capstyle="round")


def draw_D4hat(ax, cx, cy, s=0.42, color=EDGE):
    """Affine D4: a central node with four legs (the McKay star)."""
    centre = (cx, cy)
    legs = [(cx, cy + s), (cx, cy - s), (cx - s, cy), (cx + s, cy)]
    for L in legs:
        edge(ax, centre, L, color=color)
    affine_node(ax, centre, r=0.10, fill=BLUE, edge=color)        # the doublet
    affine_node(ax, legs[0], doubled=True, edge=color)            # affine node (trivial singlet)
    for L in legs[1:]:
        affine_node(ax, L, edge=color)


def draw_E6hat(ax, cx, cy, s=0.34, color=GOLD):
    """Affine E6: central node, three arms each of length two (seven nodes)."""
    centre = (cx, cy)
    arms = [90.0, 210.0, 330.0]
    tips = []
    for a in arms:
        th = np.radians(a)
        n1 = (cx + s * np.cos(th), cy + s * np.sin(th))
        n2 = (cx + 2 * s * np.cos(th), cy + 2 * s * np.sin(th))
        edge(ax, centre, n1, color=color)
        edge(ax, n1, n2, color=color)
        tips.append((n1, n2))
    affine_node(ax, centre, r=0.10, fill=GOLD, edge=color)
    # the affine node sits at the tip of the top arm (Kac mark 1)
    affine_node(ax, tips[0][1], doubled=True, edge=color)
    for (n1, n2) in tips:
        affine_node(ax, n1, edge=color)
    for i, (n1, n2) in enumerate(tips):
        if i != 0:
            affine_node(ax, n2, edge=color)


def box(ax, cx, cy, w, h, lines, fc=BOX, ec=EDGE, fs=10.5, weight=None, tc=INK):
    p = FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                       boxstyle="round,pad=0.02,rounding_size=0.06",
                       facecolor=fc, edgecolor=ec, lw=1.5, zorder=3)
    ax.add_patch(p)
    ax.text(cx, cy, lines, ha="center", va="center", fontsize=fs,
            color=tc, zorder=6, weight=weight, linespacing=1.35)


def arrow(ax, a, b, label="", color=EDGE, fs=9.5, rad=0.0, off=(0, 0), lw=1.6):
    ar = FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=14,
                         connectionstyle=f"arc3,rad={rad}", color=color, lw=lw, zorder=4)
    ax.add_patch(ar)
    if label:
        mx, my = (a[0] + b[0]) / 2 + off[0], (a[1] + b[1]) / 2 + off[1]
        ax.text(mx, my, label, ha="center", va="center", fontsize=fs,
                color=color, style="italic", zorder=7,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9))


def main():
    fig, ax = plt.subplots(figsize=(11.0, 6.4))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 6.4)
    ax.axis("off")

    # --- top spine: the group chain -------------------------------------
    box(ax, 1.7, 5.4, 2.6, 1.0,
        "24 roots of $D_4$\n(minimal vectors)", fs=11)
    box(ax, 5.5, 5.4, 3.0, 1.15,
        "$\\mathrm{SL}(2,3) = Q_8 \\rtimes \\mathbb{Z}_3$\nthe root group, $|G|=24$",
        fs=11, ec=ACCENT)

    arrow(ax, (3.0, 5.4), (4.0, 5.4),
          label="as unit\nquaternions", off=(0, 0.55), fs=9)

    # Q8 highlighted inside / below the root group
    box(ax, 5.5, 3.75, 2.3, 0.78,
        "$Q_8 = \\{\\pm1,\\pm\\mathbf{i},\\pm\\mathbf{j},\\pm\\mathbf{k}\\}$",
        fs=10.5, ec=ACCENT, fc="#f7ecec", tc=ACCENT)
    arrow(ax, (5.5, 4.82), (5.5, 4.16),
          label="order-8 core\n($\\mathbb{Z}_3$ = colour quotient)", off=(2.05, 0.0), fs=8.5, color=ACCENT)

    # --- left branch: strong isospin (McKay dual D4-hat) ----------------
    draw_D4hat(ax, 1.85, 2.05, s=0.42, color=BLUE)
    ax.text(1.85, 0.95, "strong isospin\nMcKay dual $\\widehat{D}_4$, doublet $\\mathbf{2}$",
            ha="center", va="center", fontsize=9.5, color=BLUE, weight="bold", linespacing=1.3)
    arrow(ax, (4.45, 3.55), (2.45, 2.45),
          label="McKay", color=BLUE, rad=0.12, fs=9, off=(0.15, 0.25))

    # --- right branch: weak SU(2)_L (left action) -----------------------
    box(ax, 8.7, 2.05, 3.0, 1.2,
        "weak $\\mathrm{SU}(2)_L$\nleft action of the roots\non the self-dual\nmicrorotation",
        fs=9.8, ec=GREEN, fc="#eef4f0", tc=GREEN)
    arrow(ax, (6.6, 3.6), (8.0, 2.55),
          label="left mult.\n(same $Q_8$)", color=GREEN, rad=-0.12, fs=9, off=(0.0, 0.35))

    # the shared-Q8 tie between the two isospins
    arrow(ax, (3.0, 1.7), (7.2, 1.7), color=ACCENT, rad=-0.28, lw=1.3)
    ax.text(5.1, 0.62, "one $Q_8$: strong reads its representations, weak its left action",
            ha="center", va="center", fontsize=9.5, color=ACCENT, style="italic")

    # --- E6: McKay dual of the full root group --------------------------
    draw_E6hat(ax, 9.6, 5.25, s=0.30, color=GOLD)
    ax.text(9.6, 4.15, "$\\widehat{E}_6$", ha="center", va="center",
            fontsize=12, color=GOLD, weight="bold")
    arrow(ax, (7.0, 5.4), (8.7, 5.35),
          label="McKay", color=GOLD, fs=9, off=(0, 0.32))

    plt.tight_layout()
    fig.savefig("d4_root_group_e6_chain.pdf", bbox_inches="tight")
    print("wrote d4_root_group_e6_chain.pdf")


if __name__ == "__main__":
    main()
