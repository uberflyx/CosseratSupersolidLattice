"""Figure: the edge dislocation's brane is a helical ramp, not a cylinder.

Draws the two surfaces side by side over one period of the compact direction.
Left, the cylinder the naive reading assumes: the line sits at the same
in-plane position on every layer, so its three-dimensional footprint is a
point-like segment of no determinate length.  Right, the ramp the stacking
actually builds: each layer is shifted in-plane by one Shockley vector
(a/6)<112>, so after the three layers that close the compact circle the brane
has slid one full <112> lattice period, sqrt(3) ell, and that slide is the
footprint length L_xi.

The horizontal axis is the in-plane <112> line direction in units of the
partial period d = ell/sqrt(3); the vertical axis is the compact coordinate in
units of the interlayer spacing d_111.

Writes edge_brane_helicoid.svg (Quarto edition) and .pdf (TeX edition).
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# Geometry, in units of the partial period d = ell/sqrt(3)
# --------------------------------------------------------------------------
N_LAYERS = 3            # L_4 = 3 d_111 closes the compact circle
SHOCKLEY = 1.0          # one registry shift per layer, = d by construction
LABELS = ["A", "B", "C", "A"]

BLUE = "#1f4e8c"
RUST = "#b0521a"
GREY = "0.55"


def draw_stack(ax, slide_per_layer, title, footprint_label=None):
    """One panel: N_LAYERS+1 close-packed layers, offset by slide_per_layer."""
    seg_half = 0.62                      # drawn half-length of the line segment
    for n in range(N_LAYERS + 1):
        y = n
        x0 = n * slide_per_layer
        # the close-packed layer, drawn as a faint rule
        ax.plot([-1.4, 4.2], [y, y], "-", color="0.88", lw=0.8, zorder=1)
        # the dislocation line piercing this layer
        ax.plot([x0 - seg_half, x0 + seg_half], [y, y], "-",
                color=BLUE, lw=3.0, solid_capstyle="round", zorder=3)
        ax.text(-1.55, y, LABELS[n], color=RUST, fontsize=9,
                ha="right", va="center")
        # the registry shift carrying this layer to the next
        if n < N_LAYERS:
            ax.annotate("", xy=(x0 + slide_per_layer, y + 1),
                        xytext=(x0, y),
                        arrowprops=dict(arrowstyle="->", color=RUST, lw=1.1,
                                        shrinkA=3, shrinkB=3))

    # the ruled surface swept between the layers
    xs = np.array([0.0, N_LAYERS * slide_per_layer])
    ax.fill_between(
        [xs[0] - seg_half, xs[0] + seg_half,
         xs[1] - seg_half, xs[1] + seg_half],
        [0, 0, N_LAYERS, N_LAYERS],
        [0, 0, N_LAYERS, N_LAYERS],
        color=BLUE, alpha=0.0)
    poly_x = [-seg_half, seg_half,
              N_LAYERS * slide_per_layer + seg_half,
              N_LAYERS * slide_per_layer - seg_half]
    poly_y = [0, 0, N_LAYERS, N_LAYERS]
    ax.fill(poly_x, poly_y, color=BLUE, alpha=0.13, zorder=2, lw=0)

    # the three-dimensional footprint: the projection onto one plane
    y_f = -0.75
    ax.plot([-seg_half, N_LAYERS * slide_per_layer + seg_half], [y_f, y_f],
            "-", color="0.25", lw=1.0, zorder=3)
    for xe in (-seg_half, N_LAYERS * slide_per_layer + seg_half):
        ax.plot([xe, xe], [y_f - 0.12, y_f + 0.12], "-", color="0.25", lw=1.0)
    if footprint_label:
        ax.text(0.5 * N_LAYERS * slide_per_layer, y_f - 0.30,
                footprint_label, fontsize=8.5, ha="center", va="top")

    ax.set_title(title, fontsize=9.5, pad=8)
    ax.set_xlim(-1.9, 4.5)
    ax.set_ylim(-1.6, N_LAYERS + 0.7)
    ax.axis("off")


def make_figure(stem="edge_brane_helicoid"):
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.3))

    draw_stack(axes[0], 0.0,
               "cylinder: no slide,\nno determinate length",
               footprint_label="footprint undetermined")
    draw_stack(axes[1], SHOCKLEY,
               "helical ramp: one Shockley step per layer",
               footprint_label=r"$L_\xi = 3d = \sqrt{3}\,\ell$")

    axes[1].text(4.35, 1.5,
                 r"$\frac{a}{6}\langle 112\rangle$" + "\nper layer",
                 color=RUST, fontsize=8, ha="right", va="center")
    fig.text(0.5, 0.015,
             r"in-plane $\langle 112\rangle$ (units of $d = \ell/\sqrt{3}$)"
             r"$\;\longrightarrow$;  vertical: compact direction, "
             r"$L_4 = 3\,d_{111}$",
             fontsize=7.5, ha="center")

    fig.tight_layout(rect=(0, 0.05, 1, 1))
    for ext in ("svg", "pdf"):
        fig.savefig(f"{stem}.{ext}")
        print(f"wrote {stem}.{ext}")


if __name__ == "__main__":
    make_figure()
