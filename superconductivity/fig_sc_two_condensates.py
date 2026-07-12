"""
Figure: the two condensates and the hierarchy of windings.

Three panels across five orders of magnitude in scale:
  (a) One electron: a closed winding of the VACUUM condensate at the core
      scale, whose far field in the condensate is a dipole along the spin.
  (b) A Cooper pair: two windings with the same handedness and opposite
      spins, so the two dipoles cancel and the pair is dark to the
      condensate's flow channel.
  (c) A superconducting ring: the ELECTRONIC condensate's phase winds by
      2*pi*n around the ring, quantising the trapped flux at h/2e.

Output: fig_sc_two_condensates.pdf (vector, for the monograph figures/ dir).
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Annulus

fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.8))

BLUE, RED, GREEN, GREY = "#1f77b4", "#b22222", "#2e8b57", "0.45"


def draw_winding(ax, x0, y0, spin_up=True, r=0.55, colour=BLUE):
    """One closed condensate winding: core ring, circulation arrows, spin axis."""
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(x0 + r * np.cos(th), y0 + r * np.sin(th), color=colour, lw=2.4)
    sgn = 1 if spin_up else -1
    for t in np.linspace(0, 2 * np.pi, 4, endpoint=False) + 0.4:
        dx, dy = -np.sin(t) * sgn, np.cos(t) * sgn
        ax.annotate("", xy=(x0 + r * np.cos(t) + 0.14 * dx,
                            y0 + r * np.sin(t) + 0.14 * dy),
                    xytext=(x0 + r * np.cos(t), y0 + r * np.sin(t)),
                    arrowprops=dict(arrowstyle="-|>", color=colour, lw=1.6))
    ax.annotate("", xy=(x0, y0 + sgn * 1.05), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=RED, lw=2.2))
    ax.text(x0 + 0.14, y0 + sgn * 0.95, r"$\mathbf{d}\parallel\hat{\mathbf{s}}$",
            color=RED, fontsize=10, va="center")


# ---- (a) one electron in the vacuum condensate -------------------------
ax = axes[0]
draw_winding(ax, 0, 0, spin_up=True)
# dipolar far-field sketch: a few streamline arcs
for s in (1.05, 1.5, 2.0):
    th = np.linspace(0.25, np.pi - 0.25, 100)
    ax.plot(s * np.sin(th) * 0.9, s * np.cos(th) * 1.25, color=GREY,
            lw=0.8, alpha=0.7)
    ax.plot(-s * np.sin(th) * 0.9, s * np.cos(th) * 1.25, color=GREY,
            lw=0.8, alpha=0.7)
ax.text(0, -2.55, r"far flow $\sim \kappa\ell^{2}/r^{3}$: dipole only",
        ha="center", fontsize=10)
ax.set_title(r"(a) one electron $\sim \ell$ (fm)", fontsize=11)

# ---- (b) the dark pair --------------------------------------------------
ax = axes[1]
draw_winding(ax, -1.05, 0, spin_up=True)
draw_winding(ax, 1.05, 0, spin_up=False)
ax.text(0, -2.55, r"$\mathbf{d}_1+\mathbf{d}_2 \propto \hat{\mathbf{S}}_{\rm tot}"
                  r"\;\Rightarrow\; 0$ on the singlet",
        ha="center", fontsize=10)
ax.set_title(r"(b) Cooper pair $\sim \xi$ (40 nm)", fontsize=11)

# ---- (c) the superconducting ring ---------------------------------------
ax = axes[2]
ax.add_patch(Annulus((0, 0), 1.9, 0.55, color="#c9d7e6", zorder=1))
th = np.linspace(0, 2 * np.pi, 300)
ax.plot(1.62 * np.cos(th), 1.62 * np.sin(th), color=BLUE, lw=2.2, zorder=2)
for t in np.linspace(0, 2 * np.pi, 8, endpoint=False):
    ax.annotate("", xy=(1.62 * np.cos(t) - 0.16 * np.sin(t),
                        1.62 * np.sin(t) + 0.16 * np.cos(t)),
                xytext=(1.62 * np.cos(t), 1.62 * np.sin(t)),
                arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.5), zorder=3)
ax.add_patch(Circle((0, 0), 0.35, color=GREEN, alpha=0.35, zorder=2))
ax.text(0, 0, r"$\Phi = n\,\dfrac{h}{2e}$", ha="center", va="center",
        fontsize=11, color=GREEN, zorder=4)
ax.text(0, -2.55, r"electronic phase winds: $\oint\nabla\theta\cdot d\mathbf{l}"
                  r" = 2\pi n$", ha="center", fontsize=10)
ax.set_title(r"(c) ring $\sim$ mm: $10^{23}$ dark pairs, one phase",
             fontsize=11)

for ax in axes:
    ax.set_xlim(-2.9, 2.9)
    ax.set_ylim(-2.9, 2.9)
    ax.set_aspect("equal")
    ax.axis("off")

fig.tight_layout()
fig.savefig("fig_sc_two_condensates.pdf")
print("wrote fig_sc_two_condensates.pdf")
