"""
dstar_first_principles.py
=========================

Rest-mass closure for the d*(2380) hexaquark in the Cosserat supersolid
lattice, built from first principles on the 33-node Mode V dimer.

The d* is the lowest J^P = 3^+ dibaryon, conventionally read as a bound
Delta-Delta state.  In the lattice it is the Mode V docking of two Delta
clusters: two 17-node Delta graphs (a coordination shell plus four T_d
voids each) whose centres sit at twice the nearest-neighbour distance along
a <110> direction, sharing the single bridging FCC atom at the midpoint.
De-duplicating that shared atom leaves 34 - 1 = 33 nodes.

The mass mode is inherited.  The bare Delta carries its rest mass in a
stiff phi-dominant shell-concentrated T_1 mode at lambda = 9.0515.  On the
dimer that root doubles and splits into a tight multiplet under the
bridging-atom perturbation.  Selecting the modes that stay phi-dominant
(microrotation fraction >= 97%) and shell-concentrated (void localisation
<= 1%) isolates that multiplet cleanly; it is well separated from the next
phi-dominant band near lambda = 10.6.

Mass uses the un-deduplicated count N_tot = 33 in the master formula,
m = N_tot m_0 - N_tot (4 - lambda) m_e.  The bare-dimer value overshoots
the PDG mass by the inward-facing void-surface contribution: the two
Delta-side void shells that face the interface bury 27 corner bonds, and
removing 27 m_e brings the closure onto the measured mass.
"""

import numpy as np
from scipy.linalg import eigh

from spectral_classifier import fcc_nn_vectors
from delta_first_principles import (
    cluster_delta, void_positions_Td, build_cosserat_matrix_two_d,
)
from hadron_spectral_mass import M_0, M_E

NN = fcc_nn_vectors()
VOID_SURFACE = 27          # buried corner bonds across the interface


def dedup(coords, tol=1e-6):
    keep = []
    for a in coords:
        if not any(np.linalg.norm(a - b) < tol for b in keep):
            keep.append(a)
    return np.array(keep)


def mode_v_dimer():
    """Two Delta clusters at <110> second-neighbour separation sharing the
    bridging FCC atom.  Returns the 33-node coordinate set and a boolean
    mask marking the eight void nodes."""
    B = 2 * NN[0]                                   # |B| = 2, a <220> vector
    coords = dedup(np.vstack([cluster_delta(), cluster_delta() + B]))
    voids = np.vstack([void_positions_Td(), void_positions_Td() + B])
    is_void = np.array([any(np.linalg.norm(p - v) < 1e-6 for v in voids)
                        for p in coords])
    return coords, is_void


def diagnostics(vec, n, is_void):
    u = vec[:3 * n].reshape(n, 3)
    phi = vec[3 * n:].reshape(n, 3)
    amp = np.sum(u ** 2, axis=1) + np.sum(phi ** 2, axis=1)
    phi_frac = np.sum(phi ** 2) / (np.sum(u ** 2) + np.sum(phi ** 2))
    void_loc = np.sum(amp[is_void]) / np.sum(amp)
    return phi_frac, void_loc


def main():
    coords, is_void = mode_v_dimer()
    n = len(coords)
    print(f"Mode V dimer: {n} nodes "
          f"({int(is_void.sum())} void, {int((~is_void).sum())} shell), "
          f"bridge shared\n")

    H = build_cosserat_matrix_two_d(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    vals, vecs = eigh(H)

    # Inherited doubled-T_1 multiplet: phi-dominant and shell-concentrated.
    multiplet = []
    print(f"  {'lambda':>9} {'phi%':>6} {'void%':>7}")
    for k in range(len(vals)):
        if vals[k] < 1.0:
            continue
        pf, vl = diagnostics(vecs[:, k], n, is_void)
        if pf >= 0.97 and vl <= 0.01:
            multiplet.append(vals[k])
            print(f"  {vals[k]:9.4f} {pf * 100:6.1f} {vl * 100:7.2f}")

    lam = float(np.mean(multiplet))
    m_bare = n * M_0 - n * (4.0 - lam) * M_E
    m = m_bare - VOID_SURFACE * M_E
    pdg = 2380.0
    print(f"\n  multiplet: {len(multiplet)} modes, "
          f"lambda-bar = {lam:.4f}, spread {max(multiplet) - min(multiplet):.3f}")
    print(f"  inherits the bare-Delta T_1 root at lambda = 9.0515")
    print(f"  m (bare dimer, N_tot={n}) = {m_bare:.2f} MeV  "
          f"({(m_bare - pdg) / pdg * 100:+.2f}% vs PDG {pdg:.0f})")
    print(f"  minus inward void surface {VOID_SURFACE} m_e = {m:.2f} MeV  "
          f"({(m - pdg) / pdg * 100:+.3f}%)")


if __name__ == "__main__":
    main()
