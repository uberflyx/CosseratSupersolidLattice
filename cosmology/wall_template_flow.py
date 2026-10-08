#!/usr/bin/env python3
"""
wall_template_flow.py
=====================
Toy model of the peculiar-velocity field sourced by pinned fossil-wall dark
energy (cosmology chapter, "The Hubble tension" and the fossil-wall
falsifiable consequence (iv)).

The framework's dark energy is the stored strain of a chirality domain-wall
network, pinned to the lattice with cell size L_domain = 21.6 Mpc. It is
background-degenerate with a cosmological constant but perturbation-distinct:
the walls neither dilute nor comove, so their gravity is a static template
through which matter flows.

Coupling. The fossil stress fixes the trace T_kk = -3 eps pointwise, so the
walls' Newtonian source is eps + T_kk = -2 eps (XI_N = -2): they repel, with
twice the weight of their energy. All velocities below are linear in XI_N.

Box. The foam's contrast is white noise on scales above a cell, so the
template's bulk flow falls only as R^(-1/2), carried by wavelengths longer
than R. A box several times the largest sphere is needed; a 320 Mpc box
truncates those modes and makes the flow appear to die above the cell
scale. This version uses a 1280 Mpc box.

1. Point velocities, built kinematically over a Hubble time (v = a * t_H,
   an upper estimate since Hubble damping of accumulated peculiar velocity
   is neglected): a few hundred km/s rms, resolution dependent through the
   painted wall thickness.
2. Bulk flow in spheres: falls slowly with radius, roughly as R^(-1/2), and
   is far too small to source the CosmicFlows-4 excess reported at
   >~ 100 Mpc/h (Watkins et al. 2023; Whitford et al. 2023).
3. Uncorrected mock distance ladder (SNe in a 70-150 Mpc shell, 800 per
   observer, 256 observers): observer-to-observer spread of about half a
   km/s/Mpc, an order of magnitude short of the Hubble tension.

Model: periodic Poisson-Voronoi tessellation at the cell scale; the full
dark-energy density painted onto the cell faces with active weight XI_N;
perturbation Poisson equation solved by FFT on a 384^3 grid over a 1280 Mpc
box. Every statistic is pooled over N_REAL = 8 independent foam realisations
(a single box scatters by 10-20 per cent); pass another count as the first
argument. Memory about 3 GB; runtime about 1.5 min per realisation.
"""
import numpy as np
from scipy.spatial import cKDTree

# --- constants (SI) ---
G = 6.674e-11            # m^3 kg^-1 s^-2
Mpc = 3.086e22           # m
H0 = 67.4 * 1e3 / Mpc    # s^-1
t_H = 1.0 / H0           # Hubble time, the kinematic buildup window
rho_crit = 3 * H0**2 / (8 * np.pi * G)
rho_L = 0.685 * rho_crit  # dark-energy density (mass equivalent)

# Active gravitational density of the pinned walls, in units of their energy
# density. The fossil stress has T_kk = -3 eps pointwise (energy proportional
# to proper volume), so the Newtonian source eps + T_kk is -2 eps: the walls
# repel, with twice the weight of their energy.
XI_N = -2.0

L_box, N = 1280.0 * Mpc, 384
dx = L_box / N
L_cell = 21.6 * Mpc


def wall_template(rng, L_box=L_box, N=N):
    """Velocity field (three float32 arrays, m/s) and face-weight field of
    the painted foam in a periodic box of side L_box [m] on an N^3 grid; the
    face weight is returned for wall-resident sampling. The defaults are the
    large box used throughout; a smaller, finer box tests resolution."""
    dx = L_box / N
    n_seed = int(round((L_box / L_cell) ** 3))
    seeds = rng.uniform(0, L_box, size=(n_seed, 3))
    tree = cKDTree(seeds, boxsize=L_box)
    ax = (np.arange(N) + 0.5) * dx
    Y, Z = np.meshgrid(ax, ax, indexing="ij")
    yz = np.stack([Y.ravel(), Z.ravel()], axis=1)
    lab = np.empty((N, N, N), dtype=np.int32)
    for i in range(N):
        lab[i] = tree.query(np.column_stack([np.full(N * N, ax[i]), yz]),
                            workers=-1)[1].reshape(N, N)
    # face painting: each face between differently labelled voxels carries
    # equal surface density, split half to each neighbour
    face_w = np.zeros((N, N, N), dtype=np.float32)
    for a in range(3):
        f = (lab != np.roll(lab, -1, axis=a)).astype(np.float32)
        face_w += 0.5 * f + 0.5 * np.roll(f, 1, axis=a)
    del lab
    dk = np.fft.rfftn((XI_N * rho_L * (face_w / face_w.mean() - 1.0)).astype(np.float32))
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    kz = 2 * np.pi * np.fft.rfftfreq(N, d=dx)
    K2 = k1[:, None, None]**2 + k1[None, :, None]**2 + kz[None, None, :]**2
    K2[0, 0, 0] = 1.0
    phik = (-4 * np.pi * G * dk / K2).astype(np.complex64)
    phik[0, 0, 0] = 0.0
    del dk, K2
    v = []
    for Ki in (k1[:, None, None], k1[None, :, None], kz[None, None, :]):
        v.append((np.fft.irfftn(-1j * Ki * phik, s=(N, N, N), axes=(0, 1, 2))
                  * t_H).astype(np.float32))
    return v, face_w


def bulk_flows(v, rng, R_list, n_obs=64):
    """Bulk flow |<v>| [km/s] in spheres of each radius R [Mpc] about n_obs
    random grid points; returns an (n_obs, len(R_list)) array."""
    obs_idx = rng.integers(0, N, size=(n_obs, 3))
    out = np.empty((n_obs, len(R_list)))
    for j, R_mpc in enumerate(R_list):
        nR = int(R_mpc * Mpc / dx)
        off = np.arange(-nR, nR + 1)
        OX, OY, OZ = np.meshgrid(off, off, off, indexing="ij")
        inside = (OX**2 + OY**2 + OZ**2) * dx**2 <= (R_mpc * Mpc) ** 2
        ox, oy, oz = OX[inside], OY[inside], OZ[inside]
        for i, o in enumerate(obs_idx):
            ix, iy, iz = (o[0] + ox) % N, (o[1] + oy) % N, (o[2] + oz) % N
            out[i, j] = np.linalg.norm([v[c][ix, iy, iz].mean() for c in range(3)])
    return out / 1e3


def mock_ladder(v, rng, n_obs=256, n_sn=800):
    """Uncorrected H0 bias [km/s/Mpc] for n_obs random observers, each
    fitting n_sn supernovae spread over a 70-150 Mpc shell."""
    dH0 = np.empty(n_obs)
    for k in range(n_obs):
        o = rng.uniform(0, L_box, 3)
        oi = tuple((o // dx).astype(int) % N)
        vo = np.array([v[0][oi], v[1][oi], v[2][oi]])
        u = rng.normal(size=(n_sn, 3))
        u /= np.linalg.norm(u, axis=1)[:, None]
        r = (70 + 80 * rng.random(n_sn)) * Mpc
        x = (o + u * r[:, None]) % L_box
        xi = (x // dx).astype(int) % N
        vs = np.stack([v[c][xi[:, 0], xi[:, 1], xi[:, 2]] for c in range(3)], axis=1)
        dH0[k] = np.mean(np.einsum("ij,ij->i", vs - vo, u) / r)
    return dH0 * Mpc / 1e3


R_LIST = np.array([30, 50, 75, 100, 150, 200, 300])
N_REAL = 8               # independent foam realisations (seeds 42, 43, ...)


if __name__ == "__main__":
    import sys
    n_real = int(sys.argv[1]) if len(sys.argv) > 1 else N_REAL
    print(f"grid {N}^3, box {L_box/Mpc:.0f} Mpc, dx = {dx/Mpc:.2f} Mpc, "
          f"{n_real} foam realisations")
    rms, med, B_all, dH0_all = [], [], [], []
    for s in range(n_real):
        rng = np.random.default_rng(42 + s)
        v, _ = wall_template(rng)
        rms.append(np.sqrt(np.mean(v[0].astype(np.float64)**2 + v[1]**2 + v[2]**2)) / 1e3)
        B = bulk_flows(v, rng, R_LIST)
        med.append(np.median(B, axis=0))
        B_all.append(B)
        dH0_all.append(mock_ladder(v, rng))
        del v
        print(f"  realisation {s}: rms {rms[-1]:.0f} km/s, bulk medians "
              f"{np.round(med[-1]).astype(int).tolist()} km/s, "
              f"ladder spread {dH0_all[-1].std():.3f} km/s/Mpc")
    med, B_all, dH0 = np.array(med), np.concatenate(B_all), np.concatenate(dH0_all)

    print(f"\ntemplate speed: rms = {np.mean(rms):.0f} +- {np.std(rms):.0f} km/s")
    print("\nbulk flow |<v>| in spheres (64 observers per realisation, pooled):")
    print("  R [Mpc]   pooled median   16-84th pct   realisation medians: mean +- sd")
    for j, R_mpc in enumerate(R_LIST):
        print(f"  {R_mpc:5d}      {np.median(B_all[:, j]):5.0f} km/s"
              f"      {np.percentile(B_all[:, j], 16):3.0f}-{np.percentile(B_all[:, j], 84):3.0f}"
              f"          {med[:, j].mean():5.0f} +- {med[:, j].std():3.0f}")
    print("\nmock ladder (uncorrected), SN shell 70-150 Mpc, 800 SNe/observer,")
    print(f"256 observers per realisation ({len(dH0)} in all):")
    print(f"  mean bias  = {dH0.mean():+.3f} km/s/Mpc  (repulsive walls, XI_N = -2)")
    print(f"  rms spread = {dH0.std():.3f} km/s/Mpc"
          f"  (per realisation {np.min([d.std() for d in dH0_all]):.2f}-"
          f"{np.max([d.std() for d in dH0_all]):.2f})")
    print(f"  16-84th pct = {np.percentile(dH0,16):+.3f} .. {np.percentile(dH0,84):+.3f}")
