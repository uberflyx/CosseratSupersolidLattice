#!/usr/bin/env python3
"""
wall_template_flow.py
=====================
Toy model of the peculiar-velocity field sourced by pinned fossil-wall dark
energy (cosmology chapter, "The Hubble tension" and the fossil-wall
falsifiable consequence (iv)).

The framework's dark energy is the stored strain of a chirality domain-wall
network, pinned to the lattice with cell size L_domain = 21.4 Mpc. It is
background-degenerate with a cosmological constant but perturbation-distinct:
the walls neither dilute nor comove, so their gravity is a static template
through which matter flows. This script sizes the template's observable
consequences:

1. Point velocities. Building velocity kinematically over a Hubble time
   (v = a * t0, an upper estimate since Hubble damping of accumulated
   peculiar velocity is neglected) gives ~250 km/s rms.

2. Bulk flow vs scale. The template's coherence is confined below the cell
   scale: median |<v>| falls from ~87 km/s in R = 30 Mpc spheres to ~11 km/s
   at R = 150 Mpc. It therefore CANNOT source the CosmicFlows-4 bulk-flow
   excess reported at >~ 100 Mpc/h (Watkins et al. 2023; Whitford et al.
   2023), whose defining feature is the opposite scale dependence. What it
   predicts instead is a web-locked, non-comoving component of a few tens of
   km/s on cell scales, a fingerprint no smooth dark energy produces.

3. Distance-ladder imprint. An uncorrected mock ladder (SNe in a 70-150 Mpc
   shell, 800 per observer, 256 observers) shows an observer-to-observer
   spread of +-0.2 km/s/Mpc and a mean bias of +0.03 km/s/Mpc for the
   attractive coupling sign (the repulsive sign mirrors it). The local-side
   wall channel into the Hubble tension is closed at the few-tenths level.

Model: periodic Voronoi tessellation at the cell scale; the full dark-energy
density painted onto the Voronoi boundaries; perturbation Poisson equation
solved by FFT on a 256^3 grid over a 320 Mpc box. The sign of the wall's
local gravitational coupling (attractive energy vs anisotropic-stress
corrections) is an open question in the framework; it flips the sign of the
mean ladder bias and nothing else.
"""
import numpy as np
from scipy.spatial import cKDTree

rng = np.random.default_rng(42)

# --- constants (SI) ---
G = 6.674e-11            # m^3 kg^-1 s^-2
Mpc = 3.086e22           # m
H0 = 67.4 * 1e3 / Mpc    # s^-1
t_H = 1.0 / H0           # Hubble time, the kinematic buildup window
rho_crit = 3 * H0**2 / (8 * np.pi * G)
rho_L = 0.685 * rho_crit  # dark-energy density (mass equivalent)

# --- periodic Voronoi wall network ---
L_box, N = 320.0 * Mpc, 256
dx = L_box / N
L_cell = 21.4 * Mpc                       # fossil grain scale L_domain
n_seed = int(round((L_box / L_cell) ** 3))
print(f"grid {N}^3, box {L_box/Mpc:.0f} Mpc, dx = {dx/Mpc:.2f} Mpc, "
      f"cells ~ {n_seed}")

seeds = rng.uniform(0, L_box, size=(n_seed, 3))
tree = cKDTree(seeds, boxsize=L_box)
ax1d = (np.arange(N) + 0.5) * dx
X, Y, Z = np.meshgrid(ax1d, ax1d, ax1d, indexing="ij")
pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
_, lab = tree.query(pts, workers=-1)
lab = lab.reshape(N, N, N).astype(np.int32)
del pts, X, Y, Z

wall = np.zeros((N, N, N), dtype=bool)
for axis in range(3):
    wall |= (lab != np.roll(lab, 1, axis=axis))
print(f"wall cell fraction = {wall.mean():.4f}")
del lab

# --- perturbation gravity: all DE mass on walls, mean subtracted ---
drho = np.where(wall, rho_L / wall.mean(), 0.0) - rho_L
del wall
k1 = 2 * np.pi * np.fft.fftfreq(N, d=dx)
KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
K2 = KX**2 + KY**2 + KZ**2
K2[0, 0, 0] = 1.0
phi_k = -4 * np.pi * G * np.fft.fftn(drho) / K2   # attractive sign
phi_k[0, 0, 0] = 0.0
v = [np.real(np.fft.ifftn(-1j * Ki * phi_k)) * t_H for Ki in (KX, KY, KZ)]
del drho, phi_k, KX, KY, KZ, K2

v_mag = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
print(f"template speed: rms = {np.sqrt((v_mag**2).mean())/1e3:.0f} km/s, "
      f"90th pct = {np.percentile(v_mag, 90)/1e3:.0f} km/s")
del v_mag

# --- 1. bulk flow in spheres ---
print("\nbulk flow |<v>| in spheres (median over 64 observers):")
obs_idx = rng.integers(0, N, size=(64, 3))
for R_mpc in (30, 50, 75, 100, 150):
    nR = int(R_mpc * Mpc / dx)
    off = np.arange(-nR, nR + 1)
    OX, OY, OZ = np.meshgrid(off, off, off, indexing="ij")
    inside = (OX**2 + OY**2 + OZ**2) * dx**2 <= (R_mpc * Mpc) ** 2
    ox, oy, oz = OX[inside], OY[inside], OZ[inside]
    Bs = []
    for o in obs_idx:
        ix, iy, iz = (o[0] + ox) % N, (o[1] + oy) % N, (o[2] + oz) % N
        Bs.append(np.linalg.norm([v[0][ix, iy, iz].mean(),
                                  v[1][ix, iy, iz].mean(),
                                  v[2][ix, iy, iz].mean()]))
    Bs = np.array(Bs) / 1e3
    print(f"  R = {R_mpc:4d} Mpc: median B = {np.median(Bs):5.0f} km/s "
          f"(16-84th pct: {np.percentile(Bs,16):.0f}-{np.percentile(Bs,84):.0f})")

# --- 2. uncorrected mock distance ladder ---
print("\nmock ladder (uncorrected), SN shell 70-150 Mpc, 800 SNe/observer:")
n_obs, n_sn = 256, 800
dH0 = []
for _ in range(n_obs):
    o = rng.uniform(0, L_box, 3)
    oi = tuple((o // dx).astype(int) % N)
    vo = np.array([v[0][oi], v[1][oi], v[2][oi]])
    u = rng.normal(size=(n_sn, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    r = (70 + 80 * rng.random(n_sn)) * Mpc
    x = (o + u * r[:, None]) % L_box
    xi = (x // dx).astype(int) % N
    vs = np.stack([v[0][xi[:, 0], xi[:, 1], xi[:, 2]],
                   v[1][xi[:, 0], xi[:, 1], xi[:, 2]],
                   v[2][xi[:, 0], xi[:, 1], xi[:, 2]]], axis=1)
    dH0.append(np.mean(np.einsum("ij,ij->i", vs - vo, u) / r))
dH0 = np.array(dH0) * Mpc / 1e3
print(f"  mean bias  = {dH0.mean():+.3f} km/s/Mpc  (attractive sign; flips for repulsive)")
print(f"  rms spread = {dH0.std():.3f} km/s/Mpc")
print(f"  16-84th pct = {np.percentile(dH0,16):+.3f} .. {np.percentile(dH0,84):+.3f}")
print("\nconclusion: the local-side wall channel is closed at the few-tenths level;")
print("the surviving prediction is the web-locked tens-of-km/s velocity fingerprint.")
