#!/usr/bin/env python3
"""
birefringence_walls.py
======================
Cosmic birefringence from the fossil chirality-wall network (Confrontation
chapter, "Cosmic birefringence"). Two results, both parameter-light:

1. Coherent addition vs the endpoint theorem. A pseudoscalar (axion) rotation
   telescopes over alternating chirality domains to an endpoint bounded by
   theta_ch. A twist wall is instead a screw-dislocation grid whose handedness
   does NOT survive as a coherent sum: the twist handedness alternates wall to
   wall (see birefringence_sign_audit.py). The accumulating rotation is the
   EVEN coupling carried by the globally handed Shockley-partial fault ribbons
   (see birefringence_ribbon_coupling.py): beta = c1 * theta_ch * N_walls with
   c1 = (8 w / b) g and sign set by the global A->B->C stacking. This script
   supplies the count N_walls and the sub-Poisson anisotropy, which are
   independent of that mechanism question.

2. The chord-corrected crossing count and the sub-Poisson anisotropy. A
   straight line through a Poisson-Voronoi foam crosses walls at rate
   1/(f_chord * L_cell) with f_chord = 0.694 (measured here), giving N ~ 940
   over D_LSS = 14 Gpc and beta = c1 * 0.46 deg. The foam is sub-Poisson
   (Var(N)/N ~ 0.24), so sigma_beta/beta = 1/sqrt(N_eff) ~ 1.6%, half the
   random-domain value and independent of the normalisation c1.
"""
import numpy as np
from scipy.spatial import cKDTree

alpha = 1/137.035999177
theta_ch = alpha**2/(2*np.pi)
print(f"theta_ch = {theta_ch:.4e} rad = {np.degrees(theta_ch):.4e} deg")

# ---- 1. coherent addition vs telescoping ----
N = 940
sign = np.where(np.arange(N+1) % 2 == 0, 1, -1)     # alternating chirality domains
field = sign * theta_ch                              # pseudoscalar field value per domain
beta_axion = 0.5*(field[1:] - field[:-1]).sum()      # telescopes to endpoints
beta_grid = theta_ch * N                             # single-handed grid: adds (c1 = 1)
print("\n1. mechanism:")
print(f"   pseudoscalar telescopes to {np.degrees(beta_axion):.3e} deg (<= theta_ch)")
print(f"   screw grid adds to c1*theta_ch*N = {np.degrees(beta_grid):.3f} deg (c1=1, N={N})")

# ---- 2. chord factor from a ray-traced Poisson-Voronoi foam ----
rng = np.random.default_rng(1)
L_cell = 21.4
L_box = 12*L_cell
n_seed = int(round((L_box/L_cell)**3))
seeds = rng.uniform(0, L_box, size=(n_seed, 3))
tree = cKDTree(seeds, boxsize=L_box)
step, Lpath, nray = 0.25, 4000.0, 400
pos = rng.uniform(0, L_box, size=(nray, 3))
d = np.tile(np.array([1, 0, 0.0]), (nray, 1))
last = tree.query(pos % L_box, workers=-1)[1]
cross = np.zeros(nray); dist = 0.0
while dist < Lpath:
    pos = pos + step*d
    lab = tree.query(pos % L_box, workers=-1)[1]
    cross += (lab != last); last = lab; dist += step
f_chord = (1/(cross/Lpath).mean())/L_cell
D_LSS = 14000.0
N_walls = D_LSS/(f_chord*L_cell)
print("\n2. geometry:")
print(f"   mean chord / cell size f_chord = {f_chord:.3f}")
print(f"   N_walls = D_LSS/(f_chord*L_cell) = {N_walls:.0f}")
print(f"   beta(c1=1) = {np.degrees(theta_ch*N_walls):.3f} deg")

# ---- 3. sub-Poisson scatter from transverse ray bundle ----
Npix = 64
patch = np.deg2rad(8.0)
ang = (np.arange(Npix)-Npix/2)*patch/Npix
ax, ay = np.meshgrid(ang, ang, indexing="ij")
dirs = np.stack([ax.ravel(), ay.ravel(), np.ones(ax.size)], 1)
dirs /= np.linalg.norm(dirs, axis=1)[:, None]
nr = dirs.shape[0]
pos = np.zeros((nr, 3)) + rng.uniform(0, L_box, 3)
tile = lambda: cKDTree(rng.uniform(0, L_box, (n_seed, 3)), boxsize=L_box)
tr = tile(); last = tr.query(pos % L_box, workers=-1)[1]
counts = np.zeros(nr); dist = td = 0.0
while dist < D_LSS:
    pos = pos + 3.0*dirs; td += 3.0
    lab = tr.query(pos % L_box, workers=-1)[1]
    counts += (lab != last); last = lab; dist += 3.0
    if td >= L_box:
        tr = tile(); pos = pos + rng.uniform(0, L_box, 3)
        last = tr.query(pos % L_box, workers=-1)[1]; td = 0.0
vr = counts.var()/counts.mean()
Neff = counts.mean()**2/counts.var()
print("\n3. anisotropy (c1-independent):")
print(f"   Var(N)/<N> = {vr:.3f}  (1 = Poisson; <1 = sub-Poisson foam)")
print(f"   sigma_beta/beta = 1/sqrt(N_eff) = {1/np.sqrt(Neff):.4f}")
print(f"   vs random-domain Poisson 1/sqrt(N) = {1/np.sqrt(counts.mean()):.4f}")

# ---- 4. incidence-weighting family (part of the open c1 calculation) ----
# at a crossing, cos(theta) is flux-weighted: p(c) = 2c on [0,1]. The
# per-crossing rotation may be topological (w=1), projected (w=|cos|), or
# path-length (w=1/|cos|, regularised by wall thickness). The SIGNED
# projection telescopes (frame orientation is a state function) and is
# excluded; the survivors bracket beta at unit normalisation.
c = np.linspace(1e-3, 1, 200000); p = 2*c
tz = np.trapezoid
print("\n4. incidence-weighting family (c1 = 1):")
for name, w in [("topological  w=1      ", np.ones_like(c)),
                ("projected    w=|cos|  ", c),
                ("path-length  w=1/|cos|", 1/np.clip(c, 0.05, None))]:
    wm = tz(w*p, c); w2 = tz(w**2*p, c); varw = w2 - wm**2
    beta_w = np.degrees(theta_ch)*(1.455*wm/L_cell)*D_LSS
    Nrw = 1.455*D_LSS/L_cell
    s_foam = np.sqrt(varw/(Nrw*wm**2) + 0.24/Nrw)
    s_rand = np.sqrt(varw/(Nrw*wm**2) + 1.00/Nrw)
    print(f"   {name}: beta = {beta_w:4.2f} deg, sigma/beta foam {100*s_foam:.1f}% "
          f"vs random-domain {100*s_rand:.1f}%")
