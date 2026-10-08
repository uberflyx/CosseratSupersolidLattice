"""
sn_wall_residual.py
===================
The pinned-wall velocity template as a low-redshift supernova systematic
(Confrontation chapter, dark-energy discussion).

The DESI preference for evolving dark energy is driven by the low-redshift
supernova magnitudes (Efstathiou 2025; Vincenzi et al. 2025 in rebuttal).
The fossil wall network contributes a static, NON-COMOVING velocity template
that comoving peculiar-velocity reconstructions cannot remove. If supernova
hosts and the observer sit on the walls, the sampled template does not
average to zero at low redshift.

The template is built by wall_template_flow.py (same directory): a
1280 Mpc periodic Poisson-Voronoi foam carrying the full dark-energy density
with the active weight XI_N = -2 that the fossil stress fixes, velocities
built kinematically over a Hubble time (an upper estimate). Residuals are
linear in the coupling.

Per redshift bin: the coherent sky-averaged line-of-sight velocity
difference between wall-resident hosts and a wall-resident observer, turned
into a magnitude offset delta_mu = (5/ln10) <dv_los>/(cz); the rms of that
offset over observers is reported. Wall residence is an assumption: the
walls repel matter, and whether galaxies end up on them is open.

Pooled over N_REAL = 8 foam realisations (first argument overrides).
Run from this directory. Memory about 3 GB.
"""
import sys
import numpy as np
from scipy.spatial import cKDTree
from wall_template_flow import wall_template, L_box, dx, H0, N_REAL

c_ms = 2.998e8
zbins = np.array([0.015, 0.025, 0.04, 0.06, 0.09, 0.13])


def sky_mean_residual(v, face_w, rng, n_obs=64, n_sn=600):
    """Coherent sky-mean line-of-sight velocity difference [m/s] between
    wall-resident hosts and a wall-resident observer, per redshift bin."""
    wi = np.argwhere(face_w > 0)
    wall_pos = (wi + 0.5) * dx
    wtree = cKDTree(wall_pos, boxsize=L_box)
    res = np.zeros((n_obs, len(zbins)))
    for o in range(n_obs):
        oi = wi[rng.integers(len(wi))]           # observer on a wall
        opos = (oi + 0.5) * dx
        vo = np.array([v[c][tuple(oi)] for c in range(3)])
        for bz, zc in enumerate(zbins):
            r = zc * c_ms / H0                     # below half the box for z <= 0.13
            u = rng.normal(size=(n_sn, 3))
            u /= np.linalg.norm(u, axis=1)[:, None]
            _, j = wtree.query((opos + u * r) % L_box, workers=-1)   # hosts on walls
            hv = np.stack([v[c][wi[j][:, 0], wi[j][:, 1], wi[j][:, 2]] for c in range(3)], axis=1)
            d = wall_pos[j] - opos
            d -= L_box * np.round(d / L_box)
            rhat = d / np.linalg.norm(d, axis=1)[:, None]
            res[o, bz] = np.einsum("ij,ij->i", hv - vo, rhat).mean()
    return res


n_real = int(sys.argv[1]) if len(sys.argv) > 1 else N_REAL
RES = []
for s in range(n_real):
    rng = np.random.default_rng(11 + 100 * s)
    v, face_w = wall_template(rng)
    RES.append(sky_mean_residual(v, face_w, rng))
    del v, face_w
    dmu_s = (5 / np.log(10)) * RES[-1] / (zbins * c_ms)
    print(f"  realisation {s}: delta_mu rms = {np.round(dmu_s.std(axis=0), 4).tolist()} mag")
dmu_real = np.array([((5 / np.log(10)) * R / (zbins * c_ms)).std(axis=0) for R in RES])
res = np.concatenate(RES)

print(f"\nper-observer coherent sky-mean residual, repulsive walls (XI_N = -2),")
print(f"{len(res)} wall-resident observers over {n_real} realisations:")
print("z_bin   <dv_los> km/s (mean +- rms)   delta_mu rms [mag]   (per-realisation sd)")
for bz, zc in enumerate(zbins):
    m = res[:, bz] / 1e3
    dmu = (5 / np.log(10)) * (res[:, bz] / (zc * c_ms))
    print(f" {zc:5.3f}   {m.mean():+7.1f} +- {m.std():6.1f}        {dmu.std():.4f}"
          f"               {dmu_real[:, bz].std():.4f}")
