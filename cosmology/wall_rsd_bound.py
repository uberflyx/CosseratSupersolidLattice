"""
wall_rsd_bound.py
=================
Consistency of the pinned-wall velocity template with redshift-space
distortions, and its (null) effect on the BAO+CMB neutrino-mass bound
(Confrontation chapter, neutrino-mass and dark-energy discussions).

The template is built by wall_template_flow.py (same directory): a
1280 Mpc periodic Poisson-Voronoi foam carrying the full dark-energy density
with the active weight XI_N = -2 that the fossil stress fixes, velocities
built kinematically over a Hubble time (an upper estimate).

Because the template is coherent within each cell, nearby tracers share
their motion and the PAIRWISE velocity (the observable redshift-space
distortions constrain) is far smaller than the point velocities. For
wall-resident tracers it is an outflow at 10-30 Mpc separations, compared
here with a ~30 km/s residual window on the ~300 km/s matter infall. The
estimate neglects Hubble damping and assumes the tracers stay on the walls,
which repulsive walls themselves work against. At large separations the
pairwise velocity falls to a few km/s, so the displacement at the 150 Mpc
acoustic scale stays far below DESI precision and incoherent in sign: the
BAO+CMB neutrino-mass bound is untouched by the framework's own template.

The painted walls are one grid cell thick, so the pairwise velocity at
10-15 Mpc depends on resolution; a 320 Mpc box on a 1.67 Mpc grid checks it.
Both are averaged over N_REAL = 8 foam realisations (first argument
overrides). Run from this directory. Memory about 3 GB.
"""
import sys
import numpy as np
from scipy.spatial import cKDTree
from wall_template_flow import wall_template, L_box, dx, Mpc, t_H, N_REAL

rbins = np.array([5, 10, 15, 22, 30, 40, 60, 100]) * Mpc


def pairwise_velocity(v, face_w, rng, L_box=L_box, dx=dx, n_tr=60000, n_bins=None,
                      max_pairs=300000):
    """Mean radial pairwise velocity v12 [km/s] of wall-resident tracers in
    the first n_bins separation bins (all by default); positive is outflow."""
    wi = np.argwhere(face_w > 0)
    sel = wi[rng.integers(len(wi), size=n_tr)]
    pos = (sel + 0.5) * dx
    vel = np.stack([v[c][sel[:, 0], sel[:, 1], sel[:, 2]] for c in range(3)], axis=1)
    ptree = cKDTree(pos, boxsize=L_box)
    v12 = []
    for i in range(n_bins or len(rbins) - 1):
        pairs = ptree.query_pairs(rbins[i + 1], output_type="ndarray")
        d = pos[pairs[:, 1]] - pos[pairs[:, 0]]
        d -= L_box * np.round(d / L_box)
        rr = np.linalg.norm(d, axis=1)
        idx = np.where(rr >= rbins[i])[0]
        if len(idx) > max_pairs:
            idx = rng.choice(idx, max_pairs, replace=False)
        rhat = d[idx] / rr[idx, None]
        dv = np.einsum("ij,ij->i", vel[pairs[idx, 1]] - vel[pairs[idx, 0]], rhat)
        v12.append(dv.mean() / 1e3)
    return np.array(v12)


n_real = int(sys.argv[1]) if len(sys.argv) > 1 else N_REAL
V12 = []
for s in range(n_real):
    rng = np.random.default_rng(5 + 100 * s)
    v, face_w = wall_template(rng)
    V12.append(pairwise_velocity(v, face_w, rng))
    del v, face_w
    print(f"  realisation {s}: v12 = {np.round(V12[-1], 1).tolist()} km/s")
V12 = np.array(V12)
v12, sd = V12.mean(axis=0), V12.std(axis=0)

# Resolution check: the painted walls are one grid cell thick, so the
# pairwise velocity at 10-15 Mpc depends on the grid. Repeat on a 1.67 Mpc
# grid in a 320 Mpc box (too small for the bulk flow, adequate below 40 Mpc).
L_f, N_f = 320.0 * Mpc, 192
V12f = []
for s in range(n_real):
    rng = np.random.default_rng(7 + 100 * s)
    v, face_w = wall_template(rng, L_box=L_f, N=N_f)
    V12f.append(pairwise_velocity(v, face_w, rng, L_box=L_f, dx=L_f / N_f, n_tr=30000, n_bins=5))
    del v, face_w
V12f = np.array(V12f)

print(f"\ntemplate pairwise radial velocity v12(r) for wall tracers, "
      f"mean +- sd over {n_real} realisations:")
print(" r [Mpc]   v12 [km/s]  (positive = outflow)")
for i in range(len(rbins) - 1):
    print(f"  {rbins[i]/Mpc:3.0f}-{rbins[i+1]/Mpc:3.0f}   {v12[i]:+7.1f} +- {sd[i]:4.1f}")

print(f"\nresolution check, {L_f/Mpc:.0f} Mpc box on a {L_f/N_f/Mpc:.2f} Mpc grid (r < 40 Mpc):")
for i in range(5):
    print(f"  {rbins[i]/Mpc:3.0f}-{rbins[i+1]/Mpc:3.0f}   {V12f[:, i].mean():+7.1f} +- {V12f[:, i].std():4.1f}")

# LCDM matter pairwise infall at 10-25 Mpc ~ 250-350 km/s; growth-rate data
# allow a residual of about 10 per cent
allowed = 0.10 * 300.0
vmax = np.abs(v12[1:4]).max()
vmax_f = np.abs(V12f.mean(axis=0)[1:4]).max()
print(f"\nmax |v12| (10-30 Mpc)       = {vmax:.0f} km/s (large box), {vmax_f:.0f} km/s (fine grid)")
print(f"RSD-allowed extra pairwise  ~ {allowed:.0f} km/s (10% of ~300)")
print(f"=> allowed fraction of the derived coupling: {allowed/vmax_f:.2f}-{allowed/vmax:.2f}")
disp = 0.5 * v12[-1] * 1e3 * t_H / Mpc
print(f"\nBAO-scale displacement from the 60-100 Mpc pairwise velocity: "
      f"{disp:.3f} Mpc = {100*disp/147:.3f}% of the sound horizon")
