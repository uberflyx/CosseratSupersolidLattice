#!/usr/bin/env python3
"""
wall_foam_lensing.py
====================
Lensing, CMB-temperature and void-stack signals of the pinned fossil-wall
network, set against Planck LCDM (cosmology chapter, the walls' gravity and
the dark-energy falsifiable consequences).

Physics. In the static weak field ds^2 = -(1+2Phi/c^2)c^2dt^2 +
(1-2Psi/c^2)dx^2, a source with energy density eps and stress T_ij obeys
    lap Psi       = 4 pi G eps / c^2
    lap Phi       = 4 pi G (eps + T_kk) / c^2      (moves slow matter)
    lap (Phi+Psi) = 4 pi G (2 eps + T_kk) / c^2    (bends light, ISW)
The fossil walls carry T_kk = -3 eps, so their dynamical source is -2 eps
(XI_N = -2) and their lensing source -eps (ETA_L = -1/2 of dust's 2 eps).

The wall pattern is a periodic Poisson-Voronoi foam, one cell per L^3 with
L = 21.6 Mpc, carrying the whole dark-energy density Omega_w = 0.685 on its
faces. Faces are painted on the voxel staircase and the field is normalised
to its mean, so the staircase's 3/2 area overcount drops out of the contrast.

Outputs (printed):
  1. foam contrast spectrum P_w(k): white-noise level P_w(0) ~ 0.11 L^3 and
     the cell-scale peak near kL ~ 6;
  2. dynamical-potential power of walls relative to matter at z = 0;
  3. convergence power of the walls relative to matter (Limber) for sources
     at z = 0.6, 1 and the CMB, with the pattern read either as fixed in
     proper coordinates or as comoving;
  4. ISW temperature power of the walls relative to Planck TT, for the
     comoving reading in which the potentials grow as a^2;
  5. the stacked excess surface density DeltaSigma(R) of the walls around
     foam cells, against a Hamaus-Sutter-Wandelt matter void profile with
     the DES-SV best-fit shape (delta_c = -0.60, r_s = 1.05 R_v) at the same
     radius.

Requires numpy, scipy, camb. Runtime a few minutes on two cores.
"""
import numpy as np
from scipy.spatial import cKDTree
import camb

H = 0.674
L_CELL = 21.6        # Mpc
OM_W = 0.685         # walls carry all the dark energy
XI_N = -2.0          # dynamical source / eps
ETA_L = -0.5         # lensing source / (2 eps)


def foam(box_L, n_grid, seed):
    """labels and painted contrast of a periodic Poisson-Voronoi foam (L = 1)"""
    rng = np.random.default_rng(seed)
    h = box_L / n_grid
    seeds = rng.uniform(0, box_L, size=(int(round(box_L**3)), 3))
    tree = cKDTree(seeds, boxsize=box_L)
    ax = (np.arange(n_grid) + 0.5) * h
    Y, Z = np.meshgrid(ax, ax, indexing="ij")
    yz = np.stack([Y.ravel(), Z.ravel()], axis=1)
    lab = np.empty((n_grid,) * 3, dtype=np.int32)
    for i in range(n_grid):
        lab[i] = tree.query(np.column_stack([np.full(n_grid**2, ax[i]), yz]),
                            workers=-1)[1].reshape(n_grid, n_grid)
    rho = np.zeros((n_grid,) * 3, dtype=np.float32)
    for a in range(3):
        face = (lab != np.roll(lab, -1, axis=a)).astype(np.float32)
        rho += 0.5 * face + 0.5 * np.roll(face, 1, axis=a)
    return lab, rho / rho.mean() - 1.0, h, ax


def spectrum(delta, box_L, h):
    n = delta.shape[0]
    P3 = np.abs(np.fft.rfftn(delta) * h**3) ** 2 / box_L**3
    k1 = 2 * np.pi * np.fft.fftfreq(n, d=h)
    kz = 2 * np.pi * np.fft.rfftfreq(n, d=h)
    K = np.sqrt(k1[:, None, None]**2 + k1[None, :, None]**2 + kz[None, None, :]**2)
    w = np.full(K.shape, 2.0)
    w[..., 0] = 1.0
    w[..., -1] = 1.0
    kf = 2 * np.pi / box_L
    idx = np.rint(K / kf).astype(int).ravel()
    num = np.bincount(idx, (P3 * w).ravel())
    den = np.bincount(idx, w.ravel())
    kk = np.bincount(idx, (K * w).ravel())
    m = (den > 0) & (np.arange(len(den)) > 0) & (np.arange(len(den)) < n // 2)
    return kk[m] / den[m], num[m] / den[m]


print("1. foam contrast spectrum (units L = 21.6 Mpc)")
runs = {}
for box_L, n_grid in ((48.0, 384), (16.0, 256)):
    Ps = []
    for seed in (1, 2):
        _, d, h, _ = foam(box_L, n_grid, seed)
        k, P = spectrum(d, box_L, h)
        Ps.append(P)
        del d
    runs[box_L] = (k, np.mean(Ps, axis=0))
kb, pb = runs[48.0]
kf, pf = runs[16.0]
kL = np.concatenate([kb[kb < 1.6], kf[(kf >= 1.6) & (kf < 8 * np.pi)]])
PL3 = np.concatenate([pb[kb < 1.6], pf[(kf >= 1.6) & (kf < 8 * np.pi)]])
white = float(pb[kb < 0.6].mean())
print(f"   white-noise level P_w(k->0) = {white:.3f} L^3")
print(f"   P_w(kL = 6.3) = {np.interp(6.3, kL, PL3):.4f} L^3")


def P_w(k):
    """foam contrast power [Mpc^3] at proper wavenumber k [1/Mpc]"""
    x = np.atleast_1d(k * L_CELL)
    out = np.interp(x, kL, PL3, left=white)
    hi = x > kL[-1]
    out[hi] = PL3[-1] * (kL[-1] / x[hi]) ** 2      # thin sheets: P ~ k^-2
    return out * L_CELL**3


pars = camb.set_params(H0=100 * H, ombh2=0.0224, omch2=0.120, ns=0.965,
                       As=2.1e-9, tau=0.054, mnu=0.06, lmax=2500,
                       NonLinear=camb.model.NonLinear_both)
res = camb.get_results(pars)
OM_M = (0.0224 + 0.120 + 0.06 / 93.14) / H**2
PK = camb.get_matter_power_interpolator(pars, nonlinear=True, hubble_units=False,
                                         k_hunit=False, kmax=30.0, zmax=1200)
H0C = 100 * H / 299792.458   # H0/c [1/Mpc]

print("\n2. dynamical-potential power, walls/matter, z = 0")
for k in (0.02, 0.05, 0.1, 0.2, 0.29):
    r = (XI_N * OM_W / OM_M) ** 2 * P_w(k)[0] / PK.P(0.0, k)
    print(f"   k = {k:4.2f}/Mpc  ratio = {r:.2f}")


def cl_ratio(zs, ells, reading):
    chis = res.comoving_radial_distance(zs)
    chi = np.linspace(1.0, chis - 1.0, 1500)
    z = res.redshift_at_comoving_radial_distance(chi)
    a = 1 / (1 + z)
    g = chi * (chis - chi) / chis
    out = []
    for ell in ells:
        k = (ell + 0.5) / chi
        Pw = a**-3 * P_w(k / a) if reading == "proper" else P_w(k)
        cm = np.trapezoid((g / a) ** 2 * PK.P(z, k, grid=False) / chi**2, chi)
        cw = np.trapezoid((ETA_L * OM_W / OM_M) ** 2 * (g * a**2) ** 2 * Pw / chi**2, chi)
        out.append(cw / cm)
    return out


ells = (30, 100, 300, 1000, 3000)
print("\n3. convergence power, walls/matter")
for zs in (0.6, 1.0, 1100.0):
    for reading in ("proper", "comoving"):
        r = cl_ratio(zs, ells, reading)
        print(f"   z_s = {zs:6g} {reading:8s} " +
              "  ".join(f"l={e}:{x:.4f}" for e, x in zip(ells, r)))


def cl_isw(ells_):
    chi = np.linspace(1.0, res.comoving_radial_distance(20.0), 3000)
    z = res.redshift_at_comoving_radial_distance(chi)
    a = 1 / (1 + z)
    Hz = res.hubble_parameter(z) / 299792.458
    out = []
    for ell in ells_:
        k = (ell + 0.5) / chi
        Pphi = 9 * H0C**4 * (ETA_L * OM_W) ** 2 * a**4 * P_w(k) / k**4
        out.append(np.trapezoid((2 * a * Hz) ** 2 * Pphi / chi**2, chi))
    return np.array(out)


TT = res.get_lensed_scalar_cls(CMB_unit="muK", raw_cl=True)[:, 0]
print("\n4. ISW power of the walls / Planck TT (comoving reading, Phi+Psi ~ a^2)")
ells_t = np.array([2, 5, 10, 30, 100, 300, 1000])
for e, c in zip(ells_t, cl_isw(ells_t)):
    print(f"   l = {e:4d}  ratio = {c * 2.7255e6**2 / TT[e]:.2e}")

print("\n5. stacked DeltaSigma around foam cells [M_sun/pc^2]")
lab, delta, h, ax = foam(16.0, 256, 7)
n_grid, box_L = 256, 16.0
nseed = int(round(box_L**3))
flat = lab.ravel()
vol = np.bincount(flat, minlength=nseed) * h**3
cen = np.zeros((nseed, 3))
for d_, coord in enumerate(np.meshgrid(ax, ax, ax, indexing="ij")):
    th = 2 * np.pi * coord.ravel() / box_L
    c = np.bincount(flat, np.cos(th), minlength=nseed)
    s = np.bincount(flat, np.sin(th), minlength=nseed)
    cen[:, d_] = (np.arctan2(s, c) % (2 * np.pi)) * box_L / (2 * np.pi)
del lab, flat
Req = (3 * vol / (4 * np.pi)) ** (1 / 3)
bins = np.linspace(0, 3.0, 61)
num = np.zeros(60)
den = np.zeros(60)
rng = np.random.default_rng(7)
offs = np.arange(-int(3.2 * Req.max() / h) - 1, int(3.2 * Req.max() / h) + 2)
OX, OY, OZ = np.meshgrid(offs, offs, offs, indexing="ij")
for j in rng.choice(nseed, size=600, replace=False):
    ci = cen[j] / h - 0.5
    b = np.round(ci).astype(int)
    r = np.sqrt(((b[0] + OX - ci[0]) ** 2 + (b[1] + OY - ci[1]) ** 2
                 + (b[2] + OZ - ci[2]) ** 2)) * h / Req[j]
    m = r < 3.0
    vals = delta[(b[0] + OX[m]) % n_grid, (b[1] + OY[m]) % n_grid, (b[2] + OZ[m]) % n_grid]
    ib = np.minimum((r[m] / 0.05).astype(int), 59)
    num += np.bincount(ib, vals, minlength=60)
    den += np.bincount(ib, minlength=60)
rc = 0.5 * (bins[1:] + bins[:-1])
prof = num / den


def delta_sigma(rr, d3, R):
    S = []
    for RR in R:
        s_ = np.linspace(0, np.sqrt(max(rr[-1] ** 2 - RR**2, 0)), 400)
        S.append(2 * np.trapezoid(np.interp(np.sqrt(RR**2 + s_**2), rr, d3), s_))
    S = np.array(S)
    Sbar = np.array([2 * np.trapezoid(np.interp(np.linspace(0, x, 200), R, S)
                                      * np.linspace(0, x, 200), np.linspace(0, x, 200)) / x**2
                     for x in R])
    return Sbar - S


def hsw(r, dc=-0.60, rs=1.05, al=2.1, be=9.1):
    return dc * (1 - (r / rs) ** al) / (1 + r**be)


R = np.linspace(0.05, 2.5, 50)
Rv = Req.mean() * L_CELL
rho_c = 2.775e11 * H**2                        # M_sun/Mpc^3
dS_w = ETA_L * OM_W * rho_c * Rv * delta_sigma(rc, prof, R) / 1e12
rr = np.linspace(0, 3.0, 600)
dS_m = OM_M * rho_c * Rv * delta_sigma(rr, hsw(rr), R) / 1e12
print(f"   mean cell radius R_v = {Rv:.1f} Mpc")
for x, a_, b_ in zip(R[::4], dS_w[::4], dS_m[::4]):
    print(f"   R/R_v = {x:4.2f}  walls {a_:+.3f}  matter(HSW) {b_:+.3f}")
iw = np.argmax(np.abs(dS_w))
print(f"   peak walls {dS_w[iw]:+.3f} at R/R_v = {R[iw]:.2f};"
      f" peak matter {dS_m[np.argmax(np.abs(dS_m))]:+.3f}")
