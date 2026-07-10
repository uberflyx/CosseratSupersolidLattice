"""
sn_wall_residual.py
===================
The pinned-wall velocity template as a low-redshift supernova systematic
(Confrontation chapter, dark-energy discussion).

The DESI preference for evolving dark energy is driven by the low-redshift
supernova magnitudes (Efstathiou 2025; Vincenzi et al. 2025 in rebuttal).
The fossil wall network contributes an irreducible component to that budget:
a static, NON-COMOVING velocity template that comoving peculiar-velocity
reconstructions cannot remove, sampled at wall-resident positions because
supernova hosts (galaxies) and the observer both sit on the web.

Result: a coherent, sky-averaged magnitude offset per observer of rms
0.020 mag at z = 0.015, 0.010 at z = 0.025, falling to 0.004 by z = 0.04
(kinematic upper estimate; sign depends on the observer's cell). Up to half
the disputed ~0.04 mag offset at the lowest redshifts, decaying faster with
z than the anomaly; a derived floor to the systematic budget, with a
distinguishing signature: the residuals correlate with the local web and are
epoch-independent.

Model: periodic Voronoi foam, DE on walls, perturbation Poisson -> velocity
template v = a*t_H (upper estimate). Observers ON walls; SN hosts ON walls
(galaxies live on the web). Per z-bin: coherent sky-averaged line-of-sight
velocity difference -> delta_mu(z) = (5/ln10) * <dv_los>/(cz).
Report per-observer coherent offset (mean over sky) distribution across
observers, plus the raw scatter, for both coupling signs.
"""
import numpy as np
from scipy.spatial import cKDTree

rng = np.random.default_rng(11)
G, Mpc = 6.674e-11, 3.086e22
H0 = 67.4e3/Mpc; t_H = 1/H0; c_ms = 2.998e8
rho_L = 0.685 * 3*H0**2/(8*np.pi*G)

L_box, N = 320.0*Mpc, 192
dx = L_box/N
L_cell = 21.4*Mpc
n_seed = int(round((L_box/L_cell)**3))
seeds = rng.uniform(0, L_box, (n_seed,3))
tree = cKDTree(seeds, boxsize=L_box)
ax1 = (np.arange(N)+0.5)*dx
X,Y,Z = np.meshgrid(ax1,ax1,ax1,indexing="ij")
lab = tree.query(np.stack([X.ravel(),Y.ravel(),Z.ravel()],1), workers=-1)[1].reshape(N,N,N).astype(np.int32)
del X,Y,Z
wall = np.zeros((N,N,N),bool)
for a_ in range(3): wall |= (lab != np.roll(lab,1,axis=a_))
del lab
print(f"wall fraction {wall.mean():.3f}")

drho = np.where(wall, rho_L/wall.mean(), 0.0) - rho_L
k1 = 2*np.pi*np.fft.fftfreq(N, d=dx)
KX,KY,KZ = np.meshgrid(k1,k1,k1,indexing="ij")
K2 = KX**2+KY**2+KZ**2; K2[0,0,0]=1
phik = -4*np.pi*G*np.fft.fftn(drho)/K2; phik[0,0,0]=0
v = [np.real(np.fft.ifftn(-1j*Ki*phik))*t_H for Ki in (KX,KY,KZ)]
del drho,phik,KX,KY,KZ,K2
print(f"rms point speed {np.sqrt((v[0]**2+v[1]**2+v[2]**2).mean())/1e3:.0f} km/s (upper estimate)")

# wall-cell coordinate list for host sampling
wi = np.argwhere(wall)                      # indices of wall cells
wall_pos = (wi+0.5)*dx
wtree = cKDTree(wall_pos, boxsize=L_box)
del wall

def vel_at(idx):
    return np.array([v[0][tuple(idx)], v[1][tuple(idx)], v[2][tuple(idx)]])

zbins = np.array([0.015,0.025,0.04,0.06,0.09,0.13])
n_obs, n_sn = 40, 600
res = {s: np.zeros((n_obs,len(zbins))) for s in (+1,-1)}
for o in range(n_obs):
    oi = wi[rng.integers(len(wi))]          # observer ON a wall
    opos = (oi+0.5)*dx
    vo = vel_at(oi)
    for bz, zc in enumerate(zbins):
        r = zc*c_ms/H0
        if r > L_box*0.75: r = L_box*0.75   # cap inside toy box (wraps)
        u = rng.normal(size=(n_sn,3)); u /= np.linalg.norm(u,axis=1)[:,None]
        target = (opos + u*r) % L_box
        # snap SN host to nearest wall cell (hosts live on the web)
        _, j = wtree.query(target, workers=-1)
        hp = wall_pos[j]
        hv = np.stack([v[0][wi[j][:,0],wi[j][:,1],wi[j][:,2]],
                       v[1][wi[j][:,0],wi[j][:,1],wi[j][:,2]],
                       v[2][wi[j][:,0],wi[j][:,1],wi[j][:,2]]],1)
        # line of sight from observer to host (minimum image)
        d = hp - opos; d -= L_box*np.round(d/L_box)
        rr = np.linalg.norm(d,axis=1); rhat = d/rr[:,None]
        dv = np.einsum("ij,ij->i", hv - vo, rhat)      # m/s
        res[+1][o,bz] = dv.mean()
        res[-1][o,bz] = -dv.mean()
print("\nper-observer COHERENT sky-mean residual, attractive sign:")
print("z_bin   <dv_los> km/s (mean+-rms over observers)   delta_mu [mag] rms")
for bz, zc in enumerate(zbins):
    m = res[+1][:,bz]/1e3
    dmu = (5/np.log(10))*(res[+1][:,bz]/ (zc*c_ms))
    print(f" {zc:5.3f}   {m.mean():+7.1f} +- {m.std():6.1f}        {np.abs(dmu).mean():.4f} (rms {dmu.std():.4f})")
np.savez("/tmp/sn_wall.npz", zbins=zbins, res=res[+1])
