#!/usr/bin/env python3
"""
cosserat_dw.py - Cosserat Supersolid Debye-Waller Analysis
============================================================
Companion script to "The Cosserat Supersolid" monograph (M. A. Cox, 2025).
Reproduces all numerical results for Lorentz.

Run: python3 cosserat_dw.py
"""
import numpy as np
from scipy.optimize import brentq

# --- Fundamental constants (SI) ---
c, hbar, m_e = 2.99792458e8, 1.054571817e-34, 9.1093837015e-31
eV, MeV = 1.602176634e-19, 1.602176634e-13
alpha = 1.0/137.035999177

# --- Derived lattice parameters (Monograph Sec. 2) ---
m0 = m_e/alpha; m0_MeV = m0*c**2/MeV; ell = hbar/(m0*c)

# --- Reciprocal lattices: FCC slice vs the in-plane D_4* sublattice ----------
# The photon's elastic spatial-to-spatial Bragg vector is NOT the standalone FCC
# minimum.  The four-momentum transfer must lie in D_4* AND stay in the spatial
# subspace (G.e_c = 0 with e_c = (1,1,1,1)/2).  That leaves the six half-integer
# vectors with two +1/2 and two -1/2 components, of coordinate length 1, i.e.
#   |G_min| = 2 pi sqrt(2) / ell = 8.886/ell,
# larger than the FCC |G_111| = pi sqrt(6)/ell = 7.695/ell by exactly 2/sqrt(3).
# G_111 is still the right vector for the compact-integrated density the photon
# sees (Mathieu band structure, Leggett bound); it is the wrong one for Bragg.
G111_ell = np.pi*np.sqrt(6); G111_sq_ell2 = G111_ell**2      # =6pi^2, FCC slice
Gmin_ell = 2*np.pi*np.sqrt(2); Gmin_sq_ell2 = Gmin_ell**2    # =8pi^2, in-plane D_4*

# --- FCC primitive cell -------------------------------------------------------
# One node per primitive cell of volume V_p = ell^3/sqrt(2) for nearest-neighbour
# distance ell.  Omitting the sqrt(2) (i.e. assuming a simple-cubic cell) inflates
# <u^2> by 2^(1/6) = 12%, because k_D and the mode-density weight both carry it.
Vp = 1.0/np.sqrt(2)
kD_ell = (6*np.pi**2/Vp)**(1/3)                              # = 4.375/ell

# --- PN core half-width and its form factor ----------------------------------
# w = (pi/4) d with d = ell/sqrt(3) the misfit period, so w = 0.4534 ell.  NOT
# (pi/4) ell: that is w/d expressed as though d were ell, and it overstates the
# core by sqrt(3), collapsing |F|^2 by two decades.
w_ell = (np.pi/4)/np.sqrt(3)
FF_G = np.exp(-2*Gmin_ell*w_ell)                             # at the Bragg vector

# --- Zone boundary and bounds ---
E_ZB = np.pi*np.sqrt(2)*m0_MeV                               # = |G_min| hbar c/2 = 311 MeV
fermi_dv = 7e-19; grb_dt_bound = 0.86; grb_D_c = 4e17

# --- Cosserat parameters (Born-Huang, Monograph Ch. 5) ---
N2 = 1/np.pi
mu_d = (1-N2)/(1+N2); kp_d = 1-mu_d; gm_d = mu_d; J_d = 0.1
cL = np.sqrt(3*mu_d); omega0 = np.sqrt(2*kp_d/J_d)

print("="*70)
print("COSSERAT SUPERSOLID DEBYE-WALLER ANALYSIS")
print("="*70)
print(f"\nm0={m0_MeV:.2f} MeV, ell={ell:.3e} m, N2=1/pi={N2:.4f}")
print(f"cL/c={cL:.4f}, omega0={omega0:.3f} c/ell, omegaD={kD_ell:.3f} c/ell")
print(f"|F(G_min)|^2 = {FF_G:.2e}   (in-plane D_4*, |G_min|={Gmin_ell:.4f}/ell, w={w_ell:.4f} ell)")

def transverse_dispersion(k):
    A = (mu_d+kp_d)*k**2; B = gm_d*k**2+2*kp_d
    disc = max((B+J_d*A)**2-4*J_d*(A*B-kp_d**2*k**2), 0)
    w2_ac = ((B+J_d*A)-np.sqrt(disc))/(2*J_d)
    w2_op = ((B+J_d*A)+np.sqrt(disc))/(2*J_d)
    for w2 in [w2_ac, w2_op]:
        if k < 1e-12: yield (w2, 1.0 if w2 < 1 else 0.0)
        else:
            R2 = (A-w2)**2/(kp_d**2*k**2)
            yield (w2, 1.0/(1.0+R2*J_d))

def compute_u2(f_s=0.0, n=2000):
    ks = np.linspace(1e-8, kD_ell, n); dk = ks[1]-ks[0]
    pf = Vp/(2*np.pi**2); h_ms = 1/(1-f_s)
    u2_LA = u2_TA = u2_TO = 0.0
    for k in ks:
        u2_LA += pf*h_ms/(2*cL*k)*k**2*dk
        br = list(transverse_dispersion(k))
        u2_TA += 2*pf*h_ms*br[0][1]/(2*np.sqrt(max(br[0][0],1e-30)))*k**2*dk
        u2_TO += 2*pf*h_ms*br[1][1]/(2*np.sqrt(max(br[1][0],1e-30)))*k**2*dk
    return u2_LA, u2_TA, u2_TO

# --- Simple Debye vs Cosserat ---
u2_debye = 3*Vp**(1/3)*(6*np.pi**2)**(2/3)/(8*np.pi**2)   # = 0.514 ell^2
u2_L, u2_Ta, u2_To = compute_u2(0)
u2_coss = u2_L+u2_Ta+u2_To
print(f"\nDebye (3 branch): <u^2> = {u2_debye:.4f} ell^2")
print(f"Cosserat (6 branch): <u^2> = {u2_coss:.4f} ell^2 ({(1-u2_coss/u2_debye)*100:.1f}% reduction)")
print(f"  LA={u2_L:.4f}, TA_ac={u2_Ta:.4f}, TO_opt={u2_To:.4f}")

# --- NCRI scan ---
ratio2 = (E_ZB/10000)**2
print(f"\n{'='*70}\nNCRI SCAN\n{'='*70}")
print(f"{'fs':>5} {'<u2>/l2':>8} {'Lind':>6} {'2W':>6} {'e-2W':>10} {'S_tot':>10} {'dv/c@10GeV':>12} {'stat':>5}")
print("-"*72)
for fs in [0,.2,.4,.5,.55,.58,.60,.62,.65,.70]:
    u2 = sum(compute_u2(fs)); L = np.sqrt(u2)
    tW = u2*Gmin_sq_ell2/3; DW = np.exp(-tW); S = FF_G*DW; dv = S*ratio2
    print(f"{fs:5.2f} {u2:8.3f} {L:6.3f} {tW:6.1f} {DW:10.2e} {S:10.2e} {dv:12.2e} {'PASS' if dv<fermi_dv else 'FAIL':>5}")

# --- GRB time delay -----------------------------------------------------------
# Uses the monograph's own dispersion ansatz and cosmological delay integral
# rather than a single-energy estimate:
#     dv_g/c = n_shell S (E E_BZ)^2 / (E^4 + E_BZ^4)
#     dt     = int_0^z [dv(E1(1+z')) - dv(E2(1+z'))] / [(1+z') H(z')] dz'
# Both pieces matter.  n_shell = 6 counts the in-plane D_4* first shell, and the
# redshift is not optional: below the zone boundary the anomaly grows as E^2, so
# the source-frame energies (1+z = 1.903) carry several times the observed value.
n_shell = 6
z_grb, E1_obs, E2_obs = 0.903, 31e3, 100.0        # MeV; the Fermi-LAT photon pair
H0_inv_s = 3.0857e19/67.4                          # 1/H0 in seconds
Om, OL = 0.315, 0.685


def shape(E):
    return (E*E_ZB)**2/(E**4 + E_ZB**4)


def grb_delay(S, n=4000):
    zs = np.linspace(0.0, z_grb, n); dz = zs[1]-zs[0]; tot = 0.0
    for zp in zs:
        Hs = np.sqrt(Om*(1+zp)**3 + OL)
        dv1 = n_shell*S*shape(E1_obs*(1+zp))
        dv2 = n_shell*S*shape(E2_obs*(1+zp))
        tot += (dv1-dv2)/((1+zp)*Hs)*dz
    return abs(tot)*H0_inv_s


print(f"\n{'='*70}\nGRB 090510 TIME DELAY\n{'='*70}")
print(f"{'fs':>5} {'S_tot':>10} {'dt(s)':>8} {'stat':>5}")
print("-"*35)
for fs in [.55,.58,.60,.62,.65,.70]:
    u2 = sum(compute_u2(fs)); S = FF_G*np.exp(-u2*Gmin_sq_ell2/3)
    dt = grb_delay(S)
    print(f"{fs:5.2f} {S:10.2e} {dt:10.3g} {'PASS' if dt<grb_dt_bound else 'FAIL':>5}")

# --- Critical f_s ---
def S_of(fs): return FF_G*np.exp(-sum(compute_u2(fs))*Gmin_sq_ell2/3)
def dv_res(fs): return n_shell*S_of(fs)*shape(1e4)-fermi_dv
def dt_res(fs): return grb_delay(S_of(fs))-grb_dt_bound
fs_v = brentq(dv_res,.3,.8); fs_t = brentq(dt_res,.4,.8)
print(f"\nCritical fs (velocity): {fs_v:.4f}")
print(f"Critical fs (GRB dt):  {fs_t:.4f} <-- binding")

# --- Self-consistency ---
print(f"\n{'='*70}\nSELF-CONSISTENCY: NN WAVEFUNCTION OVERLAP\n{'='*70}")
for fs in [0, .5, fs_t, .65, .80]:
    u2 = sum(compute_u2(fs)); sig2 = u2/3; S = np.exp(-1/(4*sig2))
    print(f"  fs={fs:.3f}: <u2>={u2:.3f}, overlap={S:.3f}, Z*S={12*S:.1f}")

# --- Reference table at f_s = 0.60 ---
fs_r = 0.60; u2_r = sum(compute_u2(fs_r))
tW_r = u2_r*Gmin_sq_ell2/3; DW_r = np.exp(-tW_r); S_r = FF_G*DW_r
print(f"\n{'='*70}\nREFERENCE TABLE (fs={fs_r})\n{'='*70}")
print(f"  <u2>_eff    = {u2_r:.3f} ell^2")
print(f"  Lindemann   = {np.sqrt(u2_r):.3f}")
print(f"  2W          = {tW_r:.1f}")
print(f"  e^(-2W)     = {DW_r:.2e}")
print(f"  |F(G)|^2    = {FF_G:.2e}")
print(f"  S_total     = {S_r:.2e}")
print(f"  Zone gap    = {2*np.sqrt(S_r)*E_ZB*1e6:.1f} ueV")
print(f"  dv/c(10GeV) = {n_shell*S_r*shape(1e4):.1e}  (bound: 7e-19)")
print(f"  dv/c(100MeV)= {n_shell*S_r*shape(100.0):.1e}")
print(f"  GRB dt      = {grb_delay(S_r):.3g} s  (bound: 0.86 s)")
