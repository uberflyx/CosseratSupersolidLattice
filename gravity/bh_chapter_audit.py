"""Independent audit of the black-holes chapter numbers.

Recomputes every headline quantity from CODATA constants, without importing
any framework script, so that script errors cannot propagate into the check.
Each block prints: claim in the tex, independent value, verdict.
"""

import numpy as np
from math import pi, log, sqrt, exp

# CODATA 2018/2022
hbar = 1.054571817e-34
h = 6.62607015e-34
c = 2.99792458e8
G = 6.67430e-11
kB = 1.380649e-23
me = 9.1093837015e-31
alpha = 7.2973525693e-3
sigT = 6.6524587321e-29
mp = 1.67262192369e-27
mN = 1.67492749804e-27  # neutron; nucleon ~1.674e-27
Msun = 1.98841e30
MeV = 1.602176634e-13
lP = sqrt(hbar * G / c**3)
MP = sqrt(hbar * c / G)

ell = 2.8179403262e-15          # r_e
m0 = me / alpha                  # node mass
m0c2 = m0 * c**2
aFCC = sqrt(2) * ell
rho_cr = 4 * m0 / aFCC**3        # crystal density

def check(name, claimed, computed, tol=0.05):
    ok = abs(computed / claimed - 1) < tol
    print(f"{'OK ' if ok else '***'} {name}: claimed {claimed:.4g}, computed {computed:.4g}"
          f"  (ratio {computed/claimed:.4f})")
    return ok

print("=" * 78)
print("A. CONSTANTS AND BASIC SCALES")
print("=" * 78)
check("m0 c^2 [MeV]", 70.0, m0c2 / MeV, 0.01)
check("rho_crystal [kg/m^3]", 7.9e15, rho_cr, 0.02)
check("lP [m]", 1.6e-35, lP, 0.02)
check("(ell/lP)^2", 3.045e40, (ell / lP)**2, 0.01)
alphaG = (lP / ell)**2
check("1/(4 alpha_G) [nats/cell]", 7.6e39, 1 / (4 * alphaG), 0.02)

print()
print("=" * 78)
print("B. ENTROPY SECTION (solar mass)")
print("=" * 78)
M = Msun
rs = 2 * G * M / c**2
kappa = c**2 / (2 * rs)
TH = hbar * c**3 / (8 * pi * G * M * kB)
A = 4 * pi * rs**2
S_mass = 4 * pi * G * M**2 / (hbar * c)
S_area = A / (4 * lP**2)
check("r_s [m]", 2954, rs, 0.01)
check("kappa [m/s^2]", 1.52e13, kappa, 0.01)
check("T_H [K]", 6.17e-8, TH, 0.01)
check("A [m^2]", 1.10e8, A, 0.01)
check("S/kB (mass form)", 1.049e77, S_mass, 0.01)
check("S/kB (area form)", 1.049e77, S_area, 0.01)
print(f"    identity check: mass vs area form agree to {abs(S_mass/S_area-1):.2e}")
check("A/ell^2 horizon cells", 1.4e37, A / ell**2, 0.03)
check("cells x nats/cell = S", 1.05e77, (A / ell**2) / (4 * alphaG), 0.01)

print()
print("=" * 78)
print("C. SCRAMBLING TABLE  (tex table values in [brackets])")
print("=" * 78)
for Ms, S_tab, ts_tab, ratio_tab in [(1, 7.0e66, 3.0e-3, 308),
                                     (30, 6.3e69, 95e-3, 321),
                                     (4e6, 1.1e80, 1.5e4, 369)]:
    M = Ms * Msun
    rs = 2 * G * M / c**2
    S = 4 * pi * G * M**2 / (hbar * c)
    ts = (2 * rs / c) * log(S)
    print(f"M = {Ms:g} Msun: S = {S:.3g}  [tex {S_tab:.2g}]  "
          f"WRONG by {S_tab/S:.2g}x ; t_s = {ts*1e3:.1f} ms [tex {ts_tab*1e3:.0f} ms]; "
          f"2lnS = {2*log(S):.0f} [tex {ratio_tab}]")
S_wrong = 4 * pi * G**2 * Msun**2 / (hbar * c)
print(f"    diagnosis: script's 4*pi*G^2*M^2/(hbar c) at 1 Msun = {S_wrong:.3g}"
      f"  -> reproduces the wrong tex value 7.0e66 exactly")

print()
print("=" * 78)
print("D. HAGEDORN SKIN")
print("=" * 78)
def skin(Ms):
    M = Ms * Msun
    rs = 2 * G * M / c**2
    TH = hbar * c**3 / (8 * pi * G * M * kB)
    N = sqrt(8 * pi / log(2)) * M / MP
    xi3 = 2 * log(N) * kB * TH / m0c2
    xi = xi3**(1/3)
    return 2 * rs * xi, xi, log(N)
rho1, xi1, lnN1 = skin(1)
check("ln N (1 Msun)", 89, lnN1, 0.02)
check("xi_* (1 Msun)", 2.4e-6, xi1, 0.03)
check("rho_* 1 Msun [m]", 1.4e-2, rho1, 0.05)
check("rho_* 30 Msun [m]", 0.14, skin(30)[0], 0.05)
check("rho_* SgrA* 4.3e6 [m]", 400, skin(4.3e6)[0], 0.15)
check("rho_* M87* 6.5e9 [m]", 5e4, skin(6.5e9)[0], 0.15)
# without entropy factor
M = Msun; rs = 2*G*M/c**2; TH = hbar*c**3/(8*pi*G*M*kB)
xi_no = (kB * TH / m0c2)**(1/3)
check("no-entropy depth [m]", 2.5e-3, 2 * rs * xi_no, 0.06)
check("eps0 x10 depth [m]", 6.5e-3, rho1 * 10**(-1/3), 0.06)
check("full-mixing depth factor", 0.18, (2 * lnN1)**(-1/3), 0.03)
check("deconfinement (solar)", 0.777, ((lnN1 - log(2*lnN1)) / (2*lnN1))**(1/3), 0.01)
zeta = ell / xi1
check("healing length at skin edge [m]", 1e-9, zeta, 0.20)

print()
print("=" * 78)
print("E. LINE BUDGETS AND PACKING")
print("=" * 78)
check("sqrt(8 pi/ln2)", 6.02, sqrt(8 * pi / log(2)), 0.005)
n30 = sqrt(8 * pi / log(2)) * 30 * Msun / MP
check("n (30 Msun)", 1.7e40, n30, 0.05)
nb = 30 * Msun / mN
check("baryon lines 30 Msun", 3.6e58, nb, 0.02)
check("ln Omega baryon tangle", 4e116, (log(2)/2) * nb**2, 0.15)
check("horizon exponent 30 Msun", 9.4e79, 4 * pi * (30 * Msun / MP)**2, 0.02)
spacing = (8 * pi * log(2))**0.25 * sqrt(30 * Msun / MP) * lP
check("puncture spacing 30 Msun [m]", 1.7e-15, spacing, 0.03)
Mstar = (1 / sqrt(8 * pi * log(2))) * (ell / lP)**2 * MP
check("M_* [Msun]", 80, Mstar / Msun, 0.02)
# geometric winding cap ~10 Msun: cap n=A/2 ell^2, need ln k <= ... floor when
# max winding r_s/ell reached: ln k = ln2 (M*/M)^2 and k_max ~ r_s/ell
def winding_floor():
    from scipy.optimize import brentq
    f = lambda Ms: log(2) * (Mstar / (Ms * Msun))**2 - log(2 * G * Ms * Msun / c**2 / ell)
    return brentq(f, 0.1, 79)
try:
    check("winding floor [Msun]", 10, winding_floor(), 0.5)
except Exception as e:
    print("    winding floor check skipped:", e)

print()
print("=" * 78)
print("F. ENDPOINT")
print("=" * 78)
Mmin = ell * c**2 / (2 * G)
check("M_min [kg]", 1.9e12, Mmin, 0.02)
check("kB T_star [MeV]", 5.57, (m0c2 / (4 * pi)) / MeV, 0.01)
check("T_star [K]", 6.5e10, m0c2 / (4 * pi * kB), 0.02)
check("Theta_D [K]", 2.55e12, pi * m0c2 / kB, 0.01)
check("Theta_D/T_star", 39.48, 4 * pi**2, 0.001)
check("exp(-2 m0c2/kB T_star)", 1e-11, exp(-8 * pi), 0.3)
# evaporation of floor-mass hole
t_ev = 5120 * pi * G**2 * Mmin**3 / (hbar * c**4)
tH = 4.35e17
frac = tH / (3 * t_ev)
check("fraction shed in Hubble time", 2.5e-4, frac, 0.1)
L = hbar * c**6 / (15360 * pi * G**2 * Mmin**2)
check("luminosity at floor [W]", 1e8, L, 0.3)
check("omega_max [Hz]", 3e23, pi * c / ell, 0.15)
check("relic entropy pi(ell/lP)^2", 9.5e40, pi * (ell / lP)**2, 0.02)
check("Page mass ~5e11 kg evaporates in tH", 5e11,
      (hbar * c**4 * tH / (5120 * pi * G**2))**(1/3), 0.5)

print()
print("=" * 78)
print("G. PAGE CURVE / FOKKER-PLANCK")
print("=" * 78)
check("Page fraction", 0.646, 1 - 1 / (2 * sqrt(2)), 0.001)
check("Zurek fraction", 0.568, 1 - (4/7)**1.5, 0.002)
check("(Mmin/Msun)^2", 1e-36, (Mmin / Msun)**2, 0.3)

print()
print("=" * 78)
print("H. VORTEX RESERVOIR / KNOT CHEMISTRY")
print("=" * 78)
kap = h / m0
check("circulation quantum [m^2/s]", 5.3e-6, kap, 0.02)
# Kerr 10 Msun a*=0.98
M = 10 * Msun; astar = 0.98
rg = G * M / c**2; rp = rg * (1 + sqrt(1 - astar**2)); ag = astar * rg
Om = astar * rg * c / (rp**2 + ag**2)
check("Omega_H [rad/s] (near-extremal 10 Msun)", 1e4, Om, 0.2)
nv = 2 * Om / kap
check("n_v [1/m^2]", 4e9, nv, 0.25)
check("N lines", 3e18, nv * pi * rp**2, 0.15)
check("b = sqrt(kappa/2 Om) [m]", 16e-6, sqrt(kap / (2 * Om)), 0.15)
check("E1 = (8 pi^2/5) m0c2 [GeV]", 1.107, (8 * pi**2 / 5) * m0c2 / MeV / 1e3, 0.005)
check("T_geom = m0c2/sqrt6 [MeV]", 28.6, m0c2 / sqrt(6) / MeV, 0.01)
check("2 E1 [GeV]", 2.21, 2 * (8 * pi**2 / 5) * m0c2 / MeV / 1e3, 0.01)
# bead ground radius: E(R)=E1 with rho_s = f_s m0/ell^3.  Need f_s; test f_s=1:
rho_s = m0 / ell**3
def Ering(R, xi=ell):
    return 0.5 * rho_s * kap**2 * R * (log(8 * R / xi) - 2)
from scipy.optimize import brentq
R0 = brentq(lambda R: Ering(R) - (8 * pi**2 / 5) * m0c2, 1.01 * ell, 100 * ell)
print(f"    bead R0 with f_s = 1: {R0/ell:.3f} ell   [tex: 1.68 ell]")
# try to infer f_s from the tex claim
for fs in (0.1, 0.2, 0.3, 0.5, 0.7, 1.0):
    rs_ = fs * m0 / ell**3
    E = lambda R: 0.5 * rs_ * kap**2 * R * (log(8 * R / ell) - 2)
    try:
        R0f = brentq(lambda R: E(R) - (8*pi**2/5)*m0c2, 1.001*ell, 1e4*ell)
        print(f"      f_s = {fs}: R0 = {R0f/ell:.3f} ell")
    except Exception:
        print(f"      f_s = {fs}: no root in range")
# ring E+Pc minimum
def EplusPc(R, fs=1.0):
    rs_ = fs * m0 / ell**3
    return 0.5*rs_*kap**2*R*(log(8*R/ell)-2) + rs_*kap*pi*R**2*c
from scipy.optimize import minimize_scalar
res = minimize_scalar(lambda R: EplusPc(R), bounds=(1.001*ell, 50*ell), method="bounded")
print(f"    min(E+Pc), f_s = 1: {res.fun/MeV/1e3:.2f} GeV at R = {res.x/ell:.2f} ell "
      f"[tex: 4.2 GeV]")
# ground ring impulse
Pg = rho_s * kap * pi * R0**2
print(f"    ground-ring impulse f_s = 1: {Pg*c/MeV/1e3:.2f} GeV/c  [tex: ~3 GeV/c]")

# knot masses
b16 = 16e-6
pref = rho_s * kap**2 / (4 * pi) * log(b16 / (2 * ell))
M_tref = pref * 32.74 * b16 / 2
print(f"    trefoil mass f_s = 1: {M_tref:.1f} J  [tex: 58 J]"
      f"  (fig-eight {pref*42.09*b16/2:.0f} J vs tex 75 J)")
check("trefoil J -> GeV", 3.6e11, 58 / MeV / 1e3, 0.02)

print()
print("=" * 78)
print("I. CHANNEL HORIZONS AND TIMING")
print("=" * 78)
check("second-sound trap r_s/x", 13, 3.65**2, 0.03)
D = 410 * 3.086e22
lead = (1 - 1/3.65) * D / c
check("GW150914 lead [yr]", 1e9, lead / 3.156e7, 0.35)
t_cross = 8.8e26 / (1e20 * c)
print(f"*** first sound crossing of observable universe: {t_cross*1e3:.0f} ms "
      f"(diameter) / {t_cross*1e3/2:.0f} ms (radius)  [tex: 'about a millisecond']")

print()
print("=" * 78)
print("J. INTERIOR GEOMETRY")
print("=" * 78)
check("Lambda_vac [1/m^2]", 1e-10, 8 * pi * G * rho_cr / c**2, 0.5)
check("r_a(vac) [m]", 1.43e5, sqrt(3 * c**2 / (8 * pi * G * rho_cr)), 0.02)
rhoDE = 6.6e-27
check("r_a(DE) [m]", 1.6e26, sqrt(3 * c**2 / (8 * pi * G * rhoDE)), 0.03)
check("r_a(DE) [Gpc]", 5, sqrt(3 * c**2 / (8 * pi * G * rhoDE)) / 3.086e25, 0.05)
rs1 = 2 * G * Msun / c**2
print(f"    (rs/ra)^4 solar: {(rs1/1.56e26)**4:.2g}  [tex: ~1e-92]  (is 1.3e-91)")
rho0 = 3 * c**2 / (8 * pi * G * rs1**2)
check("rho_0 uniform interior [kg/m^3]", 1.8e19, rho0, 0.05)
check("rho_0/rho_crystal", 2300, rho0 / rho_cr, 0.05)

print()
print("=" * 78)
print("K. POCKETS")
print("=" * 78)
rho_fl = 3 * rho_cr
check("rho_fluid [kg/m^3]", 2.4e16, rho_fl, 0.02)
Mmin_p = sqrt(3 * c**6 / (32 * pi * G**3 * rho_fl))
check("pocket floor [Msun]", 28, Mmin_p / Msun, 0.03)
check("pocket floor R_S [km]", 83, 2 * G * Mmin_p / c**2 / 1e3, 0.03)
check("upper floor (rho_crystal) [Msun]", 48, sqrt(3) * Mmin_p / Msun, 0.03)
# table rows
for Ms, RS_t, Rp_t in [(1e2, 2.95e5, 1.30e5), (1e6, 2.95e9, 2.80e6),
                       (1e9, 2.95e12, 2.80e7), (1e10, 2.95e13, 6.04e7)]:
    M = Ms * Msun
    RS = 2 * G * M / c**2
    Rp = (3 * M / (4 * pi * rho_fl))**(1/3)
    print(f"    M = {Ms:g}: R_S {RS:.3g} [tex {RS_t:.3g}], "
          f"R_pocket {Rp:.3g} [tex {Rp_t:.3g}]")
# heavy seed region size 1e5 Msun ~ 1300 km
Rp5 = (3 * 1e5 * Msun / (4 * pi * rho_fl))**(1/3)
check("1e5 Msun fluid region diameter [km]", 1300, 2 * Rp5 / 1e3, 0.10)
# Salpeter
tS = 0.1 / 0.9 * sigT * c / (4 * pi * G * mp)
check("Salpeter time e=0.1 [Myr]", 50, tS / 3.156e13, 0.03)
check("Edd lum timescale [Myr]", 450, sigT * c / (4 * pi * G * mp) / 3.156e13, 0.03)
check("e-folds in 750 Myr", 15, 750 / 50, 0.01)
check("growth 28 Msun x e^15 [Msun]", 1e8, 28 * exp(15), 0.15)
check("M_typical [1e9 Msun]", 1, 4.5e5 * 22**3 / 5 / 1e9, 0.1)

print()
print("=" * 78)
print("L. MERGER / OBSERVATIONAL")
print("=" * 78)
mu0 = rho_cr * c**2
rs30 = 2 * G * 30 * Msun / c**2
EGB = 3 * mu0 * ell * rs30**2
Mc2 = 30 * Msun * c**2
print(f"    E_GB = 3 mu0 ell rs^2 = {EGB:.3g} J = {EGB/Mc2:.2g} Mc^2  [tex ~1e-20 OK]")
mid = 6 * lP**2 / (ell * rs30)
print(f"*** middle identity 6 lP^2/(ell rs) = {mid:.3g}, vs actual ratio "
      f"{EGB/Mc2:.3g}: WRONG by {EGB/Mc2/mid:.2g}x")
print(f"    correct identity: E_GB/Mc^2 = 6 G rho ell rs/c^2 = "
      f"{6*G*rho_cr*ell*rs30/c**2:.3g}")
check("r_s/ell", 1e18, rs30 / ell, 3.2)

print()
print("=" * 78)
print("M. KERR ENTROPY CORRECTION CLAIMS")
print("=" * 78)
TH_kerr = hbar * c**3 / (4*pi*G*10*Msun*kB) * sqrt(1-0.98**2) / (1+sqrt(1-0.98**2))
gap_ratio = (8*pi**2/5)*m0c2 / (kB * TH_kerr)
print(f"    bead gap / kB T_H (10 Msun, a*=0.98) = {gap_ratio:.2g}  [tex: 6e21]")
A_kerr = 4*pi*((G*10*Msun/c**2*(1+sqrt(1-0.98**2)))**2 + (0.98*G*10*Msun/c**2)**2)
S_kerr = A_kerr / (4 * lP**2)
print(f"    S_BH (Kerr 10 Msun) = {S_kerr:.3g};  6e26/S = {6e26/S_kerr:.2g} [tex ~1e-52]")
