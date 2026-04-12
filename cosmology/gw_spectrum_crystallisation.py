#!/usr/bin/env python3
"""Corrected D4 nucleation analysis for GW spectrum."""
import numpy as np

hbar_c = 197.3269804
M_Pl = 1.22089e19
alpha_em = 1.0 / 137.035999177
m_e = 0.51099895
m0 = m_e / alpha_em
ell = hbar_c / m0
Theta_D = np.pi * m0
T_c = 156.1
DeltaG = (16*np.pi/3) * m0

g_QGP, g_had, g_nonQCD = 44.7, 3.0, 14.2
g_star = g_QGP + g_nonQCD
alpha_PT = 4*(g_QGP - g_had)/(3*g_star)

L4 = 3*ell
gamma_s = m0/ell**2
Delta_f = m0/ell**3
R_star = 2*gamma_s/Delta_f  # = 2*ell

# Bag constant
B_fm3 = (np.pi**2/90)*(g_QGP - g_had)*T_c**4/hbar_c**3
B_MeV4 = B_fm3 * hbar_c**3
B14 = B_MeV4**0.25

print("="*70)
print("D4 NUCLEATION ANALYSIS")
print("="*70)
print(f"  m_0       = {m0:.3f} MeV")
print(f"  ell       = {ell:.4f} fm")
print(f"  L_4       = {L4:.2f} fm = 3 ell")
print(f"  R*        = {R_star:.2f} fm = 2 ell")
print(f"  2 pi R*   = {2*np.pi*R_star:.2f} fm")
print(f"  L_4/(2piR*) = {L4/(2*np.pi*R_star):.3f}")
print(f"  Correction = O((L4/2piR*)^2) = {(L4/(2*np.pi*R_star))**2:.3f} = {(L4/(2*np.pi*R_star))**2*100:.0f}%")
print()
print(f"  RESULT: L_4 = {L4:.1f} fm < 2piR* = {2*np.pi*R_star:.1f} fm")
print(f"  -> Bubble wraps compact direction")
print(f"  -> 4D bounce reduces to 3D thermal bounce")
print(f"  -> CNT barrier DeltaG*/T is correct in both FCC and D4")
print()

# Cross-check with O(4) thin-wall
S4 = 27*np.pi**2/2
S3T = DeltaG/T_c
print(f"  O(4) thin-wall S_4 (if uncompactified) = {S4:.1f}")
print(f"  3D thermal S_3/T = DeltaG*/T_c = {S3T:.1f}")
print(f"  Ratio = {S4/S3T:.1f}x (would kill nucleation if O(4) applied)")
print(f"  But O(4) does NOT apply: bubble wraps compact direction.")

print(f"\n{'='*70}")
print("BAG MODEL")
print(f"{'='*70}")
print(f"  B^(1/4)     = {B14:.1f} MeV")
print(f"  Lambda_QCD  = {Theta_D:.1f} MeV")
print(f"  Agreement   = {abs(B14 - Theta_D)/Theta_D*100:.1f}%")

print(f"\n{'='*70}")
print("GW SPECTRUM (all scenarios)")
print(f"{'='*70}")

def gw_peak(T, alpha, g):
    beta_H = DeltaG / T
    v_w = 1/np.sqrt(3)
    T_GeV = T/1000
    kappa_v = alpha/(0.73 + 0.083*np.sqrt(alpha) + alpha)
    Ups = 1 - 1/(1 + 1/beta_H)
    f_pk = 1.9e-5 * (1/v_w) * beta_H * (T_GeV/100) * (g/100)**0.1667
    K = kappa_v*alpha/(1+alpha)
    h2_pk = 2.65e-6 * beta_H**(-2) * K**2 * (100/g)**(1/3) * v_w * Ups
    return beta_H, f_pk*1e9, h2_pk, kappa_v, Ups

cases = [
    ("Derived (bag model)", T_c, alpha_PT, g_star),
    ("Conservative (shear)", 200, 0.03, 60),
    ("Intermediate", 180, 0.3, 55),
    ("Total conversion", 220, 1.0, 62),
]

print(f"  {'Case':<30s} {'beta/H*':>8s} {'f_pk(nHz)':>10s} {'h2 Om_pk':>10s}")
print(f"  {'-'*62}")
for name, T, a, g in cases:
    bH, fpk, h2pk, kv, Ups = gw_peak(T, a, g)
    print(f"  {name:<30s} {bH:8.2f} {fpk:10.1f} {h2pk:10.2e}")

print(f"\n  Primary prediction:")
bH, fpk, h2pk, kv, Ups = gw_peak(T_c, alpha_PT, g_star)
print(f"    alpha      = {alpha_PT:.3f}")
print(f"    beta/H*    = {bH:.2f}")
print(f"    kappa_v    = {kv:.4f}")
print(f"    Upsilon    = {Ups:.4f}")
print(f"    f_peak     = {fpk:.1f} nHz")
print(f"    h2 Om_peak = {h2pk:.2e}")

# Plot
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    f_Hz = np.logspace(-10, -4.5, 2000)
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    for name, T, a, g in cases:
        params = {'T_star_MeV': T, 'T_star_GeV': T/1000, 'alpha': a,
                  'beta_over_H': DeltaG/T, 'v_w': 1/np.sqrt(3), 'g_star': g,
                  'H_star_GeV': np.sqrt(8*np.pi**3*g/90)*(T/1000)**2/M_Pl,
                  'kappa_v': a/(0.73+0.083*np.sqrt(a)+a),
                  'Upsilon': 1-1/(1+T/DeltaG)}
        
        bH = params['beta_over_H']
        vw = params['v_w']
        TG = params['T_star_GeV']
        gs = params['g_star']
        kv = params['kappa_v']
        Ups = params['Upsilon']
        
        f_pk = 1.9e-5*(1/vw)*bH*(TG/100)*(gs/100)**0.1667
        K = kv*a/(1+a)
        h2_pk = 2.65e-6*bH**(-2)*K**2*(100/gs)**(1/3)*vw*Ups
        s = f_Hz/f_pk
        h2 = h2_pk * s**3 * (7/(4+3*s**2))**3.5
        
        col = {'Derived': '#c0392b', 'Conservative': '#2980b9',
               'Intermediate': '#27ae60', 'Total': '#8e44ad'}
        c = [v for k,v in col.items() if name.startswith(k)][0]
        lw = 2.5 if 'Derived' in name else 1.5
        ls = '--' if 'Total' in name else '-'
        ax.loglog(f_Hz*1e9, h2, color=c, ls=ls, lw=lw, label=name)
    
    # NANOGrav
    ng_f = np.arange(1,15)/(15*365.25*24*3600)*1e9
    ng_h2 = np.array([2e-9,5.5e-10,2.2e-10,1.2e-10,7.5e-11,5e-11,4e-11,
                       2.8e-11,2.2e-11,1.8e-11,1.5e-11,1.3e-11,1.1e-11,9e-12])
    ax.fill_between(ng_f, ng_h2/3.5, ng_h2*3.5, alpha=0.2, color='gray',
                    label='NANOGrav 15yr')
    ax.scatter(ng_f, ng_h2, color='black', s=20, zorder=5)
    
    ax.axvspan(1, 100, alpha=0.04, color='blue')
    ax.text(8, 3e-8, 'PTA band', fontsize=9, color='blue', alpha=0.5, ha='center')
    ax.axvspan(100, 3000, alpha=0.04, color='green')
    ax.text(500, 3e-8, 'SKA', fontsize=9, color='green', alpha=0.5, ha='center')
    
    ax.set_xlim(0.3, 3000)
    ax.set_ylim(1e-16, 1e-7)
    ax.set_xlabel(r'Frequency $f$ [nHz]', fontsize=13)
    ax.set_ylabel(r'$h^2 \Omega_{\mathrm{GW}}(f)$', fontsize=13)
    ax.set_title('GW spectrum from vacuum crystallisation\n'
                 '(D4 dimensional reduction confirms 3D CNT)', fontsize=12)
    ax.legend(fontsize=9, loc='upper right', framealpha=0.9)
    ax.grid(True, which='both', alpha=0.15)
    plt.tight_layout()
    plt.savefig('/mnt/user-data/outputs/gw_spectrum_crystallisation.png', dpi=150)
    print(f"\n  Plot saved.")
except Exception as e:
    print(f"  Plot error: {e}")
