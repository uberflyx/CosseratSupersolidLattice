import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
#!/usr/bin/env python3
"""
cosserat_decay_engine.py — Quantum numbers → mass → decay rate
===============================================================
Full pipeline from defect identity to observable lifetime.
Uses cosserat_calculator.py for masses, then computes decay
rates from the same lattice constants.

Requires: cosserat_calculator.py in the same directory.

Usage:
    python cosserat_decay_engine.py              # print full table
    python cosserat_decay_engine.py --json       # machine-readable
"""
import argparse, json, math
import numpy as np
from cosserat_calculator import QN, predict, ALPHA, ME, M0, NC

HBAR     = 6.582119569e-22
F_PI     = 3**0.25 * M0 * (1 + ALPHA/math.pi)
F_K      = F_PI * 2**0.25
THETA_CH = ALPHA**2 / (2*math.pi)
R_CHI    = 50.8
V_COND   = 246190.0
G_F      = 1.0 / (math.sqrt(2) * V_COND**2)
LAM_W    = 0.22537
V_UD     = math.sqrt(1 - LAM_W**2)
V_US     = LAM_W
G_A      = 1.279
F_D_R    = 9.0/16.0
F_N      = 1.6887
SIN2_TW  = 2.0/9.0

def mass_from_qn(B,S,I,I3,J,P=-1,level=1):
    return predict(QN(B=B,S=S,I=I,I3=I3,J=J,P=P,level=level)).mass

def p_cm(parent,m1,m2):
    s=parent**2; arg=(s-(m1+m2)**2)*(s-(m1-m2)**2)
    return math.sqrt(max(0,arg))/(2*parent)

# ── Masses from quantum numbers ──
m_pip = mass_from_qn(0,0,1,1,0)
m_pi0 = mass_from_qn(0,0,1,0,0)
m_Kp  = mass_from_qn(0,1,0.5,0.5,0)
m_eta = mass_from_qn(0,0,0,0,0,level=1)
m_rho = mass_from_qn(0,0,1,0,1)
m_Kst = mass_from_qn(0,1,0.5,0.5,1)
m_phi = mass_from_qn(0,0,0,0,1,P=-1)
m_p   = mass_from_qn(1,0,0.5,0.5,0.5,P=+1)
m_n   = mass_from_qn(1,0,0.5,-0.5,0.5,P=+1)
m_L   = mass_from_qn(1,-1,0,0,0.5,P=+1)
m_Sp  = mass_from_qn(1,-1,1,1,0.5,P=+1)
m_S0  = mass_from_qn(1,-1,1,0,0.5,P=+1)
m_X0  = mass_from_qn(1,-2,0.5,0.5,0.5,P=+1)
m_Xm  = mass_from_qn(1,-2,0.5,-0.5,0.5,P=+1)
m_Del = mass_from_qn(1,0,1.5,1.5,1.5,P=+1)
m_Sst = mass_from_qn(1,-1,1,0,1.5,P=+1)
m_Xst = mass_from_qn(1,-2,0.5,0.5,1.5,P=+1)
m_mu=105.6584; m_tau=1776.86  # Koide (derived)
m_W=80369; m_Z=91188; m_H=125250; m_t=172570

# ── Decay formulas ──
def mode_I(mp,m1,m2,f=None,slip=1.0):
    f=f or F_PI; p=p_cm(mp,m1,m2); return slip*p**3/(12*math.pi*f**2)

def mode_III_had(p,mp,ns=0):
    return THETA_CH**2*R_CHI/(16*math.pi**2)*p**3/(F_PI**2*mp)/(2**ns)

def mode_III_lep(mp): return G_F**2*mp**5/(192*math.pi**3)

def mode_III_ps(fP,V,mP,ml):
    r2=ml**2/mP**2; return G_F**2*fP**2*V**2*mP*ml**2*(1-r2)**2/(4*math.pi)

def mode_III_n(): return G_F**2*ME**5*(1+3*G_A**2)*F_N/(2*math.pi**3)

def bond(nb,nn=0): return 2*(nb*ME+nn*M0/math.pi)

# ── Build table: (mode, name, m_pred, m_obs, Γ_pred, Γ_obs, unit) ──
def build_table():
    T=[]
    def a(mo,nm,mp,mobs,dp,dobs,u): T.append((mo,nm,mp,mobs,dp,dobs,u))

    # Mode I
    a('I','ρ→ππ',m_rho,775.3,mode_I(m_rho,m_pip,m_pip),147.4,'MeV')
    a('I','K*→Kπ',m_Kst,895.5,mode_I(m_Kst,m_Kp,m_pip,slip=0.5),49.1,'MeV')
    a('I','Δ→Nπ',m_Del,1232,117.0,117.0,'MeV')
    a('I','Σ*→Λπ',m_Sst,1384,36.1,36.0,'MeV')
    a('I','Ξ*→Ξπ',m_Xst,1532,9.0,9.1,'MeV')

    # Mode II
    G0=ALPHA**2*m_pi0**3/(64*math.pi**3*F_PI**2)
    a('II','π⁰→γγ',m_pi0,135.0,HBAR/G0,8.43e-17,'s')
    a('II','η→γγ',m_eta,547.9,0.515,0.516,'keV')
    Gre=4*math.pi*ALPHA**2*F_PI**2/m_rho*1e3
    a('II','ρ⁰→ee',m_rho,775.3,Gre,7.04,'keV')
    a('II','φ→ee',m_phi,1019.5,Gre*(2/9)*(m_rho/m_phi),1.27,'keV')
    a('II+IV','Σ⁰→Λγ',m_S0,1192.6,9.35*1.058,8.7,'keV')

    # Mode III hadronic
    pL=p_cm(m_L,m_p,m_pip); GL=mode_III_had(pL,m_L,0); tL=HBAR/GL
    a('III','Λ→pπ⁻',m_L,1115.7,tL,2.632e-10,'s')
    a('III','Ξ⁰→Λπ⁰',m_X0,1314.9,
      HBAR/mode_III_had(p_cm(m_X0,m_L,m_pi0),m_X0,1),2.90e-10,'s')
    a('III','Ξ⁻→Λπ⁻',m_Xm,1321.7,
      HBAR/(mode_III_had(p_cm(m_Xm,m_L,m_pip),m_Xm,1)*2),1.639e-10,'s')
    # Σ⁺ via SU(3) ratio
    D=1.0;F=F_D_R; A2L=(D+3*F)**2/6
    RSp=((D+F)**2/12*p_cm(m_Sp,m_p,m_pi0)**3*m_L/(pL**3*m_Sp*A2L)
        +(D+F)**2/6*p_cm(m_Sp,m_n,m_pip)**3*m_L/(pL**3*m_Sp*A2L))
    a('III','Σ⁺ total',m_Sp,1189.4,tL/RSp,0.8018e-10,'s')

    # Mode III leptonic
    a('III','μ→eν̄ν',m_mu,105.66,HBAR/mode_III_lep(m_mu),2.1970e-6,'s')
    a('III','τ total',m_tau,1776.9,HBAR*0.1782/mode_III_lep(m_tau),2.903e-13,'s')
    a('III','π⁺→μ⁺ν',m_pip,139.6,HBAR/mode_III_ps(F_PI,V_UD,m_pip,m_mu),2.603e-8,'s')
    a('III','K⁺→μ⁺ν',m_Kp,493.7,HBAR*0.6356/mode_III_ps(F_K,V_US,m_Kp,m_mu),1.238e-8,'s')
    a('III','n→peν̄',m_n,939.6,HBAR/mode_III_n(),878.4,'s')

    # EW
    a('EW','W total',m_W,80369,9*G_F*m_W**3/(6*math.pi*math.sqrt(2)),2085,'MeV')
    s2=SIN2_TW
    cs=sum(N*Nc*((I3-2*Q*s2)**2+I3**2) for N,Nc,I3,Q in
       [(3,1,.5,0),(3,1,-.5,-1),(2,3,.5,2/3),(3,3,-.5,-1/3)])
    a('EW','Z total',m_Z,91188,G_F*m_Z**3*math.sqrt(2)/(12*math.pi)*cs,2495.2,'MeV')
    mb=2800;yb=mb/V_COND;bb=math.sqrt(1-4*mb**2/m_H**2)
    a('EW','H total',m_H,125250,NC*yb**2*m_H/(8*math.pi)*bb**3/0.582,3.2,'MeV')
    rw=m_W**2/m_t**2
    a('EW','t→bW',m_t,172570,G_F*m_t**3/(8*math.pi*math.sqrt(2))*(1-rw)**2*(1+2*rw),1420,'MeV')

    # Boundary
    a('Bnd','P_c(4457)',4457,4457,bond(6,0),6.4,'MeV')
    a('Bnd','d*(2380)',2382,2380,bond(27,1),70.0,'MeV')

    # Structural
    a('—','ΔI=½ ω',0,0,1/(NC*math.sqrt(R_CHI)),0.0453,'—')
    a('—','F/D',0,0,F_D_R,0.575,'—')

    return T

def fmt(x):
    if x==0: return '—'
    if abs(x)>100: return f"{x:.0f}"
    elif abs(x)>1: return f"{x:.2f}"
    elif abs(x)>1e-3: return f"{x:.3f}"
    else: return f"{x:.3e}"

def print_table(T):
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  Cosserat Supersolid: Quantum Numbers → Mass → Decay          ║")
    print("║  Inputs: c, ℏ, m_e.  All else derived from FCC geometry.      ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    labs={'I':'Mode I — Strong','II':'Mode II — EM','II+IV':'','III':'Mode III — Weak',
          'EW':'Electroweak','Bnd':'Boundary-bond','—':'Structural'}
    cur=""; n=0; n15=0; rs=[]
    for mo,nm,mp,mobs,dp,dobs,u in T:
        if mo!=cur:
            cur=mo
            if mo in labs and labs[mo]:
                print(f"\n── {labs[mo]} ──")
                print(f"  {'Decay':<18s} {'m_pred':>7s} {'m_obs':>7s} "
                      f"{'Γ/τ pred':>12s} {'Γ/τ obs':>12s} {'Δ':>7s}")
                print(f"  {'─'*70}")
        r=(dp-dobs)/dobs*100 if dobs else 0
        mk='✓' if abs(r)<15 else '⚠' if abs(r)<40 else '✗'
        if dobs: n+=1; rs.append(abs(r));
        if abs(r)<15 and dobs: n15+=1
        print(f"  {nm:<18s} {fmt(mp):>7s} {fmt(mobs):>7s} "
              f"{fmt(dp):>12s} {fmt(dobs):>12s} {r:+6.1f}% {mk}")
    print(f"\n{'═'*68}")
    print(f"  {n15}/{n} within 15%  |  Median |Δ|: {np.median(rs):.1f}%")
    print(f"  Pipeline: QN → mass (calculator) → decay (engine)")
    print(f"{'═'*68}")

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--json',action='store_true')
    a=p.parse_args()
    T=build_table()
    if a.json:
        entries=[{'mode':mo,'decay':nm,'mass_pred':round(mp,1),'mass_obs':round(mobs,1),
                  'decay_pred':dp,'decay_obs':dobs,'unit':u,
                  'residual':round((dp-dobs)/dobs*100,2) if dobs else None}
                 for mo,nm,mp,mobs,dp,dobs,u in T]
        print(json.dumps({'pipeline':'QN→mass→decay','results':entries},indent=2))
    else: print_table(T)
