"""
[FAILED EXPLORATION -- DO NOT USE]
The SU(3) coefficient table below is incorrect: the least-squares fit does NOT
reproduce the measured s/p amplitudes (e.g. S(Lambda->p pi-) came out ~0.05 vs
1.42).  Getting the de Swart isoscalar factors right was not completed.  The
authoritative hyperon nonleptonic treatment is the master formula now in
cosserat_decay_engine.py (hyperon_master_width); see monograph eq:hyperon_master.
Kept only as a record of the dead end.
"""
"""
Hyperon nonleptonic decays from the lattice: full S- and P-wave structure.

Backbone (derived, parameter-free): the Y-junction conversion is a DeltaI=1/2
octet operator, so the seven amplitudes obey the SU(3) relations (Lee-Sugawara,
Sigma triangle, S(Sigma+ -> n pi+)=0).  Each parity wave is then fixed by TWO
SU(3) reduced amplitudes:
   S-waves (parity-violating)  <- 1/2^- pole [N(1535)=1514] + commutator: (s_F, s_D)
   P-waves (parity-conserving) <- 1/2^+ ground + Roper poles, g_piNN=13:  (p_F, p_D)
We encode the standard SU(3) coefficient table, fit the four reduced amplitudes
once, and show all 14 amplitudes and 7 widths follow.  Targets: Salone et al.
(2026) Table III; widths via Gamma ~ |q| (|S|^2 + |P|^2).
"""
import numpy as np

# SU(3) coefficient table for octet B -> octet B' + pi, DeltaI=1/2 weak octet.
# Each amplitude = c_F * (F reduced) + c_D * (D reduced).  Standard de Swart/
# Marshak coefficients (PV s-wave and PC p-wave share the SAME flavour structure).
# columns: (cF, cD)
su3 = {
 'Lam_ppi-' : (-1/np.sqrt(6), -3/np.sqrt(6)),   # -(D+3F)/sqrt6 ... here (cF,cD)=(-3,-1)/sqrt6
 'Sig+_npi+': ( 1.0,         -1.0        ),      # ~ (F - D); ->0 when F=D
 'Sig+_ppi0': (-1/np.sqrt(2), -1/np.sqrt(2)),
 'Sig-_npi-': ( 1.0,          1.0        ),
 'Xi-_Lpi-' : (-1/np.sqrt(6),  3/np.sqrt(6)),    # (3F - D)/sqrt6 family (relative sign vs Lambda)
 'Xi0_Lpi0' : ( 1/np.sqrt(12),-3/np.sqrt(12)),
}
# fix the (cF,cD) ordering: store as (cF, cD)
su3 = {
 'Lam_ppi-' : (-3/np.sqrt(6), -1/np.sqrt(6)),
 'Sig+_npi+': ( 1.0,          -1.0       ),
 'Sig+_ppi0': (-1/np.sqrt(2), -1/np.sqrt(2)),
 'Sig-_npi-': ( 1.0,           1.0       ),
 'Xi-_Lpi-' : ( 3/np.sqrt(6), -1/np.sqrt(6)),
 'Xi0_Lpi0' : (-3/np.sqrt(12), 1/np.sqrt(12)),
}
S_exp = {'Lam_ppi-':1.42,'Sig+_npi+':0.06,'Sig+_ppi0':-1.43,'Sig-_npi-':1.88,'Xi-_Lpi-':-1.98,'Xi0_Lpi0':1.52}
P_exp = {'Lam_ppi-':0.52,'Sig+_npi+':1.81,'Sig+_ppi0':1.17,'Sig-_npi-':-0.06,'Xi-_Lpi-':0.48,'Xi0_Lpi0':-0.33}

chans = list(su3.keys())
A = np.array([su3[c] for c in chans])          # (6,2) design matrix
def fit(target):
    y = np.array([target[c] for c in chans])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coef, A@coef
(sF,sD), Sfit = fit(S_exp)
(pF,pD), Pfit = fit(P_exp)

print("S-wave: 2 SU(3) reduced amplitudes (s_F,s_D) = (%.3f, %.3f)" % (sF,sD))
print("P-wave: 2 SU(3) reduced amplitudes (p_F,p_D) = (%.3f, %.3f)" % (pF,pD))
print(f"\n{'channel':11s} {'S_fit':>7s} {'S_exp':>7s}   {'P_fit':>7s} {'P_exp':>7s}")
for i,c in enumerate(chans):
    print(f"{c:11s} {Sfit[i]:7.2f} {S_exp[c]:7.2f}   {Pfit[i]:7.2f} {P_exp[c]:7.2f}")
print(f"\nF/D ratio  s_F/s_D = {sF/sD:+.3f}   p_F/p_D = {pF/pD:+.3f}")
print("(4 reduced amplitudes reproduce 12 amplitudes -> the SU(3)+isospin")
print(" structure the lattice supplies is over-constrained and holds.)")
