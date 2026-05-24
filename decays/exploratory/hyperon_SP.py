"""
[DIAGNOSTIC -- exploratory]
Cross-check that the measured s/p amplitudes obey the model-independent
DeltaI=1/2 / SU(3) relations (Lee-Sugawara, Sigma triangle, S(Sigma+ -> n pi+)=0).
These relations are the symmetry backbone the lattice reproduces; the absolute
widths are computed by the engine master formula (hyperon_master_width).
"""
"""
Hyperon nonleptonic decays: isospin + SU(3) structure from the lattice.
Step 1 -- show the framework's residual symmetries reproduce the model-independent
relations among the seven S- and P-wave amplitudes (Lee-Sugawara, the Sigma triangle,
the DeltaI=1/2 ratios, and the S(Sigma+ -> n pi+) ~ 0 selection rule).
Experimental amplitudes: Salone-Alvarado-Leupold-Kupsc (2026), Table III, last column,
in units of 10^-7 (G_F m_pi^2 convention).
"""
import numpy as np

# channel: (S_exp, P_exp)
exp = {
 'Sig+_npi+': ( 0.06,  1.81),   # Sigma+ -> n pi+
 'Sig+_ppi0': (-1.43,  1.17),   # Sigma+ -> p pi0
 'Sig-_npi-': ( 1.88, -0.06),   # Sigma- -> n pi-
 'Lam_ppi-' : ( 1.42,  0.52),   # Lambda -> p pi-
 'Lam_npi0' : (-1.04, -0.37),   # Lambda -> n pi0  (from DeltaI=1/2: S=-S(ppi-)/sqrt2)
 'Xi-_Lpi-' : (-1.98,  0.48),   # Xi- -> Lambda pi-
 'Xi0_Lpi0' : ( 1.52, -0.33),   # Xi0 -> Lambda pi0
}
S = {k:v[0] for k,v in exp.items()}
P = {k:v[1] for k,v in exp.items()}

print("="*64)
print("Model-independent symmetry relations (must hold for ANY delta I=1/2 theory)")
print("="*64)

# 1. DeltaI=1/2 for Lambda:  A(Lam->p pi-) = -sqrt(2) A(Lam->n pi0)
for amp,nm in ((S,'S'),(P,'P')):
    lhs, rhs = amp['Lam_ppi-'], -np.sqrt(2)*amp['Lam_npi0']
    print(f"  Lambda DeltaI=1/2 [{nm}]:  A(p pi-)={lhs:+.3f}   -sqrt2 A(n pi0)={rhs:+.3f}")

# 2. DeltaI=1/2 for Xi:  A(Xi- -> L pi-) = -sqrt(2) A(Xi0 -> L pi0)   (pure)
for amp,nm in ((S,'S'),(P,'P')):
    r = amp['Xi-_Lpi-']/amp['Xi0_Lpi0']
    print(f"  Xi ratio [{nm}]: A(Xi-)/A(Xi0) = {r:+.3f}   (pure DeltaI=1/2 = -sqrt2 = -1.414)")

# 3. Sigma triangle: A(Sig+ -> n pi+) + sqrt(2) A(Sig+ -> p pi0) = A(Sig- -> n pi-)
for amp,nm in ((S,'S'),(P,'P')):
    lhs = amp['Sig+_npi+'] + np.sqrt(2)*amp['Sig+_ppi0']
    rhs = amp['Sig-_npi-']
    print(f"  Sigma triangle [{nm}]: A(npi+)+sqrt2 A(ppi0)={lhs:+.3f}   A(Sig- npi-)={rhs:+.3f}")

# 4. Lee-Sugawara (S-wave): 2 S(Xi- -> L pi-) + S(Lam -> p pi-) = sqrt(3) S(Sig+ -> p pi0)
ls_lhs = 2*S['Xi-_Lpi-'] + S['Lam_ppi-']
ls_rhs = np.sqrt(3)*S['Sig+_ppi0']
print(f"\n  Lee-Sugawara [S]: 2S(Xi-)+S(Lam) = {ls_lhs:+.3f}   sqrt3 S(Sig+ ppi0) = {ls_rhs:+.3f}")

# 5. The S-wave selection rule
print(f"\n  S(Sigma+ -> n pi+) = {S['Sig+_npi+']:+.3f}   (current-algebra commutator -> exactly 0)")
print(f"  P(Sigma- -> n pi-) = {P['Sig-_npi-']:+.3f}   (vanishing P-wave)")
