"""Independent audit of the spectral-mass-formula chapter: master-formula
arithmetic and PDG confrontations.

Recomputes every quoted closure m = N*m0 - N*(4 - lambda)*m_e from CODATA
constants alone (m0 = m_e/alpha), and every quoted PDG deviation against
PDG 2024 values read from the project's PDG booklet. No framework scripts
are imported. Eigenvalues are taken as quoted here; their independent
recomputation from the cluster geometry is the second audit stage.
"""
from math import sqrt

# CODATA 2022
me = 0.51099895069        # MeV
alpha = 7.2973525643e-3
m0 = me / alpha
print(f"m0 = {m0:.4f} MeV   (tex: 70.025)\n")

# PDG 2024 (from project pdg.pdf booklet text)
PDG = dict(pi_pm=139.57039, pi0=134.9768, eta=547.862, etap=957.78,
           K_pm=493.677, Kstar_pm=891.67, rho=775.26, omega=782.66,
           phi=1019.461, f2=1275.4, a2=1318.2, f2p=1517.3, K2star=1427.3,
           eta1295=1294.0, p=938.27208816, n=939.5654205)
pi_iso = (2*PDG['pi_pm'] + PDG['pi0'])/3
nuc_iso = (PDG['p'] + PDG['n'])/2

def m(N, lam): return N*m0 - N*(4.0 - lam)*me

rows = [
 # (state, N, lambda, tex mass, PDG value, tex-quoted deviation note)
 ("pi (iso)",        2, 2.0,            138.007, pi_iso,   "-0.024%"),
 ("eta",             8, 0.8229,         547.21,  PDG['eta'], "-0.12%"),
 ("eta'",           14, 1.000,          958.89,  PDG['etap'], "+0.12%"),
 ("K (pm)",          7, 5.0,            493.75,  PDG['K_pm'], "+0.016%"),
 ("K*(892)",        13, 1.6626,         894.80,  PDG['Kstar_pm'], "+0.35%"),
 ("rho(770)",       11, 4.891,          775.29,  PDG['rho'], ""),
 ("omega(782)",     11, 6.241,          782.88,  PDG['omega'], ""),
 ("phi(1020)",      14, 8.441,          1012.13, PDG['phi'], ""),
 ("eta(1295)",      18, (11+sqrt(17))/2,1293.21, PDG['eta1295'], "-0.06%"),
 ("f2(1270)",       18, 5.580,          1274.99, PDG['f2'], ""),
 ("a2(1320)",       18, 10.480,         1320.06, PDG['a2'], "+0.14%"),
 ("f2 Eg (pred)",   18, 4.844,          1268.22, None, "forward"),
 ("f2'(1525)",      21, 8.554,          1519.39, PDG['f2p'], ""),
 ("nucleon iso",    13, 8.303,          938.913, nuc_iso, "6-7 keV residual"),
 ("p alt-conv",     13, 8.066,          937.335, nuc_iso, "miss 1.6 MeV"),
]
print(f"{'state':14s} {'m_indep':>9s} {'m_tex':>9s} {'d_tex':>7s} {'PDG':>9s} {'dev%':>7s}")
for name,N,lam,mtex,pdg,note in rows:
    mi = m(N,lam)
    dev = (mi-pdg)/pdg*100 if pdg else float('nan')
    flag = "" if abs(mi-mtex) < 0.05 else "  <-- MISMATCH vs tex"
    print(f"{name:14s} {mi:9.3f} {mtex:9.3f} {mi-mtex:7.3f} {(pdg or 0):9.3f} {dev:7.3f}{flag}  {note}")

print("\n-- analytic eigen-identities --")
lo,hi = (7-sqrt(29))/2,(7+sqrt(29))/2
print(f"A2g pair (7+-sqrt29)/2 = {lo:.4f}, {hi:.4f}   (tex 0.807, 6.193)  sum {lo+hi}, prod {lo*hi:.1f}")
a,b = (11-sqrt(17))/2,(11+sqrt(17))/2
print(f"A1u pair (11-+sqrt17)/2 = {a:.4f}, {b:.4f}  (tex 3.4384, 7.5616)")
sa,sb = (9-sqrt(17))/2,(9+sqrt(17))/2
print(f"shell basis (9-+sqrt17)/2: trace {sa+sb:.0f} det {sa*sb:.0f}   (tex: trace 9, det 16)")
print(f"hex cap 7*m0 = {7*m0:.2f}  (tex 490.18);  m_K = 7*(m0+me) = {7*(m0+me):.3f}")
print(f"m_phi - m_K* = {m(14,8.441)-m(13,1.6626):.2f} MeV  (tex: 117)")
print(f"d* 34*m0 = {34*m0:.1f} MeV  (PDG d*(2380): 2380)")
print(f"pion split PDG: {PDG['pi_pm']-PDG['pi0']:.2f} MeV (tex 4.59)")
print(f"phi vs PDG: {(m(14,8.441)-PDG['phi'])/PDG['phi']*100:+.2f}%")

print("\n-- coupling-selection and width checks (added after physics audit) --")
kc_mu = 2/(3.141592653589793 - 2)
print(f"m_K at alpha_Cos=1: {7*(m0+me):.3f}; at kc/mu={kc_mu:.4f}: {7*(m0+kc_mu*me):.3f}  (PDG 493.677+-0.015)")
fpi = 3**0.25*m0*(1 + alpha/3.141592653589793)
def Gam_rho(mrho, mpi=139.57039):
    p = (mrho**2/4 - mpi**2)**0.5
    return p**3/(12*3.141592653589793*fpi**2)
print(f"Gamma_rho: spectral mass 775.29 -> {Gam_rho(775.29):.2f} MeV (tex 147.0, obs 147.4+-0.8)")
print(f"           superseded 775.7    -> {Gam_rho(775.7):.2f} MeV (was quoted 147.3)")
g = 775.29/(2**0.5*fpi); print(f"g = {g:.3f} (tex 5.935); g_orp = {3*g**2/(8*3.141592653589793**2*fpi)*1000:.2f} GeV^-1 (tex 14.49)")
print(f"omega-rho split: spectral 11*1.350*me = {11*1.350*me:.2f} MeV (PDG 7.40); superseded 12*me = {12*me:.2f}")
