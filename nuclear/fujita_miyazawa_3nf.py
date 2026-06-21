"""
The two-pion-exchange three-nucleon force, parameter-free, and where it binds.

Cosserat Supersolid Lattice framework. The deuteron showed that two-pion exchange between two
nucleons is just the docking (the intermediate-range central attraction) plus a tensor that
opposes the one-pion force, so it adds nothing new there. The one thing two-pion exchange does
that has NO two-body analogue is the three-nucleon force: a pion goes from nucleon 1 to nucleon
2, the kick lifts nucleon 2 into a Delta(1232), and the Delta fires a second pion on to nucleon
3. This is the Fujita-Miyazawa force (Prog. Theor. Phys. 17, 360, 1957). It cannot appear in the
deuteron at all, and it is exactly what the two-nucleon force needs to bind 3H, 3He and 4He,
which it otherwise underbinds by ~1 MeV.

This script shows the force is parameter-free in the framework. In chiral effective field theory
the two-pion-3NF strength is carried by the dimension-two low-energy constants c_3 and c_4, whose
dominant part is fixed by the Delta-isobar through resonance (Delta) saturation:

  c_3 = -2 c_4 = -4 h_A^2 / (9 Delta),     Delta = m_Delta - m_N        (EHM, RMP 81, 1773, Eq.2.47)

Every quantity on the right is derived in the framework:
  - h_A, the pi-N-Delta axial coupling, is h_A = 3 g_A/(2 sqrt2) from SU(4)/SU(6) (large-Nc),
    the same spin-flavour structure that fixes the nucleon moment ratio mu_n/mu_p = -2/3;
  - Delta = m_Delta - m_N uses the framework's spectral masses, m_Delta = 1234.31 MeV (lambda=9.05
    void-extended cluster) and m_N = 938.918 MeV.
So c_3 and c_4 carry no fitted input. They come out at the standard Delta-saturation value, which
is the dominant piece of the empirical c_3 ~ -3 to -5.6 GeV^-1 (the rest is sub-leading non-Delta
physics). With realistic three-body wavefunctions this force supplies ~0.9 MeV of attraction in
the triton, about 10% of its 8.48 MeV binding, which is precisely the amount the two-nucleon force
misses. The framework's geometric counterpart is the triple-shared atom of Sec. nuclear_extension:
a lattice node touched by all three nucleon shells at once, with no two-body analogue.

Author: Mitch Cox.  Run: python fujita_miyazawa_3nf.py
"""

import numpy as np

# ---- framework-derived inputs (CODATA / framework spectral masses), nothing tuned ----
M0      = 0.51099895/7.2973525693e-3      # m_e/alpha = 70.025 MeV
G_A     = 1.279                            # nucleon axial coupling (derived, Sec. FD_ratio etc.)
F_PI    = 3**0.25*M0                        # 92.16 MeV
M_PI    = 2*M0                              # 140.05 MeV
M_N     = 938.918                          # nucleon isoscalar mass (derived)
M_DELTA = 1234.31                          # Delta(1232) spectral mass (derived, lambda=9.05)
HBARC   = 197.3269804

DELTA_GeV = (M_DELTA - M_N)/1000.0          # Delta-nucleon splitting in GeV
H_A       = 3*G_A/(2*np.sqrt(2.0))          # pi-N-Delta coupling, SU(4)/SU(6) (large-Nc)

# ---- Delta-saturation of the two-pion-3NF low-energy constants (EHM Eq. 2.47) ----
c3 = -4*H_A**2/(9*DELTA_GeV)
c4 = -c3/2.0

print("=== Two-pion-exchange three-nucleon force: parameter-free coupling ===")
print(f"  g_A   = {G_A}          (derived)")
print(f"  f_pi  = {F_PI:.2f} MeV    (3^1/4 m0)")
print(f"  m_pi  = {M_PI:.2f} MeV   (2 m0)")
print(f"  m_Delta - m_N = {M_DELTA:.1f} - {M_N:.1f} = {DELTA_GeV*1000:.1f} MeV   (both spectral)")
print(f"  h_A = 3 g_A/(2 sqrt2) = {H_A:.4f}   [SU(4)/SU(6), same rule as mu_n/mu_p = -2/3]")
print()
print("  Delta-saturation  c_3 = -2 c_4 = -4 h_A^2/(9 Delta):")
print(f"    c_3 = {c3:+.2f} GeV^-1    (empirical c_3 ~ -3.2 to -5.6; Delta saturates ~2/3 of it)")
print(f"    c_4 = {c4:+.2f} GeV^-1    (empirical c_4 ~ +2.0 to +3.4)")
print()
print("  The framework reproduces the standard Delta-saturation value with nothing fitted.")
print("  This force has no two-body shadow, so it did nothing in the deuteron; in the triton")
print("  it supplies the ~0.9 MeV (~10% of 8.48 MeV) that the two-nucleon force underbinds.")
print("  Geometric counterpart in the framework: the triple-shared atom (Sec. nuclear_extension).")

# ---- scale of the force: the pion range that sets where it acts ------------------------
r_pi = HBARC/M_PI
print(f"\n  Each of the two pion legs has range hbar c/m_pi = {r_pi:.2f} fm, so the force is")
print(f"  intermediate-range and strongest where three nucleons crowd inside ~1.5 fm, exactly")
print(f"  the close-packed geometry of the tritium and helium-4 compounds.")
