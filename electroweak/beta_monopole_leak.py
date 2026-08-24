#!/usr/bin/env python3
"""Predicted second-sound energy leak per beta decay, and its confrontation
with the tritium Q-value comparison.

E_2nd = (1/2) C_11 eps^2 V  with C_11 = 3 mubar, V = N ell^3, mubar ell^3 = m_0 c^2
      = (3/2) eps^2 N m_0 c^2

eps = theta_ch = alpha^2/(2 pi)   the chiral monopole projection of a beta flip
N   = 13                          nodes in a nucleon coordination cluster

The leak is INDEPENDENT of Q, so it is a fixed endpoint offset, not a fractional
one, and it is absorbed by a free endpoint parameter -- hence invisible to the
neutrino-mass fit and visible only in the Q comparison.
"""
import numpy as np

alpha  = 1/137.035999
th_ch  = alpha**2/(2*np.pi)
m0c2   = 70.04e6          # eV, node rest energy
N      = 13               # nodes per nucleon cluster

def leak(eps=th_ch, N=N):
    return 1.5*eps**2*N*m0c2

print(f"theta_ch = {th_ch:.4e}")
print(f"normalisation check, eps=1 N=1 : {1.5*1.0**2*1*m0c2/1e6:.1f} MeV  (corpus: 105 MeV)")
print(f"PREDICTED LEAK  = {leak()*1e3:.1f} meV per beta decay, independent of Q\n")

Q_pt, dQ_pt = 18592.071, 0.022     # Medina Restrepo & Myers, PRL 131, 243002 (2023)
Q_k,  dQ_k  = 18591.49,  0.50      # KATRIN
print(f"Penning trap Q = {Q_pt:.3f} +- {dQ_pt:.3f} eV")
print(f"KATRIN       Q = {Q_k:.2f} +- {dQ_k:.2f} eV")
print(f"difference     = {Q_pt-Q_k:.3f} eV = {(Q_pt-Q_k)/dQ_k:.2f} sigma, KATRIN LOW (correct sign)")
print(f"prediction is {dQ_k/leak():.1f}x below current sensitivity\n")

print("sensitivity of the prediction to the reaction volume:")
for n, lab in [(1,"one lattice cell"), (13,"one nucleon cluster"), (39,"whole triton")]:
    print(f"   N={n:<3d} {lab:<22s} {leak(N=n)*1e3:7.1f} meV")
