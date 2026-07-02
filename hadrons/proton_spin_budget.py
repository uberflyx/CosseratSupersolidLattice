"""
The proton spin budget from the Cosserat rolling lock.

CLAIM. The fraction of the proton's spin carried by quark spins (the
singlet axial charge a0 = Delta Sigma of polarised DIS) equals the
Cosserat coupling number N^2 = 1/pi = 0.3183, the framework's fraction
of elastic energy carried by rotation rather than translation
(sec:rolling_constraint). The remaining 1 - 1/pi = 0.6817 is the locked
condensate-circulation share, which the Standard Model books as quark
orbital angular momentum plus gluon angular momentum.

DERIVATION SKELETON (semiclassical; full write-up in the monograph).
  1. Every quark is a locked pair: a partial dislocation (crystal sector)
     bound to a condensate phase vortex (circulation sector). The medium
     therefore holds a hadron's angular momentum in two reservoirs:
     the microrotation field phi (intrinsic, per-node spin: in the
     dictionary this is the quark spinor content) and the moment of the
     displacement/condensate flow (orbital: r x rho u_dot + r x rho_s v_s).
  2. Polarised DIS measures the axial current, which couples to the
     spinor, i.e. to the microrotation reservoir only. Decay amplitudes
     and current matrix elements are EDGE quantities in the framework's
     NLO taxonomy, exactly the class where the rotational fraction N^2
     survives (monograph, results overview).
  3. The rolling constraint is holonomic: node rotation and translation
     are phase-locked at the contact, so both reservoirs are driven at
     the defect's one Compton frequency. For harmonic reservoirs at a
     common frequency, the action (= angular momentum for a rotational
     mode) partitions as the energy partitions: J_i = E_i / omega.
  4. The energy partition is the Cosserat number: rotational share N^2,
     translational/flow share 1 - N^2. Hence
        Delta Sigma = N^2 = 1/pi,   L_q + J_g = (1/2)(1 - 1/pi).
  5. The electron shows no analogous crisis because it is an elementary
     screw: a probe scatters off the whole locked object, and there is
     no internal partition to expose.

BONUS PREDICTION. With the octet axial charge a8 from hyperon beta
decays, Delta s = (a0 - a8)/3 follows with no new input.

Experimental anchors (verified July 2026):
  HERMES (Q^2 = 5 GeV^2, MSbar):        a0 = 0.330 +/- 0.025
  HERMES+COMPASS+JLab combination:      a0 = 0.33  +/- 0.014
  COMPASS band:                         0.26 < a0 < 0.36
  Lattice chiQCD (MSbar, 2 GeV):        a0 = 0.405 +/- 0.045
  Ellis-Jaffe (Delta s = 0) baseline:   a0 = 0.59   (the 'crisis' gap)
  Octet charge a8: 0.585 +/- 0.025 (SU(3) fits); 0.675 (older analyses)
"""

import numpy as np

N2 = 1.0 / np.pi                     # Cosserat coupling number (derived)

print("-" * 72)
print("PROTON SPIN BUDGET FROM THE ROLLING LOCK")
print("-" * 72)
print(f"Prediction  Delta Sigma = N^2 = 1/pi = {N2:.4f}")
print(f"Ji budget   (1/2)DeltaSigma = {N2/2:.4f}  |  L_q + J_g = {(1-N2)/2:.4f}")
print(f"Circulation share of proton spin      = 1 - 1/pi = {1-N2:.4f}")
print("-" * 72)

anchors = [("HERMES (Q2=5, MSbar)",        0.330, 0.025),
           ("HERMES+COMPASS+JLab comb.",   0.330, 0.014),
           ("COMPASS band midpoint",       0.310, 0.050),
           ("Lattice chiQCD (MSbar 2GeV)", 0.405, 0.045),
           ("Ellis-Jaffe baseline",        0.590, 0.0)]
print(f"{'anchor':34s}{'value':>8s}{'pull':>10s}")
for name, v, s in anchors:
    pull = (v - N2)/s if s else float('inf')
    print(f"{name:34s}{v:8.3f}{pull:9.2f}s")
print("-" * 72)

# Strange-sea polarisation from the octet charge, Delta s = (a0 - a8)/3.
for a8, sa8, tag in [(0.585, 0.025, "SU(3) fit"), (0.675, 0.020, "older analyses")]:
    ds  = (N2 - a8)/3
    sds = sa8/3
    print(f"Delta s  (a8 = {a8:.3f}, {tag:14s}): {ds:+.4f} +/- {sds:.4f}")
print("Global-fit comparison: Delta s ~ -0.08 to -0.12 (inclusive DIS route);")
print("SIDIS at measured x prefers less negative, a known open tension.")
print("-" * 72)
print("Scale note: a0 runs only through the two-loop anomaly, so a")
print("scale-free geometric value is compatible; the framework number")
print("carries no Q^2. EIC will test the orbital share L_q + J_g = 0.341")
print("directly through generalised parton distributions (Ji sum rule).")
