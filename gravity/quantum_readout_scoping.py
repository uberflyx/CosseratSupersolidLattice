"""Quantum channels for detecting the covariance wave: honest numbers.

The wave IS a travelling squeeze of zero-point statistics, so three quantum
detection strategies suggest themselves. This script prices each.

[1] DIRECT VARIANCE CHANNEL. Measure the vacuum variance of a probe mode and
    look for its modulation at the wave frequency. The wave delivers a squeeze
    parameter r ~ h/2 ~ 5e-22. Resolving a fractional variance change eps with
    N independent samples needs N ~ 2/eps^2; the shortfall against any
    realistic sampling budget is computed below. Conclusion: hopeless directly;
    coherent interferometry (LIGO) is the optimal amplifier of the SAME
    squeezing operation, because a large coherent field converts a quadrature
    squeeze into a first-moment phase shift read at shot noise.

[2] QUANTA PER EVENT. The occupation of the radiation field at detection is
    astronomically classical: the number of quanta crossing a LIGO-scale
    aperture during a merger is ~1e36. Statistics, not clicks, is the only
    quantum observable, which is the Parikh-Wilczek-Zahariade noise-kernel
    programme.

[3] THE FRAMEWORK'S OWN PREDICTION. The emission vertex is pair creation
    (two-mode squeezing of the transverse branch), so the radiation should
    arrive in SQUEEZED-PAIR statistics, not coherent-graviton statistics.
    PWZ showed the noise a gravitational field imprints on geodesic deviation
    depends on its quantum state (vacuum / coherent / squeezed); the framework
    picks the squeezed column a priori. Single-quantum proposals (Tobar et al.
    2024) would face a two-quantum continuum edge rather than a graviton line;
    stimulated absorption is classically equivalent either way (Carney,
    Domcke, Rodd 2024), so the discriminator lives in spontaneous-exchange
    statistics.
"""

import numpy as np

G = 6.674e-11; c = 2.998e8; hbar = 1.055e-34

print("=" * 72)
print("[1] direct variance channel")
print("=" * 72)
h = 1e-21
r = h / 2
eps = 2 * r                        # fractional variance modulation ~ 2r
N_needed = 2 / eps**2
for name, rate, T in (("optical sampling, 1 yr", 1e14, 3.15e7),
                      ("optical sampling, 100 yr", 1e14, 3.15e9)):
    N_avail = rate * T
    print(f"    {name}: N_avail = {N_avail:.1e},  N_needed = {N_needed:.1e},"
          f"  shortfall = {np.sqrt(N_needed/N_avail):.1e} in amplitude")
print("    -> ten orders short even with a century of optical-rate sampling.")
print("       A coherent carrier of n photons boosts the same squeeze to a")
print("       first-moment signal ~ sqrt(n) larger: that amplifier is LIGO.")

print()
print("=" * 72)
print("[2] quanta per event")
print("=" * 72)
f = 100.0
hdot = 2 * np.pi * f * h
F = c**3 / (16 * np.pi * G) * hdot**2          # GW energy flux [W/m^2]
A, T = 16.0, 0.2                                # aperture ~ (4 m)^2, 0.2 s burst
Nq = F * A * T / (hbar * 2 * np.pi * f)
print(f"    flux = {F:.2e} W/m^2;  quanta through {A:.0f} m^2 in {T} s: {Nq:.1e}")
print("    -> occupation ~ 1e29 through this aperture alone: deep classical")
print("       regime; only field STATISTICS")
print("       (noise kernels, correlations) carry quantum information.")

print()
print("=" * 72)
print("[3] the discriminating statistics")
print("=" * 72)
print("    emission vertex = pair creation -> two-mode squeezed radiation.")
print("    PWZ noise kernel: squeezed states enhance the induced noise over")
print("    coherent states (exponentially in the squeeze parameter), so the")
print("    framework predicts the ENHANCED column, falsifiable in principle.")
print("    Single-quantum sensing (Tobar 2024): continuum edge, not a line;")
print("    stimulated response classical either way (Carney-Domcke-Rodd).")
print("    Static-sector entanglement (BMV): the constraint channel is carried")
print("    by a quantised medium -> entanglement expected; a decoherence null")
print("    would falsify. A positive alone proves less than advertised: this")
print("    framework entangles through the constraint sector with no graviton")
print("    pole at all, a live counterexample to 'entanglement implies graviton'.")
