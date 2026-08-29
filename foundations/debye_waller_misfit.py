#!/usr/bin/env python3
"""The ghost-lattice reading: the misfit harmonics are Debye-Waller suppressed,
and the fade rate stops being a free function.

THE STRUCTURAL POINT (Mitch's question, 2026-08-29)
A supersolid is a condensate with a periodic density modulation, not a set of
objects joined by bonds. Nothing is attached to anything, so "at what distance
does a bond cease to exist" is not a well-posed question, and the sharp
truncation radius r_c that d4_wall_microscopy.py showed to be a knife edge was
a parameter of a description that does not apply. The framework already says
so: at the Mott point a node is caged collectively by twelve neighbours' zero-
point pressure, the well is flat-bottomed, and the effective potential is the
BARE pair potential averaged over the zero-point fluctuations
(sec:gamma4_scp), with the self-consistent amplitude sigma = 0.354 ell.

WHAT THAT DOES TO THE MISFIT POTENTIAL
Averaging a periodic potential over Gaussian node smearing is a convolution,
and convolution in real space is multiplication in Fourier space. For a
corrugation of period d with harmonics V_m, the dressed harmonics are

    V_m(dressed) = V_m(bare) * exp(-m^2 k_1^2 sigma_1^2 / 2),   k_1 = 2 pi/d,

which is exactly a Debye-Waller factor: the same suppression that weakens
high-order Bragg peaks in a warm crystal, here produced by zero-point rather
than thermal motion. The slope sum that misfit_slope_sum.py showed to be the
ONLY property of the gamma-surface the core width can feel is therefore

    S = 1 + sum_{m>=2} m (V_m/V_1)_bare exp(-k_1^2 sigma_1^2 (m^2 - 1)/2).

Three consequences, and they replace the open fade law:
 1. S -> 1 exponentially in sigma. Smearing washes the overtones out and
    leaves a pure sinusoid, which is the Frenkel form the alpha paper assumed
    and whose fidelity it flagged as the source of the 0.1 per cent residual.
    The residual is therefore not a modelling defect: it is the leftover
    second harmonic.
 2. The fade of a contact is no longer a free function. There is no contact to
    fade. The single scale is sigma, which the Mott self-consistency fixes.
 3. Because the suppression is exponential in sigma^2, S - 1 is a very sharp
    probe of the smearing amplitude.

THE NUMBER
With the per-axis projection sigma_1 = sigma/sqrt3 = 0.20438 ell and the
partial-hop period d = ell/sqrt3, the second-harmonic suppression is
5.986e-4, and the on-file alpha demand S_alpha = 1.001067 (committed in
misfit_slope_sum.py BEFORE this reading existed) requires a bare second
harmonic of

    (V_2/V_1)_bare = 0.891,

i.e. a bare corrugation whose second harmonic is comparable to its first,
which is what a hard-contact geometric landscape looks like. NOT ADOPTED as a
result: it is a consistency statement, since the bare surface has not been
computed independently. What makes it worth recording is that no scale was
tuned to get there, sigma comes from the Koide/Mott sector and d from the
lattice, and the answer landed within 12 per cent of unity rather than orders
away.

WHAT WOULD FALSIFY IT
S is exponentially sensitive to sigma_1: a 1e-3 target on S pins sigma_1 to
6.3 per cent, and alpha's own band pins it to 4e-10 ell. So an independent
computation of the bare (V_2/V_1) turns this into a sharp test of the Mott
amplitude, and vice versa. If the bare second harmonic comes out at 0.1 or at
5, this reading fails.
"""
import numpy as np

D_HOP = 1/np.sqrt(3)
K1 = 2*np.pi/D_HOP
SIGMA_MOTT = 0.354
SIGMA_AXIS = SIGMA_MOTT/np.sqrt(3)
S_ALPHA = 1.001067          # from foundations/misfit_slope_sum.py, on file first

def dw(m, sig=SIGMA_AXIS):
    """Debye-Waller suppression of harmonic m relative to m = 1."""
    return np.exp(-0.5*K1**2*sig**2*(m*m - 1))

def slope_sum(ratios, sig=SIGMA_AXIS):
    """S from bare harmonic ratios [(m, V_m/V_1), ...] with m >= 2."""
    return 1.0 + sum(m*r*dw(m, sig) for m, r in ratios)

def main():
    print(f"d = {D_HOP:.5f} ell,  k1 = {K1:.4f} /ell,  sigma_1 = {SIGMA_AXIS:.5f} ell")
    print(f"Debye-Waller: DW(2)/DW(1) = {dw(2):.4e},  DW(3)/DW(1) = {dw(3):.2e}")
    v2 = (S_ALPHA - 1)/(2*dw(2))
    print(f"\nS_alpha = {S_ALPHA} demands bare (V_2/V_1) = {v2:.4f}")
    print("  hard-contact landscapes have V_2/V_1 of order unity; this sits")
    print("  within 12 per cent of it, with no scale tuned to arrange that.")
    print(f"\nsensitivity: sigma_1 sweep at fixed bare V_2/V_1 = {v2:.3f}")
    for s in [0.18, 0.19, SIGMA_AXIS, 0.21, 0.22]:
        print(f"   sigma_1 = {s:.5f}  ->  S = {slope_sum([(2, v2)], s):.6f}")
    print("\nS is exponentially sensitive to the smearing, so an independent")
    print("bare-surface computation turns this into a sharp test either way.")

if __name__ == "__main__":
    main()
