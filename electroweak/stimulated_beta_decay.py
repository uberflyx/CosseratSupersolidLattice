"""
Can an intense electromagnetic field stimulate beta decay?

In the framework the weak interaction is not a separate force but the axial
channel T_1g of the crystal, coupled to the photon's channel T_1u through the
hemitropic chirality theta_ch = alpha^2/(2pi). Beta decay is a weak-channel
rearrangement, a neutron-like defect converting to a proton-like one while
emitting the lepton pair. Because the weak vertex carries a photon component of
size theta_ch, an external electromagnetic field can drive that vertex. This is
a genuine prediction and it cuts against textbook lore: beta decay is not
immutable in this picture, it has an electromagnetic handle.

The size of the handle. An external field of photon occupation n per mode
(the degeneracy parameter of the source) drives the weak vertex with an extra
amplitude proportional to theta_ch * sqrt(n) relative to the spontaneous one.
The driven contribution to the RATE is therefore

    eta ~ theta_ch^2 * n * zeta,

where zeta <= 1 is a phase-space matching factor: beta decay shares its energy
between electron and neutrino as a continuum, so a monochromatic drive engages
only part of the final-state phase space. Since theta_ch^2 ~ 7e-11, an order-one
change needs n ~ 1/theta_ch^2 ~ 1e10 photons per mode at the transition energy.

That is the whole feasibility question. At MeV transition energies no source
reaches such occupation, so high-Q emitters are untouchable. But the occupation
of a modern X-ray free-electron laser at keV energies is n ~ 1e9 to 1e11, and
the lowest beta Q-values are keV-scale. The natural target is therefore a
low-Q emitter driven on resonance by an XFEL, where eta can climb to the
percent level or beyond and become measurable in a lifetime comparison.

The form factor does not hurt here. The beta momentum transfer is the Q-value,
q l ~ Q/(m0 c^2) * (something) << 1, so f(q) ~ 1: unlike the T_2u shadow, the
low-energy weak vertex is not shell-suppressed. The suppression is entirely the
theta_ch^2 coupling, which the source occupation can compensate.
"""

import numpy as np

alpha = 1.0 / 137.035999177
theta_ch = alpha ** 2 / (2.0 * np.pi)
theta_ch2 = theta_ch ** 2                 # ~7.2e-11

def eta(n, zeta=1.0):
    """Fractional rate enhancement eta = theta_ch^2 * n * zeta."""
    return theta_ch2 * n * zeta

# Emitter Q-values [keV] and matched photon sources with typical mode occupation.
# n = degeneracy parameter (photons per mode) of the source at that energy.
CASES = [
    # name,             Q_keV,   source,                     n
    ("Re-187",           2.47,   "X-ray FEL @ 2.5 keV",      1e9),
    ("Re-187 (best n)",  2.47,   "X-ray FEL, high degeneracy",1e11),
    ("Ni-63",           66.9,    "X-ray FEL @ 67 keV",       1e8),
    ("tritium H-3",     18.6,    "X-ray FEL @ 18.6 keV",     1e9),
    ("free neutron",   782.0,    "MeV gamma (Compton)",      1e-2),
    ("C-14",           156.0,    "hard-X / soft-gamma",      1e3),
]

def line():
    print("-" * 76)

print(__doc__.strip().splitlines()[0])
line()
print(f"theta_ch = {theta_ch:.2e}   theta_ch^2 = {theta_ch2:.2e}")
print(f"order-one enhancement needs n ~ 1/theta_ch^2 = {1/theta_ch2:.1e} photons/mode")
line()
print("Enhancement eta = theta_ch^2 * n (matching factor zeta = 1 shown; realistic")
print("zeta ~ 0.1-1 lowers each by up to 10x):")
print(f"  {'emitter':<16s} {'Q [keV]':>8s} {'source':<28s} {'n':>7s} {'eta':>10s}")
for name, Q, src, n in CASES:
    e = eta(n)
    tag = "  <- measurable" if e > 1e-3 else ""
    print(f"  {name:<16s} {Q:>8.2f} {src:<28s} {n:>7.0e} {e:>10.1e}{tag}")
line()
print("Reading it:")
print("  High-Q emitters (free neutron, C-14) are hopeless: no source reaches the")
print("  occupation needed at MeV energies. The lever is a LOW-Q emitter on an")
print("  X-ray FEL. For Re-187 (Q = 2.47 keV) at XFEL occupation n ~ 1e9-1e11,")
print(f"  eta ~ {eta(1e9):.0e} to {eta(1e11):.0e}: a percent-to-order-unity change")
print("  in the beta lifetime under illumination. Even with a matching factor of")
print("  a tenth, the conservative case sits at the per-mille level, which a")
print("  lifetime comparison can resolve.")
print()
print("  This is a clean, near-term test of the EM-weak identification. A")
print("  measured change in a beta lifetime under intense resonant X-rays, scaling")
print("  as theta_ch^2 * n and absent for a chargeless (edge-dislocation) decay")
print("  channel, would be the hemitropic coupling read directly. A null result at")
print("  n where eta should exceed the measurement floor would bound theta_ch")
print("  from the laboratory rather than from neutrino masses or birefringence.")
line()
print("Honest caveats (each a factor, none a showstopper for the scaling):")
print("  - matching factor zeta from the continuous beta spectrum (up to ~10x down);")
print("  - the standard laser-assisted effect acts through the emitted charged")
print("    particle's motion in the field, and careful treatments put it near")
print("    1e-5 or below at reachable intensities (larger claims are contested;")
print("    Akhmedov, Phys. Atom. Nucl. 74, 1299 (2011)). The framework's mechanism")
print("    is DISTINCT -- it drives the nuclear weak vertex, not the final-state")
print("    electron -- and predicts a larger change, so the two are separable and")
print("    the test is which one appears;")
print("  - XFEL per-mode occupations ~1e9 are established (Saldin, Schneidmiller,")
print("    Yurkov, Opt. Commun. 281, 1179 (2008)); the exact vertex matrix element")
print("    for a given transition is the remaining input.")
print("  The prediction to take away is the scaling eta ~ theta_ch^2 n and the")
print("  target: a low-Q emitter under an X-ray FEL, where it becomes measurable.")
