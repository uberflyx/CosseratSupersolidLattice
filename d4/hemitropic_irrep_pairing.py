"""
Which force channels the hemitropic (chiral) coupling connects.

The forces section assigns each force of nature to one O_h irrep and finds five
channels that carry no single-node force. The weak interaction is recovered
nonetheless, because its channel T_1g mixes at second order with the photon's
T_1u through the hemitropic cross-term W_ch, whose strength is the chirality
theta_ch = alpha^2/(2pi). This script asks the group-theory question that
generalises that step: the hemitropic coupling is a pseudoscalar, so it carries
the parity-odd singlet A_1u of O_h, and it connects two channels Gamma_a and
Gamma_b exactly when their product contains A_1u. Tensoring by A_1u flips a
representation's parity and preserves its orbital letter, so the coupling pairs
each channel with its opposite-parity partner. The pairs are computed below from
the character table.

The consequence: the "silent" channels are not truly decoupled. Each mixes, at
order theta_ch, with a channel that does carry a force -- E_u with the tensor
channel E_g, T_2u with the strong channel T_2g, A_1u (the dilaton) with gravity
A_1g, and A_2g (baryon winding) with the eta-prime A_2u. Only the weak pairing
T_1g <-> T_1u produces a new long-range force, because only there is the partner
(electromagnetism) massless. The other partners are massive and short-ranged
(strong, tensor) or suppressed to alpha^19 (gravity), so the silent channels are
the parity-odd shadows of those forces, not new forces.
"""

import numpy as np

# O_h classes: E, 8C3, 6C2, 6C4, 3C2', i, 6S4, 8S6, 3sigma_h, 6sigma_d
SIZES = np.array([1, 8, 6, 6, 3, 1, 6, 8, 3, 6])
ORDER = 48

CHAR = {
    'A1g': [1, 1, 1, 1, 1,  1, 1, 1, 1, 1],
    'A2g': [1, 1, -1, -1, 1,  1, -1, 1, 1, -1],
    'Eg':  [2, -1, 0, 0, 2,  2, 0, -1, 2, 0],
    'T1g': [3, 0, -1, 1, -1,  3, 1, 0, -1, -1],
    'T2g': [3, 0, 1, -1, -1,  3, -1, 0, -1, 1],
    'A1u': [1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
    'A2u': [1, 1, -1, -1, 1, -1, 1, -1, -1, 1],
    'Eu':  [2, -1, 0, 0, 2, -2, 0, 1, -2, 0],
    'T1u': [3, 0, -1, 1, -1, -3, -1, 0, 1, 1],
    'T2u': [3, 0, 1, -1, -1, -3, 1, 0, 1, -1],
}
CHAR = {k: np.array(v) for k, v in CHAR.items()}
IRREPS = list(CHAR)

FORCE = {'A1g': 'gravity', 'T2g': 'strong', 'Eg': 'tensor (f2)', 'T1u': 'electromagnetism',
         'A2u': "eta' / strong-CP", 'T1g': 'weak (2nd order)', 'A2g': 'baryon-number winding',
         'A1u': 'dilaton 0- (silent)', 'Eu': 'silent', 'T2u': 'silent'}

def inner(a, b):
    return int(round(np.sum(SIZES * CHAR[a] * CHAR[b]) / ORDER))

def mult_A1u(a, b):
    """Multiplicity of A1u in a (x) b = <chi_a chi_b, chi_A1u>."""
    return int(round(np.sum(SIZES * CHAR[a] * CHAR[b] * CHAR['A1u']) / ORDER))

def line():
    print("-" * 72)

# Sanity: orthonormal characters.
assert all(inner(x, x) == 1 for x in IRREPS)
assert inner('Eg', 'T2u') == 0 and inner('A1g', 'A2g') == 0

print(__doc__.strip().splitlines()[0])
line()
print("Hemitropic pairing (channels linked when A1u is in their product):")
print(f"  {'gerade channel':>28s}        {'partner (parity-flip)':<24s}")
gerade = ['A1g', 'A2g', 'Eg', 'T1g', 'T2g']
LONG = {'T1u': 'MASSLESS -> genuine force'}
for a in gerade:
    partners = [b for b in IRREPS if mult_A1u(a, b) == 1]
    assert len(partners) == 1, f"{a}: expected one partner, got {partners}"
    b = partners[0]
    note = LONG.get(b, 'massive/short-range or alpha^19 -> P-odd shadow')
    print(f"  {a:>4s} [{FORCE[a]:<22s}] <-> {b:<4s} [{FORCE[b]:<20s}]  {note}")
line()
print("Reading it:")
print("  Every 'silent' channel pairs with a carrying one at order theta_ch:")
print("    E_u  with the tensor channel E_g,")
print("    T_2u with the strong channel T_2g,")
print("    A_1u (dilaton) with gravity A_1g,")
print("    A_2g (baryon winding) with the eta-prime A_2u.")
print("  Only T_1g <-> T_1u makes a new long-range force, because only")
print("  electromagnetism is massless. So there is still no sixth force, but the")
print("  silent channels are the parity-odd shadows of the strong, tensor, and")
print("  gravitational channels, carried by the same chirality that runs the weak")
print("  interaction -- not channels that couple to nothing.")
print()
print("  Open: whether the strongest shadow, T_2u on the strong channel at")
print("  theta_ch ~ 8.5e-6, leaves a measurable trace of parity violation in the")
print("  strong sector is a matrix-element question, and a second place the")
print("  chirality parameter could be read.")
