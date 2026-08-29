#!/usr/bin/env python3
"""The inter-family hopping t: from a data bound to a mechanics-anchored range.

THE OBJECT
Four {111} plane families each carry the Z_3 stacking ring, and the gauge-
invariant inter-family coupling, after Frank's rule, is delta H = 3 t I: a
pure identity shift (paper, projection and protection theorem). So t moves no
angle at any value; it shifts all three masses together, the solar splitting
through the difference, and the data bound it at |t| < 0.040 meV. The open
item asks what the cross-slip mechanics itself gives.

WHY THE NAIVE BARRIER ARGUMENT SAYS NOTHING
The hop is a cross-slip event: the dissociated edge must constrict its
stacking-fault ribbon to change glide plane. In a metal that constriction is
thermally activated and rare. Here the defect tunnels, and the WKB exponent
for a barrier of height Delta_f crossed over one interplanar spacing by a
defect of inertia m_1 = 5 meV is

    S/hbar = sqrt(2 m_1 Delta_f) d_111 / (hbar c)  ~  6e-6,

utterly unsuppressed: a 5 meV object does not feel a barrier as a wall.
The paper's coherence argument used exactly this. The suppression of t is
therefore amplitude-based, not barrier-based: the hop must proceed through a
virtual constricted intermediate, and its smallness is the smallness of a
second-order process, not of an exponent.

THE FAULT GAP, CALIBRATED INTERNALLY
The constricted intermediate is a fault-channel excitation. Its gap Delta_f
is fixed by the sector's own constants: the flavour coupling C = g^2/Delta_f
with C = 30.30 meV, and the vertex g^2 = K_w m_D calibrated by the
rotational-channel consistency check (edge_moment_overlap.py):
K_w = 6.628 meV against the gap m_D = 130.9 MeV reproduces the ring's K_w
with overlap 0.62. The same g^2 then gives

    Delta_f = g^2 / C = K_w m_D / C = 28.6 MeV = m_D / 4.57,

the framework's quantitative form of "the fault is anomalously cheap": the
fault intermediate lies 4.6 times below the rotational one, and the pure
number 4.57 = C/K_w = (sqrt(3)+tan(delta))/(2 tan(delta)) ... = the mirror
ratio already derived. This number is an input the mirror-constants item
(open_items) has been waiting for, and it is flagged there.

THE RANGE FOR t, AND WHAT IT PROMOTES
The gate vertex g_x connecting the twelve-mode sectors through the
constricted state is not yet derived, so |t| is bracketed by the two
structures the mechanics allows:
  upper (elastic gate):    |t| = h * (g / Delta_f)     one vertex, one gap
  lower (second order):    |t| = h^2 / Delta_f         the hop itself twice
The two ends tell different stories, and the honest close of the item is
the map between them rather than a single verdict. At the lower end t is
invisible at any foreseeable precision. At the upper end 6 t (m2 - m1) is
0.6 per cent of the solar splitting, 0.45 of the experimental error, which
is exactly the order the paper's next-order analysis asks a correction to
supply. The sign is physics, not convention: a bonding gate (symmetric
combination lowered) gives t < 0 and pushes the solar splitting further
below the data; an antibonding gate gives t > 0 and closes most of the
-0.66 sigma residual. The gate vertex is thereby promoted from nuisance to
candidate next-order physics, and deriving its sign and magnitude becomes
the sharpest cheap test the twelve-mode sector offers.
"""

import numpy as np

# ---- sector constants [meV unless noted]
H      = 14.4          # intra-family hopping
KW     = 6.628         # screw-channel on-site energy
C_EDGE = 30.30         # edge-channel rank-one coupling
M_D    = 130.9e6       # rotational gap [meV]
M1     = 5.03          # lightest mass
M2     = 10.00
M3     = 50.5
DM21_MEV2 = 74.71      # solar splitting [meV^2]
T_BOUND   = 0.040      # data bound on |t| [meV]
D111_ELL  = np.sqrt(2.0/3.0)   # interplanar spacing [ell]
DM31_MEV2 = 2523.0             # atmospheric splitting [meV^2]

def main():
    m0_mev = 0.51099895 * 137.035999177          # node mass [MeV]
    # WKB non-suppression
    g2 = KW * M_D                                 # meV^2
    delta_f = g2 / C_EDGE                         # meV
    s_wkb = np.sqrt(2 * (M1/ (m0_mev*1e9)) * (delta_f/(m0_mev*1e9))) * D111_ELL
    print("=" * 66)
    print("Inter-family hopping t from the cross-slip mechanics")
    print("=" * 66)
    print(f"  fault gap Delta_f = K_w m_D / C     : {delta_f/1e6:8.1f} MeV")
    print(f"                    = m_D / (C/K_w)   : m_D / {C_EDGE/KW:.3f}")
    print(f"  WKB exponent for m_1 over Delta_f   : {s_wkb:.1e}   (no barrier suppression)")
    print()
    t_upper = H * np.sqrt(g2) / delta_f
    t_lower = H**2 / delta_f
    print(f"  t upper (elastic gate, h g/Delta_f) : {t_upper:.2e} meV")
    print(f"  t lower (second order, h^2/Delta_f) : {t_lower:.2e} meV")
    print(f"  data bound                          : {T_BOUND:.2e} meV")
    print(f"  margin to bound                     : {T_BOUND/t_upper:8.0f}x  to {T_BOUND/t_lower:.0f}x")
    print()
    print("residual shifts (identity: m_k -> m_k + 3t), both signs shown:")
    for name, tmag in [("upper", t_upper), ("lower", t_lower)]:
        for sgn, lab in [(+1, "antibonding t>0"), (-1, "bonding t<0")]:
            t = sgn * tmag
            d21 = 6 * t * (M2 - M1) / DM21_MEV2
            d31 = 6 * t * (M3 - M1) / DM31_MEV2
            print(f"  {name:5s} {lab:16s}: dDm21/Dm21 = {d21:+.1e} "
                  f"({d21/0.013:+.2f} exp sigma), dDm31/Dm31 = {d31:+.1e} "
                  f"({d31/0.0084:+.2f}), d(sum) = {9*t:+.2e} meV")
    print()
    print("current residuals for orientation: Dm2_21 at -0.66 sigma (low),")
    print("Dm2_31 at +0.59 sigma (high), ratio R at -0.87 sigma.")
    print()
    print("CONCLUSION. The mechanics brackets |t| between 7e-6 and 1.5e-2 meV,")
    print("inside the 0.040 meV data bound everywhere, so the protection")
    print("theorem's consistency check closes: no angle moves, and the sum")
    print("moves by at most 0.13 meV against a 1.4 meV budget. But the upper")
    print("end is not invisible: an antibonding gate at full strength supplies")
    print("+0.45 sigma on the solar splitting, most of the current residual,")
    print("in the direction and at the order the next-order analysis requests.")
    print("The gate vertex, its sign first, is promoted to the sector's")
    print("cheapest live derivation target.")

if __name__ == "__main__":
    main()
