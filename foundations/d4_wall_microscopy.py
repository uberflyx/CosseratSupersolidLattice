#!/usr/bin/env python3
"""D4 wall microscopy: what reproduces exactly, and why r_c cannot be pinned.

v2 step 1 of the pi/4 programme (open_items 2026-08-29). The task was to pin
the truncated-shifted Morse cutoff r_c on the monograph's printed rigid fault
microscopy before computing any slope sum. It does not pin, and the reason is
worth more than the number would have been.

WHAT REPRODUCES EXACTLY (the geometry and the potential are certified)
The D4 {111} wall has 24 bonds per node, 12 crossing. From a layer-0 node the
crossing set is 9 bonds, and the partial slip b = d t_hat (d = 1/sqrt3,
t_hat = (1,1,-2,0)/sqrt6) maps their lengths to
    2 -> sqrt(2/3) = 0.81650   (compressed to the interplanar spacing)
    4 -> sqrt(5/3) = 1.29099   (stretched)
    2 -> 1                      (unchanged)
    1 -> sqrt2 = 1.41421        (stretched out of the contact shell)
exactly the 2 + 4 + remainder the monograph prints. With the gravitational
sector's pinned Morse (a = 7/3, D = 9/98, V''(1) = k_n = 1) the two strained
classes store +0.02620 and +0.02230 k_n per bond, matching the printed digits,
and their column sum is 0.14170 against the printed 0.1416. Geometry, slip
vector, potential and bookkeeping are therefore all confirmed independently.

WHAT DOES NOT PIN, AND WHY THAT IS THE RESULT
The printed remainder, -0.0051 k_n, is reproduced only at r_c = 1.29706, and
that value is a knife edge: gamma_rigid runs from 0.1408 to 0.1203 as r_c goes
from 1.292 to 1.316, i.e. a 15 per cent swing across a 2 per cent change in
cutoff, and r_c = 1.29706 lands on no lattice length (nearest structural
candidate 13/10, 0.23 per cent away; the shell it must sit between is
sqrt(5/3) = 1.291 to sqrt2 = 1.414). Restricting to nearest-neighbour contacts
alone makes the remainder provably non-negative (it equals D + V(r_c) >= 0),
so the printed negative remainder requires further shells and a sharp cutoff
placed to a part in a thousand.

The slope sum inherits the fragility and worse. Across the admissible window
the central-channel S swings 0.26 -> 2.97 with no plateau, and the harmonic
ratios change sign. A sharp cutoff sitting between two shells creates a
spurious contact event exactly where the misfit harmonics are measured.

WHERE THAT POINTS
The monograph's own lossless argument already forbids the sharp cutoff: a
contact that broke while still loaded would radiate a snap, so the stiffness
must fade continuously to zero as a contact opens, and a hard cut-off produces
"spurious contact events". This computation is the central-channel evidence
for the same conclusion, and it shows the fade is needed for cohesion as well
as for the tangential springs. Substituting a smooth fade of width 0.10
stabilises the harmonics (V_2/V_1 falls from order 0.2 to order 0.01 and S
lands near 0.9-1.1 for both shell placements), which is the behaviour a
physical misfit line should have.

CONSEQUENCE FOR THE PROGRAMME: r_c is not an input the printed microscopy can
supply, and should not be fitted to it. The fade width and the fade law are
the same single free function the monograph identifies as the handover law,
now shown to control the central channel too. The pi/4 target
S_alpha = 1.001067 remains untouched and unconsulted in any choice above.
"""
import numpy as np

A, D = 7/3, 9/98
S2 = 1/np.sqrt(2)
W = np.array([1, 1, 1, 0])/np.sqrt(3)
T = np.array([1, 1, -2, 0])/np.sqrt(6)
LAY, DHOP = 1/np.sqrt(6), 1/np.sqrt(3)

def V(r):
    return D*(np.exp(-2*A*(r-1)) - 2*np.exp(-A*(r-1)))

def star():
    out = []
    for i in range(4):
        for j in range(i+1, 4):
            for si in (-1, 1):
                for sj in (-1, 1):
                    v = np.zeros(4); v[i] = si*S2; v[j] = sj*S2
                    out.append((int(round((v@W)/LAY)), v))
    return out

def crossing_classes():
    from collections import Counter
    return Counter(round(np.linalg.norm(v + DHOP*T), 5) for m, v in star() if m >= 1)

def main():
    print("crossing-bond image under the partial slip:", dict(crossing_classes()))
    for r, n in [(np.sqrt(2/3), 2), (np.sqrt(5/3), 4)]:
        print(f"  {n} bonds -> r = {r:.5f}, dV = {V(r)-V(1):+.5f} k_n each")
    core = 2*(V(np.sqrt(2/3))-V(1)) + 4*(V(np.sqrt(5/3))-V(1))
    print(f"  column core sum = {core:.5f} k_n   (monograph prints 0.1416)")
    print(f"  printed total 0.1365 implies remainder -0.0051, reachable only at")
    print(f"  r_c = 1.29706, a knife edge: see the docstring for the sensitivity.")

if __name__ == "__main__":
    main()
