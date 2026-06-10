#!/usr/bin/env python3
"""
isospin_isotensor_count.py
==========================
Isospin mass splittings from mediator-occupancy counting on the FCC lattice.

The two model-independent observables
-------------------------------------
Within a multiplet the splitting carries a part linear in I_3 (the
screw-versus-edge character cost) and a part bilinear in the quark charges
(the Coulomb-type term on charge-active links). Two combinations in the
light spectrum cancel the linear part exactly and expose the bilinear term
alone:

  (1) the Sigma curvature   Sigma+ + Sigma- - 2 Sigma0   (pure isotensor)
  (2) the pion splitting    pi+- - pi0  (strong part negligible; the 2022
      lattice QCD+QED value 4.534(60) MeV confirms it is essentially EM)

Both are integers in units of the electron mass:

  Sigma curvature = 3.004 m_e        pion splitting = 8.989 m_e

In the framework's currency, where one unit charge product across one bond
costs exactly 1 m_e because e^2/(4 pi eps0 r_e) = m_e c^2, these are link
counts. This script derives the 3 from mediator occupancy at the Y-junction
and tests the consequences across the octet.

The occupancy rule
------------------
A nearest-neighbour pair interacts through its direct bond plus the
common-nearest-neighbour sites on the bond's perpendicular bisector (four of
them in FCC). A mediator site that is already occupied by another charge or
by the junction core is not a free mediator: its bonds are the direct bonds
of OTHER pairs, and counting them again double-counts. The rule:

  N_CB(pair) = 1 (direct) + number of UNOCCUPIED common-NN mediators.

  Isolated cell pair (meson):    4 mediators free          -> N_CB = 5
  Y-junction pair (baryon):      centre + third arm occupy
                                 2 of the 4 mediators      -> N_CB = 3

The 5 is the count that fixes m_d - m_u = 5 m_e; the 3 is the Sigma
curvature. One rule, both integers, no fitted quantity.

Consequences tested below
-------------------------
With the bilinear coefficient fixed at C = 3 m_e per unit q_i q_j for every
octet baryon (same junction geometry), the linear u->d character cost per
multiplet is extracted from the observed splittings. It climbs a ladder of
about 4 m_e per strange quark, so the strange climb lives in the LINEAR
term, not in the bilinear coefficient. Coleman-Glashow then holds because
the ladder is linear in strangeness while C is universal.

Inputs: CODATA m_e, PDG 2024 masses and magnetic moments, FCC geometry.
No fitted parameters anywhere.

Accompanies: M. A. Cox, "The Cosserat Supersolid" (2026).
https://doi.org/10.5281/zenodo.18636501

Author: Mitchell A. Cox, University of the Witwatersrand
"""

import numpy as np
from itertools import combinations, product

# ---------------------------------------------------------------- constants
M_E = 0.51099895069          # electron rest energy [MeV], CODATA 2022

# PDG 2024 masses [MeV]
M = {
    "p":   938.27208943, "n":   939.56542194,
    "Sig+": 1189.37,     "Sig0": 1192.642,    "Sig-": 1197.449,
    "Xi0": 1314.86,      "Xi-":  1321.71,
    "pi+": 139.57039,    "pi0":  134.9768,
    "K0":  497.611,      "K+":   493.677,
}
SIG_ERR = {"Sig+": 0.07, "Sig0": 0.024, "Sig-": 0.030}

# PDG magnetic moments [nuclear magnetons]; Sig0 from the quark-model
# average of its charged partners (unmeasured).
MU = {"p": 2.79284734, "n": -1.91304276,
      "Sig+": 2.458, "Sig-": -1.160, "Sig0": 0.5 * (2.458 - 1.160),
      "Xi0": -1.250, "Xi-": -0.6507}

Q_U, Q_D = 2.0 / 3.0, -1.0 / 3.0


# ------------------------------------------------------- FCC geometry tools
def is_fcc(p, tol=1e-9):
    """FCC site test in a/2 integer coordinates: integer triple, even sum."""
    p = np.asarray(p, dtype=float)
    if np.any(np.abs(p - np.round(p)) > tol):
        return False
    return int(round(p.sum())) % 2 == 0


def common_nn(p1, p2, span=3):
    """All FCC sites at nearest-neighbour distance (sqrt 2 in a/2 units)
    from both p1 and p2, excluding the endpoints themselves."""
    out = []
    for xyz in product(range(-span, span + 1), repeat=3):
        p = np.array(xyz, dtype=float)
        if not is_fcc(p):
            continue
        if np.array_equal(p, p1) or np.array_equal(p, p2):
            continue
        if (abs(np.linalg.norm(p - p1) - np.sqrt(2)) < 1e-9 and
                abs(np.linalg.norm(p - p2) - np.sqrt(2)) < 1e-9):
            out.append(p)
    return out


def n_cb(pair, occupied):
    """Charge-active link count of a pair under the occupancy rule:
    1 direct bond + the common-NN mediators NOT in `occupied`."""
    meds = common_nn(*pair)
    free = [m for m in meds
            if not any(np.array_equal(m, o) for o in occupied)]
    return 1 + len(free), meds, free


# --------------------------------------------------- charge bilinear helper
def pair_bilinear(charges):
    """sum_{i<j} q_i q_j for a charge list."""
    return sum(qa * qb for qa, qb in combinations(charges, 2))


def report_occupancy():
    print("=" * 72)
    print("1. The occupancy rule on the two canonical geometries")
    print("=" * 72)

    # Isolated cell pair: a meson's two nodes, nothing else occupied.
    a, b = np.array([0., 0., 0.]), np.array([1., 1., 0.])
    n, meds, free = n_cb((a, b), occupied=[])
    print(f"\nIsolated cell pair {tuple(int(x) for x in a)}-{tuple(int(x) for x in b)}:")
    print(f"  common-NN mediators: {len(meds)}, occupied: 0, free: {len(free)}")
    print(f"  N_CB = 1 + {len(free)} = {n}   "
          f"(the count behind m_d - m_u = 5 m_e)")

    # Y-junction: three arms mutually NN, centre at the origin.
    centre = np.array([0., 0., 0.])
    arms = [np.array([1., 1., 0.]), np.array([1., 0., 1.]),
            np.array([0., 1., 1.])]
    print("\nY-junction (centre at origin, arms mutually NN):")
    counts = []
    for (i, j) in combinations(range(3), 2):
        k = 3 - i - j
        occ = [centre, arms[k]]
        n, meds, free = n_cb((arms[i], arms[j]), occupied=occ)
        tags = []
        for m in meds:
            if np.array_equal(m, centre):
                tags.append(f"{tuple(int(x) for x in m)} = junction centre [occupied]")
            elif np.array_equal(m, arms[k]):
                tags.append(f"{tuple(int(x) for x in m)} = third arm [occupied]")
            else:
                tags.append(f"{tuple(int(x) for x in m)} free")
        print(f"  pair {i+1}-{j+1}: mediators -> " + "; ".join(tags))
        print(f"           N_CB = 1 + {len(free)} = {n}")
        counts.append(n)
    assert all(c == counts[0] for c in counts)
    print(f"\n  Every junction pair: N_CB = {counts[0]}  "
          f"(the Sigma curvature, see below)")
    return counts[0]


def report_isotensor():
    print("\n" + "=" * 72)
    print("2. The two model-independent isotensor observables")
    print("=" * 72)
    d2 = M["Sig+"] + M["Sig-"] - 2 * M["Sig0"]
    err = np.sqrt(SIG_ERR["Sig+"]**2 + SIG_ERR["Sig-"]**2
                  + (2 * SIG_ERR["Sig0"])**2)
    dpi = M["pi+"] - M["pi0"]
    print(f"\n  Sigma curvature  S+ + S- - 2 S0 = {d2:.4f} +- {err:.3f} MeV"
          f"  =  {d2 / M_E:.4f} +- {err / M_E:.3f} m_e")
    print(f"  3 m_e = {3 * M_E:.4f} MeV  (deviation {100*(d2/(3*M_E)-1):+.2f}%,"
          f" inside 1 sigma)")
    print(f"\n  Pion splitting   pi+- - pi0     = {dpi:.4f} MeV"
          f"            =  {dpi / M_E:.4f} m_e")
    print(f"  9 m_e = {9 * M_E:.4f} MeV  (deviation {100*(dpi/(9*M_E)-1):+.2f}%)")
    print("\n  Both cancel the linear I_3 term exactly: the Sigma curvature by")
    print("  construction, the pion because its strong part is O((m_d-m_u)^2)")
    print("  through pi0-eta mixing, ~0.1 MeV. The pion is spin 0, so no")
    print("  magnetic self-energy enters; the 9 is clean of moments.")
    return d2


def report_octet_ladder(C_units):
    print("\n" + "=" * 72)
    print("3. Octet decomposition with the universal bilinear C = "
          f"{C_units} m_e")
    print("=" * 72)
    C = C_units * M_E

    # Observed splittings and their bilinear content.
    obs = {
        "n - p":        (M["n"] - M["p"],
                         pair_bilinear([Q_U, Q_D, Q_D])
                         - pair_bilinear([Q_U, Q_U, Q_D]), 1),
        "Sig- - Sig+":  (M["Sig-"] - M["Sig+"],
                         pair_bilinear([Q_D, Q_D, Q_D])
                         - pair_bilinear([Q_U, Q_U, Q_D]), 2),
        "Xi-  - Xi0":   (M["Xi-"] - M["Xi0"],
                         pair_bilinear([Q_D, Q_D, Q_D])
                         - pair_bilinear([Q_U, Q_D, Q_D]), 1),
    }
    # NOTE on charges: the strange quark carries the same -1/3 as the down,
    # so within each multiplet the bilinear difference only needs the light
    # charges swapped; writing s as -1/3 above is exact for the differences.

    print(f"\n  {'splitting':14}{'obs [MeV]':>11}{'d(sum qq)':>11}"
          f"{'bilinear':>10}{'linear':>9}{'per swap':>10}")
    linear = {}
    for name, (dm, dqq, n_swap) in obs.items():
        bil = C * dqq
        lin = dm - bil
        linear[name] = lin / n_swap
        print(f"  {name:14}{dm:11.4f}{dqq:11.4f}{bil:10.4f}"
              f"{lin:9.4f}{lin / n_swap / M_E:10.4f}")
    dn = linear["n - p"] / M_E
    ds = linear["Sig- - Sig+"] / M_E
    dx = linear["Xi-  - Xi0"] / M_E
    print(f"\n  Linear u->d cost ladder [m_e]:  N: {dn:.3f}   "
          f"Sigma: {ds:.3f}   Xi: {dx:.3f}")
    print(f"  Ladder steps: {ds - dn:.3f} and {dx - ds:.3f} m_e per "
          f"strange quark")
    print("  -> the strange climb sits in the LINEAR character term,")
    print("     about 4 m_e per hex cap; the bilinear stays flat at 3.")

    cg = (M["Sig-"] - M["Sig+"]) - (M["n"] - M["p"]) - (M["Xi-"] - M["Xi0"])
    print(f"\n  Coleman-Glashow residual: {cg:+.3f} MeV. With C universal it")
    print("  equals the ladder's deviation from exact linearity in s.")
    return dn, ds, dx


def report_magnetic_sensitivity():
    print("\n" + "=" * 72)
    print("4. Magnetic-moment sensitivity (proton-calibrated)")
    print("=" * 72)
    # Field energy of a moment mu over the PN core radius: the companion
    # script nucleon_magnetic_selfenergy.py gives U_M(p) ~ 0.18-0.21 MeV.
    # Calibrate k = U_M / mu^2 on the proton and apply to the hyperons.
    for UMp in (0.18, 0.21):
        k = UMp / MU["p"] ** 2
        mag2nd = k * (MU["Sig+"]**2 + MU["Sig-"]**2 - 2 * MU["Sig0"]**2)
        print(f"\n  U_M(p) = {UMp:.2f} MeV  ->  k = {k:.4f} MeV/mu_N^2")
        print(f"  magnetic part of the Sigma curvature: {mag2nd:+.3f} MeV"
              f" ({mag2nd / M_E:+.2f} m_e)")
        print(f"  Coulomb-only reading of the curvature: "
              f"{(1.535 - mag2nd) / M_E:.2f} m_e")
    print("\n  Read on the lattice, electric and magnetic are one shear-sector")
    print("  field; if the 1 m_e per link covers the whole EM flux of the")
    print("  pair, the raw 3.004 is the right comparison and the moment terms")
    print("  are not separate. If they are separate, the Coulomb-only value")
    print("  drifts to ~2.7 m_e and the integer is partly accidental. The")
    print("  pion 9 (spin 0) is immune to this ambiguity either way.")


def report_pion():
    print("\n" + "=" * 72)
    print("5. The pion 9: what link count the bilinear rule implies")
    print("=" * 72)
    # pi+ = (u dbar): node charges (2/3, 1/3). pi0 = (u ubar - d dbar)/sqrt2:
    # average of (2/3,-2/3) and (1/3,-1/3). Self terms sum q^2 are EQUAL
    # (5/9 both), so they cancel in the difference and only the bilinear
    # survives -- the same structure as the Sigma curvature.
    bil_pip = (2/3) * (1/3)
    bil_pi0 = 0.5 * ((2/3) * (-2/3) + (1/3) * (-1/3))
    dbil = bil_pip - bil_pi0
    print(f"\n  d(q1 q2) between pi+- and pi0: {dbil:.4f}  (= 1/2 exactly)")
    print(f"  observed 9 m_e  ->  implied link count N = 9 / (1/2) = 18")
    print("  18 = 2 N_c^2: the cell pair's full node-channel count, the same")
    print("  2 N_c^2 that divides the quark base mass m_0/(2 N_c^2) and that")
    print("  sets the charm antibonding multiplicity. Candidate, not derived:")
    print("  the meson bilinear runs over every node-channel slot of the")
    print("  pair, where the baryon pair runs over its 3 free links.")
    print("  Open: derive 18 from the pair's gauge-mode geometry.")


def report_falsifiable():
    print("\n" + "=" * 72)
    print("6. Falsifiable checks the rule makes")
    print("=" * 72)
    # Decuplet Sigma*(1385): same Y-junction geometry -> same curvature 3 m_e.
    Sp, S0, Sm = 1382.83, 1383.7, 1387.2
    eSp, eS0, eSm = 0.34, 1.0, 0.5
    d2 = Sp + Sm - 2 * S0
    err = np.sqrt(eSp**2 + eSm**2 + (2 * eS0)**2)
    print(f"\n  Sigma*(1385) curvature: {d2:.2f} +- {err:.2f} MeV "
          f"= {d2 / M_E:.1f} +- {err / M_E:.1f} m_e")
    print("  Prediction: 3 m_e = 1.53 MeV if the decuplet junction shares the")
    print("  octet geometry. Current errors cover it; a sharper Sigma*0 mass")
    print("  is a direct test.")
    print("\n  Kaon: K0 - K+ mixes one u->d swap with the charged kaon's EM.")
    dK = M["K0"] - M["K+"]
    print(f"  K0 - K+- = {dK:.3f} MeV = {dK / M_E:.2f} m_e; the swap cost on")
    print("  the kaon's hex-cap pair is the open quantity it isolates once")
    print("  the EM piece is counted by the same rule.")




def report_audit():
    """Structural audit: the claims the monograph footnote makes, executed.

    1. Ledger disjointness: the fifteen bonds counted by the three junction
       pairs (3 direct + 12 through free mediators) are mutually disjoint and
       disjoint from the core's three arm bonds.
    2. Route identity: the blocked route through the third quark uses exactly
       the other two pairs' direct bonds; the blocked route through the centre
       uses exactly the core's arm bonds.
    3. Exteriority: every free mediator lies outside the 13-node cluster.
    4. Invariance: a rotated and translated junction returns the same counts.
    5. Cap collision: the hex-cap placement of the companion classifier
       (spectral_classifier.py) on the (-1,-1,-1) inactive direction lands one
       layer above the quark face and would occupy all six free mediators,
       one per pair; the other three inactive directions touch none. The
       light-strange count C_ls therefore depends on cap placement, and is
       left open in the monograph.
    """
    print("\n" + "=" * 72)
    print("7. Structural audit (the footnote's claims, executed)")
    print("=" * 72)
    centre = np.array([0., 0., 0.])
    arms = [np.array([1., 1., 0.]), np.array([1., 0., 1.]),
            np.array([0., 1., 1.])]
    core = {frozenset([tuple(centre), tuple(a)]) for a in arms}

    def route_bonds(i, j, m):
        return [frozenset([tuple(i), tuple(m)]), frozenset([tuple(m), tuple(j)])]

    used, ok = set(), True
    for (i, j) in combinations(range(3), 2):
        k = 3 - i - j
        meds = common_nn(arms[i], arms[j])
        free = [m for m in meds if not (np.array_equal(m, centre)
                                        or np.array_equal(m, arms[k]))]
        counted = [frozenset([tuple(arms[i]), tuple(arms[j])])]
        for m in free:
            counted += route_bonds(arms[i], arms[j], m)
        for b in counted:
            ok &= b not in used
            used.add(b)
        via_arm = set(route_bonds(arms[i], arms[j], arms[k]))
        directs = {frozenset([tuple(arms[i]), tuple(arms[k])]),
                   frozenset([tuple(arms[j]), tuple(arms[k])])}
        via_ctr = set(route_bonds(arms[i], arms[j], centre))
        assert via_arm == directs, "third-arm route is not the other directs"
        assert via_ctr <= core, "centre route is not the core's arm bonds"
        for m in free:
            assert np.linalg.norm(m) > np.sqrt(2) + 1e-9, "free mediator on cluster"
    assert ok and not (used & core), "ledger double counts a bond"
    print(f"  ledger: {len(used)} counted bonds, mutually disjoint and "
          f"disjoint from the {len(core)} core bonds  [PASS]")
    print("  blocked routes are identically other ledger lines  [PASS]")
    print("  all free mediators lie outside the 13-node cluster  [PASS]")

    R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], float)
    T = np.array([2., 0., 0.])
    c2, a2 = R @ centre + T, [R @ a + T for a in arms]
    for (i, j) in combinations(range(3), 2):
        k = 3 - i - j
        meds = common_nn(a2[i], a2[j])
        free = [m for m in meds if not (np.allclose(m, c2)
                                        or np.allclose(m, a2[k]))]
        assert len(meds) == 4 and len(free) == 2
    print("  rotated and translated junction: same counts  [PASS]")

    # Hex-cap collision, in a/2 units: 3 FCC sites on the {111} plane at
    # projection -2a/sqrt(3) along the chosen inactive direction, nearest
    # the axis (the construction of spectral_classifier.py).
    def cap(direction):
        d = np.array(direction, float) / np.sqrt(3)
        target = -4.0 / np.sqrt(3)
        sites = [np.array(x, float) for x in
                 product(range(-5, 6), repeat=3) if is_fcc(x)]
        on = [p for p in sites if abs(p @ d - target) < 1e-6]
        on.sort(key=lambda p: np.linalg.norm(p - (p @ d) * d))
        return on[:3]

    freeset = {tuple(int(x) for x in m)
               for (i, j) in combinations(range(3), 2)
               for m in common_nn(arms[i], arms[j])
               if not (np.array_equal(m, centre)
                       or np.array_equal(m, arms[3 - i - j]))}
    print("  cap collision with the six free mediators, by inactive direction:")
    for d in [(-1, -1, -1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]:
        hit = {tuple(int(x) for x in p) for p in cap(d)} & freeset
        print(f"    cap on {d}: blocks {sorted(hit) if hit else 'none'}")
    print("  -> C_ls depends on cap placement; open in the monograph.")


def main():
    n_junction = report_occupancy()
    report_isotensor()
    report_octet_ladder(C_units=n_junction)
    report_magnetic_sensitivity()
    report_pion()
    report_falsifiable()
    report_audit()


if __name__ == "__main__":
    main()
