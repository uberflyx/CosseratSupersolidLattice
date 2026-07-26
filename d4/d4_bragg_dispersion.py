"""
d4_bragg_dispersion.py
======================

Velocity anomaly of a spatial photon from elastic Bragg scattering off the
four-dimensional D_4 vacuum lattice.

Physics
-------
The vacuum lattice of the framework is the D_4 root lattice in 4D, whose
spatial slice is FCC. A spatial photon Bragg-scattering elastically transfers
a four-momentum G that must satisfy two conditions:

    (1) G in D_4* (the reciprocal of the underlying 4D lattice), and
    (2) G perpendicular to the compact direction e_c = (1,1,1,1)/2,
        so that the scattered photon remains in the spatial sector.

The set of allowed G is therefore the in-plane sublattice of D_4*, which
consists of two cosets:

    integer  Z^4 with sum = 0  (smallest: 12 vectors at coord length sqrt(2))
    centred  (Z + 1/2)^4 with sum = 0  (smallest: 6 vectors at coord length 1)

Mapping coord units to physical lengths: the D_4 minimal vectors have coord
length sqrt(2) and physical length ell (the FCC nearest-neighbour distance,
identified with r_e via the bootstrap m_0 c^2 = hbar c / ell). The reciprocal
coord unit is therefore 2*pi*sqrt(2) / ell, and the smallest in-plane D_4*
vectors have physical magnitude

    |G_min| = 2*pi*sqrt(2) / ell  (six vectors).

The Bragg amplitude at each shell is suppressed by the product of three
geometric factors:

    |F(G)|^2     = exp(-2 |G| w)         Peierls-Nabarro form factor
    exp(-2W_ac)  = exp(-|G|^2 <u^2>/3)   acoustic Debye-Waller
    NCRI factor  = exp(-2W_ac f_s/(1-f_s))  supersolid enhancement

with w = 0.783 ell the PN core width, <u^2> = 0.49 ell^2 the Cosserat-corrected
zero-point displacement, and f_s the superfluid fraction. The bootstrap
self-consistency condition fixes f_s = 4/5.

The velocity anomaly at probe energy E is summed over shells with the
nearly-free-wave perturbative formula

    delta_v_g / c ~ sum_shells [ n_shell * S(G) * (E_BZ_shell / E)^2 ]

valid for E far enough from any zone boundary, where E_BZ_shell = hbar c |G|/2.

Output
------
The script prints a geometry summary, the bootstrap-headline numbers, a
sensitivity sweep over f_s, the GRB-saturation value of f_s, and the FCC
standalone reference (kept as a comparison only; FCC is the right object for
the spatial-density questions that arise elsewhere in the framework, but
not for elastic spatial-to-spatial Bragg).

Companion to the Lorentz invariance derivation in the monograph.
"""

from __future__ import annotations

import math
from itertools import product
from typing import Dict, List, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Framework parameters (all derived; no fits)
# ---------------------------------------------------------------------------

HBAR_C_MEV_FM = 197.3269788          # hbar*c in MeV*fm (CODATA)
ELL_FM = 2.82                        # FCC nearest-neighbour spacing = r_e
M0_C2_MEV = HBAR_C_MEV_FM / ELL_FM   # bootstrap node mass-energy ~70 MeV

W_OVER_ELL = 0.783                   # PN core half-width in units of ell
U2_OVER_ELL2 = 0.49                  # <u^2>/ell^2, Cosserat-corrected Debye

F_S_BOOTSTRAP = 4.0 / 5.0            # bootstrap fixed-point superfluid fraction
F_S_CONSERVATIVE = 0.65              # conservative sensitivity reference

# GRB 090510 reference: photon at observed 10 GeV, Fermi-LAT bound 0.86 s.
# The chapter's perturbative integral over the FRW cosmology was performed
# numerically; we reuse its linear scaling delta_t ~ delta_v/c here.
GRB_REF_DELAY_S = 0.86               # observational upper bound (Fermi-LAT)
GRB_REF_DELTA_V = 5.20e-21           # delta_v/c at 10 GeV that saturates 0.86s

# Probe energies for routine evaluation (in MeV)
E_GAMMA_RAY_MEV = 1.0e4              # 10 GeV (Fermi-LAT band)
E_NEUTRINO_MEV = 6.0e6               # 6 PeV (IceCube Glashow resonance)


# ---------------------------------------------------------------------------
# In-plane D_4* sublattice geometry
# ---------------------------------------------------------------------------

def in_plane_d4_star_shells(max_norm_coord: float = 2.5
                            ) -> Dict[float, List[np.ndarray]]:
    """Enumerate in-plane D_4* lattice points, grouped by coordinate norm.

    Parameters
    ----------
    max_norm_coord
        Maximum coordinate norm to enumerate. Coord units are such that the
        D_4 minimal lattice vectors have length sqrt(2); the smallest in-plane
        D_4* vector has coord length 1.

    Returns
    -------
    dict
        Keys are coord norms (rounded to 6 dp); values are lists of 4-vectors.

    Notes
    -----
    D_4* has two cosets relative to Z^4: the integer points themselves and
    the centred points (Z + 1/2)^4. The in-plane condition G.e_c = 0 selects:

        integer:  sum(G_i) = 0
        centred:  G_i = h_i + 1/2 with sum(h_i) = -2

    The sets are disjoint (integer vs half-integer components), so no
    inter-coset deduplication is required.
    """
    shells: Dict[float, List[np.ndarray]] = {}
    N = int(math.ceil(max_norm_coord)) + 1

    # Coset 1: integer Z^4 with sum = 0 (excluding the zero vector)
    for v in product(range(-N, N + 1), repeat=4):
        if sum(v) != 0:
            continue
        if all(x == 0 for x in v):
            continue
        norm = math.sqrt(sum(x * x for x in v))
        if norm > max_norm_coord:
            continue
        shells.setdefault(round(norm, 6), []).append(np.array(v, dtype=float))

    # Coset 2: centred (Z + 1/2)^4 with sum = 0  =>  h in Z^4 with sum = -2
    for h in product(range(-N - 1, N + 1), repeat=4):
        if sum(h) != -2:
            continue
        v = np.array(h, dtype=float) + 0.5
        norm = float(np.linalg.norm(v))
        if norm > max_norm_coord:
            continue
        shells.setdefault(round(norm, 6), []).append(v)

    return shells


def shell_summary(shells: Dict[float, List[np.ndarray]],
                  n_shells: int = 5) -> List[Tuple[float, int]]:
    """Return [(coord_norm, multiplicity)] for the first n_shells shells."""
    return [(k, len(shells[k])) for k in sorted(shells.keys())[:n_shells]]


# ---------------------------------------------------------------------------
# Bragg suppression and dispersion
# ---------------------------------------------------------------------------

def bragg_suppression(g_norm_coord: float,
                      f_s: float,
                      w_over_ell: float = W_OVER_ELL,
                      u2_over_ell2: float = U2_OVER_ELL2
                      ) -> Dict[str, float]:
    """Return all suppression factors for a single Bragg shell.

    Parameters
    ----------
    g_norm_coord
        Reciprocal-vector coord norm (1 for the smallest in-plane D_4* shell).
    f_s
        Superfluid fraction; must satisfy 0 <= f_s < 1.
    w_over_ell, u2_over_ell2
        PN core width and zero-point displacement in units of ell.

    Returns
    -------
    dict
        Keys: G_ell (= |G| ell), F2 (= |F(G)|^2), DW_ac, DW_total, DW_factor,
        S_G (= F2 * DW_factor).
    """
    if not 0.0 <= f_s < 1.0:
        raise ValueError(f"f_s must be in [0, 1); got {f_s}")

    g_ell = g_norm_coord * 2.0 * math.pi * math.sqrt(2.0)
    f_squared = math.exp(-2.0 * g_ell * w_over_ell)
    dw_ac = (g_ell ** 2) * u2_over_ell2 / 3.0
    dw_total = dw_ac / (1.0 - f_s)
    dw_factor = math.exp(-dw_total)

    return {
        "G_ell": g_ell,
        "F2": f_squared,
        "DW_ac": dw_ac,
        "DW_total": dw_total,
        "DW_factor": dw_factor,
        "S_G": f_squared * dw_factor,
    }


def velocity_anomaly_d4(shells: Dict[float, List[np.ndarray]],
                        f_s: float,
                        e_probe_mev: float,
                        m0_c2_mev: float = M0_C2_MEV
                        ) -> Tuple[float, List[Dict]]:
    """Total spatial-photon velocity anomaly summed over in-plane D_4* shells.

    Implements the perturbative second-order Bloch formula

        delta_v_g/c = sum_shells n_shell * S(G) * (E_BZ_shell / E)^2

    valid for E far from every shell's zone boundary E_BZ_shell.

    Returns
    -------
    delta_v_total : float
    contributions : list of per-shell dicts
    """
    contributions: List[Dict] = []
    delta_v_total = 0.0

    for norm in sorted(shells.keys()):
        n_in_shell = len(shells[norm])
        sup = bragg_suppression(norm, f_s)
        e_bz_mev = m0_c2_mev * sup["G_ell"] / 2.0
        if e_probe_mev > e_bz_mev:
            contribution = n_in_shell * sup["S_G"] * (e_bz_mev / e_probe_mev) ** 2
        else:
            contribution = 0.0
        delta_v_total += contribution
        contributions.append(
            {"norm": norm, "n_shell": n_in_shell, "G_ell": sup["G_ell"],
             "F2": sup["F2"], "DW_factor": sup["DW_factor"], "S_G": sup["S_G"],
             "E_BZ_MeV": e_bz_mev, "delta_v": contribution}
        )

    return delta_v_total, contributions


def velocity_anomaly_fcc_reference(f_s: float,
                                   e_probe_mev: float,
                                   m0_c2_mev: float = M0_C2_MEV
                                   ) -> Dict[str, float]:
    """Standalone-FCC reference: 8 vectors at |G_111| = pi sqrt(6)/ell.

    This is the calculation one obtains by treating the spatial slice as an
    isolated 3D FCC crystal and ignoring the four-dimensional kinematic
    constraint. It is included for comparison only; for elastic spatial-to-
    spatial Bragg, the in-plane D_4* result is the correct one.
    """
    g_ell_fcc = math.pi * math.sqrt(6.0)
    n_fcc = 8
    f2 = math.exp(-2.0 * g_ell_fcc * W_OVER_ELL)
    dw_ac = (g_ell_fcc ** 2) * U2_OVER_ELL2 / 3.0
    dw_factor = math.exp(-dw_ac / (1.0 - f_s))
    s_g = f2 * dw_factor
    e_bz_mev = m0_c2_mev * g_ell_fcc / 2.0
    delta_v = (n_fcc * s_g * (e_bz_mev / e_probe_mev) ** 2
               if e_probe_mev > e_bz_mev else 0.0)
    return {"G_ell": g_ell_fcc, "n": n_fcc, "F2": f2, "DW_factor": dw_factor,
            "S_G": s_g, "E_BZ_MeV": e_bz_mev, "delta_v": delta_v}


def grb_time_delay(delta_v_at_10gev: float) -> float:
    """Estimated GRB 090510 differential delay (seconds) under linear scaling."""
    return GRB_REF_DELAY_S * (delta_v_at_10gev / GRB_REF_DELTA_V)


def saturation_fs(shells: Dict[float, List[np.ndarray]],
                  e_probe_mev: float = E_GAMMA_RAY_MEV
                  ) -> float:
    """Find f_s at which the predicted GRB delay just reaches 0.86 s.

    Uses a binary search since delta_v(f_s) is monotonic in f_s.
    """
    f_lo, f_hi = 0.01, 0.99
    target = GRB_REF_DELTA_V
    for _ in range(80):
        f_mid = 0.5 * (f_lo + f_hi)
        dv, _ = velocity_anomaly_d4(shells, f_mid, e_probe_mev)
        if dv > target:
            f_lo = f_mid    # too much violation; need larger f_s
        else:
            f_hi = f_mid
    return 0.5 * (f_lo + f_hi)


def sensitivity_sweep(shells: Dict[float, List[np.ndarray]],
                      f_s_values,
                      e_probe_mev: float = E_GAMMA_RAY_MEV
                      ) -> List[Dict[str, float]]:
    """Compute delta_v/c and GRB delay across a sweep of f_s values."""
    rows = []
    for f_s in f_s_values:
        dv, contribs = velocity_anomaly_d4(shells, f_s, e_probe_mev)
        s_first = contribs[0]["S_G"] if contribs else float("nan")
        rows.append({"f_s": f_s, "S_first": s_first, "delta_v": dv,
                     "grb_delay_s": grb_time_delay(dv)})
    return rows


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def _hr(char: str = "=", n: int = 78) -> str:
    return char * n


def report() -> None:
    print(_hr())
    print(" D_4 BRAGG DISPERSION")
    print(" Velocity anomaly of a spatial photon, in-plane D_4* sublattice")
    print(_hr())
    print()
    print(f"  Lattice spacing ell    = {ELL_FM} fm  (= r_e)")
    print(f"  Node mass m_0 c^2      = {M0_C2_MEV:.3f} MeV  "
          f"(bootstrap: hbar c / ell)")
    print(f"  PN core w/ell          = {W_OVER_ELL}")
    print(f"  <u^2>/ell^2 (Cosserat) = {U2_OVER_ELL2}")
    print(f"  Probe energy           = {E_GAMMA_RAY_MEV / 1000:.1f} GeV "
          f"(Fermi-LAT band)")
    print()

    shells = in_plane_d4_star_shells(max_norm_coord=2.5)
    print(_hr("-"))
    print("  In-plane D_4* shells (coord norm, multiplicity):")
    for norm, mult in shell_summary(shells, n_shells=5):
        g_ell = norm * 2.0 * math.pi * math.sqrt(2.0)
        e_bz = M0_C2_MEV * g_ell / 2.0
        print(f"    norm = {norm:6.4f}   n = {mult:3d}   "
              f"|G| ell = {g_ell:7.4f}   E_BZ = {e_bz:6.1f} MeV")
    print()

    # ---- Bootstrap headline ------------------------------------------------
    f_s = F_S_BOOTSTRAP
    dv_d4, contribs = velocity_anomaly_d4(shells, f_s, E_GAMMA_RAY_MEV)
    print(_hr("-"))
    print(f"  BOOTSTRAP HEADLINE (f_s = {f_s} = 4/5)")
    print(_hr("-"))
    header = (f"    {'shell':<8}{'|G|*ell':>10}{'n':>4}"
              f"{'|F|^2':>13}{'DW':>13}{'S(G)':>13}"
              f"{'E_BZ':>9}{'delta_v':>14}")
    print(header)
    for c in contribs[:3]:
        print(f"    {c['norm']:<8.4f}{c['G_ell']:>10.4f}{c['n_shell']:>4d}"
              f"{c['F2']:>13.3e}{c['DW_factor']:>13.3e}{c['S_G']:>13.3e}"
              f"{c['E_BZ_MeV']:>9.1f}{c['delta_v']:>14.3e}")
    print(f"\n  Total delta_v/c at 10 GeV     = {dv_d4:.3e}")
    print(f"  Implied GRB 090510 delay      = {grb_time_delay(dv_d4):.3e} s")
    print(f"  Fermi-LAT bound on this delay = {GRB_REF_DELAY_S} s")
    margin = GRB_REF_DELAY_S / max(grb_time_delay(dv_d4), 1e-300)
    print(f"  Margin below bound            = {margin:.2e}x")
    print()

    # ---- IceCube cross-check at 6 PeV --------------------------------------
    dv_icecube, _ = velocity_anomaly_d4(shells, f_s, E_NEUTRINO_MEV)
    print(f"  IceCube Glashow at 6 PeV (bootstrap f_s):")
    print(f"    delta_v/c = {dv_icecube:.3e}")
    print()

    # ---- Sensitivity sweep -------------------------------------------------
    print(_hr("-"))
    print("  SENSITIVITY: delta_v/c at 10 GeV as a function of f_s")
    print(_hr("-"))
    sweep = sensitivity_sweep(
        shells,
        f_s_values=[0.50, 0.53, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85],
    )
    print(f"    {'f_s':>5}{'S_first_shell':>18}{'delta_v/c':>15}"
          f"{'GRB delay (s)':>18}{'verdict':>16}")
    for r in sweep:
        if r["grb_delay_s"] > GRB_REF_DELAY_S:
            verdict = "FAIL"
        elif r["grb_delay_s"] > 0.5 * GRB_REF_DELAY_S:
            verdict = "saturating"
        elif r["grb_delay_s"] > 1e-6:
            verdict = f"safe ({GRB_REF_DELAY_S / r['grb_delay_s']:.0e}x)"
        else:
            verdict = "undetectable"
        print(f"    {r['f_s']:>5.2f}{r['S_first']:>18.3e}{r['delta_v']:>15.3e}"
              f"{r['grb_delay_s']:>18.3e}{verdict:>16}")
    print()

    f_sat = saturation_fs(shells)
    print(f"  GRB-saturation f_s (delay = {GRB_REF_DELAY_S} s): {f_sat:.4f}")
    print(f"  Bootstrap f_s = {F_S_BOOTSTRAP:.4f} sits "
          f"{(F_S_BOOTSTRAP - f_sat) * 100:+.1f} percentage points above this.")
    print()

    # ---- FCC standalone reference (comparison only) ------------------------
    print(_hr("-"))
    print("  STANDALONE-FCC REFERENCE (comparison only)")
    print("  (Treats the spatial slice as a 3D FCC crystal in isolation.")
    print("   Eight vectors at |G_111| = pi sqrt(6)/ell. The corresponding")
    print("   scattering channel is kinematically forbidden in D_4 because it")
    print("   would transfer compact-direction momentum.)")
    print(_hr("-"))
    for f_s_label, f_s_val in [("conservative", F_S_CONSERVATIVE),
                                ("bootstrap   ", F_S_BOOTSTRAP)]:
        fcc = velocity_anomaly_fcc_reference(f_s_val, E_GAMMA_RAY_MEV)
        dv_d4_at_fs, _ = velocity_anomaly_d4(shells, f_s_val, E_GAMMA_RAY_MEV)
        print(f"    f_s = {f_s_val:.2f} ({f_s_label}):  "
              f"FCC delta_v/c = {fcc['delta_v']:.3e},  "
              f"D_4 delta_v/c = {dv_d4_at_fs:.3e},  "
              f"ratio D_4/FCC = {dv_d4_at_fs / fcc['delta_v']:.3e}")
    print()
    print(_hr())


if __name__ == "__main__":
    report()
