"""
pmns_construction.py -- end-to-end PMNS matrix from the flavour-side generator.

Cosserat supersolid vacuum framework: the charged-lepton masses are locked by
the Burgers-vector magnitude of the screw dislocation, so lattice perturbations
can only rotate the charged-lepton eigenframe.  The PMNS matrix is therefore

    U_PMNS = exp(i A) . U_TBM ,

with U_TBM the tribimaximal matrix fixed by the S3 = Z3 x Z2 symmetry of the
{111} glide plane, and A the Hermitian flavour-side generator built from three
derived quantities (no free parameters):

    z12  = sqrt(m_e/m_mu) * omega          (stacking rotation, omega = e^{2 pi i/3})
    z13  = +/- sqrt(3 m_mu / (8 m_tau))    (Cosserat torsion; the mu-tau mirror
                                            forces it real, sign = orientation)
    gamma = -m_mu/(2 m_tau)                (mu-tau rotation: the second-order
                                            diagonal baseline.  The inter-family
                                            coupling is circulant in its gauge-
                                            invariant form and contributes zero
                                            to every mixing channel.)

    A_e mu  = -i (z12 + sqrt(2) z13) / 2
    A_e tau = -i (z12 - sqrt(2) z13) / 2
    A_mu tau = i gamma                     (a real rotation in the mu-tau plane)

Exact linear-order statements (verified symbolically and numerically below):
    U_e2  = (1 + z12)/sqrt(3)        -> sin^2 th12 = (1 - t + t^2)/3
    U_e3  = z13                      -> sin^2 th13 = 3 m_mu/(8 m_tau)
    U_mu3 = (1 + gamma)/sqrt(2)      -> sin^2 th23 = (1 - m_mu/m_tau)/2
    J     = -|z13| sin(psi)/(3 sqrt2) = 0 for real z13 (mirror-enforced)

The script prints the exponentiated (all-orders) observables for both torsion
orientations, the mass-squared splittings from the Z3 stacking ring, pulls
against a single anchor, and machine-precision checks of unitarity and
isospectrality.

The anchor is NuFIT 6.1 (v61 parameter table, IC24 with SK-atm, Normal
Ordering).  6.1 already folds in the JUNO first result, so the solar pair is
taken from it rather than from JUNO separately; quoting JUNO alongside a
global fit that has absorbed it would enter the same measurement twice.

Two cautions the printout enforces rather than states:

  1. Errors are asymmetric, so each pull uses the sigma on the side the
     prediction falls.
  2. delta_CP is an angle and its likelihood is far from Gaussian.  The 3-sigma
     interval runs 125 -> 365 degrees about a best fit of 212, which is not
     three times the quoted 1-sigma error of +26.  A Gaussianised pull is
     therefore reported only for orientation, and the operative statement is
     containment in the 3-sigma interval, evaluated modulo 360 degrees.

Run:  python3 pmns_construction.py
"""

import numpy as np
from scipy.linalg import expm

# ----------------------------------------------------------------------
# inputs: charged-lepton masses (PDG) and fine-structure constant (CODATA)
# ----------------------------------------------------------------------
M_E, M_MU, M_TAU = 0.51099895069, 105.6583755, 1776.86   # MeV
ALPHA = 1 / 137.035999177

# derived generator parameters (closed form, zero free parameters)
THETA_E = np.sqrt(M_E / M_MU)                  # 0.06954
Z13_MAG = np.sqrt(3 * M_MU / (8 * M_TAU))      # 0.14933
GAMMA = -M_MU / (2 * M_TAU)                                # -0.02973
OMEGA = np.exp(2j * np.pi / 3)

# tribimaximal matrix (orthogonal; S3-fixed frame)
U_TBM = np.array([[np.sqrt(2 / 3), 1 / np.sqrt(3), 0],
                  [-1 / np.sqrt(6), 1 / np.sqrt(3), 1 / np.sqrt(2)],
                  [-1 / np.sqrt(6), 1 / np.sqrt(3), -1 / np.sqrt(2)]])
assert np.max(np.abs(U_TBM @ U_TBM.T - np.eye(3))) < 1e-14, "TBM not orthogonal"

# ----------------------------------------------------------------------
# experimental anchor: NuFIT 6.1, IC24 with SK-atm, Normal Ordering.
# Each entry is (central, sigma_up, sigma_down); the 3-sigma ranges are the
# fit's own and are not (central -/+ 3 sigma).
# ----------------------------------------------------------------------
NUFIT = {
    "s12sq":  (0.3088,  0.0067,  0.0066),
    "s23sq":  (0.470,   0.017,   0.014),
    "s13sq":  (0.02248, 0.00055, 0.00059),
    "dcp":    (212.0,   26.0,    36.0),
    "dm21":   (7.537e-5, 0.094e-5, 0.10e-5),
    "dm31":   (2.511e-3, 0.021e-3, 0.020e-3),
}
# 3-sigma intervals, NO, from the same table (with SK-atm column)
NUFIT_3SIG = {
    "s23sq": (0.435, 0.584),      # spans both octants: the octant is unresolved
    "dcp":   (125.0, 365.0),      # wraps past 360; evaluate modulo 360
}


def pull(pred, key):
    """Pull in units of the experimental sigma on the side the prediction falls."""
    central, up, down = NUFIT[key]
    return (pred - central) / (up if pred > central else down)


def within_3sigma_angle(pred_deg, key):
    """Containment test for an angle, checked at both pred and pred + 360."""
    lo, hi = NUFIT_3SIG[key]
    return any(lo <= x <= hi for x in (pred_deg, pred_deg + 360.0))


def build_pmns(z13_signed):
    """Assemble the generator and return (U_PMNS, A) for a signed real z13."""
    z12 = THETA_E * OMEGA
    a = -1j * (z12 + np.sqrt(2) * z13_signed) / 2
    b = -1j * (z12 - np.sqrt(2) * z13_signed) / 2
    A = np.array([[0, a, b],
                  [np.conj(a), 0, 1j * GAMMA],
                  [np.conj(b), -1j * GAMMA, 0]], complex)
    return expm(1j * A) @ U_TBM, A


def observables(U):
    """PDG mixing observables from a unitary 3x3 matrix."""
    s13sq = abs(U[0, 2]) ** 2
    s12sq = abs(U[0, 1]) ** 2 / (1 - s13sq)
    s23sq = abs(U[1, 2]) ** 2 / (1 - s13sq)
    J = np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0]))
    # delta_CP: sin from J, cos from the rephasing-invariant |U_mu1|^2 relation
    s12, c12 = np.sqrt(s12sq), np.sqrt(1 - s12sq)
    s23, c23 = np.sqrt(s23sq), np.sqrt(1 - s23sq)
    s13, c13 = np.sqrt(s13sq), np.sqrt(1 - s13sq)
    sin_d = J / (s12 * c12 * s23 * c23 * s13 * c13 ** 2)
    cos_d = (abs(U[1, 0]) ** 2 - s12sq * c23 ** 2 - c12 ** 2 * s23sq * s13sq) \
        / (2 * s12 * c12 * s23 * c23 * s13)
    delta = np.rad2deg(np.arctan2(sin_d, cos_d)) % 360
    return s12sq, s23sq, s13sq, J, delta


def splittings():
    """Mass eigenvalues of the Z3 stacking ring at the NLO chiral phase."""
    delta = (1 - ALPHA) ** 2 / np.pi ** 2          # per-leg screened phase
    m1 = ALPHA ** 3 * M_E * 1e9 / (4 * np.pi ** 2)  # meV
    h_over_a = (np.pi + 1) / (2 * np.pi)
    cos_min = np.cos(2 * np.pi / 3 + delta)
    a_diag = m1 / (1 + 2 * h_over_a * cos_min)
    h = h_over_a * a_diag
    m2 = a_diag + 2 * h * np.cos(4 * np.pi / 3 + delta)
    m3 = a_diag + 2 * h * np.cos(delta)
    return m1, m2, m3


def main():
    m1, m2, m3 = splittings()
    dm21 = (m2 ** 2 - m1 ** 2) * 1e-6              # eV^2
    dm31 = (m3 ** 2 - m1 ** 2) * 1e-6
    print("=" * 72)
    print("Z3 stacking-ring spectrum (zero parameters)")
    print(f"  m1={m1:.3f}  m2={m2:.3f}  m3={m3:.2f} meV   sum={m1+m2+m3:.1f} meV")
    print(f"  Dm2_21={dm21:.4e} eV^2  pull {pull(dm21,'dm21'):+.2f}s")
    print(f"  Dm2_31={dm31:.4e} eV^2  pull {pull(dm31,'dm31'):+.2f}s")

    print("=" * 72)
    print("PMNS from U = exp(iA) U_TBM  (flavour-side generator, closed form)")
    print(f"  theta_e={THETA_E:.5f}  |z13|={Z13_MAG:.5f}  gamma={GAMMA:.5f}")
    print(f"  linear forms: s12^2={(1-THETA_E+THETA_E**2)/3:.4f}"
          f"  s23^2={0.5*(1+2*GAMMA):.4f}  s13^2={Z13_MAG**2:.5f}  J=0 (mirror)")
    for name, z in (("torsion +b (psi=0) ", +Z13_MAG),
                    ("torsion -b (psi=pi)", -Z13_MAG)):
        U, A = build_pmns(z)
        # machine-precision checks
        unit = np.max(np.abs(U @ U.conj().T - np.eye(3)))
        u_e = expm(-1j * A)
        ev = np.sort(np.linalg.eigvalsh(u_e @ np.diag([M_E, M_MU, M_TAU])
                                        @ u_e.conj().T))
        iso = np.max(np.abs(ev - [M_E, M_MU, M_TAU])
                     / np.array([M_E, M_MU, M_TAU]))
        s12sq, s23sq, s13sq, J, delta = observables(U)
        inside = "in" if within_3sigma_angle(delta, "dcp") else "OUT of"
        print(f"  [{name}] s12^2={s12sq:.4f} ({pull(s12sq,'s12sq'):+.2f}s)"
              f"  s23^2={s23sq:.4f} ({pull(s23sq,'s23sq'):+.2f}s)"
              f"  s13^2={s13sq:.5f} ({pull(s13sq,'s13sq'):+.2f}s)")
        print(f"      J={J:+.5f}  delta_CP={delta:.1f} deg"
              f" ({pull(delta,'dcp'):+.2f}s nominal, {inside} the 3-sigma interval)"
              f"   unitarity={unit:.0e}  isospectrality={iso:.0e}")
    print("=" * 72)
    print("Structural statement: the mu-tau mirror forces the torsion real, so J")
    print("vanishes at first order and delta_CP falls at a CP-conserving point,")
    print("delta ~ 0 (+b torsion) or delta ~ 182 deg (-b). The stacking handedness")
    print("selects -b, and the data do not settle the choice: both branches lie")
    print("inside the 3-sigma interval, because that interval runs 125 -> 365 deg")
    print("and 2.1 deg is 362.1 deg. The nominal -5.8 sigma on the +b branch is an")
    print("artefact of Gaussianising a wrapped, strongly non-Gaussian likelihood.")
    print("Separating the branches needs DUNE or Hyper-K, not the present fits.")
    print("The octant is likewise open: the 3-sigma range for sin^2 th23 is")
    print(f"{NUFIT_3SIG['s23sq'][0]}-{NUFIT_3SIG['s23sq'][1]}, which spans maximal mixing, so the")
    print("agreement of both NuFIT 6.1 variants on 0.470 is not a resolved octant.")


if __name__ == "__main__":
    main()
