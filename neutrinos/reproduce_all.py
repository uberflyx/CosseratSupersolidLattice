"""Reproduce every neutrino-sector prediction from four measured inputs.

Inputs: the three charged-lepton masses and the fine-structure constant.
Nothing else enters.  The script computes the neutrino mass spectrum, both
mass-squared splittings, the mass sum, all three mixing angles, the CP
phase, and the Jarlskog invariant, then prints each against its sharpest
published anchor (mid-2026).  It closes with finite-difference input
sensitivities, so the reader can see how each prediction moves when an
input moves.

Chain (all derived in the monograph):
  theta_ch = alpha^2/(2 pi)            lattice chirality angle
  m1   = alpha^3 m_e/(4 pi^2)          lightest mass (theta_ch^2 x parent scale)
  delta = (1-alpha)^2/pi^2             per-leg screened chiral phase
  rho  = h/A = (pi+1)/(2 pi)           hopping ratio (two-channel overlap)
  m_k  = A [1 + 2 rho cos(2pi/3 + delta + 2pi k/3)]
  U    = exp(iA_gen) U_TBM             flavour-row generator, closed form:
         z12 = sqrt(m_e/m_mu) omega,   omega = exp(2i pi/3)
         |z13| = sqrt(3 m_mu/(8 m_tau)),  psi = pi
         gamma = 2 t_cross |z13| - m_mu/(2 m_tau),  t_cross = 1/(2 sqrt 3)
"""

import numpy as np
from scipy.linalg import expm

# ----- the four inputs (PDG / CODATA) --------------------------------------
M_E, M_MU, M_TAU = 0.51099895069, 105.6583755, 1776.86   # MeV
ALPHA = 1.0 / 137.035999177

OMEGA = np.exp(2j * np.pi / 3)
S2, S3, S6 = np.sqrt(2), np.sqrt(3), np.sqrt(6)
U_TBM = np.array([[2/S6, 1/S3, 0], [-1/S6, 1/S3, 1/S2], [-1/S6, 1/S3, -1/S2]])


def predict(m_e=M_E, m_mu=M_MU, m_tau=M_TAU, alpha=ALPHA):
    """Full prediction set from the four inputs.  Returns a dict."""
    # --- spectrum -----------------------------------------------------------
    m1 = alpha**3 * m_e * 1e9 / (4 * np.pi**2)        # meV
    delta = (1 - alpha)**2 / np.pi**2                 # chiral phase [rad]
    rho = (np.pi + 1) / (2 * np.pi)                   # h/A
    phi = 2 * np.pi / 3 + delta
    A = m1 / (1 + 2 * rho * np.cos(phi))              # meV
    masses = np.sort([A * (1 + 2 * rho * np.cos(phi + 2*np.pi*k/3))
                      for k in range(3)])             # m1, m2, m3 [meV]

    # --- mixing -------------------------------------------------------------
    theta_e = np.sqrt(m_e / m_mu)
    z12 = theta_e * OMEGA
    z13 = np.sqrt(3 * m_mu / (8 * m_tau)) * np.exp(1j * np.pi)   # psi = pi
    gamma = np.sqrt(m_mu / (8 * m_tau)) - m_mu / (2 * m_tau)
    a = -0.5j * (z12 + S2 * z13)
    b = -0.5j * (z12 - S2 * z13)
    A_gen = np.array([[0, a, b],
                      [np.conj(a), 0, 1j * gamma],
                      [np.conj(b), -1j * gamma, 0]])
    U = expm(1j * A_gen) @ U_TBM

    s13 = abs(U[0, 2])**2
    s12 = abs(U[0, 1])**2 / (1 - s13)
    s23 = abs(U[1, 2])**2 / (1 - s13)
    J = np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0]))
    denom = (np.sqrt(s12*(1-s12)) * np.sqrt(s23*(1-s23))
             * np.sqrt(s13) * (1 - s13))
    dcp = 180.0 - np.degrees(np.arcsin(np.clip(J / denom, -1, 1)))

    m1_, m2_, m3_ = masses
    return {
        "m1 [meV]": m1_, "m2 [meV]": m2_, "m3 [meV]": m3_,
        "Sum m_nu [meV]": masses.sum(),
        "Dm2_21 [1e-5 eV^2]": (m2_**2 - m1_**2) * 1e-1,
        "Dm2_31 [1e-3 eV^2]": (m3_**2 - m1_**2) * 1e-3,
        "sin2_th12": s12, "sin2_th23": s23, "sin2_th13": s13,
        "delta_CP [deg]": dcp, "J": J, "m_bb [meV]": 0.0,   # Dirac
    }


ANCHORS = {   # (best fit, sigma+, sigma-, source)
    "sin2_th12":          (0.3092, 0.0087, 0.0087, "JUNO 2025"),
    "Dm2_21 [1e-5 eV^2]": (7.50,   0.12,   0.12,   "JUNO 2025"),
    "sin2_th13":          (0.02195, 0.00054, 0.00058, "NuFit 6.0"),
    "sin2_th23":          (0.56,   0.03,   0.05,   "T2K+NOvA"),
    "delta_CP [deg]":     (177.0,  19.0,   20.0,   "NuFit 6.0 (NO)"),
    "Dm2_31 [1e-3 eV^2]": (2.534,  0.025,  0.023,  "NuFit 6.0"),
}

if __name__ == "__main__":
    p = predict()
    print("Predictions from (m_e, m_mu, m_tau, alpha) alone:\n")
    chi2 = 0.0
    for key, val in p.items():
        line = f"  {key:<20} = {val:10.5g}"
        if key in ANCHORS:
            best, sp, sm, src = ANCHORS[key]
            pull = (val - best) / (sp if val > best else sm)
            chi2 += pull**2
            line += f"   vs {best:g} ({src}): pull {pull:+.2f} sigma"
        print(line)
    print(f"\n  joint chi^2 = {chi2:.2f} for 6 observables, 0 parameters")

    # --- input sensitivities: d(obs)/d(ln input), i.e. % shift per % --------
    print("\nSensitivities (% change in observable per 1% change in input):")
    obs = ["Dm2_21 [1e-5 eV^2]", "Dm2_31 [1e-3 eV^2]", "sin2_th12",
           "sin2_th23", "sin2_th13", "delta_CP [deg]", "Sum m_nu [meV]"]
    inputs = {"m_e": "m_e", "m_mu": "m_mu", "m_tau": "m_tau", "alpha": "alpha"}
    header = f"  {'':<20}" + "".join(f"{k:>9}" for k in inputs)
    print(header)
    base = predict()
    for o in obs:
        row = f"  {o:<20}"
        for k in inputs:
            kw = {k: dict(m_e=M_E, m_mu=M_MU, m_tau=M_TAU, alpha=ALPHA)[k] * 1.01}
            shifted = predict(**kw)
            row += f"{100*(shifted[o]/base[o] - 1):>9.2f}"
        print(row)
