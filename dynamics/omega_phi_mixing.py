"""One mixing angle for the omega and the phi leptonic widths.

THE PROBLEM.  Ideal mixing puts the omega and the phi in pure quark states,
omega = (u ubar + d dbar)/sqrt2 and phi = s sbar, so their photon couplings follow from
quark charges alone: relative to the rho they carry 1/9 and 2/9.  Scaled off the
measured rho, that puts the omega 21% high and the phi 6% low.  The signs are opposite,
which no rescaling can fix and which is the signature of a rotation rather than an error
of size.

THE FIX.  The physical states are rotated away from ideal by a small angle delta, so

    f_omega / f_omega(ideal) = cos d - sqrt2 sin d
    f_phi   / f_phi(ideal)   = cos d + sin d / sqrt2

to first order, with the same d in both and opposite effect, which is exactly the
pattern the two residuals show.  Gamma_ee goes as f^2.

Two measurements, one parameter, so the fit can fail; it does not.  The two isoscalars
agree on the angle within their errors and the joint fit brings both below one standard
deviation.  The resulting nonet angle also lands on the standard phenomenological value,
which is a check the fit had no way to arrange.

SCOPE.  The isoscalars are scaled off the MEASURED rho here rather than off the
framework's own f_pi.  That is deliberate and it is not free: see vector_normalisation.py,
where the difference between the two routes is shown to be the rho's own leptonic excess
and is used to establish that the excess belongs to the whole light vector sector.  The
angle itself is fitted, not derived; deriving it needs the off-diagonal element of the
isoscalar mass matrix, which the framework does not yet supply.
"""
import numpy as np
from scipy import optimize

GAMMA_RHO_EE = 7.04          # keV, PDG constrained fit
M_RHO = 775.26               # MeV

# name: (measured Gamma_ee [keV], error, mass [MeV], charge weight relative to the rho)
ISO = {"omega": (7.38e-5 * 8.68e3, 0.22e-5 * 8.68e3, 782.66, 1 / 9),
       "phi": (2.979e-4 * 4.249e3, 0.033e-4 * 4.249e3, 1019.461, 2 / 9)}

R = {"omega": lambda d: np.cos(d) - np.sqrt(2) * np.sin(d),
     "phi": lambda d: np.cos(d) + np.sin(d) / np.sqrt(2)}


def ideal(name):
    """Ideal-mixing leptonic width, scaled off the measured rho by charge weight."""
    _, _, m, w = ISO[name]
    return w * GAMMA_RHO_EE * M_RHO / m


print("=== ideal mixing, scaled off the measured rho ===")
for k in ISO:
    o, e, m, w = ISO[k]
    print("  %-6s weight %.4f  ideal %6.4f keV   measured %6.4f +- %.4f   %+5.1f%%  %5.1f sigma"
          % (k, w, ideal(k), o, e, 100 * (ideal(k) / o - 1), (ideal(k) - o) / e))
print("  Opposite signs, so this is a rotation and not a mis-scaling.\n")

print("=== the angle each isoscalar asks for on its own ===")
for k in ISO:
    o, e, _, _ = ISO[k]
    d = optimize.brentq(lambda x: ideal(k) * R[k](x)**2 - o, -0.3, 0.3)
    dp = optimize.brentq(lambda x: ideal(k) * R[k](x)**2 - (o + e), -0.3, 0.3)
    print("  %-6s delta = %.2f +- %.2f degrees" % (k, np.degrees(d), abs(np.degrees(dp - d))))
print("  They agree within their errors, having been given no chance to.\n")


def chi2(d):
    return sum(((ideal(k) * R[k](d)**2 - ISO[k][0]) / ISO[k][1])**2 for k in ISO)


fit = optimize.minimize_scalar(chi2, bounds=(-0.2, 0.2), method="bounded")
d_fit = fit.x
print("=== joint fit ===")
print("  delta = %.2f degrees, chi2 = %.2f on one degree of freedom"
      % (np.degrees(d_fit), fit.fun))
for k in ISO:
    o, e, _, _ = ISO[k]
    print("  %-6s %5.1f sigma ideal  ->  %.2f sigma mixed" % (k, (ideal(k) - o) / e,
                                                              (ideal(k) * R[k](d_fit)**2 - o) / e))

theta_ideal = np.arctan(1 / np.sqrt(2))
print("\n=== the check the fit could not arrange ===")
print("  nonet angle = theta_ideal + delta = %.2f + %.2f = %.2f degrees"
      % (np.degrees(theta_ideal), np.degrees(d_fit), np.degrees(theta_ideal + d_fit)))
print("  against the standard phenomenological value near 39 degrees.  Nothing in the")
print("  two leptonic widths knows about that number.")
