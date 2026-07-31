"""Site covariance against bond covariance: the last factor in the readout chain.

The readout kernel (readout_coefficient.py) prices an anisotropic BOND-RELATIVE
covariance q_ij = <Du_i Du_j>, built from the difference in displacement across
a bond, Du = u(R) - u(0). The gravitational wave modulates the SITE covariance
<u_i u_j> instead. The two are different quantities: a bond stretches only where
its two ends move differently, so neighbours that jitter in step leave their
bond alone however violently each end moves. This script computes the conversion
between them.

Translational invariance makes the relation exact:

    q_ij(R) = 2 [ <u_i u_j> - <u_i(0) u_j(R)> ].

A passing wave acts on the displacement field as a small traceless linear map
u -> (1 + eps) u, and to first order both the variance and the correlation
respond through the same contraction, so the conversion collapses to one number

    C(R) = dq_ij / d<u_i u_j> = 2 [1 - g(R)],
    g(R) = <u_i(0) u_j(R)> / <u_i u_j>,

with g the normalised displacement correlation at the bond vector. Its limits
say what C measures. Long-wavelength modes move neighbours together (g -> 1, no
bond stretches); zone-scale modes leave neighbours independent (g -> 0, each
bond feels the site amplitude twice, once from each end).

Four computations follow.

[1] THE SINC-SQUARED IDENTITY, symbolically. In the isotropic Debye model the
    framework already uses for <u^2> = 0.514 l^2, the polarisation sum over the
    three degenerate acoustic branches is the identity, so the correlation
    tensor is isotropic and the angular average of exp(i k.R) is sinc(kR).
    Weighting each mode by its zero-point amplitude 1/2w ~ 1/k gives a closed
    form, g(R) = sinc^2(k_D R / 2), verified here with sympy rather than
    asserted. At the nearest-neighbour separation k_D l / 2 = 2.188, so
    g(l) = 0.139: nodes a bond apart are already almost independent.

[2] THE WEIGHTED INTEGRAL. The perturbation does not carry the cloud's own
    weighting. The squeeze reaches the displacement field only through the
    transverse (u, phi) hybridisation, and the gravity-carrying branch loses its
    displacement content as kappa_eff ~ (k l)^2 dies away at long wavelength. The
    readable part of the squeeze therefore sits at short wavelength, precisely
    where neighbours are least correlated. Reweighting the same integral by that
    branch's displacement content gives g_W(l) = 0.070 at the physical
    microinertia j = l^2, hence C = 1.86 and a total gain
    Lambda_site = Lambda_read * C = 24.7.

[3] THE PER-BOND KERNEL CHECK. C multiplies Lambda_read as a scalar only if
    every bond that survives the anisotropic kernel sees the same conversion.
    The one-leg rule says the crossing bonds contribute nothing there, so the
    survivors are all in-slice nearest neighbours at R = l. Running the kernel
    with per-bond conversions rather than a common factor confirms it: the gain
    moves by under two per cent.

[4] THE ANISOTROPY SWEEP. The Debye model is isotropic by construction. A real
    face-centred cubic spectrum is not, and the cleanest bracket on that error
    is to strip the model to its transverse branches alone, which makes the
    correlation depend on the bond's orientation relative to the polarisation.
    Sweeping that spread together with the microinertia over a factor of two
    either side of j = l^2 brackets the gain.

The hybridisation model is the one covariance_readout.py uses for the transfer
ratio, reproduced here so the two scripts cannot drift apart.
"""

import numpy as np
import sympy as sp
from itertools import permutations, product

# --- framework constants ------------------------------------------------------
MU = 1.0                            # shear modulus, units of itself
KC = 2.0 / (np.pi - 2.0)            # kappa_c / mu, from N^2 = kappa_c/[2(mu+kappa_c)] = 1/pi
RHO = 1.0                           # mass density
C = 1.0                             # shear-wave speed, sqrt(mu/rho)
L = 1.0                             # lattice spacing l, the unit of length here

G3 = -7.0                           # gamma_3 = V''' d / V''    (gravity sector)
G4 = 343.0 / 9.0                    # gamma_4 = V'''' d^2 / V''  (Morse consistency)

# Debye radius: one wavevector per primitive cell, V_p = l^3 / sqrt(2) for FCC
K_D = (6.0 * np.pi**2 * np.sqrt(2.0))**(1.0 / 3.0) / L


# =============================================================================
# [1] the sinc-squared identity
# =============================================================================
print("=" * 74)
print("[1] THE UNWEIGHTED CORRELATION: g(R) = sinc^2(k_D R / 2)")
print("=" * 74)

k, R, kd = sp.symbols('k R k_D', positive=True)
# zero-point weight 1/2w ~ 1/k against the k^2 measure; angular average -> sinc(kR)
num = sp.integrate(k**2 * k**-1 * sp.sin(k * R) / (k * R), (k, 0, kd))
den = sp.integrate(k**2 * k**-1, (k, 0, kd))
g_sym = sp.simplify(num / den)
closed = 2 * (1 - sp.cos(kd * R)) / (kd * R)**2
sinc2 = (sp.sin(kd * R / 2) / (kd * R / 2))**2
print(f"    integral    g(R) = {sp.simplify(g_sym)}")
print(f"    equals 2[1-cos(k_D R)]/(k_D R)^2 : "
      f"{sp.simplify(g_sym - closed) == 0}")
print(f"    equals sinc^2(k_D R / 2)         : "
      f"{sp.simplify(sp.expand_trig(closed - sinc2)) == 0}")

g_unweighted = float(sinc2.subs({kd: K_D, R: L}))
print(f"\n    k_D l   = {K_D * L:.4f}      (k_D l / 2 = {K_D * L / 2:.4f})")
print(f"    g(l)    = {g_unweighted:.4f}   -> C = 2[1 - g] = {2 * (1 - g_unweighted):.3f}")
print("    Nodes a bond apart are already almost independent: this is what a")
print("    lattice deep into its own zero-point motion should look like.")


# =============================================================================
# [2] the weighted integral
# =============================================================================
print()
print("=" * 74)
print("[2] THE WEIGHTED CORRELATION: what the squeeze actually reaches")
print("=" * 74)


def gravity_branch(k_val, j_over_l2):
    """Displacement content and frequency of the gravity-carrying branch.

    Per wavevector the transverse sector is a 2x2 problem in the mass-weighted
    coordinates (sqrt(rho) u, sqrt(J) phi), identical to covariance_readout.py:

        H(k) = [[(mu + kc_eff) k^2 / rho,        g(k)              ],
                [g(k),                    (gam k^2 + 4 kc_eff) / J ]],
        g(k) = 2 kc_eff k / sqrt(rho J),   kc_eff = KC (k l)^2,

    with the bending modulus gam tuned to keep both branches gapless at c as
    k -> 0. The eigenvector of the phi-dominant branch splits into u- and
    phi-content; its u-content squared is the fraction of that branch's squeeze
    that lands in the displacement cloud, which is the part light can read.

    Returns (u-content squared, angular frequency).
    """
    J = RHO * j_over_l2 * L**2
    kce = KC * (k_val * L)**2
    gam = C**2 * J - 4 * kce * L**2 / np.pi**2
    a = (MU + kce) * k_val**2 / RHO
    b = (max(gam, 0.05 * C**2 * J) * k_val**2 + 4 * kce) / J
    off = 2 * kce * k_val / np.sqrt(RHO * J)
    w2, V = np.linalg.eigh(np.array([[a, off], [off, b]]))
    br = int(abs(V[1, 0]) < abs(V[1, 1]))       # the phi-dominant branch
    return V[0, br]**2, np.sqrt(w2[br])


def g_weighted(R_val, j_over_l2, n_k=4000):
    """Normalised correlation with each mode weighted by the gravity branch's
    displacement content times its zero-point amplitude."""
    ks = np.linspace(1e-4, K_D, n_k)
    cu2, om = np.array([gravity_branch(kk, j_over_l2) for kk in ks]).T
    wt = ks**2 * cu2 / (2.0 * om)
    return np.trapezoid(wt * np.sinc(ks * R_val / np.pi), ks) / np.trapezoid(wt, ks)


g_W = g_weighted(L, 1.0)
CONV = 2.0 * (1.0 - g_W)
LAM_READ = 239.0 / 18.0
print(f"    unweighted cloud   g(l)   = {g_unweighted:.4f}")
print(f"    weighted squeeze   g_W(l) = {g_W:.4f}   (at the physical j = l^2)")
print(f"    conversion         C      = 2[1 - g_W] = {CONV:.3f}")
print(f"    total gain         Lambda_site = Lambda_read * C = "
      f"{LAM_READ:.4f} * {CONV:.3f} = {LAM_READ * CONV:.2f}")
print("    The squeeze is pushed to short wavelength because the gravity branch")
print("    loses its displacement content as kappa_eff ~ (k l)^2 dies away, and")
print("    short wavelength is where neighbours are least correlated. So the")
print("    conversion rises rather than falls: a squeeze written on the sites")
print("    reaches the bonds nearly in full.")

print(f"\n    {'j / l^2':>9} {'g_W(l)':>9} {'C':>8} {'Lambda_site':>13}")
for j_val in (0.1, 0.5, 1.0, 2.0, 10.0):
    gw = g_weighted(L, j_val)
    print(f"    {j_val:>9} {gw:>9.4f} {2*(1-gw):>8.3f} {LAM_READ*2*(1-gw):>13.2f}")


# =============================================================================
# [3] per-bond kernel check
# =============================================================================
print()
print("=" * 74)
print("[3] PER-BOND KERNEL CHECK: is C a scalar multiplier?")
print("=" * 74)


def d4_bonds():
    """The 24 nearest-neighbour bonds of the D4 lattice, unit length."""
    v = set()
    for p in permutations(range(4), 2):
        for s1, s2 in product((1, -1), repeat=2):
            w = [0, 0, 0, 0]
            w[p[0]], w[p[1]] = s1, s2
            v.add(tuple(w))
    return np.array(sorted(v), float) * np.sqrt(0.5)


BONDS = d4_bonds()


def mu_eff(bonds, e4, z4, q3, conv=None):
    """Probe modulus for polarisation e4, propagation z4, under a SITE
    covariance q3 converted to bond-relative covariance by conv (per bond).

    This is the kernel of readout_coefficient.py with the conversion pulled
    inside the bond loop, so a per-bond C can be tested against a common one.
    """
    tot = 0.0
    for i, n in enumerate(bonds):
        ns = n[:3]                                  # spatial leg
        q = q3 * (1.0 if conv is None else conv[i])
        nqn = ns @ q @ ns
        dr = (np.trace(q) - nqn) / 2.0              # transverse fluctuations lengthen
        dr2 = nqn                                   # longitudinal mean-square stretch
        v2eff = 1.0 + G3 * dr + 0.5 * G4 * dr2      # dressed stiffness, units of V''
        teff = dr + 0.5 * G3 * dr2                  # fluctuation-induced tension
        ne, nz = n @ e4, n @ z4
        tot += v2eff * ne**2 * nz**2 + teff * nz**2 * (1.0 - ne**2)
    return tot


EX = np.array([1.0, 0, 0, 0])
EY = np.array([0, 1.0, 0, 0])
EZ = np.array([0, 0, 1.0, 0])
Q0 = 1e-6                                           # linear-response probe amplitude


def gain(conv=None):
    """(c_x - c_y)/c per unit site-covariance amplitude, for an h_+ pattern."""
    mu0 = mu_eff(BONDS, EX, EZ, np.zeros((3, 3)), conv)
    qp = Q0 * np.diag([1.0, -1.0, 0.0])
    return (mu_eff(BONDS, EX, EZ, qp, conv) - mu_eff(BONDS, EY, EZ, qp, conv)) \
        / (2.0 * mu0 * Q0)


def correlation_grid(transverse_only, n=180):
    """Normalised correlation g_xx(R) on a grid over the Debye sphere.

    With three degenerate branches the polarisation sum is the identity and the
    result is the isotropic sinc^2 above. Keeping only the transverse branches
    inserts the projector (delta_ij - khat_i khat_j), which makes g depend on
    the bond's orientation relative to the polarisation axis: a derived
    anisotropy rather than an assumed one.
    """
    ax = np.linspace(-K_D, K_D, n)
    KX, KY, KZ = np.meshgrid(ax, ax, ax, indexing='ij')
    kmag = np.sqrt(KX**2 + KY**2 + KZ**2)
    inside = (kmag < K_D) & (kmag > 1e-9)
    kx, ky, kz, km = KX[inside], KY[inside], KZ[inside], kmag[inside]
    weight = 1.0 / (2.0 * km)                       # zero-point amplitude at speed c
    proj = (1.0 - (kx / km)**2) if transverse_only else np.ones_like(km)

    def g_of(Rvec):
        phase = np.cos(kx * Rvec[0] + ky * Rvec[1] + kz * Rvec[2])
        return float((weight * proj * phase).sum() / (weight * proj).sum())

    return g_of


g_iso = correlation_grid(False)
g_tra = correlation_grid(True)
legs = np.array([np.linalg.norm(b[:3]) for b in BONDS])
in_slice = legs > 0.99                              # spatial leg of full length l

g_shell_iso = np.array([g_iso(b[:3]) for b in BONDS])
g_shell_tra = np.array([g_tra(b[:3]) for b in BONDS])
print(f"    degenerate 3-branch, in-slice bonds: g_xx = "
      f"{g_shell_iso[in_slice].min():.4f} to {g_shell_iso[in_slice].max():.4f}"
      f"   (analytic {g_unweighted:.4f})")
print(f"    transverse-only,     in-slice bonds: g_xx = "
      f"{g_shell_tra[in_slice].min():.4f} to {g_shell_tra[in_slice].max():.4f}"
      f"   (mean {g_shell_tra[in_slice].mean():.4f})")

# rescale the transverse-only spread onto the weighted mean, then go per bond
spread = g_shell_tra / g_shell_tra[in_slice].mean()
conv_common = np.full(len(BONDS), CONV)
conv_perbond = 2.0 * (1.0 - g_W * spread)

lam_bare = gain(None)
lam_common = gain(conv_common)
lam_perbond = gain(conv_perbond)
print(f"\n    Lambda_read  (bond covariance, no conversion) = {lam_bare:.4f}"
      f"   [239/18 = {LAM_READ:.4f}]")
print(f"    Lambda_site  (common C = {CONV:.2f})              = {lam_common:.2f}")
print(f"    Lambda_site  (per-bond C, {conv_perbond.min():.2f} to "
      f"{conv_perbond.max():.2f})       = {lam_perbond:.2f}")
print(f"    difference = {100 * (lam_perbond / lam_common - 1):+.1f} per cent")
print("    The crossing bonds drop out of the anisotropic kernel by the one-leg")
print("    rule, so every surviving bond is an in-slice nearest neighbour at")
print("    R = l, and cubic symmetry makes the twelve equivalent. The conversion")
print("    multiplies Lambda_read as a scalar.")


# =============================================================================
# [4] anisotropy and mode-weighting sweep
# =============================================================================
print()
print("=" * 74)
print("[4] SWEEP: how far can the gain move?")
print("=" * 74)
print(f"    {'j / l^2':>9} {'isotropic':>11} {'per-bond aniso':>16}")
sweep = []
for j_val in (0.5, 0.75, 1.0, 1.5, 2.0):
    gw = g_weighted(L, j_val)
    iso = gain(np.full(len(BONDS), 2.0 * (1.0 - gw)))
    ani = gain(2.0 * (1.0 - gw * spread))
    sweep += [iso, ani]
    print(f"    {j_val:>9} {iso:>11.2f} {ani:>16.2f}")
print(f"\n    bracket: {min(sweep):.1f} to {max(sweep):.1f}"
      f"   (central value {LAM_READ * CONV:.1f})")
print("    The last factor in the radiative normalisation is a computed number")
print("    with a bracket of about fifteen per cent, not an unknown of order one.")
