"""
The two-nucleon potential V(r), derived from lattice mechanics.

Cosserat Supersolid Lattice framework. This script carries out the derivation of
the central nucleon-nucleon potential, with no quantity tuned to the answer. Every
input is fixed elsewhere in the monograph.

  PER-BOND STRENGTH (analytic, not fitted).  The defect-defect interaction is
  V_nm(r) = b_n . G(r) . b_m with G the static Cosserat Green's function. For two
  charge-active (screw) arms it is the Coulomb form V = alpha hbar c / r. At one
  lattice spacing,
        eps_bond = alpha hbar c / ell = alpha m0 c^2 = m_e c^2 = 0.511 MeV,
  using the bootstrap hbar c / ell = m0 c^2 and alpha = m_e / m0. So the m_e in the
  bond-counting rule is DERIVED, the Coulomb energy of unit charges one cell apart.

  DOCKING ATTRACTION (well).  Two coordination shells dock by sharing charge-active
  surface bonds; the number shared grows as the two charge clouds overlap. Modelling
  each cloud as the Peierls-Nabarro squared-Lorentzian density rho ~ (w^2 + r^2)^-2,
  the cloud-cloud overlap is itself a squared-Lorentzian of width 2w (analytic). It
  is normalised so the deuteron's four-bond dock at r = ell gives -4 m_e:
        V_dock(r) = -4 m_e * [ (4w^2 + ell^2) / (4w^2 + r^2) ]^2.
  This is the CENTRAL attraction, the bond count spread over separation.

  PAULI HARD CORE.  Each nucleon carries a dense core of radius R_conf = ell/2pi.
  When two cores overlap (r < 2 R_conf) exclusion forbids the shared volume; the
  cost is M times the overlapped fraction Omega(r) (lens volume of two spheres).
  M = kappa_ex * E_conf, E_conf = 2 sqrt(sigma) = 4 pi m0 ~ 880 MeV, and kappa_ex is
  the one order-unity knob, the high-density stiffness that sets the neutron-star
  maximum mass. Default kappa_ex = 1.

  ONE-PION TAIL (new in this version).  The longest-range piece is one-pion exchange.
  The pion is the two-node stacking fault, m_pi = 2 m0, with decay constant
  f_pi = 3^{1/4} m0, both lattice-fixed. Its range hbar c / m_pi = 1.41 fm exceeds
  the couple-stress range 0.75 fm, so the pion sets the outer edge of V(r). Its
  central (spin-spin) part is the piece a single-channel S-wave can carry:
        V_ope^C(r) = (f^2/4pi)(m_pi/3)(tau1.tau2)(sigma1.sigma2) e^{-x}/x ,  x = r m_pi/hbar c,
  with the pion-nucleon coupling itself lattice-derived through Goldberger-Treiman,
        f^2/4pi = (m_pi g_A / 2 f_pi)^2 / 4pi = g_A^2 / (4 pi sqrt 3) ,
  using m_pi = 2 m0 and f_pi = 3^{1/4} m0; g_A = 1.279 is the derived axial coupling.
  For 3S1 (T=0) and 1S0 (T=1) alike, (tau1.tau2)(sigma1.sigma2) = -3, so the central
  OPE is attractive in both channels.

  WHY THE TAIL MATTERS.  The central docking well alone is deep (-18 MeV) but narrow
  (~1.5 fm), and a square-well estimate needs V0 R^2 > pi^2 hbar^2 / 8 mu ~ 102
  MeV fm^2 to bind a two-nucleon state (mu = M_N/2). The docking gives only ~40 and
  does not bind. The one-pion tail reaches to ~3 fm and widens the well; this script
  tests whether the derived tail closes the gap and binds the deuteron; the singlet
  1S0 channel is the same central force without the tensor, which leaves it just
  unbound (the dineutron virtual state).

  ISOSPIN.  The central force is the SAME in both channels. The deuteron (3S1, T=0)
  and the dineutron (1S0, T=1) share the same docking pocket. In the singlet the two
  spins are opposite, so the two nucleons are distinguishable and may share a docking
  node exactly as the proton and neutron do. What Pauli forbids is the like-nucleon
  TRIPLET 3S1 (same spin and same flavour), which confines two neutrons to the singlet
  1S0. The singlet is the channel where the tensor switches off (it needs S=1). So the
  deuteron binds on the tensor, while the dineutron, the same deep pocket without the
  tensor, just fails to bind: the 1S0 near-threshold virtual state, a_s = -23.7 fm
  (computed in nn_scattering_length.py). The hard core is common to both.

  APPROXIMATION, STATED PLAINLY.  This is a single-channel S-wave calculation. The
  couple-stress contact and the one-pion TENSOR force are S_12 operators; they act on
  3S1 only through 3S1-3D1 coupling and are absent from a single-channel central
  solve. The number below is therefore the central-force binding; the tensor pieces,
  which Sec. tensor_deuteron shows dominate the deuteron's D-state, are the next
  refinement (a coupled-channel solve). The couple-stress range is drawn on the
  figure for orientation but is NOT added to the central well, since adding a tensor
  force as if it were central would spuriously over-bind.

The overlap identity used for the dock:
  rho(r) = 1/(r^2+w^2)^2 has 3D transform rho_hat(k) = (pi^2/w) e^{-w k}. The
  cloud-cloud overlap is the inverse transform of |rho_hat|^2, a squared-Lorentzian
  of width 2w; the normalisation cancels in V_dock(r)/V_dock(ell).

Author: Mitch Cox.  Run: python nn_potential_derived.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from numba import njit
from scipy.linalg import eigh_tridiagonal

mpl.rcParams.update({"font.size": 10, "axes.linewidth": 0.8,
                     "figure.dpi": 140, "savefig.dpi": 140, "font.family": "serif"})

# ---- Inputs, all fixed elsewhere (CODATA / PDG / framework) -------------------
HBARC  = 197.3269804      # MeV fm
M_E    = 0.51099895       # MeV, electron mass = screw-screw Coulomb at one spacing
ALPHA  = 7.2973525693e-3  # fine-structure constant
M0     = M_E / ALPHA      # MeV ~ 70.025, node rest energy
ELL    = HBARC / M0       # fm ~ 2.8179, lattice constant = classical electron radius
R_CONF = ELL / (2*np.pi)  # fm ~ 0.4485, hard-core radius (junction confinement)
RP     = 0.8409           # fm, proton charge radius (PDG), the cloud width scale
D_SAT  = 1.80             # fm, saturation inter-nucleon spacing
M_N    = 938.918          # MeV, isoscalar nucleon mass (framework)

SQRT_SIGMA = 2*np.pi*M0           # MeV ~ 440, string tension^(1/2)
E_CONF     = 2*SQRT_SIGMA         # MeV ~ 880, core confinement energy = 4 pi m0
KAPPA_EX   = 1.0                  # the one order-unity EOS knob (high-density stiffness)
M_CORE     = KAPPA_EX*E_CONF      # hard-core strength

N_BOND  = 4                       # deuteron Mode-P charge-active face bonds
W_CLOUD = RP                      # squared-Lorentzian cloud width ~ charge radius

# One-pion tail (all lattice-derived) -----------------------------------------
M_PI   = 2*M0                     # MeV ~ 140, two-node stacking fault
F_PI   = 3**0.25 * M0            # MeV ~ 92.2, pion decay constant
M_D    = 4*M0/np.sqrt(np.pi-2)    # MeV ~ 262, Cosserat rotational gap (couple-stress)
XI     = HBARC/M_D                # fm ~ 0.753, couple-stress Yukawa range
G_C    = 1.0/(np.pi-1.0)          # ~ 0.4669, rotational-channel vertex coupling
Z1     = 12                       # FCC nearest-neighbour coordination
G_A    = 1.279                    # derived axial coupling (Sec. nucleon_axial)
F2_4PI = G_A**2/(4*np.pi*np.sqrt(3.0))   # = (m_pi g_A/2 f_pi)^2/4pi, pion-nucleon coupling
LAM_PI = HBARC/M_PI               # fm ~ 1.409, one-pion range
TAU_SIGMA = -3.0                  # (tau1.tau2)(sigma1.sigma2), same for 3S1 and 1S0
OPE_C  = F2_4PI*(M_PI/3.0)*TAU_SIGMA      # MeV, central OPE coefficient (attractive)
A_COUPLE = G_C**2 * Z1 * M0       # MeV ~ 183, couple-stress depth scale (TENSOR; not central)

print("=== Derived scales (nothing tuned) ===")
print(f"  eps_bond = alpha hbar c / ell = {ALPHA*HBARC/ELL:.5f} MeV   (= m_e = {M_E:.5f})")
print(f"  ell = {ELL:.4f} fm   R_conf = {R_CONF:.4f} fm   2 R_conf = {2*R_CONF:.4f} fm")
print(f"  M_core = kappa_ex * 4 pi m0 = {M_CORE:.1f} MeV   (kappa_ex = {KAPPA_EX})")
print(f"  one-pion: m_pi = {M_PI:.2f} MeV   f_pi = {F_PI:.2f} MeV   range = {LAM_PI:.3f} fm")
print(f"  f^2/4pi = g_A^2/(4 pi sqrt3) = {F2_4PI:.4f}  (g_A = {G_A})   OPE_C = {OPE_C:.3f} MeV")
print(f"  couple-stress: m_D = {M_D:.1f} MeV   xi = {XI:.3f} fm   g_c = {G_C:.4f}   (tensor)")


# ---- The pieces --------------------------------------------------------------
@njit(cache=True)
def core_overlap(r, a):
    """Normalised lens volume of two spheres radius a at separation r."""
    if r >= 2*a:
        return 0.0
    lens = (np.pi/12.0)*(2*a - r)**2*(4*a + r)
    return lens/((4.0/3.0)*np.pi*a**3)

@njit(cache=True)
def v_dock(r, w, ell, depth):
    """Central docking attraction, squared-Lorentzian normalised to -depth at ell."""
    num = (4*w*w + ell*ell)**2
    den = (4*w*w + r*r)**2
    return -depth*num/den

@njit(cache=True)
def v_ope_c(r, coeff, lam):
    """Central (spin-spin) one-pion-exchange Yukawa, coeff * e^{-x}/x, x = r/lam."""
    x = r/lam
    if x < 1e-9:
        x = 1e-9
    return coeff*np.exp(-x)/x

@njit(cache=True)
def v_couple(r, amp, xi):
    """Couple-stress Yukawa (TENSOR). Plotted for range only; not in central solve."""
    return -amp*(xi/r)*np.exp(-r/xi)

@njit(cache=True)
def v_central(r, w, ell, a, m_core, dock_depth, ope_coeff, lam, ope_on):
    V = m_core*core_overlap(r, a) + v_dock(r, w, ell, dock_depth)
    if ope_on > 0.5:
        V += v_ope_c(r, ope_coeff, lam)
    return V


# ---- Bound-state solver (radial matrix eigensolver, single-channel S-wave) ----
def binding_matrix(w, ell, a, m_core, dock_depth, ope_coeff, lam, ope_on,
                   mu, hbarc, rmax=40.0, h=0.02):
    """Lowest S-wave eigenvalue of -hbar^2/2mu u'' + V(r) u = E u on a radial grid,
    u(0)=u(rmax)=0. Returns E (MeV); E<0 is a bound state. Validated against square
    wells (35 MeV/2.1 fm -> -2.6 MeV; 100 MeV/1.0 fm -> unbound)."""
    N = int(rmax/h)
    rr = (np.arange(N) + 1)*h
    V = np.array([v_central(ri, w, ell, a, m_core, dock_depth, ope_coeff, lam, ope_on)
                  for ri in rr])
    t = hbarc*hbarc/(2.0*mu*h*h)
    diag = 2.0*t + V
    offd = -t*np.ones(N - 1)
    ev = eigh_tridiagonal(diag, offd, select='i', select_range=(0, 0))[0]
    return ev[0]   # lowest eigenvalue; E<0 is a true bound state

mu_red = M_N/2.0
depth_T0 = N_BOND*M_E      # full dock, -4 m_e at r = ell (common to BOTH channels)

# Lowest S-wave eigenvalue as the derived pieces switch on (E0<0 => bound)
E_dock = binding_matrix(W_CLOUD, ELL, R_CONF, M_CORE, depth_T0,   0.0, LAM_PI, 0.0, mu_red, HBARC)
E_full = binding_matrix(W_CLOUD, ELL, R_CONF, M_CORE, depth_T0, OPE_C, LAM_PI, 1.0, mu_red, HBARC)

def E_scaled(lam):
    """Lowest eigenvalue with the central attraction (dock + OPE) scaled by lam."""
    return binding_matrix(W_CLOUD, ELL, R_CONF, M_CORE, depth_T0*lam, OPE_C*lam,
                          LAM_PI, 1.0, mu_red, HBARC)

# threshold (E0 = 0) and full-binding (E0 = -2.225) scale factors, by bisection
def bisect_lam(target):
    lo, hi = 1.0, 8.0
    for _ in range(60):
        mid = 0.5*(lo + hi)
        if -E_scaled(mid) < target: lo = mid
        else: hi = mid
    return mid
lam_thresh = bisect_lam(0.0)
lam_deut   = bisect_lam(2.2246)

r = np.linspace(0.25, 4.0, 1600)
V_T0 = np.array([v_central(ri, W_CLOUD, ELL, R_CONF, M_CORE, depth_T0, OPE_C, LAM_PI, 1.0) for ri in r])
imin = np.argmin(V_T0)

print("\n=== Central T=0 well (docking + one-pion tail + hard core) ===")
print(f"  well minimum: V = {V_T0[imin]:.2f} MeV at r = {r[imin]:.3f} fm")
print(f"  V at 1.0 fm = {v_central(1.0,W_CLOUD,ELL,R_CONF,M_CORE,depth_T0,OPE_C,LAM_PI,1.0):.2f} MeV"
      f"   V at 1.8 fm = {v_central(D_SAT,W_CLOUD,ELL,R_CONF,M_CORE,depth_T0,OPE_C,LAM_PI,1.0):.2f} MeV")

def tag(E): return f"bound by {-E:.3f} MeV" if E < 0 else "unbound (lowest state at box floor)"
print("\n=== Single-channel S-wave result (central force only) ===")
print(f"  docking only              : {tag(E_dock)}")
print(f"  docking + one-pion tail    : {tag(E_full)}")
print(f"  (this central force is common to 3S1 and 1S0; the 1S0 scattering length is")
print(f"   computed in nn_scattering_length.py)")
print(f"\n  The derived central force sits right at the deuteron threshold: it crosses to")
print(f"  binding only when the attraction is scaled by lambda = {lam_thresh:.2f}, and reaches the")
print(f"  measured 2.225 MeV at lambda = {lam_deut:.2f}. So the lattice central force supplies")
print(f"  ~{100/lam_deut:.0f}% of the depth the deuteron needs. The remaining factor ~{lam_deut:.1f} is the")
print(f"  TENSOR force (one-pion + couple-stress, the S_12 operator) acting through 3S1-3D1")
print(f"  coupling, which Sec. tensor_deuteron shows dominates the deuteron's binding and")
print(f"  D-state. This is the textbook result that the deuteron is unbound without the")
print(f"  tensor force; computing it directly is the coupled-channel step that follows.")


# ---- Figure ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.6, 5.0))
ax.axhline(0, color="0.6", lw=0.8, zorder=1)
ax.axvspan(0.25, 2*R_CONF, color="0.90", zorder=0)
ax.text(0.33, 165, "cores overlap\n$r<2R_{\\rm conf}$", fontsize=8, color="0.30",
        ha="left", va="top")
# couple-stress range band (tensor; orientation only)
ax.axvspan(2*R_CONF, 2*R_CONF+XI, color="#fdf2e0", zorder=0)
ax.text(2*R_CONF+0.02, -55, "couple-stress\nrange $\\xi$ (tensor)", fontsize=7.4,
        color="#7f6000", ha="left", va="bottom")

ax.plot(r, V_T0, color="#1f4e79", lw=2.4, zorder=4,
        label=r"central force, common to $^3S_1$ and $^1S_0$")
ax.text(0.985, 0.80,
        "same pocket in both channels;\n"
        r"$^3S_1$ (deuteron) binds via the tensor," "\n"
        r"$^1S_0$ (dineutron) has no tensor:" "\n"
        r"near-threshold virtual state, $a_s=-23.7$ fm",
        transform=ax.transAxes, fontsize=8.0, color="#1f4e79", ha="right", va="top")

for x, lab in [(R_CONF, r"$R_{\rm conf}=\ell/2\pi$"),
               (D_SAT,  r"$d_{\rm sat}$"),
               (LAM_PI, r"$\hbar c/m_\pi$"),
               (ELL,    r"$\ell=r_e$")]:
    ax.axvline(x, color="0.55", lw=0.8, ls=":", zorder=2)
    ax.text(x, -58, lab, rotation=90, fontsize=7.4, color="0.40", ha="right", va="bottom")

ax.annotate(f"central well min\n$V={V_T0[imin]:.0f}$ MeV at {r[imin]:.2f} fm",
            xy=(r[imin], V_T0[imin]), xytext=(1.62, -36), fontsize=8.2, color="#1f4e79",
            ha="left", va="center",
            arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=1.0))
ax.annotate("one-pion tail widens the well\n$\\hbar c/m_\\pi=1.41$ fm, $f^2/4\\pi=g_A^2/4\\pi\\sqrt{3}$",
            xy=(2.2, v_central(2.2,W_CLOUD,ELL,R_CONF,M_CORE,depth_T0,OPE_C,LAM_PI,1.0)),
            xytext=(2.35, -30), fontsize=8.2, color="#1f4e79", ha="left", va="center",
            arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=1.0))
ax.annotate("Pauli hard core $M=\\kappa_{\\rm ex}\\,4\\pi m_0$",
            xy=(0.60, 105), xytext=(0.98, 120), fontsize=8.4, color="#7f6000",
            ha="left", va="center",
            arrowprops=dict(arrowstyle="->", color="#7f6000", lw=1.0))

ax.set_xlabel(r"nucleon$-$nucleon separation  $r$  [fm]")
ax.set_ylabel(r"central potential  $V(r)$  [MeV]")
ax.set_xlim(0.25, 4.0)
ax.set_ylim(-60, 190)
ax.legend(loc="upper right", frameon=False, fontsize=8.4, bbox_to_anchor=(1.0, 0.99))
ax.set_title("The two-nucleon central potential, derived from the lattice", fontsize=11)
fig.tight_layout()
fig.savefig("nn_potential_derived.pdf", bbox_inches="tight")
print("\nwrote nn_potential_derived.pdf")
