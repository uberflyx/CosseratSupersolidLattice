"""
Chiral two-pion exchange for the deuteron, parameter-free, and what it does (and does not) add.

Cosserat Supersolid Lattice framework. The coupled-channel solve (nn_deuteron_coupled.py)
binds the deuteron to 91% with two derived tensors, one-pion and couple-stress, leaving ~0.2
MeV. The natural question is whether two-pion exchange supplies it. This script computes the
leading (NLO) chiral two-pion-exchange potential parameter-free, from the same lattice-derived
inputs that fix the one-pion force (g_A, f_pi, m_pi), and tests it in the deuteron solve.

  INPUTS (framework-derived, nothing tuned):
    m_pi = 2 m0 = 140.05 MeV,  f_pi = 3^{1/4} m0 = 92.16 MeV,  g_A = 1.279.

  NLO TWO-PION AMPLITUDE (Holt, Kaiser, Weise, Phys. Rept. / arXiv:1304.6350, Eq. 22;
  momentum space, verified against the primary literature):
    V_2pi(q) = (1/384 pi^2 f_pi^4){4 m_pi^2(1+4g_A^2-5g_A^4) + q^2(1+10g_A^2-23g_A^4)
                                    - 48 g_A^4 m_pi^4/(4 m_pi^2+q^2)} L(q) tau1.tau2
             + (3 g_A^4/64 pi^2 f_pi^4) L(q) (sig1.sig2 q^2 - sig1.q sig2.q),
    L(q) = (sqrt(4 m_pi^2+q^2)/q) ln[(q+sqrt(4 m_pi^2+q^2))/2 m_pi].
  At this order there are only three structures: an ISOVECTOR CENTRAL (W_C), an ISOSCALAR
  SPIN-SPIN (V_S), and an ISOSCALAR TENSOR (V_T). Crucially there is NO isoscalar central,
  which is what protects against double-counting the framework's docking.

  COORDINATE SPACE (spectral / dispersive transform, derived here explicitly).
  Continuing q -> i mu gives the mass spectrum Im L(i mu) = (pi/2) sqrt(mu^2-4 m_pi^2)/mu for
  mu > 2 m_pi, and any central amplitude with spectral function eta(mu) transforms as
    V(r) = 1/(2 pi^2 r) int_{2 m_pi}^inf mu eta(mu) e^{-mu r} dmu.
  The Fourier transform of L(q) itself is the closed form L~(r) = m_pi/(2 pi r^2) K_1(2 m_pi r).
  Differentiating it for the spin-spin (-2/3 A nabla^2) and tensor (1/3 A (d2-d1/r)) operators,
  with A = 3 g_A^4/(64 pi^2 f_pi^4), gives the local potentials as single spectral integrals:

    V_T(r) = g_A^4/(256 pi^3 f_pi^4 r) int sqrt(mu^2-4m_pi^2) mu^2 [1+3/(mu r)+3/(mu r)^2] e^{-mu r} dmu
    V_S(r) = -g_A^4/(128 pi^3 f_pi^4 r) int sqrt(mu^2-4m_pi^2) mu^2 e^{-mu r} dmu
    W_C(r) = 1/(1536 pi^3 f_pi^4 r) int P(-mu^2) sqrt(mu^2-4m_pi^2) e^{-mu r} dmu,
      P(-mu^2) = 4 m_pi^2(1+4g_A^2-5g_A^4) - mu^2(1+10g_A^2-23g_A^4) - 48 g_A^4 m_pi^4/(4 m_pi^2-mu^2).
  (The 1/(4 m_pi^2-mu^2) pole is integrable: it cancels one power of the sqrt at threshold.)

  THE DEUTERON (T=0, S=1) sees tau1.tau2 = -3 and sig1.sig2 = +1, so the two-pion
  contribution is central = -3 W_C + V_S and tensor = V_T.

  WHAT IT SHOWS.  The two-pion CENTRAL is strongly attractive (-72 MeV at 1 fm), four times
  the docking's -16 MeV: it is the same intermediate-range attraction the framework already
  carries as the docking (node-sharing), here in bare, unrenormalised form. Adding it double-
  counts. The two-pion TENSOR (+14.5 MeV at 1 fm) has the OPPOSITE sign to the one-pion tensor
  (-82 MeV), so it fights the binding: added alone it drops the deuteron from 91% to 39%. The
  parameter-free long-range two-pion therefore does NOT supply the missing 0.2 MeV by addition.
  Its central piece is the docking; its tensor opposes the pion. The deuteron's last tenth of
  an MeV is a short-range balance, the threshold position of the central attraction against the
  hard core, set in the framework by the static bond count. That balance is the open piece.

Author: Mitch Cox.  Run: python nn_tpe_chiral.py
"""

import numpy as np
from scipy import integrate
from scipy.sparse import diags, bmat
from scipy.sparse.linalg import eigsh

# ---- derived inputs (CODATA / framework), nothing tuned ----------------------
HBARC = 197.3269804
M0    = 0.51099895/7.2973525693e-3      # m_e/alpha = 70.025 MeV
ELL   = HBARC/M0
R_CONF= ELL/(2*np.pi)
RP    = 0.8409
M_N   = 938.918
HB2M  = HBARC*HBARC/(2.0*(M_N/2.0))
M_CORE= 4*np.pi*M0
NB    = 4
W     = RP
MPI   = 2*M0
FPI   = 3**0.25*M0
GA    = 1.279
F2    = GA**2/(4*np.pi*np.sqrt(3.0))
LAM   = HBARC/MPI
OPE_C = F2*(MPI/3.0)*(-3.0)
OPE_T = F2*(MPI/3.0)*(-3.0)
RCUT  = 2*R_CONF
MD    = 4*M0/np.sqrt(np.pi-2.0)
XI    = HBARC/MD
GC    = 1.0/(np.pi-1.0)
COUPLE_T = -(GC**2/(4*np.pi))*(MD/3.0)
M2    = 4*MPI**2

# ---- NLO two-pion coordinate-space potentials (single spectral integrals) -----
def _root(u):
    return np.sqrt(u*u - M2)

def V_T_tpe(rfm):
    """isoscalar tensor (positive: opposes the one-pion tensor in the deuteron)."""
    if rfm < 0.5:                      # cut off inside the cores, like the one-pion tensor
        return 0.0
    r = rfm/HBARC
    pref = GA**4/(256*np.pi**3*FPI**4*r)
    f = lambda u: _root(u)*u*u*(1.0+3.0/(u*r)+3.0/(u*r)**2)*np.exp(-u*r)
    return pref*integrate.quad(f, 2*MPI, 100*MPI, limit=300)[0]

def V_S_tpe(rfm):
    """isoscalar spin-spin (attractive in the S=1 deuteron)."""
    if rfm < 0.5:
        return 0.0
    r = rfm/HBARC
    pref = -GA**4/(128*np.pi**3*FPI**4*r)
    f = lambda u: _root(u)*u*u*np.exp(-u*r)
    return pref*integrate.quad(f, 2*MPI, 100*MPI, limit=300)[0]

def W_C_tpe(rfm):
    """isovector central (attractive in the deuteron through tau1.tau2 = -3)."""
    if rfm < 0.5:
        return 0.0
    r = rfm/HBARC
    a = 1+4*GA**2-5*GA**4
    b = 1+10*GA**2-23*GA**4
    pref = 1.0/(1536*np.pi**3*FPI**4*r)
    def f(u):
        return (4*MPI**2*a - u*u*b - 48*GA**4*MPI**4/(M2-u*u))*_root(u)*np.exp(-u*r)
    return pref*integrate.quad(f, 2*MPI, 100*MPI, limit=400)[0]

def ope_tensor(rfm):
    x = rfm/LAM
    return OPE_T*(1.0+3.0/x+3.0/x**2)*np.exp(-x)/x

def v_dock(rfm):
    return -NB*0.51099895*((4*W*W+ELL*ELL)**2)/((4*W*W+rfm*rfm)**2)

# ---- table: the two-pion force against the one-pion tensor and the docking ----
print("=== NLO two-pion exchange, parameter-free (g_A, f_pi, m_pi all derived) ===")
print(f"  m_pi = {MPI:.2f} MeV   f_pi = {FPI:.2f} MeV   g_A = {GA}")
print("\n  r[fm]   V_T(2pi)  V_S(2pi)  W_C(2pi) | deuteron 2pi: central(-3W_C+V_S)  tensor(V_T)"
      " | OPE_T   docking")
for rfm in (0.9, 1.0, 1.2, 1.5, 2.0):
    vt, vs, wc = V_T_tpe(rfm), V_S_tpe(rfm), W_C_tpe(rfm)
    print(f"  {rfm:4.2f}  {vt:8.2f}  {vs:8.2f}  {wc:8.2f} |     central ={-3*wc+vs:8.2f}"
          f"   tensor ={vt:7.2f} | {ope_tensor(rfm):7.1f} {v_dock(rfm):7.1f}")
print("  (MeV. Two-pion central is the docking in bare form; two-pion tensor opposes OPE.)")

# ---- deuteron solve: the net effect of dropping the long-range two-pion in -----
def omega(r, a):
    return np.where(r < 2*a, (2*a-r)**2*(4*a+r)/(16*a**3), 0.0)
def yrad(r):
    x = r/LAM; return np.exp(-x)/x
def trad(r, lam):
    x = r/lam; return (1.0+3.0/x+3.0/x**2)*np.exp(-x)/x
def reg(r, wr):
    return 1.0/(1.0+np.exp(-(r-RCUT)/wr))

def deuteron(tpe_T=False, tpe_S=False, wr=0.12, rmax=25.0, h=0.02):
    N = int(rmax/h); r = (np.arange(N)+1)*h; f = reg(r, wr)
    VTt = np.array([V_T_tpe(ri) for ri in r])
    VSs = np.array([V_S_tpe(ri) for ri in r])
    V_C = M_CORE*omega(r, R_CONF) + v_dock(r) + f*OPE_C*yrad(r)
    if tpe_S:
        V_C = V_C + f*VSs
    V_T = f*OPE_T*trad(r, LAM) + f*COUPLE_T*trad(r, XI)
    if tpe_T:
        V_T = V_T + f*VTt
    t = HB2M/(h*h); cent = HB2M*6.0/(r*r); off = -t*np.ones(N-1)
    A = diags([off, 2*t+V_C, off], [-1, 0, 1])
    D = diags([off, 2*t+cent+V_C-2.0*V_T, off], [-1, 0, 1])
    B = diags([np.sqrt(8.0)*V_T], [0])
    H = bmat([[A, B], [B, D]], format='csr')
    vals, vecs = eigsh(H, k=2, which='SA')
    i0 = int(np.argmin(vals)); E0 = vals[i0]; v = vecs[:, i0]
    u, w = v[:N], v[N:]; nrm = np.sum(u*u+w*w)
    return E0, 100*np.sum(w*w)/nrm

print("\n=== Deuteron with the parameter-free long-range two-pion dropped in ===")
for tag, tt, ts in (("OPE + couple-stress (current)        ", False, False),
                    ("  + two-pion tensor                  ", True,  False),
                    ("  + two-pion tensor + spin-spin       ", True,  True)):
    E, PD = deuteron(tpe_T=tt, tpe_S=ts)
    print(f"  {tag}: E = {E:7.3f} MeV ({100*(-E)/2.2246:3.0f}%)   P_D = {PD:5.2f}%")
print("  measured                              : E = -2.225 MeV (100%)")
print("\n  The two-pion central is the docking (already in the solve); its tensor opposes the")
print("  one-pion tensor and reduces the binding. Parameter-free two-pion does not add the")
print("  missing 0.2 MeV: the residual is the short-range threshold balance, not a missing tail.")
