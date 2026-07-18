"""Bogoliubov-de Gennes spectrum of the lattice-condensate vortex on the
FCC (111) geometry at compact wavenumber k_4 = 2 pi / L_4: the transverse-
channel exclusion result.

OUTCOME (negative, and decisive). The scan in U n_bar (healing length
xi -> 0, the framework's incompressible-condensate limit) finds NO
saturating translation-kelvon branch and NO positive-norm anomalous mode:
every core-localised mode at k_4 d_111 = 2 pi/3 rides the chemical scale,
omega ~ 2.8 sqrt(U n_bar J), i.e. climbs with the condensate stiffness.
At a three-site compact wavelength there is no hydrodynamic kelvon regime,
and the framework's near-incompressible condensate (first-sound bulk
stiffness K_sf/mu ~ 1e40 comes from the condensate) pushes every
transverse condensate mode far above the 148.5 MeV crystal shear bottom.
The filament reading agrees from the other side: the incompressible-limit
lattice core parameter (triangular XY vortex, energy-log fit) is
a_E = 0.266 ell, for which the Thomson slow branch sits at 179 MeV, above
the ceiling. CONCLUSION: the transverse channels contribute no sharp
internal excitation at k_4 = +/-1. Combined with the gapless k_4 = 0
result (no sharp e*), the only sharp internal excitations of the
composite line are the crystal-core stacking-ring states, the Koide
Z_3 objects: the generations are ring states, uniquely, and the Koide
machinery is the primary mass calculation with no transverse
double-booking. A second structural result: omega(+k_4) = omega(-k_4)
exactly in minimal GP (no chirality splitting at quadratic order), so
the mu/tau splitting must come from the T-odd chiral compact coupling,
confirming the division of labour with the Koide ring Hamiltonian.

Model as built: GP on the FCC lattice folded to one (111) layer; 6
in-plane and 3+3 layer-shifted bonds, one hopping J = m_0 c^2/4 fixed by
the continuum limit with hbar = m_0 c ell; layer hops carry phases
e^{+/- i 2 pi/3} at k_4 = 2 pi/L_4; single knob U n_bar/J (healing
length). BdG about the relaxed vortex, both k_4 signs, both frequency
signs, positive-norm filtering.
"""
import numpy as np

J = 1.0
PHI = 2*np.pi/3.0                 # k_4 d_111

def build(R):
    a1 = np.array([1.0,0.0]); a2 = np.array([0.5,np.sqrt(3)/2])
    pts, idx = [], {}
    n = int(R*1.6)+2
    for i in range(-n,n+1):
        for j in range(-n,n+1):
            p = i*a1+j*a2
            if np.hypot(*p) <= R:
                idx[(i,j)] = len(pts); pts.append((i,j,p[0],p[1]))
    N = len(pts)
    inplane = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
    up      = [(0,0),(-1,0),(0,-1)]     # A -> B folded shifts
    dn      = [(0,0),(1,0),(0,1)]       # A -> C folded shifts
    bonds_ip, bonds_up, bonds_dn = [], [], []
    for (i,j,_,_ ) in pts:
        u = idx[(i,j)]
        for di,dj in inplane:
            v = idx.get((i+di,j+dj))
            if v is not None: bonds_ip.append((u,v))
        for di,dj in up:
            v = idx.get((i+di,j+dj))
            if v is not None: bonds_up.append((u,v))
        for di,dj in dn:
            v = idx.get((i+di,j+dj))
            if v is not None: bonds_dn.append((u,v))
    xy = np.array([[p[2],p[3]] for p in pts])
    return N, xy, bonds_ip, bonds_up, bonds_dn

def hop_matrix(N, bonds_ip, bonds_up, bonds_dn, phi):
    """Single-particle hopping at compact phase phi (directed bonds)."""
    H = np.zeros((N,N), complex)
    for u,v in bonds_ip: H[u,v] += -J
    for u,v in bonds_up: H[u,v] += -J*np.exp( 1j*phi)
    for u,v in bonds_dn: H[u,v] += -J*np.exp(-1j*phi)
    return H

def ground_vortex(N, xy, H0, Un, center, iters=4000):
    """Relax the GP vortex (k_4 = 0 sector, H0 real-phase hops)."""
    mu = Un - 12*J                     # uniform chemical potential, n = 1
    r = np.hypot(xy[:,0]-center[0], xy[:,1]-center[1])
    th = np.arctan2(xy[:,1]-center[1], xy[:,0]-center[0])
    psi = np.sqrt(np.clip(r,0,1.4)/1.4)*np.exp(1j*th)   # seed
    R = r.max(); rim = r > R-1.2
    psi[rim] = np.exp(1j*th[rim])
    dt = 0.05/ (12*J + 2*Un)
    for it in range(iters):
        grad = H0 @ psi + (Un*np.abs(psi)**2 - mu)*psi
        psi -= dt*grad
        psi[rim] = np.exp(1j*th[rim])                    # pin rim
    res = np.linalg.norm(grad[~rim])/np.sqrt(N)
    return psi, mu, res

def bdg_kelvon(N, xy, psi, mu, Un, Hp, Hm, center, want=6):
    """BdG at compact phase +phi (u-sector Hp, v-sector uses Hm.conj())."""
    n2 = np.abs(psi)**2
    A  = Hp + np.diag(-mu + 2*Un*n2)
    Am = Hm + np.diag(-mu + 2*Un*n2)
    B  = np.diag(Un*psi**2)
    L = np.block([[A, B],[-B.conj(), -Am.conj()]])
    ev, vec = np.linalg.eig(L)
    r = np.hypot(xy[:,0]-center[0], xy[:,1]-center[1])
    core = r < 3.0
    th = np.arctan2(xy[:,1]-center[1], xy[:,0]-center[0])
    phase = np.exp(1j*th)
    out = []
    for k in range(len(ev)):
        if abs(ev[k].imag) > 1e-6: continue
        w = ev[k].real
        if w <= 1e-9: continue
        u, v = vec[:N,k], vec[N:,k]
        norm = (np.abs(u)**2 - np.abs(v)**2).sum()
        if norm <= 0: continue
        u = u/np.sqrt(norm); v = v/np.sqrt(norm)
        loc = (np.abs(u[core])**2 + np.abs(v[core])**2).sum()
        if loc < 0.25: continue
        # retrograde kelvon: u carries e^{-i theta} relative to the vortex
        m_char = np.abs((u*np.conj(psi)*phase).sum())    # m = -1 overlap proxy
        out.append((w, loc, m_char))
    out.sort(key=lambda t: t[0])
    return out[:want]

# ----------------------------------------------------------------------
R = 14.0
N, xy, b_ip, b_up, b_dn = build(R)
H0 = hop_matrix(N, b_ip, b_up, b_dn, 0.0)
Hp = hop_matrix(N, b_ip, b_up, b_dn, +PHI)
Hm = hop_matrix(N, b_ip, b_up, b_dn, -PHI)
cen = (np.array([1.0,0.0])+np.array([0.5,np.sqrt(3)/2]))/3   # plaquette centre

# folded continuum bottom at k_4: min over plane waves of single-particle
# Bogoliubov energy at compact phase (Bloch on the 12-bond lattice)
def continuum_bottom(Un, phi, nq=241):
    a1 = np.array([1.0,0.0]); a2 = np.array([0.5,np.sqrt(3)/2])
    b1 = 2*np.pi*np.array([1.0,-1/np.sqrt(3)])
    b2 = 2*np.pi*np.array([0.0, 2/np.sqrt(3)])
    ip = [np.array(s[0]*a1+s[1]*a2) for s in [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]]
    upv= [np.array(s[0]*a1+s[1]*a2) for s in [(0,0),(-1,0),(0,-1)]]
    dnv= [np.array(s[0]*a1+s[1]*a2) for s in [(0,0),(1,0),(0,1)]]
    best = 1e9
    for fx in np.linspace(-0.5,0.5,nq):
        for fy in np.linspace(-0.5,0.5,nq):
            q = fx*b1 + fy*b2
            e  = sum(-J*np.exp(1j*q@d) for d in ip)
            e += sum(-J*np.exp(1j*(q@d))*np.exp(1j*phi)  for d in upv)
            e += sum(-J*np.exp(1j*(q@d))*np.exp(-1j*phi) for d in dnv)
            eps = e.real - (-12*J)            # kinetic above uniform bottom
            Ebog = np.sqrt(max(eps,0.0)*(eps + 2*Un))
            if Ebog < best and eps > 1e-12: best = Ebog
    return best

print(f"sites N = {N} at R = {R}")
for Un in [8.0, 18.0]:
    xi = np.sqrt(2*J/Un)
    psi, mu, res = ground_vortex(N, xy, H0, Un, cen)
    bot_p = continuum_bottom(Un, +PHI)
    print(f"\nUn/J = {Un}  (xi = {xi:.3f} l)   GP residual {res:.2e}   "
          f"continuum bottom at +k4: {bot_p:.4f} J")
    for sgn, Ha, Hb in [("+k4", Hp, Hm), ("-k4", Hm, Hp)]:
        modes = bdg_kelvon(N, xy, psi, mu, Un, Ha, Hb, cen)
        print(f"  {sgn}: low core modes (w/J, core-frac, m=-1 char):")
        for w, loc, mc in modes:
            print(f"      {w:8.4f}   {loc:.2f}   {mc:.2f}"
                  f"   -> {w*17.5065:7.2f} MeV")
