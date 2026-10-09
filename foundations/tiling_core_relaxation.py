"""Internal relaxation of four-registry pairs and the width of a screw core.

Companion to the monograph chapters "Internal order and sublattices" and "The fine
structure constant". The pair model is the conditional one of
four_registry_pair_order.py: each node carries a rigid pair of radius a, neighbours
interact through four equal central Morse contacts, and the pair rotor has inertia
multiplier eta. Nothing here fixes the vacuum's constituent forces.

The script asks how much the pairs' internal response softens (i) a long-wavelength
anti-plane shear and (ii) the slip across one {111} glide wall, and what that does to
the Cosserat Peierls-Nabarro screw core. It computes, in order:

1. The self-consistent pair state (seeded with the tetrahedral texture) and the joint
   second-derivative blocks: the angular Hessian, the centre Hessian of each bond and
   their mixed derivative, at compact momentum zero. The centre lattice may be
   dilated to remove the virial stress of the fixed box.
2. For a centre displacement u e^{ik.X}, equal on all four registries, the energy at
   three levels of internal response: angular states fixed in laboratory axes, each
   node's state free to turn rigidly ("free rotation"), and fully relaxed. The relaxed
   energy is the Schur complement D - G^dag H^{-1} G of the Bloch blocks.
3. The relaxation of a slip across the glide wall, as a zone integral over waves
   along the wall normal with weight 1/(4 sin^2(k_n h_p/2)), where h_p = ell/sqrt6 is
   the projected plane spacing. Two profiles: a sharp step between adjacent planes,
   and the D4 wall whose interleaved half-layer takes the share of the slip that
   minimises its energy (always one half).
4. The equilibrium core width w/d of the screw (Cosserat anti-plane kernel at
   N^2 = 1/pi, gamma = mu ell^2, arctan profiles), with the kernel multiplied by the
   relaxed-to-reference ratio at every wavevector inside |K| <= kmax. The misfit is
   multiplied either by the wall's ratio (discrete reading) or by the long-wave ratio
   along the normal (the continuum convention the bare equilibrium uses for the
   microrotation).

Checks (always): the long-wave moduli reproduce the frozen and relaxed screw moduli
of the chapter (a = 0.3 ell, eta = 600); the bare wall-to-step ratio tends to the
exact 5/8 of the interleaved half-layer as a -> 0; the unrelaxed core gives
w/d = 0.83664. Checks (--checks, at the reference point): the split of the long-wave
relaxation between in-slice and compact bonds; the shift when the registries'
relative displacements also relax; the dependence on kmax; a direct evaluation of the
kernel against the interpolated map.

Units: ell = hbar = m0 = c = 1. Runtime: about a minute per parameter point; --scan
adds about forty minutes and --checks about ten on one core.
"""
import argparse
import json
import os
import sys

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import RegularGridInterpolator
from scipy.linalg import eigh
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from four_registry_pair_order import (BETA, DIRECTORS, KN, angular_grid, charge,  # noqa: E402
                                      neighbour_registry, real_harmonics, roots)

H_PLANE = 1/np.sqrt(6)                     # projected {111} plane spacing h_p at k4 = 0
LINE = np.array([1, -1, 0, 0])/np.sqrt(2)  # screw line and slip direction
NORMAL = np.array([1, 1, 1, 0])/np.sqrt(3)
GLIDE = np.r_[np.cross(LINE[:3], NORMAL[:3]), 0.]   # in-plane, normal to the line


# ------------------------------------------------------------ ordered state and blocks

def joint_blocks(a=.3, eta=600., sx=1., s4=1., n=16, L=8, tol=1e-11):
    """Hartree state and bond blocks with spatial bond components scaled by sx, compact by s4."""
    pts, w = angular_grid(n)
    Y, T = real_harmonics(pts, L)
    rs = roots()
    B = 1/(2*a*a*eta)
    R = rs/np.sqrt(2)*np.array([sx, sx, sx, s4])
    target = np.array([neighbour_registry(charge(r)) for r in rs]).T     # [site, bond]
    delta = np.zeros((len(pts), len(pts), 4))
    delta[:, :, :3] = a*(pts[None, :, :] - pts[:, None, :])
    W, F, C = [], [], []
    for Rb in R:
        x = Rb + delta
        r = np.linalg.norm(x, axis=2)
        nh = x/r[:, :, None]
        e = np.exp(-BETA*(r - 1))
        d1, d2 = (e - e*e)/BETA, 2*e*e - e            # W' and W''
        W.append((e*e - 2*e)/(2*BETA**2))
        F.append(d1[:, :, None]*nh)
        C.append((d2 - d1/r)[:, :, None, None]*nh[:, :, :, None]*nh[:, :, None, :]
                 + (d1/r)[:, :, None, None]*np.eye(4))
    W, F, C = np.array(W), np.array(F), np.array(C)

    def fields(g):
        p = w[None, :]*(g @ Y.T)**2
        v = np.zeros_like(p)
        for b in range(24):
            v += p[target[:, b]] @ W[b].T
        return p, v, float(np.mean((g*g) @ (B*T)) + KN*np.sum(p*v)/8)

    g = np.zeros((4, len(T)))
    for i, d in enumerate(DIRECTORS):                  # tetrahedral seed
        V = -20*(3*(pts @ d)**2 - 1)/2
        g[i] = eigh(np.diag(T) + Y.T @ ((w*V)[:, None]*Y), subset_by_index=[0, 0])[1][:, 0]
        g[i] *= np.sign(g[i, 0])
    for _ in range(3000):
        p, v, energy = fields(g)
        new = np.array([eigh(np.diag(B*T) + KN*Y.T @ ((w*v[i])[:, None]*Y),
                             subset_by_index=[0, 0])[1][:, 0] for i in range(4)])
        new *= np.sign(np.sum(new*g, axis=1))[:, None]
        if np.max(abs(new - g)) < tol:
            break
        step = .5
        while True:
            trial = (1 - step)*g + step*new
            trial /= np.linalg.norm(trial, axis=1)[:, None]
            if fields(trial)[2] <= energy + 2e-14 or step < 1e-6:
                break
            step *= .5
        g = trial
    p, v, energy = fields(g)
    ng = len(T) - 1
    h0 = np.zeros((4*ng, 4*ng))
    tangent = []
    for i in range(4):
        ev, vec = eigh(np.diag(B*T) + KN*Y.T @ ((w*v[i])[:, None]*Y))
        h0[i*ng:(i + 1)*ng, i*ng:(i + 1)*ng] = np.diag(2*(ev[1:] - ev[0]))
        tangent.append(2*(w*(Y @ g[i]))[:, None]*(Y @ vec[:, 1:]))
    Hb = np.empty((4, 24, ng, ng))
    Gb = np.empty((4, 24, ng, 4))
    Cb = np.empty((4, 24, 4, 4))
    force = np.empty((4, 24, 4))
    for i in range(4):
        for b in range(24):
            j = target[i, b]
            Hb[i, b] = KN*tangent[i].T @ W[b] @ tangent[j]
            Gb[i, b] = KN*np.einsum('me,mnd,n->ed', tangent[i], F[b], p[j], optimize=True)
            Cb[i, b] = KN*np.einsum('m,mndf,n->df', p[i], C[b], p[j], optimize=True)
            force[i, b] = KN*np.einsum('m,mnd,n->d', p[i], F[b], p[j], optimize=True)
    Q = np.einsum('an,ni,nj->aij', p, pts, pts) - np.eye(3)/3
    return dict(a=a, eta=eta, sx=sx, s4=s4, R=R, target=target, h0=h0, Hb=Hb, Gb=Gb, Cb=Cb,
                ng=ng, energy=energy,
                virial=np.einsum('ibc,bd->cd', force, R)/8,    # dE/d(strain) per node
                alignment=float(np.sqrt(np.mean(np.sum(Q*Q, axis=(1, 2))))))


def stress_free(a, eta):
    """Dilations (sx, s4) at which the spatial and compact virial stresses vanish."""
    s4 = 1.
    for _ in range(8):
        sx = brentq(lambda x: np.trace(joint_blocks(a, eta, x, s4)['virial'][:3, :3]),
                    1., 1. + 3*a*a, xtol=1e-6)
        s4_new = brentq(lambda y: joint_blocks(a, eta, sx, y)['virial'][3, 3],
                        .9, 1.1 + a*a, xtol=1e-6)
        if abs(s4_new - s4) < 1e-5:
            return sx, s4_new
        s4 = s4_new
    return sx, s4


# ------------------------------------------------------------ Bloch energies

class Bloch:
    """Bare, free-rotation and relaxed energies of u e^{ik.X} (u equal on all registries)."""

    def __init__(self, m):
        ng = m['ng']
        N = 4*ng
        self.m, self.R = m, m['R']
        self.MH = np.zeros((24, N, N))
        self.MD = np.zeros((24, 16, 16))
        self.MG = np.zeros((24, N, 16))
        self.D0 = np.zeros((16, 16))
        self.G0 = np.zeros((N, 16))
        for i in range(4):
            ii, ui = slice(i*ng, (i + 1)*ng), slice(4*i, 4*i + 4)
            for b in range(24):
                j = m['target'][i, b]
                jj, uj = slice(j*ng, (j + 1)*ng), slice(4*j, 4*j + 4)
                self.MH[b, ii, jj] += m['Hb'][i, b]
                self.D0[ui, ui] += m['Cb'][i, b]
                self.MD[b, ui, uj] -= m['Cb'][i, b]
                self.G0[ii, ui] -= m['Gb'][i, b]
                self.MG[b, ii, uj] += m['Gb'][i, b]
        # rigid rotations of each registry's state: the k = 0 response to a uniform rotation
        H0 = m['h0'] + self.MH.sum(0)
        S = np.einsum('ibea,bc->ieac', m['Gb'], self.R).reshape(N, 16)
        V = np.zeros((N, 12))
        for k, (p, q) in enumerate([(1, 2), (2, 0), (0, 1)]):
            Om = np.zeros((4, 4))
            Om[p, q], Om[q, p] = -1, 1
            gen = -np.linalg.solve(H0, S @ Om.ravel())
            for i in range(4):
                V[i*ng:(i + 1)*ng, 3*i + k] = gen[i*ng:(i + 1)*ng]
        # a state with no orientation at all has nothing to rotate: the free-rotation level
        # then coincides with the laboratory one
        self.V = V if np.linalg.norm(V) > 1e-8*np.sqrt(N) else None

    def blocks(self, k):
        ph = np.exp(1j*(self.R @ k))
        return (self.m['h0'] + np.tensordot(ph, self.MH, 1),
                self.D0 + np.tensordot(ph, self.MD, 1),
                self.G0 + np.tensordot(ph, self.MG, 1))

    def angular_min(self, k):
        """Lowest eigenvalue of the angular Hessian H(k); positive for a locally stable state."""
        return float(np.linalg.eigvalsh(self.blocks(k)[0])[0])

    def relaxed_with_registries(self, k, pol=LINE):
        """Bare energy and the relaxed energy with the registries' relative (and compact)
        centre displacements also free, at fixed registry-averaged spatial displacement."""
        HH, D, G = self.blocks(k)
        S = D - G.conj().T @ np.linalg.solve(HH, G)
        avg = np.zeros((3, 16))
        for i in range(4):
            avg[:, 4*i:4*i + 3] = np.eye(3)/4
        Z = np.linalg.svd(avg)[2][3:].T                  # displacements with zero average
        u = np.tile(pol, 4).astype(complex)
        b = Z.T @ S @ u
        return np.array([(u.conj() @ D @ u).real,
                         (u.conj() @ S @ u - b.conj() @ np.linalg.solve(Z.T @ S @ Z, b)).real])

    def energies(self, k, pol=LINE):
        HH, D, G = self.blocks(k)
        u = np.tile(pol, 4).astype(complex)
        Gu = G @ u
        bare = (u.conj() @ D @ u).real
        relaxed = bare - (Gu.conj() @ np.linalg.solve(HH, Gu)).real
        if self.V is None:
            free = bare
        else:
            gr = self.V.T @ Gu
            free = bare - (gr.conj() @ np.linalg.solve(self.V.T @ HH @ self.V, gr)).real
        return np.array([bare, free, relaxed])


def long_wave_moduli(bl, direction, eps=3e-4):
    """Bare, free-rotation, relaxed shear moduli (m0 c^2/ell^3) of a long anti-plane wave."""
    return bl.energies(eps*direction)*2*np.sqrt(2)/(4*eps**2)


def wall_integrals(bl, sx=1.):
    """Zone integrals A_j = <E_j w>, B_j = <E_j w cos(kappa h)> for the three levels."""
    def E(kap):
        return bl.energies(kap*NORMAL/sx)
    wgt = lambda kap: 1/(4*np.sin(kap*H_PLANE/2)**2)
    top = np.pi/H_PLANE
    A = np.array([quad(lambda q, j=j: E(q)[j]*wgt(q), 1e-9, top, limit=400, epsrel=1e-10)[0]
                  for j in range(3)])
    Bc = np.array([quad(lambda q, j=j: E(q)[j]*wgt(q)*np.cos(q*H_PLANE), 1e-9, top,
                        limit=400, epsrel=1e-10)[0] for j in range(3)])
    return A, Bc


def wall_energies(A, Bc):
    """Sharp step, and the wall with half-layer slip x minimising each level's energy.

    The profile x Theta(m >= 1) + (1 - x) Theta(m >= 2) has
    |u(kappa)|^2 = [x^2 + (1-x)^2 + 2x(1-x) cos(kappa h)] / (4 sin^2(kappa h/2)).
    The energy is (2x^2 - 2x + 1) A + 2x(1-x) B, stationary at x = 1/2 whenever A != B,
    so the half-layer sits at the midpoint at every level of internal response.
    """
    x = .5
    return A, (x*x + (1 - x)**2)*A + 2*x*(1 - x)*Bc, x


# ------------------------------------------------------------ Peierls-Nabarro core

N2 = 1/np.pi
MU, GAM = 1., 1.
KAP = 2*N2*MU/(1 - 2*N2)
MUB, MUT = MU + KAP/2, MU + KAP
D_HOP, D111 = 1/np.sqrt(3), np.sqrt(2/3)
_XG, _WG = np.polynomial.legendre.leggauss(600)


def _kraw(K2):
    """Cosserat anti-plane stiffness after eliminating the microrotation (times K^2)."""
    return K2*(GAM*MUT*K2 + 2*KAP*MUB)/(GAM*K2 + 2*KAP)


def line_kernel(k, keep, cut):
    """K(k) = k^2/(2 pi) int dky Kraw(K) keep(k, ky)/K^4; keep = 1 beyond |ky| = cut."""
    th = np.arctan(cut/k)*_XG
    ky = k*np.tan(th)
    jac = k/np.cos(th)**2*np.arctan(cut/k)
    K2 = k*k + ky*ky
    inner = np.sum(_WG*jac*_kraw(K2)*keep(np.full_like(ky, k), ky)/K2**2)
    u = (_XG + 1)/(2*cut)
    K2o = k*k + 1/u**2
    outer = 2*np.sum(_WG/(2*cut)*_kraw(K2o)/K2o**2/u**2)
    return k*k*(inner + outer)/(2*np.pi)


def core_width(keep=lambda a, b: 1., misfit=1., cut=16., n=900):
    """w/d from int K(k)/k e^{-2kw} dk = misfit * mubar/(2 d111) (arctan profiles)."""
    ks = np.geomspace(1e-5, 80, n)
    K = np.array([line_kernel(k, keep, cut) for k in ks])
    lk = np.log(ks)
    f = lambda w: np.trapezoid(K*np.exp(-2*ks*w), lk) + MUB*ks[0]/2 - misfit*MUB/(2*D111)
    return brentq(f, .2, 1.5, xtol=1e-9)/D_HOP


def relaxation_map(bl, sx=1., kmax=16., dk=.1, nth=72, energy=None):
    """Ratios relaxed/reference on a polar grid in the (glide, normal) plane.

    Returns keep-functions (1 - f) for the laboratory and free-rotation references and
    the lowest angular-Hessian eigenvalue found on a subgrid of the map."""
    energy = energy or bl.energies
    Kg = np.r_[0., np.geomspace(1e-3, .05, 6), np.arange(.1, kmax - 1e-9, dk), kmax]
    th = np.linspace(0, np.pi, nth + 1)
    dirs = [np.cos(t)*GLIDE + np.sin(t)*NORMAL for t in th]
    E = np.array([[energy(max(K, 1e-4)*d/sx) for d in dirs] for K in Kg])
    hmin = min(bl.angular_min(K*d/sx) for K in Kg[::8] for d in dirs[::6])

    def keeper(ratio):
        it = RegularGridInterpolator((Kg, th), ratio)

        def keep(a, b):
            K = np.hypot(a, b)
            t = np.mod(np.arctan2(b, a), np.pi)
            out = np.ones_like(K)
            ok = K <= kmax
            out[ok] = it(np.column_stack([K[ok], t[ok]]))
            return out
        return keep
    if E.shape[-1] == 2:                                 # bare and relaxed only
        return keeper(E[..., 1]/E[..., 0]), None, hmin
    return keeper(E[..., 2]/E[..., 0]), keeper(E[..., 2]/E[..., 1]), hmin


# ------------------------------------------------------------ one parameter point

def analyse(a, eta, sx=1., s4=1., kmax=16.):
    m = joint_blocks(a, eta, sx, s4)
    bl = Bloch(m)
    out = {'a': a, 'eta': eta, 'sx': sx, 's4': s4, 'alignment': m['alignment'],
           'virial_diag': np.diag(m['virial']).tolist()}
    for name, d in (('glide', GLIDE), ('normal', NORMAL)):
        e = long_wave_moduli(bl, d/sx)
        out['long_' + name] = {'moduli': e.tolist(), 'f_lab': 1 - e[2]/e[0], 'f_free': 1 - e[2]/e[1]}
    A, Bc = wall_integrals(bl, sx)
    step, wall, x = wall_energies(A, Bc)
    out['step'] = {'f_lab': 1 - step[2]/step[0], 'f_free': 1 - step[2]/step[1]}
    out['wall'] = {'f_lab': 1 - wall[2]/wall[0], 'f_free': 1 - wall[2]/wall[1],
                   'half_layer_slip': x, 'bare_wall_over_step': wall[0]/step[0]}
    lab, free, hmin = relaxation_map(bl, sx, kmax)
    out['angular_hessian_min'] = {'zone_centre': bl.angular_min(np.zeros(4)), 'map': hmin}
    out['w_over_d'] = {
        'lab_wall': core_width(lab, 1 - out['wall']['f_lab'], kmax),
        'free_wall': core_width(free, 1 - out['wall']['f_free'], kmax),
        'lab_step': core_width(lab, 1 - out['step']['f_lab'], kmax),
        'free_step': core_width(free, 1 - out['step']['f_free'], kmax),
        'lab_kernel_only': core_width(lab, 1., kmax),
        'lab_wall_misfit_only': core_width(misfit=1 - out['wall']['f_lab'], cut=kmax),
        # continuum convention: the misfit takes the full long-wave fraction along the normal
        'lab_longwave_misfit': core_width(lab, 1 - out['long_normal']['f_lab'], kmax),
        'free_longwave_misfit': core_width(free, 1 - out['long_normal']['f_free'], kmax)}
    return out


def profile_along_normal(a=.3, eta=600., npts=121):
    """Relaxed fractions of waves along the wall normal across the zone (for plotting)."""
    bl = Bloch(joint_blocks(a, eta))
    kap = np.linspace(1e-3, np.pi/H_PLANE, npts)
    E = np.array([bl.energies(q*NORMAL) for q in kap])
    return {'kn_hp_over_pi': (kap*H_PLANE/np.pi).tolist(),
            'f_lab': (1 - E[:, 2]/E[:, 0]).tolist(), 'f_free': (1 - E[:, 2]/E[:, 1]).tolist()}


def checks():
    out = {}
    bl = Bloch(joint_blocks(.3, 600.))
    p, q = np.array([1, 0, 0, 0.]), np.array([0, 1, -1, 0])/np.sqrt(2)
    t = np.array([0, 1, 1, 0])/np.sqrt(2)
    out['screw_moduli_p'] = (bl.energies(3e-4*p, t)*2*np.sqrt(2)/(4*9e-8)).tolist()
    out['screw_moduli_q'] = (bl.energies(3e-4*q, t)*2*np.sqrt(2)/(4*9e-8)).tolist()
    out['screw_moduli_chapter'] = {'p': [2.245284, 1.576013], 'q': [1.294132, 1.077250]}
    A, Bc = wall_integrals(Bloch(joint_blocks(.01, 600.)))
    step, wall, _ = wall_energies(A, Bc)
    out['bare_wall_over_step_small_pair'] = wall[0]/step[0]
    out['unrelaxed_core'] = core_width()
    return out


def extra_checks(a=.3, eta=600.):
    """Slow checks at the reference point (see the module docstring)."""
    out = {}
    m = joint_blocks(a, eta)
    bl = Bloch(m)
    kzone = np.pi/H_PLANE
    # which bonds carry the coupling: zero the mixed derivatives of one bond family
    compact = np.array([r[3] != 0 for r in roots()])
    for name, mask in (('in_slice_only', ~compact), ('compact_only', compact)):
        mm = dict(m)
        Gb = m['Gb'].copy()
        Gb[:, ~mask] = 0
        mm['Gb'] = Gb
        b2 = Bloch(mm)
        out['split_' + name] = [float(1 - e[2]/e[0]) for e in
                                (b2.energies(1e-4*NORMAL), b2.energies(.5*kzone*NORMAL))]
    # registries' relative displacements also relaxed (at the relaxed level only)
    wgt = lambda q: 1/(4*np.sin(q*H_PLANE/2)**2)
    Ew = lambda q: bl.relaxed_with_registries(q*NORMAL)
    A = np.array([quad(lambda q, j=j: Ew(q)[j]*wgt(q), 1e-9, kzone, limit=400, epsrel=1e-9)[0]
                  for j in range(2)])
    Bc = np.array([quad(lambda q, j=j: Ew(q)[j]*wgt(q)*np.cos(q*H_PLANE), 1e-9, kzone,
                        limit=400, epsrel=1e-9)[0] for j in range(2)])
    f_wall = 1 - (A[1] + Bc[1])/(A[0] + Bc[0])
    keep, _, _ = relaxation_map(bl, energy=bl.relaxed_with_registries)
    lab, _, _ = relaxation_map(bl)
    Aw, Bw = wall_integrals(bl)
    f0 = 1 - (Aw[2] + Bw[2])/(Aw[0] + Bw[0])
    out['registry_relative'] = {'f_wall': f_wall, 'w_lab_wall': core_width(keep, 1 - f_wall),
                                'w_lab_wall_reference': core_width(lab, 1 - f0)}
    # dependence on where the relaxation is cut off
    cuts = {}
    for km in (kzone, 16., 24., 32.):
        kp, _, _ = relaxation_map(bl, kmax=km)
        cuts['%.4f' % km] = core_width(kp, 1 - f0, km)
    out['cut_dependence_lab_wall'] = cuts
    # interpolated map against direct evaluation of the kernel
    direct = lambda x, y: np.array([1 - (lambda e: e[2]/e[0])(bl.energies(u*GLIDE + v*NORMAL))
                                    if np.hypot(u, v) <= 16. else 0. for u, v in zip(x, y)])
    out['kernel_direct_vs_map'] = {str(k): [line_kernel(k, lab, 16.),
                                            line_kernel(k, lambda x, y: 1 - direct(x, y), 16.)]
                                   for k in (.3, 1., 2., 4.)}
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--scan', action='store_true', help='full parameter scan (slow)')
    ap.add_argument('--checks', action='store_true', help='slow checks at the reference point')
    ap.add_argument('--output', default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                     'tiling_core_relaxation.json'))
    args = ap.parse_args()
    res = {'checks': checks()}
    print(json.dumps(res['checks'], indent=1), flush=True)
    res['profile'] = profile_along_normal()
    res['points'] = [analyse(.3, 600.)]
    print(json.dumps(res['points'][0], indent=1), flush=True)
    if args.scan:
        for a in (.2, .25, .3, .35, .4):
            for eta in (400., 600., 1000.):
                if (a, eta) != (.3, 600.):
                    res['points'].append(analyse(a, eta))
                    print(json.dumps(res['points'][-1]), flush=True)
        for a in (.2, .3, .4):
            sx, s4 = stress_free(a, 600.)
            res['points'].append(analyse(a, 600., sx, s4))
            print(json.dumps(res['points'][-1]), flush=True)
    if args.checks:
        res['extra_checks'] = extra_checks()
        print(json.dumps(res['extra_checks'], indent=1), flush=True)
    with open(args.output, 'w') as fh:
        json.dump(res, fh, indent=1)
