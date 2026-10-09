"""Harmonic curvature of the D4 {111} glide wall, and the kernel that belongs with it.

Companion to the monograph's appendix on the dislocation core (the D4 dictionary and the
tangential sector) and to the chapter "The fine structure constant" (the equilibrium core width).

The lattice is the D4 contact lattice in its anti-plane sector: tangent spheres of diameter ell = 1
on the 24 nearest-neighbour bonds, normal stiffness k_n = 1, tangential stiffness k_t = r k_n at the
rolling point r = 1/(3 pi - 5) (where N^2 = 3r/(1 + 5r) = 1/pi). A node carries a displacement u
along the line L = (1,-1,0,0)/sqrt2 and a rotation phi in span(T, H), with T = (1,1,-2,0)/sqrt6 the
glide direction and H = (1,1,1,0)/sqrt3 the wall normal. The relative displacement at a contact is
delta = (u_j - u_i) L - (1/2)(phi_i + phi_j) x s, with s the spatial part of the bond, and the bond
energy is (1/2)[(k_n - k_t)(s.delta)^2 + k_t |delta|^2]. The contacts alone leave the transverse
curvature modulus negative, so a rotational spring (1/2) k_phi |phi_j - phi_i|^2 on every bond can
be added to supply gamma = mu ell^2, i.e. gamma/kappa_c = (pi - 2)/2, the value the chapter's
kernel uses.

The script computes, all per node or per cell of one node in every plane normal to H:
1. The long-wave moduli: frozen (k_n + 5k_t)/2 and rotation-relaxed (k_n + 2k_t)/2 per node, so
   N^2 = 3r/(1 + 5r), the same along T and H.
2. The Frenkel reference for a wall: the planes normal to H are h_p = d111/2 apart, so a column
   holds two nodes per interplanar spacing d111 and the Frenkel curvature mubar/d111 is
   2 c_rel/d111^2 = (3/2)(k_n + 2k_t) per column. Checked independently against the slope of the
   eigenstrain kernel at long wavelength, which reads mubar A per cell off the same stack.
3. The wall stiffness against a uniform slip along L, in Frenkel units:
   frozen, one cut between adjacent planes           (4/3)(k_n+5k_t)/(k_n+2k_t) = 4 pi/(3(pi-1))
   frozen, double joint with the half-layer at s/2   5/8 of that                 = 5 pi/(6(pi-1))
   rotations relaxed as well, contacts alone or with the rotational spring
   and, finally, every plane except the two that bound the double joint left free to slide
   (the exact harmonic partition at zero wavevector). The rigid-block values are repeated with
   all three rotation components free for a slip along the line and along the partial hop t,
   which agree at this order, and the one-cut relaxation is checked against the rolling
   residual eta_0 = 0.68834 of the appendix.
4. The kernel that belongs with that partition: hold the disregistry between the two bounding
   planes at Bloch wavevector k along T, relax everything else, subtract the zero-wavevector value.
   Compared with the bulk dispersion of the same lattice and with the continuum Cosserat kernel;
   the arctangent equilibrium is solved with each kernel against the same wall. Holding thicker
   slabs rigid shows how far the split itself moves the answer.

Units: ell = k_n = 1. Runtime: about a minute. Requires numpy, scipy.
"""
import itertools

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spl
from scipy.optimize import brentq

L_HAT = np.array([1, -1, 0, 0]) / np.sqrt(2)
T_HAT = np.array([1, 1, -2, 0]) / np.sqrt(6)
H_HAT = np.array([1, 1, 1, 0]) / np.sqrt(3)
HP = 1 / np.sqrt(6)                      # spacing of the planes normal to H
D111 = 2 * HP                            # {111} interplanar spacing sqrt(2/3)
DHOP = 1 / np.sqrt(3)                    # partial hop, the misfit period
R_ROLL = 1 / (3 * np.pi - 5)
N2 = 1 / np.pi
GAMMA_OVER_KAPPA = (np.pi - 2) / 2       # gamma = mu ell^2 with kappa_c/mu = 2/(pi - 2)


def roots():
    out = []
    for i, j in itertools.combinations(range(4), 2):
        for a, b in itertools.product((-1, 1), repeat=2):
            v = np.zeros(4)
            v[i], v[j] = a, b
            out.append(v / np.sqrt(2))
    return np.array(out)


class ContactLattice:
    """Anti-plane sector of the D4 contact lattice; x = (u, phi_T, phi_H) per node."""

    def __init__(self, r=R_ROLL, kphi=0.0):
        self.r, self.kphi = r, kphi
        self.bonds = []
        for n in roots():
            s = n[:3]
            A = np.zeros((3, 3))
            A[:, 0] = L_HAT[:3]
            A[:, 1] = -0.5 * np.cross(T_HAT[:3], s)
            A[:, 2] = -0.5 * np.cross(H_HAT[:3], s)
            K = (1 - r) * np.outer(s, s) + r * np.eye(3)
            self.bonds.append((A.T @ K @ A, n))

    def bloch(self, k):
        """Per-node matrix D(k) with energy per node (1/2) x^dag D x; k = (k_T, k_H)."""
        D = np.zeros((3, 3), complex)
        for M, n in self.bonds:
            e = np.exp(1j * (k[0] * (n @ T_HAT) + k[1] * (n @ H_HAT)))
            Tb = np.diag([e - 1, 1 + e, 1 + e])
            D += 0.5 * Tb.conj().T @ M @ Tb
            if self.kphi:
                P = np.diag([0, 1, 1]) * (e - 1)
                D += 0.5 * self.kphi * P.conj().T @ P
        return D


def relaxed(D):
    return (D[0, 0] - D[0, 1:] @ np.linalg.solve(D[1:, 1:], D[1:, 0])).real


def moduli(lat, eps=1e-4):
    """Frozen and relaxed per-node shear coefficients, gap 2 kappa_c V, transverse curvature."""
    D = lat.bloch(np.array([eps, 0.0]))
    froz, rel = D[0, 0].real / eps ** 2, relaxed(D) / eps ** 2
    g0 = lat.bloch(np.zeros(2))[1:, 1:].real
    k = 0.05
    gk = lat.bloch(np.array([k, 0.0]))[1:, 1:].real
    return froz, rel, g0[0, 0], (gk[1, 1] - g0[1, 1]) / k ** 2


def kphi_for(ratio):
    """Rotational spring giving transverse gamma/kappa_c = ratio: each bond adds 3 k_phi k^2
    per node to the rotational dispersion, since the 24 unit roots give sum (n.T)^2 = 6."""
    froz, rel, gap, gT = moduli(ContactLattice())
    return (ratio * gap / 2 - gT) / 3.0


def stack_matrix(lat, ke, P):
    """Sparse K with energy per cell (1/2) X^dag K X; plane p holds (u, phi_T, phi_H)."""
    rows, cols, vals = [], [], []
    A = np.diag([-1.0, 1.0, 1.0]).astype(complex)
    for M, n in lat.bonds:
        dp = int(round((n @ H_HAT) / HP))
        e = np.exp(1j * ke * (n @ T_HAT))
        B = np.diag([e, e, e])
        blocks = {(0, 0): 0.5 * A.conj().T @ M @ A, (0, 1): 0.5 * A.conj().T @ M @ B,
                  (1, 0): 0.5 * B.conj().T @ M @ A, (1, 1): 0.5 * B.conj().T @ M @ B}
        if lat.kphi:
            Ar, Br = np.diag([0, -1.0, -1.0]), np.diag([0, e, e])
            for key, (X, Y) in {(0, 0): (Ar, Ar), (0, 1): (Ar, Br), (1, 0): (Br, Ar), (1, 1): (Br, Br)}.items():
                blocks[key] = blocks[key] + 0.5 * lat.kphi * X.conj().T @ Y
        ps = np.arange(max(0, -dp), min(P, P - dp))
        for (a, b), blk in blocks.items():
            pa = ps if a == 0 else ps + dp
            pb = ps if b == 0 else ps + dp
            for i in range(3):
                for j in range(3):
                    if blk[i, j] != 0:
                        rows.append(3 * pa + i)
                        cols.append(3 * pb + j)
                        vals.append(np.full(len(ps), blk[i, j]))
    return sps.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(3 * P, 3 * P))


def held_energy(lat, ke, held, P=800):
    """2E per cell when the u of the listed planes (offsets from the wall's lower plane) take the
    given amplitudes and every other coordinate relaxes."""
    K = stack_matrix(lat, ke, P)
    p0 = P // 2 - 1
    ih = np.array([3 * (p0 + o) for o, _ in held])
    v = np.array([a for _, a in held], complex)
    ifr = np.setdiff1d(np.arange(3 * P), ih)
    Khh = K[ih][:, ih].toarray()
    Kfh = K[ifr][:, ih]
    X = spl.splu(K[ifr][:, ifr].tocsc()).solve(Kfh @ v)
    return 2 * 0.5 * (v.conj() @ Khh @ v - (Kfh @ v).conj() @ X).real


def profile_energy(lat, profile, P=80, relax_rotations=True):
    """2E per cell for a zero-wavevector slip profile u_p (unit slip), rotations frozen or relaxed."""
    K = stack_matrix(lat, 0.0, P).toarray().real
    p0 = P // 2 - 1
    u = np.array([profile(p - p0) for p in range(P)], float)
    iu = np.arange(0, 3 * P, 3)
    ip = np.setdiff1d(np.arange(3 * P), iu)
    E2 = u @ K[np.ix_(iu, iu)] @ u
    if relax_rotations:
        b = K[np.ix_(ip, iu)] @ u
        E2 -= b @ np.linalg.solve(K[np.ix_(ip, ip)], b)
    return E2


def wall_any_direction(slip, kphi=0.0, r=R_ROLL, P=80):
    """Rigid-block walls for a uniform slip along an arbitrary in-wall direction, with all three
    rotation components free: returns 2E per cell for (one cut, double joint), rotations relaxed."""
    E3 = np.eye(3)
    bl = []
    for n in roots():
        s = n[:3]
        A = np.zeros((3, 4))
        A[:, 0] = slip[:3]
        for a in range(3):
            A[:, 1 + a] = -0.5 * np.cross(E3[a], s)
        K = (1 - r) * np.outer(s, s) + r * np.eye(3)
        bl.append((A.T @ K @ A, int(round((n @ H_HAT) / HP))))
    K = np.zeros((4 * P, 4 * P))
    for M, dp in bl:
        for p in range(max(0, -dp), min(P, P - dp)):
            q = p + dp
            S = np.zeros((4, 4 * P))
            S[0, 4 * q] += 1
            S[0, 4 * p] -= 1
            for a in range(3):
                S[1 + a, 4 * p + 1 + a] += 1
                S[1 + a, 4 * q + 1 + a] += 1
            K += 0.5 * S.T @ M @ S
            if kphi:
                Sp = np.zeros((3, 4 * P))
                for a in range(3):
                    Sp[a, 4 * q + 1 + a] += 1
                    Sp[a, 4 * p + 1 + a] -= 1
                K += 0.5 * kphi * Sp.T @ Sp
    p0 = P // 2 - 1
    iu = np.arange(0, 4 * P, 4)
    ip = np.setdiff1d(np.arange(4 * P), iu)
    out = []
    for prof in (lambda m: 1.0 if m >= 1 else 0.0, lambda m: 1.0 if m >= 2 else (0.5 if m == 1 else 0.0)):
        u = np.array([prof(p - p0) for p in range(P)], float)
        b = K[np.ix_(ip, iu)] @ u
        out.append(u @ K[np.ix_(iu, iu)] @ u - b @ np.linalg.lstsq(K[np.ix_(ip, ip)], b, rcond=None)[0])
    return out


def eigenstrain_slope(lat, ke=0.02, P=1600):
    """2 K(ke)/ke, where K = 2E/|s|^2 is the energy of a slip wave s e^{i ke x} imposed as an
    eigenstrain on the bonds that cross the cut, everything relaxed. A screw kernel tends to
    (mubar A/2)|k| at long wavelength, so this tends to mubar A per cell."""
    K = stack_matrix(lat, ke, P)
    p0 = P // 2 - 1
    F = np.zeros(3 * P, complex)
    c0 = 0.0
    A = np.diag([-1.0, 1.0, 1.0]).astype(complex)
    for M, n in lat.bonds:
        dp = int(round((n @ H_HAT) / HP))
        if dp == 0:
            continue
        dx = n @ T_HAT
        B = np.diag([np.exp(1j * ke * dx)] * 3)
        for p in range(p0 - 3, p0 + 4):
            q = p + dp
            if not (min(p, q) <= p0 < max(p, q)):
                continue
            sig = np.array([(1.0 if q > p else -1.0) * np.exp(0.5j * ke * dx), 0, 0])
            v = M @ sig
            F[3 * p:3 * p + 3] += 0.5 * A.conj().T @ v
            F[3 * q:3 * q + 3] += 0.5 * B.conj().T @ v
            c0 += 0.5 * (sig.conj() @ v).real
    X = spl.spsolve(K, F)
    return 2 * (c0 - (F.conj() @ X).real) / ke


def continuum_kernel(k, mubarA, ratio):
    q = np.sqrt(2 / ratio)
    return 0.5 * mubarA / (1 - N2) * ((1 - N2) * k + N2 * k * k / np.sqrt(k * k + q * q))


def arctan_width(Kfun, Gamma):
    """Arctangent equilibrium int_0^inf (K(k)/k) e^{-2kw} dk = Gamma/2, returned as w/d."""
    ks = np.r_[np.geomspace(1e-4, 0.3, 300), np.linspace(0.3005, 12.0, 3000)]
    Kv = np.array([Kfun(k) for k in ks])
    return brentq(lambda w: np.trapezoid(Kv / ks * np.exp(-2 * ks * w), ks) - Gamma / 2,
                  0.05, 3.0, xtol=1e-12) / DHOP


def main():
    r = R_ROLL
    print('D4 contact lattice, anti-plane sector, r = k_t/k_n = 1/(3pi-5) = %.6f' % r)
    for label, kp in (('contacts alone', 0.0), ('with rotational spring, gamma = mu ell^2',
                                                 kphi_for(GAMMA_OVER_KAPPA))):
        lat = ContactLattice(kphi=kp)
        froz, rel, gap, gT = moduli(lat)
        DH = lat.bloch(np.array([0.0, 1e-4]))
        print('\n%s (k_phi = %.5f)' % (label, kp))
        print('  per node: frozen %.6f [(1+5r)/2 = %.6f]  relaxed %.6f [(1+2r)/2 = %.6f]  along H %.6f'
              % (froz, (1 + 5 * r) / 2, rel, (1 + 2 * r) / 2, relaxed(DH) / 1e-8))
        print('  N^2 = 1 - relaxed/frozen = %.8f (1/pi = %.8f); transverse gamma/kappa_c = %.4f'
              % (1 - rel / froz, N2, gT / (gap / 2)))
        mubarA = rel / HP
        GF = mubarA / D111
        print('  Frenkel curvature per cell 2 c_rel/d111^2 = %.6f [(3/2)(1+2r) = %.6f];'
              ' eigenstrain slope gives mubar A = %.5f (c_rel/h_p = %.5f)'
              % (GF, 1.5 * (1 + 2 * r), eigenstrain_slope(lat), mubarA))
        one_cut = lambda m: 1.0 if m >= 1 else 0.0
        wall = lambda m: 1.0 if m >= 2 else (0.5 if m == 1 else 0.0)
        print('  wall stiffness / Frenkel:')
        print('    frozen, one cut                        %.5f  [4pi/(3(pi-1)) = %.5f]'
              % (profile_energy(lat, one_cut, relax_rotations=False) / GF, 4 * np.pi / (3 * (np.pi - 1))))
        print('    frozen, double joint, half-layer at 1/2 %.5f  [5pi/(6(pi-1)) = %.5f]'
              % (profile_energy(lat, wall, relax_rotations=False) / GF, 5 * np.pi / (6 * (np.pi - 1))))
        print('    one cut, rotations relaxed              %.5f' % (profile_energy(lat, one_cut) / GF))
        print('    double joint, rotations relaxed         %.5f' % (profile_energy(lat, wall) / GF))
        G0 = held_energy(lat, 0.0, [(0, -0.5), (2, 0.5)])
        print('    every other plane free to slide too     %.5f' % (G0 / GF))
        for name, sd in (('the line L', L_HAT), ('the partial hop t', T_HAT)):
            a, b = wall_any_direction(sd, kphi=kp)
            print('    slip along %-18s (all three rotations free): one cut %.5f, double joint %.5f'
                  % (name, a / GF, b / GF))
        if kp == 0.0:
            eta0 = 0.68834033
            print('    one cut against eta_0: (2 + 10 eta_0 r)/((3/2)(1 + 2r)) = %.5f'
                  % ((2 + 10 * eta0 * r) / (1.5 * (1 + 2 * r))))
            continue
        # bulk dispersion against the continuum Cosserat stiffness mu_tot (k^2 + p^2)/(k^2 + q^2)
        q2 = 2 / GAMMA_OVER_KAPPA
        print('  bulk anti-plane stiffness c(k)/k^2 against the continuum (k along T, then H):')
        for kk in (0.5, 1, 2, 3, 4):
            cont = froz * (kk * kk + (rel / froz) * q2) / (kk * kk + q2)
            vals = [relaxed(lat.bloch(kk * np.array(d))) / kk ** 2 for d in ((1, 0), (0, 1))]
            print('    k = %.1f: lattice %.4f, %.4f   continuum %.4f   ratio %.3f'
                  % (kk, vals[0], vals[1], cont, vals[0] / cont))
        # kernels and the arctangent equilibrium against the same wall
        kg = np.r_[np.geomspace(0.01, 0.3, 12), np.linspace(0.35, 12, 90)]
        Kc = lambda k: continuum_kernel(k, mubarA, GAMMA_OVER_KAPPA)
        print('  arctangent equilibrium w/d against the wall of each partition:')
        print('    continuum kernel, Frenkel wall            %.4f' % arctan_width(Kc, GF))
        for name, held in (('the two planes bounding the double joint', [(0, -0.5), (2, 0.5)]),
                           ('two rigid planes each side', [(-1, -0.5), (0, -0.5), (2, 0.5), (3, 0.5)]),
                           ('three rigid planes each side', [(-2, -0.5), (-1, -0.5), (0, -0.5),
                                                             (2, 0.5), (3, 0.5), (4, 0.5)])):
            G = held_energy(lat, 0.0, held)
            Kel = np.array([held_energy(lat, k, held) for k in kg]) - G
            ratio = Kel / kg
            Kl = lambda k, ratio=ratio: k * np.interp(k, kg, ratio)
            print('    held: %-40s wall %.4f GF; kernel/continuum at k = 1, 2: %.3f, %.3f;'
                  ' w/d continuum kernel %.4f, lattice kernel %.4f'
                  % (name, G / GF, Kl(1.0) / Kc(1.0), Kl(2.0) / Kc(2.0), arctan_width(Kc, G), arctan_width(Kl, G)))


if __name__ == '__main__':
    main()
