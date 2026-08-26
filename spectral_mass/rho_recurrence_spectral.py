"""Does the rho's first recurrence close on the spectral mass formula?

The framework carries two accounts of a radial excitation and never reconciles
them.  The Roper is the proton's 13-node shell extended by six second-neighbour
sites to N = 19, its mass closing on the master formula

    m = N m_0 - N (4 - lambda) m_e.

The vector tower is a Regge ladder, m_n^2 = m_rho^2 + n (2 pi sigma).  One is
linear in node count, the other in mass squared, and nothing connects them.

This script runs the Roper construction on the rho.  Take the 11-node crossed
fault, add the same six second-neighbour sites along the cubic axes, classify
the extended cluster's modes by D_2h irrep, and ask where the rho's B_3u mode
goes.  If the resulting mass lands near the ladder's 1348 MeV the two accounts
are one; if it cannot, the chapter is carrying two theories of one object.
"""
import sys
import os

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import composite_clusters as cc
from cosserat_classifier import build_cosserat_matrix

M_0 = 70.0253243          # MeV, node mass m_e/alpha
M_E = 0.51099895069       # MeV, CODATA electron mass
SIGMA = 193583.84         # MeV^2, derived string tension
M_RHO = 775.286           # MeV, the N = 11 closure


def mass(N, lam):
    """The master spectral formula."""
    return N * M_0 - N * (4.0 - lam) * M_E


def lam_for(N, m):
    """Invert it: which eigenvalue would a cluster of N nodes need?"""
    return 4.0 + (m - N * M_0) / (N * M_E)


def add_second_shell(cluster):
    """The six second-neighbour FCC sites along the cubic axes, as for the Roper."""
    s = np.sqrt(2.0)
    extra = np.array([[+s, 0, 0], [-s, 0, 0],
                      [0, +s, 0], [0, -s, 0],
                      [0, 0, +s], [0, 0, -s]])
    return np.vstack([cluster, extra])


def d2h_ops():
    """The eight D_2h elements, in the crossed fault's own axes."""
    ax = {'X': np.array([1., 0., 0.]),
          'Y': np.array([0., 1., 1.]) / np.sqrt(2.0),
          'Z': np.array([0., -1., 1.]) / np.sqrt(2.0)}
    def c2(n):
        return 2.0 * np.outer(n, n) - np.eye(3)
    E = np.eye(3)
    I = -np.eye(3)
    R = {'E': E, 'C2X': c2(ax['X']), 'C2Y': c2(ax['Y']), 'C2Z': c2(ax['Z'])}
    for k in list(R):
        R['i*' + k] = I @ R[k]
    return R


def rep(R, coords):
    """Lift a 3x3 rotation to the 6N-dimensional (u, phi) space.

    u is polar and phi is axial, so phi picks up det(R): a reflection reverses
    a rotation's sense where it leaves a displacement's direction alone.
    """
    n = len(coords)
    D = np.zeros((6 * n, 6 * n))
    d = round(np.linalg.det(R))
    for i, x in enumerate(coords):
        y = R @ x
        j = int(np.argmin(np.linalg.norm(coords - y, axis=1)))
        assert np.linalg.norm(coords[j] - y) < 1e-8, "cluster not invariant"
        D[3*j:3*j+3, 3*i:3*i+3] = R                        # u block, polar
        D[3*n+3*j:3*n+3*j+3, 3*n+3*i:3*n+3*i+3] = d * R    # phi block, axial
    return D


CHAR = {  # D_2h characters, order E C2Z C2Y C2X i iC2Z iC2Y iC2X
    'Ag':  (1, 1, 1, 1, 1, 1, 1, 1),
    'B1g': (1, 1, -1, -1, 1, 1, -1, -1),
    'B2g': (1, -1, 1, -1, 1, -1, 1, -1),
    'B3g': (1, -1, -1, 1, 1, -1, -1, 1),
    'Au':  (1, 1, 1, 1, -1, -1, -1, -1),
    'B1u': (1, 1, -1, -1, -1, -1, 1, 1),
    'B2u': (1, -1, 1, -1, -1, 1, -1, 1),
    'B3u': (1, -1, -1, 1, -1, 1, 1, -1),
}
ORDER = ['E', 'C2Z', 'C2Y', 'C2X', 'i*E', 'i*C2Z', 'i*C2Y', 'i*C2X']


def classify(coords):
    """Eigenvalues of the Cosserat matrix, each tagged with its D_2h irrep."""
    M = build_cosserat_matrix(coords, K_u=1.0, K_phi=1.0, alpha=1.0)
    R = d2h_ops()
    D = {k: rep(R[k], coords) for k in ORDER}
    for k, Dk in D.items():
        assert np.allclose(Dk @ M, M @ Dk, atol=1e-8), "%s is not a symmetry" % k
    vals, vecs = np.linalg.eigh(M)
    proj = {}
    for name, chi in CHAR.items():
        P = sum(c * D[k] for c, k in zip(chi, ORDER)) / 8.0
        proj[name] = P
    tags = []
    for i in range(len(vals)):
        v = vecs[:, i]
        w = {nm: float(v @ (P @ v)) for nm, P in proj.items()}
        tags.append(max(w, key=w.get))
    return vals, vecs, tags


def main():
    base, _ = cc.cluster_crossed_fault()
    ext = add_second_shell(base)

    print("=" * 74)
    print("the rho, as the chapter closes it")
    vals, vecs, tags = classify(base)
    b3u = [(v, i) for i, (v, t) in enumerate(zip(vals, tags)) if t == 'B3u']
    lam_rho = min(v for v, _ in b3u if v > 1e-6)
    print("  N = %d   lowest B_3u lambda = %.4f   m = %.2f MeV  (chapter 775.29)"
          % (len(base), lam_rho, mass(len(base), lam_rho)))
    print("  eigenvalue range %.4f to %.4f" % (vals.min(), vals.max()))
    rho_vec = vecs[:, [i for v, i in b3u if abs(v - lam_rho) < 1e-9][0]]

    print("\n" + "=" * 74)
    print("the same extension the Roper uses, six second-neighbour sites")
    vals2, vecs2, tags2 = classify(ext)
    N2 = len(ext)
    print("  N = %d   eigenvalue range %.4f to %.4f" % (N2, vals2.min(), vals2.max()))

    # embed the rho mode in the larger space, zero on the new sites
    emb = np.zeros(6 * N2)
    emb[:6 * len(base)] = rho_vec
    emb /= np.linalg.norm(emb)

    cand = [(i, vals2[i], abs(float(emb @ vecs2[:, i])))
            for i in range(len(vals2)) if tags2[i] == 'B3u' and vals2[i] > 1e-6]
    cand.sort(key=lambda r: -r[2])
    print("\n  B_3u modes by overlap with the embedded rho mode:")
    print("    %-8s %-10s %-10s %s" % ("lambda", "overlap", "m [MeV]", "vs ladder 1348.1"))
    for i, lam, ov in cand[:6]:
        m = mass(N2, lam)
        print("    %-8.4f %-10.4f %-10.2f %+.1f%%" % (lam, ov, m, 100 * (m / 1348.1 - 1)))

    print("\n" + "=" * 74)
    print("what the formula would REQUIRE")
    for label, m in (("ladder n=1", np.sqrt(M_RHO**2 + 2 * np.pi * SIGMA)),
                     ("observed rho(1450)", 1465.0),
                     ("observed rho(1700)", 1720.0)):
        print("  %-20s m = %7.1f MeV" % (label, m))
        for N in (11, 13, 17, 19, 21, 23):
            lam = lam_for(N, m)
            flag = "" if 0.0 <= lam <= vals2.max() + 2 else "   <- outside the spectrum"
            print("      N = %2d  needs lambda = %8.3f%s" % (N, lam, flag))
    print("\n  reachable mass window at N = %d, lambda in [%.2f, %.2f]:  %.1f to %.1f MeV"
          % (N2, 0.0, vals2.max(), mass(N2, 0.0), mass(N2, vals2.max())))


if __name__ == "__main__":
    main()
