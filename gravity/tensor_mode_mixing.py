#!/usr/bin/env python3
"""
tensor_mode_mixing.py -- does a strained vacuum couple the helicity-2
bond-length mode P to the lattice's gapless channels, and to which one?

Seven modes per wavevector on the FCC contact lattice: displacement (3),
microrotation (3) and one helicity-2 bond-length amplitude P.  A background
strain eps_bg modulates each bond's stiffness through the anharmonicity,
k_n -> k_n (1 + xi n.eps_bg.n) with xi = -7, which is what switches the
coupling on.  The inertia of P is not derived; it is set to m l^2 here, and the
crossing energy scales as its inverse square root.

Result (wave along z): a uniaxial background couples P+ to the longitudinal
displacement only; a shear with one leg along z couples P+ to the transverse
displacement and opens an avoided crossing with the photon, of size
dw/w_P ~ 1.9 eps_bg, while Px stays decoupled.  Helicity about the propagation
axis is conserved: P carries +-2, and only an odd-helicity background reaches
the photon's +-1.
"""
import numpy as np
from scipy.optimize import minimize_scalar
from impulse_response_audit import fcc_shell, cross_matrix, R_RATIO

XI, SHELL, KZ = -7.0, fcc_shell(), np.array([0.0, 0.0, 1.0])
MASS = np.array([1, 1, 1, 0.533, 0.533, 0.533, 1.0])   # node, J = 0.533 m l^2, I_P = m l^2 (assumed)


def blocks(k, bg, P):
    D = np.zeros((7, 7), complex)
    for R in SHELL:
        n = R / np.linalg.norm(R); kn = 1 + XI * (n @ bg @ n); kt = R_RATIO * kn
        ph = np.exp(1j * np.dot(k, R)); A = (ph - 1) * np.eye(3); B = -0.5 * (ph + 1) * cross_matrix(R)
        P1 = np.outer(n, n); Q = np.eye(3) - P1; Mn, Mu, Mp = P1 @ A, Q @ A, Q @ B
        D[:3, :3] += kn * (Mn.conj().T @ Mn) + kt * (Mu.conj().T @ Mu)
        D[:3, 3:6] += kt * (Mu.conj().T @ Mp); D[3:6, :3] += kt * (Mp.conj().T @ Mu); D[3:6, 3:6] += kt * (Mp.conj().T @ Mp)
        d = n @ P @ n; row = -kn * d * (Mn.conj().T @ n)
        D[:3, 6] += row; D[6, :3] += row.conj(); D[6, 6] += kn * d * d
    return 0.5 * D


def splitting(k, bg, P):
    Mi = np.diag(1 / np.sqrt(MASS)); H = Mi @ blocks(k * KZ, bg, P) @ Mi
    w = np.sqrt(np.clip(np.linalg.eigvalsh((H + H.conj().T) / 2), 0, None))
    wp = np.sqrt(blocks(0 * KZ, bg, P)[6, 6].real); i = np.argsort(abs(w - wp))[:2]
    return abs(w[i[0]] - w[i[1]]), wp


if __name__ == "__main__":
    Pp = np.diag([1, -1, 0.0]) / np.sqrt(2); Px = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0.0]]) / np.sqrt(2)
    backgrounds = {"uniaxial x": np.diag([1, 0, 0.0]), "shear xz": np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0.0]])}
    for bname, B in backgrounds.items():
        for pname, P in (("P+", Pp), ("Px", Px)):
            for eps in (1e-3, 1e-2):
                r = minimize_scalar(lambda k: splitting(k, eps * B, P)[0], bounds=(0.6, 1.2), method="bounded", options={"xatol": 1e-10})
                dw, wp = splitting(r.x, eps * B, P)
                print(f"{bname:10s} {pname}  eps={eps:6.0e}  dw/w_P={dw / wp:9.2e}  at k l = {r.x:.3f}")
