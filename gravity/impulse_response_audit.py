#!/usr/bin/env python3
"""
impulse_response_audit.py -- what static kernel each kind of source produces on
the micropolar contact lattice.

The gravity chapter's no-go is stated over *field* types: no linear field of the
lattice carries helicity 2, and a transverse-traceless bond-length change pays
the full shear modulus with nothing to relax it.  This script asks the same
question over *source* types instead, which is the sharper statement: given the
lattice's own dynamical matrix, enumerate every source the symmetry allows and
classify the static response kernel each one produces as k -> 0.

    pole      E(k) ~ 1/k^2   long-ranged, Coulomb or Newton
    contact   E(k) ~ const    no far field
    derivative contact  E(k) ~ k^2

The lattice is the FCC slice with the framework's rolling contact: a normal
spring k_n along each bond and a tangential spring k_t = r k_n resisting the
part of the relative displacement that the two nodes' microrotations do not
account for.  Each node carries a displacement u and a microrotation phi, so
the dynamical matrix is 6x6 at every wavevector.

Sources, by the irreducible representation they belong to:

    force            T_1u   a point force on a node
    torque           T_1g   a point couple on a node
    eigenstrain A1g  A_1g   an isotropic change of natural bond length
    eigenstrain TT   E_g + T_2g  a traceless change of natural bond length,
                            which is the bond-length tensor P_ij of the
                            gravity chapter

For each source the script reports the exponent p in E(k) ~ k^p, fitted over two
decades in k, and the residual energy that the lattice cannot relax.  A pole
(p = -2) in the transverse-traceless channel would be a tensor sector appearing
from inside the lattice; a contact (p = 0) is the gap.
"""

import numpy as np

# ---------------------------------------------------------------- the lattice
R_RATIO = 1.0 / (2.0 * np.pi - 3.0)   # k_t/k_n at the slice's rolling point
K_N = 1.0                              # normal contact stiffness, units of itself


def fcc_shell(a=1.0):
    """The twelve nearest neighbours of an FCC lattice, bond length a."""
    out = []
    for s1 in (+1, -1):
        for s2 in (+1, -1):
            out += [np.array([s1, s2, 0.0]), np.array([s1, 0.0, s2]), np.array([0.0, s1, s2])]
    seen, shell = set(), []
    for v in out:
        key = tuple(np.round(v, 9))
        if key not in seen:
            seen.add(key)
            shell.append(v / np.sqrt(2.0) * a)
    return shell


def cross_matrix(v):
    """The matrix C with C @ w = v x w."""
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


def dynamical_matrix(k, shell, k_n=K_N, r=R_RATIO):
    """The 6x6 static stiffness at wavevector k, ordered (u_x,u_y,u_z,phi_x,phi_y,phi_z).

    Bond energy, per bond:  1/2 k_n (Delta.n)^2 + 1/2 k_t |Delta_perp - thetabar x R|^2
    with Delta = u(R) - u(0) and thetabar the mean microrotation of the two ends.
    """
    k_t = r * k_n
    D = np.zeros((6, 6), dtype=complex)
    for R in shell:
        n = R / np.linalg.norm(R)
        phase = np.exp(1j * np.dot(k, R))
        A = (phase - 1.0) * np.eye(3)                 # Delta      = A u
        B = -0.5 * (phase + 1.0) * cross_matrix(R)    # -thetabar x R = B phi
        P_par = np.outer(n, n)
        P_perp = np.eye(3) - P_par
        # normal part: k_n |P_par A u|^2
        Mn = P_par @ A
        # tangential part: k_t |P_perp (A u) + B phi|^2, the lever arm R carried in B
        Mt_u = P_perp @ A
        Mt_p = P_perp @ B
        D[:3, :3] += k_n * (Mn.conj().T @ Mn) + k_t * (Mt_u.conj().T @ Mt_u)
        D[:3, 3:] += k_t * (Mt_u.conj().T @ Mt_p)
        D[3:, :3] += k_t * (Mt_p.conj().T @ Mt_u)
        D[3:, 3:] += k_t * (Mt_p.conj().T @ Mt_p)
    return 0.5 * D          # each bond shared between its two nodes


# ---------------------------------------------------------------- the sources
def source_force(k, shell, direction=np.array([1.0, 0.0, 0.0])):
    """A point force: couples directly to u, with no k dependence."""
    S = np.zeros(6, dtype=complex)
    S[:3] = direction
    return S


def source_torque(k, shell, direction=np.array([1.0, 0.0, 0.0])):
    """A point couple: couples directly to phi."""
    S = np.zeros(6, dtype=complex)
    S[3:] = direction
    return S


def source_eigenstrain(k, shell, eps_star, k_n=K_N):
    """A prescribed change of natural bond length, delta_b = n.eps*.n.

    The linear term in the energy is -k_n sum_b delta_b (Delta.n), so the
    generalised force carries one factor of the phase difference, hence one
    factor of k at long wavelength.
    """
    S = np.zeros(6, dtype=complex)
    for R in shell:
        n = R / np.linalg.norm(R)
        delta = n @ eps_star @ n
        phase = np.exp(1j * np.dot(k, R))
        S[:3] += k_n * delta * (phase.conjugate() - 1.0) * n
    return 0.5 * S


def eigenstrain_bare_energy(shell, eps_star, k_n=K_N):
    """The energy the eigenstrain would cost with the lattice held rigid."""
    tot = 0.0
    for R in shell:
        n = R / np.linalg.norm(R)
        tot += k_n * (n @ eps_star @ n) ** 2
    return 0.5 * 0.5 * tot


# ---------------------------------------------------------------- the audit
def interaction_energy(k, shell, S):
    """E(k) = 1/2 S^dagger D^+ S, the energy the medium returns to the source."""
    D = dynamical_matrix(k, shell)
    Dp = np.linalg.pinv(D, rcond=1e-12)
    return float(np.real(0.5 * S.conj() @ Dp @ S))


def slope(ks, vals):
    """Fitted exponent p in val ~ k^p."""
    good = [(k, v) for k, v in zip(ks, vals) if v > 0]
    if len(good) < 3:
        return float("nan")
    x = np.log([g[0] for g in good])
    y = np.log([g[1] for g in good])
    return float(np.polyfit(x, y, 1)[0])


def classify(p):
    if p < -1.5:
        return "pole (long-ranged)"
    if -0.5 < p < 0.5:
        return "contact"
    if 1.5 < p < 2.5:
        return "derivative contact"
    return f"other (p = {p:.2f})"


def main():
    shell = fcc_shell()
    khat = np.array([0.0, 0.0, 1.0])            # propagation along z
    ks = np.logspace(-4, -2, 9)

    eps_iso = np.eye(3) / 3.0
    eps_tt = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]) / np.sqrt(2.0)
    eps_mixed = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]) / np.sqrt(2.0)

    cases = [
        ("force            (T_1u)", lambda k: source_force(k * khat, shell), None),
        ("torque           (T_1g)", lambda k: source_torque(k * khat, shell), None),
        ("eigenstrain iso  (A_1g)", lambda k: source_eigenstrain(k * khat, shell, eps_iso), eps_iso),
        ("eigenstrain long (E_g)",  lambda k: source_eigenstrain(k * khat, shell, eps_mixed), eps_mixed),
        ("eigenstrain TT   (E_g,T_2g)", lambda k: source_eigenstrain(k * khat, shell, eps_tt), eps_tt),
    ]

    print(f"micropolar contact lattice, FCC shell, k_t/k_n = {R_RATIO:.5f}")
    print(f"{'source':28s} {'E(k) exponent':>14s}  {'kernel':22s} {'unrelaxed energy':>18s}")
    print("-" * 88)
    for name, make, eps in cases:
        vals = [interaction_energy(k * khat, shell, make(k)) for k in ks]
        p = slope(ks, vals)
        if eps is None:
            residual = ""
        else:
            bare = eigenstrain_bare_energy(shell, eps)
            relaxed = vals[0]
            residual = f"{bare - relaxed:18.6f}"
        print(f"{name:28s} {p:14.3f}  {classify(p):22s} {residual:>18s}")

    print()
    print("Reading: a pole is a long-range force, a contact is a gap.  A pole in the")
    print("transverse-traceless row would be a tensor sector inside the lattice.")


if __name__ == "__main__":
    main()
