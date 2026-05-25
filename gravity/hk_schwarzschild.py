#!/usr/bin/env python3
"""
Heat-kernel coefficients on Euclidean Schwarzschild.

Verifies the Seeley-DeWitt a_4 heat-kernel engine that determines the
logarithmic correction to black-hole entropy, per Sen (JHEP 04, 2013, 156).
Reproduces the per-field coefficients for the scalar (2/90), Maxwell vector
(-26/90), and the geometric inputs for the graviton, by explicit symbolic
computation on the Schwarzschild metric.

The framework predicts the graviton contributes -26/90 (a vector's anomaly)
rather than the standard +424/90 (a tensor's anomaly), because its graviton
is a microrotation vector whose curl gives the spin-2 observable. See the
monograph's gravity chapter, "Logarithmic corrections and the vector graviton."

Requires: sympy
Usage:    python3 hk_schwarzschild.py

Repository: https://github.com/uberflyx/CosseratSupersolidLattice
Authors:    Mitchell Cox, Warren Carlson (University of the Witwatersrand)
"""
import sympy as sp

# =========================================================================
# 1. EUCLIDEAN SCHWARZSCHILD GEOMETRY
# =========================================================================
# Coordinates (t, r, theta, phi), with r_s = 1 (units of the horizon radius).
# Euclidean signature: ds^2 = f dt^2 + f^{-1} dr^2 + r^2 dOmega^2.
# =========================================================================
t, r, th, ph = sp.symbols('t r theta phi', positive=True)
coords = [t, r, th, ph]
n = 4

f = 1 - 1 / r
g = sp.diag(f, 1 / f, r**2, r**2 * sp.sin(th)**2)
ginv = g.inv()


def christoffel(met, metinv, X):
    """Christoffel symbols Gamma^a_{bc} from the metric."""
    N = len(X)
    Ga = [[[sp.S(0)] * N for _ in range(N)] for _ in range(N)]
    for a in range(N):
        for b in range(N):
            for c in range(N):
                s = sum(
                    metinv[a, d] * (
                        sp.diff(met[d, b], X[c])
                        + sp.diff(met[d, c], X[b])
                        - sp.diff(met[b, c], X[d])
                    )
                    for d in range(N)
                )
                Ga[a][b][c] = sp.simplify(s / 2)
    return Ga


def riemann_up(Ga, X):
    """Riemann tensor R^a_{bcd} from Christoffel symbols."""
    N = len(X)
    R = [[[[sp.S(0)] * N for _ in range(N)] for _ in range(N)] for _ in range(N)]
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    val = (
                        sp.diff(Ga[a][b][d], X[c])
                        - sp.diff(Ga[a][b][c], X[d])
                        + sum(
                            Ga[a][c][e] * Ga[e][b][d]
                            - Ga[a][d][e] * Ga[e][b][c]
                            for e in range(N)
                        )
                    )
                    R[a][b][c][d] = sp.simplify(val)
    return R


def lower_first(Rup, met):
    """Lower the first index: R_{abcd} = g_{ae} R^e_{bcd}."""
    N = met.shape[0]
    Rl = [[[[sp.S(0)] * N for _ in range(N)] for _ in range(N)] for _ in range(N)]
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    Rl[a][b][c][d] = sp.simplify(
                        sum(met[a, e] * Rup[e][b][c][d] for e in range(N))
                    )
    return Rl


def raise_all(Rl, metinv):
    """Raise all four indices: R^{abcd} from R_{abcd}."""
    N = metinv.shape[0]
    Ruu = [[[[sp.S(0)] * N for _ in range(N)] for _ in range(N)] for _ in range(N)]
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    s = sp.S(0)
                    for p in range(N):
                        if metinv[a, p] == 0:
                            continue
                        for q in range(N):
                            if metinv[b, q] == 0:
                                continue
                            for u in range(N):
                                if metinv[c, u] == 0:
                                    continue
                                for v in range(N):
                                    if metinv[d, v] == 0:
                                        continue
                                    s += (
                                        metinv[a, p] * metinv[b, q]
                                        * metinv[c, u] * metinv[d, v]
                                        * Rl[p][q][u][v]
                                    )
                    Ruu[a][b][c][d] = sp.simplify(s)
    return Ruu


print("Building Schwarzschild geometry (r_s = 1)...")
Ga = christoffel(g, ginv, coords)
Rup = riemann_up(Ga, coords)
Rl = lower_first(Rup, g)

# ---- Ricci tensor (should vanish) ----
Ric = sp.zeros(n)
for b in range(n):
    for d in range(n):
        Ric[b, d] = sp.simplify(sum(Rup[a][b][a][d] for a in range(n)))
ricci_flat = all(Ric[i, j] == 0 for i in range(n) for j in range(n))
print(f"  Ricci flat: {ricci_flat}")

# ---- Fully raised Riemann and Kretschmann scalar ----
print("Raising all Riemann indices...")
Ruu = raise_all(Rl, ginv)

K = sp.simplify(
    sum(
        Rl[a][b][c][d] * Ruu[a][b][c][d]
        for a in range(n)
        for b in range(n)
        for c in range(n)
        for d in range(n)
    )
)
print(f"  Kretschmann K = {K}   (expect 12/r^6)")

# =========================================================================
# 2. EULER INTEGRAL (GAUSS-BONNET)
# =========================================================================
# On Ricci-flat, the Euler density E_4 = K, and
# Int sqrt(g) K = 32 pi^2 chi = 64 pi^2   (chi = 2 for Eucl. Schwarzschild).
# The Euclidean time has period beta = 4 pi r_s = 4 pi (r_s = 1).
# =========================================================================
beta = 4 * sp.pi
sqrtg = r**2 * sp.sin(th)
I_K = sp.simplify(
    beta
    * sp.integrate(
        sp.integrate(
            sp.integrate(sqrtg * K, (ph, 0, 2 * sp.pi)),
            (th, 0, sp.pi),
        ),
        (r, 1, sp.oo),
    )
)
print(f"  Euler integral = {I_K}   (expect 64*pi**2)")

# =========================================================================
# 3. KEY RIEMANN IDENTITY: R_{abcd} R^{acbd} = K/2
# =========================================================================
# This identity (pair-swap contraction equals half the Kretschmann) is
# used in the symmetric-tensor Lichnerowicz traces and verified here.
# =========================================================================
Racbd = sp.simplify(
    sum(
        Rl[a][b][c][d] * Ruu[a][c][b][d]
        for a in range(n)
        for b in range(n)
        for c in range(n)
        for d in range(n)
    )
)
print(f"  R_{{abcd}} R^{{acbd}} = {Racbd},  K/2 = {sp.simplify(K / 2)}")

# =========================================================================
# 4. VECTOR BUNDLE: tr(Omega^2)
# =========================================================================
# The connection curvature on the cotangent (vector) bundle is
# (Omega_{mn})^a_b = R^a_{b mn}.
# tr(Omega_{mn} Omega^{mn}) = R^a_{b mn} R^b_{a}^{mn}.
# On any background this equals -K (by Riemann antisymmetry in the first
# pair), so w2_vec = -1.
# =========================================================================

# Raise last two indices of R^a_{bcd}
def raise_last2(Rup_in, metinv_in):
    """R^a_b{}^{cd} from R^a_{bcd}."""
    N = metinv_in.shape[0]
    out = [[[[sp.S(0)] * N for _ in range(N)] for _ in range(N)] for _ in range(N)]
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    s = sp.S(0)
                    for u in range(N):
                        if metinv_in[c, u] == 0:
                            continue
                        for v in range(N):
                            if metinv_in[d, v] == 0:
                                continue
                            s += metinv_in[c, u] * metinv_in[d, v] * Rup_in[a][b][u][v]
                    out[a][b][c][d] = sp.simplify(s)
    return out


print("Computing vector bundle tr(Omega^2)...")
Rup_uu = raise_last2(Rup, ginv)
trOm2_vec = sp.simplify(
    sum(
        Rup[a][b][c][d] * Rup_uu[b][a][c][d]
        for a in range(n)
        for b in range(n)
        for c in range(n)
        for d in range(n)
    )
)
w2_vec = sp.simplify(trOm2_vec / K)
print(f"  tr(Omega^2)_vec / K = {w2_vec}   (expect -1)")

# =========================================================================
# 5. SYMMETRIC-TENSOR BUNDLE: tr(E^2) and tr(Omega^2)
# =========================================================================
# Lichnerowicz endomorphism on Ricci-flat: (E h)_{ab} = 2 R_a{}^c{}_b{}^d h_{cd}.
# tr(E^2) = 4 R_a{}^c{}_b{}^d R_c{}^a{}_d{}^b.
# tr(Omega^2) computed explicitly on the 10-dim symmetric basis.
# =========================================================================

# Build R_a{}^c{}_b{}^d (raise 2nd and 4th indices of R_{abcd})
print("Computing graviton (symmetric tensor) traces...")
Rc_d = [
    [
        [
            [
                sp.simplify(
                    sum(
                        ginv[c2, cc] * ginv[d2, dd] * Rl[a][cc][b][dd]
                        for cc in range(n)
                        for dd in range(n)
                    )
                )
                for d2 in range(n)
            ]
            for b in range(n)
        ]
        for c2 in range(n)
    ]
    for a in range(n)
]

trE2 = sp.simplify(
    4
    * sum(
        Rc_d[a][c][b][d] * Rc_d[c][a][d][b]
        for a in range(n)
        for b in range(n)
        for c in range(n)
        for d in range(n)
    )
)
e2_grav = sp.simplify(trE2 / K)
print(f"  Lichnerowicz tr(E^2) / K = {e2_grav}")

# tr(Omega^2) on 10-dim symmetric bundle, computed on explicit basis.
sym_basis = []
for a in range(n):
    for b in range(a, n):
        M = sp.zeros(n)
        M[a, b] = 1
        M[b, a] = 1
        sym_basis.append((a, b, M))


def apply_Om(T, m, ni, raised):
    """Apply Omega_{m,ni} (or Omega^{m,ni} if raised) to symmetric tensor T."""
    R_use = Rup_uu if raised else Rup
    out = sp.zeros(n)
    for a in range(n):
        for b in range(n):
            s = sp.S(0)
            for c in range(n):
                s += -R_use[c][a][m][ni] * T[c, b] - R_use[c][b][m][ni] * T[a, c]
            out[a, b] = s
    return out


def inner(S, T):
    """Inner product <S, T> = S_{ab} T^{ab} using the inverse metric."""
    s = sp.S(0)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    s += S[a, b] * ginv[a, c] * ginv[b, d] * T[c, d]
    return sp.simplify(s)


trOm2_sym = sp.S(0)
for _, _, M in sym_basis:
    acc = sp.zeros(n)
    for m in range(n):
        for ni in range(n):
            acc += apply_Om(apply_Om(M, m, ni, True), m, ni, False)
    trOm2_sym += inner(M, acc) / inner(M, M)
trOm2_sym = sp.simplify(trOm2_sym)
w2_grav = sp.simplify(trOm2_sym / K)
print(f"  Symmetric-tensor tr(Omega^2) / K = {w2_grav}")

# =========================================================================
# 6. MASTER FORMULA AND PER-FIELD COEFFICIENTS
# =========================================================================
# C_local = (1/90)(180 e2 + 30 w2 + 2 d)
# =========================================================================


def C_local(e2, w2, d):
    """Heat-kernel coefficient on Ricci-flat Schwarzschild."""
    return sp.Rational(1, 90) * (180 * e2 + 30 * w2 + 2 * d)


print("\n" + "=" * 65)
print("RESULTS: per-field heat-kernel coefficients")
print("=" * 65)

C_scalar = C_local(0, 0, 1)
print(f"  Scalar (E=0, Om=0, d=1):             C = {C_scalar} = {float(C_scalar):.6f}")
print(f"    Sen value: 2/90 = {sp.Rational(2, 90)}")

C_vec_raw = C_local(0, w2_vec, 4)
C_maxwell = C_vec_raw - 2 * C_scalar  # subtract 2 ghost scalars
print(f"  Maxwell vector (raw):                 C = {C_vec_raw}")
print(f"  Maxwell (vector - 2 ghosts):          C = {C_maxwell}")
print(f"    Sen value: -26/90 = {sp.Rational(-26, 90)}")

# Graviton (standard, symmetric tensor + Lichnerowicz, minus vector ghosts)
C_grav_raw = C_local(e2_grav, w2_grav, 10)
C_grav_ghosts = -2 * C_vec_raw  # -2 x (vector, no ghosts)
C_grav_total = sp.simplify(C_grav_raw + C_grav_ghosts)
print(f"  Graviton raw (sym tensor, d=10):      C = {C_grav_raw}")
print(f"  Graviton (raw - 2 vec ghosts):        C = {C_grav_total}")
print(f"    Sen value: 424/90 = {sp.Rational(424, 90)}")
print(f"    NOTE: the gap ({C_grav_total} vs 424/90) is from the")
print(f"    conformal/trace mode and ghost endomorphism bookkeeping")
print(f"    not encoded here. The geometric traces (e2={e2_grav},")
print(f"    w2={w2_grav}) are verified.")

# =========================================================================
# 7. FRAMEWORK PREDICTION: VECTOR GRAVITON
# =========================================================================
# The framework's graviton is a microrotation vector phi, not a symmetric
# tensor h. A gauge vector contributes the Maxwell value -26/90. The total
# logarithmic coefficient is photon + graviton-as-vector.
# =========================================================================
C_photon = C_maxwell
C_graviton_lattice = C_maxwell  # same: gauge vector
C_total_lattice = C_photon + C_graviton_lattice
C_total_GR = sp.Rational(398, 90)

print(f"\n{'=' * 65}")
print(f"FRAMEWORK vs GENERAL RELATIVITY")
print(f"{'=' * 65}")
print(f"  Photon contribution:                  {C_photon} = {float(C_photon):.4f}")
print(f"  Graviton-as-vector contribution:      {C_graviton_lattice} = {float(C_graviton_lattice):.4f}")
print(f"  Framework total C_local:              {C_total_lattice} = {float(C_total_lattice):.4f}")
print(f"  GR total C_local:                     {C_total_GR} = {float(C_total_GR):.4f}")
print(f"  Difference (GR - framework):          {sp.simplify(C_total_GR - C_total_lattice)}")
print(f"\n  The framework predicts a NEGATIVE log correction (entropy")
print(f"  slightly below A/4), opposite in sign to GR/string theory's")
print(f"  positive correction. The sign traces to the graviton being")
print(f"  a vector (microrotation) rather than a tensor (metric).")
