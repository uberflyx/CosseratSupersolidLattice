"""First harmonic of a Peierls potential in a lattice: the misfit's part and the strain energy's.

Companion to the monograph chapter "The fine structure constant" (the equilibrium core width).
The chapter identifies the coupling with the first Fourier harmonic of the misfit energy summed
over the lattice rows, divided by the continuum misfit energy, and evaluates it on the arctangent
core, where it equals exp(-2 pi w/d). This script asks four questions of that identification.

1. Rigid arctangent core, misfit and strain energy both summed on rows.
   A core of half-width w (rows at spacing d = 1) is displaced by X from a row. The elastic energy
   is a lattice quadratic form (1/2N) sum_theta A(theta)|FFT(Delta)|^2, so A(theta) is the stiffness
   per row of a slip wave with phase theta per row, theta in the Brillouin zone. The first
   harmonic of the misfit energy is +(Gamma w/pi) e^{-2 pi w} cos(2 pi X). The strain energy carries
   its own first harmonic, of opposite sign, to leading order
       -(1/pi) e^{-2 pi w} cos(2 pi X) int_0^pi A(theta) [ 1/(theta (2 pi - theta))
                                        - e^{-2 theta w}/(theta (2 pi + theta)) ] dtheta,
   because the kink's content beyond the zone boundary aliases into the zone and beats against its
   content inside it (by Poisson summation the sampled profile at phase theta is the sum of the
   continuum transform at theta + 2 pi m; the neighbouring pairs m = (-1, 0) and (0, 1) carry the
   first harmonic, and the next pairs are down by a further e^{-2 pi w}). The closed
   form is checked against direct lattice sums for three kernels: the continuum kernel |theta|/2 cut
   at the zone boundary (at Cauchy balance the ratio of the two harmonics tends to ln 2 for a wide
   core), |sin(theta/2)|, and the exact kernel of the rectangular lattice of part 2.

2. Exact Peierls ratio of an anti-plane screw in a rectangular nearest-neighbour lattice.
   Spacing d along the glide line and h across it, shear modulus mu, sinusoidal misfit on the bonds
   that cross the glide plane with curvature Gamma = mu/h (the classical Peierls-Nabarro choice,
   whose continuum solution is the arctangent with w = h/2). The two half-lattices are eliminated
   exactly, leaving a ring of rows with kernel A(theta) = (mu d/2h)(1 - lambda)/lambda, where
   lambda + 1/lambda = 2 + 4 (h/d)^2 sin^2(theta/2). The barrier is the energy difference between
   the kink centred on a row and between rows, each relaxed with its mirror symmetry imposed, for a
   kink-antikink pair on the ring. --check2d relaxes the full two-dimensional lattice instead and
   compares the profile.

3. Free minimiser of the chapter's continuum functional.
   Cosserat anti-plane kernel at N^2 = 1/pi with gamma = mu ell^2 (q = 1.87/ell), Frenkel misfit of
   curvature Gamma = mubar/d111, period d = ell/sqrt3. The profile is free on a fine periodic grid
   (kink-antikink pair). Reported: the first-harmonic ratio in the chapter's reading (the Fourier
   coefficient of the misfit energy density at k1 = 2 pi/d over its mean) and in the slip-density
   reading, its energy against the best arctangent pair on the same grid, the nearest complex
   singularity of the core from the decay rate of the slip density's Fourier transform, and the
   shape against the arctangent family (peak, half-width at half maximum, far tail). The same
   construction with the Cauchy kernel is the control: it must return the arctangent of width d111/2.

4. The same functional sampled on its own rows (spacing d).
   The misfit is read only at the rows and the strain energy is discretised in three standard ways:
   the kernel cut at the zone boundary (the continuum energy of the band-limited interpolant), the
   kernel evaluated at the lattice wavenumber (2/d) sin(kd/2), and the continuum energy of the slip
   interpolated linearly between rows. For each, the relaxed row-centred and between-row cores give
   the full barrier, and each core translated rigidly through its band-limited interpolant gives the
   misfit's harmonic in the chapter's reading. In the first scheme the strain energy of a translated
   core does not change, so the two fixed-core values bracket the full barrier.

Units: d = 1 in parts 1, 2 and 4 with mu = 1 (mubar = 1 in part 4); ell = 1 and mubar = 1 in part 3.
Runtime: parts 1-2 under a minute; part 3 about four minutes; part 4 about two; --check2d a minute more.
Requires numpy, scipy.
"""
import argparse

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize


# ---------------------------------------------------------------------------------------------
# Part 1 and 2: rows on a ring, kink-antikink pair with mirror symmetry imposed
# ---------------------------------------------------------------------------------------------
def symmetry_map(N, centre):
    """Map a base vector b to the full ring Delta_n = cst_n + sgn_n * b[idx_n].

    The kink sits at row 0 ('row', Delta_0 = 1/2) or between rows -1 and 0 ('bond'); the
    antikink sits half a ring away with the same centring. Mirror symmetry about the kink,
    Delta -> 1 - Delta, and the reflection taking the kink to the antikink fix the ring from
    a quarter of it. idx = -1 marks a held value."""
    h, q = N // 2, N // 4
    idx = np.full(N, -1)
    sgn = np.zeros(N)
    cst = np.zeros(N)
    if centre == 'row':
        for n in range(1, q + 1):
            idx[n] = idx[h - n] = n - 1
            sgn[n] = sgn[h - n] = 1.0
        cst[0] = cst[h] = 0.5
    else:
        for n in range(q):
            idx[n] = idx[h - 1 - n] = n
            sgn[n] = sgn[h - 1 - n] = 1.0
    for n in range(h, N):                    # Delta_{n + N/2} = 1 - Delta_n
        idx[n] = idx[n - h]
        sgn[n] = -sgn[n - h]
        cst[n] = 1.0 - cst[n - h]
    return idx, sgn, cst, q


def arctan_base(N, centre, w):
    """Arctangent profile 1/2 + arctan(x/w)/pi sampled on the base rows."""
    q = N // 4
    x = np.arange(1, q + 1) if centre == 'row' else np.arange(q) + 0.5
    return 0.5 + np.arctan(x / w) / np.pi


def full_ring(b, N, centre):
    idx, sgn, cst, _ = symmetry_map(N, centre)
    D = cst.copy()
    free = idx >= 0
    D[free] += sgn[free] * b[idx[free]]
    return D


def frenkel(Gamma):
    """gamma(Delta) = (Gamma/4 pi^2)(1 - cos 2 pi Delta), curvature Gamma at Delta = 0 (d = 1)."""
    c = Gamma / (4 * np.pi ** 2)
    return (lambda D: c * (1 - np.cos(2 * np.pi * D)),
            lambda D: c * 2 * np.pi * np.sin(2 * np.pi * D))


def elastic(A, D):
    """(1/2N) sum_theta A(theta) |FFT(Delta)|^2 and its gradient."""
    F = np.fft.fft(D)
    return 0.5 * np.sum(A * np.abs(F) ** 2) / len(D), np.real(np.fft.ifft(A * F))


def relax(A, gam, dgam, N, centre, b0, gtol=1e-13):
    idx, sgn, cst, nb = symmetry_map(N, centre)
    free = idx >= 0

    def f(b):
        D = cst.copy()
        D[free] += sgn[free] * b[idx[free]]
        Eel, g = elastic(A, D)
        g = g + dgam(D)
        gb = np.bincount(idx[free], weights=sgn[free] * g[free], minlength=nb)
        return Eel + np.sum(gam(D)), gb

    res = minimize(f, b0, jac=True, method='L-BFGS-B',
                   options={'maxiter': 400000, 'gtol': gtol, 'ftol': 1e-18, 'maxcor': 80})
    return res, full_ring(res.x, N, centre)


def rect_kernel(N, hd, mu=1.0):
    """Exact kernel per row of two rectangular half-lattices (spacing 1 along the line, hd
    across), joined only through the disregistry Delta = u(upper row) - u(lower row)."""
    th = 2 * np.pi * np.fft.fftfreq(N)
    s = 2 * hd * np.abs(np.sin(th / 2))
    lam = 1 + s * s / 2 - s * np.sqrt(1 + s * s / 4)
    return (mu / (2 * hd)) * (1 - lam) / lam


def rect_kernel_fn(hd, mu=1.0):
    def A(t):
        s = 2 * hd * np.sin(t / 2)
        lam = 1 + s * s / 2 - s * np.sqrt(1 + s * s / 4)
        return (mu / (2 * hd)) * (1 - lam) / lam
    return A


def elastic_harmonic_coefficient(Afun, w):
    """Closed form: E_el^(1) = coefficient * e^{-2 pi w} cos(2 pi X) for a rigid arctangent core."""
    f = lambda t: Afun(t) * (-1.0 / (t * (2 * np.pi - t)) + np.exp(-2 * t * w) / (t * (2 * np.pi + t)))
    return quad(f, 0, np.pi, limit=400)[0] / np.pi


def part1(N=2 ** 15):
    print('Part 1. Rigid arctangent core on rows: first harmonics of the misfit and strain energies')
    th = 2 * np.pi * np.fft.fftfreq(N)
    hd = np.sqrt(2.0)
    kernels = [('continuum |theta|/2, zone-cut', 0.5 * np.abs(th), lambda t: 0.5 * t),
               ('|sin(theta/2)|', np.abs(np.sin(th / 2)), lambda t: np.sin(t / 2)),
               ('rectangular lattice, h/d = sqrt2', rect_kernel(N, hd), rect_kernel_fn(hd))]
    gam, _ = frenkel(1.0)
    worst = 0.0
    for name, A, Af in kernels:
        for w in (0.6, 0.8, 1.0, 1.5):
            Eel, Emis = {}, {}
            for c in ('row', 'bond'):
                D = full_ring(arctan_base(N, c, w), N, c)
                Eel[c] = elastic(A, D)[0]
                Emis[c] = np.sum(gam(D))
            # two kinks switch centring, and cos(0) - cos(pi) = 2
            c_el = (Eel['row'] - Eel['bond']) / 4
            c_mis = (Emis['row'] - Emis['bond']) / 4
            e = np.exp(-2 * np.pi * w)
            pred_el = elastic_harmonic_coefficient(Af, w) * e
            pred_mis = 1.0 * w / np.pi * e
            worst = max(worst, abs(c_el / pred_el - 1), abs(c_mis / pred_mis - 1))
            print('  %-34s w=%.2f  misfit %+.5e (closed %+.5e)  strain %+.5e (closed %+.5e)'
                  '  strain/misfit at Cauchy balance %.4f'
                  % (name, w, c_mis, pred_mis, c_el, pred_el, -2 * np.pi * elastic_harmonic_coefficient(Af, w)))
    print('  largest relative miss of the closed forms: %.1e (next alias pairs and the third harmonic,'
          ' since the row-minus-between-row difference cancels even harmonics)' % worst)
    big = -2 * np.pi * elastic_harmonic_coefficient(lambda t: 0.5 * t, 50.0)
    print('  wide-core limit, continuum kernel: strain/misfit = %.5f, ln 2 = %.5f' % (big, np.log(2)))


def part2(hd=np.sqrt(2.0), check2d=False):
    print('\nPart 2. Exact Peierls ratio, rectangular anti-plane lattice, h/d = %.4f' % hd)
    gam, dgam = frenkel(1.0 / hd)                       # Gamma = mu/h
    w0 = hd / 2                                         # continuum Peierls-Nabarro half-width
    for N in (1024, 4096, 16384):
        A = rect_kernel(N, hd)
        out = {}
        for c in ('row', 'bond'):
            res, D = relax(A, gam, dgam, N, c, arctan_base(N, c, w0))
            out[c] = (res.fun, np.sum(gam(D)), res.success, D)
        dE = 0.5 * (out['row'][0] - out['bond'][0])     # barrier of one kink
        Emis = 0.25 * (out['row'][1] + out['bond'][1])  # mean misfit energy of one kink
        alpha = dE / (4 * Emis)
        print('  N=%6d  barrier %.6e  misfit energy %.6f (continuum %.6f)  ratio %.6e  W = %.5f'
              '  [continuum e^{-2 pi w} = %.6e, w/d = %.5f]  converged %s'
              % (N, dE, Emis, w0 / hd / (2 * np.pi), alpha, -np.log(alpha) / (2 * np.pi),
                 np.exp(-2 * np.pi * w0), w0, out['row'][2] and out['bond'][2]))
    print('  misfit energy as an arctangent half-width: w/d = %.4f' % (2 * np.pi * Emis * hd))
    D = out['bond'][3]
    print('  between-row kink, rows -2..1:', np.round(np.r_[D[-2:], D[:2]], 6))
    # rigid arctangent of the continuum width: the two harmonics separately
    for w in (w0, 2 * np.pi * Emis * hd):
        Ee, Em = {}, {}
        for c in ('row', 'bond'):
            Dr = full_ring(arctan_base(N, c, w), N, c)
            Ee[c] = elastic(A, Dr)[0]
            Em[c] = np.sum(gam(Dr))
        print('  rigid arctangent w/d = %.4f: misfit harmonic %+.4e, strain harmonic %+.4e'
              % (w, (Em['row'] - Em['bond']) / 4, (Ee['row'] - Ee['bond']) / 4))
    if check2d:
        check_two_dimensional(hd, out['bond'][3])


def check_two_dimensional(hd, D1, N=256, M=128):
    """Relax the full lattice (rows -M..M-1, gap bonds between rows -1 and 0 sinusoidal,
    free outer rows) for the between-row kink and compare the disregistry with the reduction."""
    kx, kz = hd, 1.0 / hd                                # mu h/d and mu d/h with mu = d = 1
    gam, dgam = frenkel(1.0 / hd)

    def energy(U):
        U = U.reshape(2 * M, N)
        E = 0.5 * kx * np.sum((np.roll(U, -1, axis=1) - U) ** 2)
        dz = U[1:] - U[:-1]
        mask = np.ones(2 * M - 1, bool)
        mask[M - 1] = False
        E += 0.5 * kz * np.sum(dz[mask] ** 2)
        D = U[M] - U[M - 1]
        E += np.sum(gam(D))
        g = np.zeros_like(U)
        gx = kx * (U - np.roll(U, -1, axis=1))
        g += gx - np.roll(gx, 1, axis=1)
        t = kz * dz
        t[M - 1] = 0.0
        g[1:] += t
        g[:-1] -= t
        gd = dgam(D)
        g[M] += gd
        g[M - 1] -= gd
        return E, g.ravel()

    D0 = full_ring(arctan_base(N, 'bond', 1.0), N, 'bond')
    U0 = np.zeros((2 * M, N))
    U0[M:] = D0 / 2
    U0[:M] = -D0 / 2
    res = minimize(energy, U0.ravel(), jac=True, method='L-BFGS-B',
                   options={'maxiter': 200000, 'gtol': 1e-11, 'ftol': 1e-18, 'maxcor': 50})
    U = res.x.reshape(2 * M, N)
    D = U[M] - U[M - 1]
    print('  two-dimensional check (N=%d rows, %d layers per side): rows -2..1' % (N, M),
          np.round(np.r_[D[-2:], D[:2]], 6), ' reduction', np.round(np.r_[D1[-2:], D1[:2]], 6))


# ---------------------------------------------------------------------------------------------
# Part 3: free minimiser of the continuum functional
# ---------------------------------------------------------------------------------------------
D111 = np.sqrt(2.0 / 3.0)
DHOP = 1.0 / np.sqrt(3.0)
N2 = 1.0 / np.pi
Q = np.sqrt(2.0 / ((np.pi - 2.0) / 2.0))                 # q^2 = 2 kappa_c/gamma with gamma = mu ell^2


def cosserat_kernel(k):
    """Glide-line kernel in units of mubar: (1/2)[|k| + (N^2/(1 - N^2)) k^2/sqrt(k^2 + q^2)]."""
    k = np.abs(k)
    return 0.5 * (k + N2 / (1 - N2) * k * k / np.sqrt(k * k + Q * Q))


def cauchy_kernel(k):
    return 0.5 * np.abs(k)


def pair_energy(N, dx, kernel):
    """Energy of a kink-antikink pair on the periodic grid (continuum functional, Frenkel misfit)
    and its gradient, with the grid, the kernel on it and the misfit amplitude."""
    Lb = N * dx
    x = (np.arange(N) - N // 4) * dx
    k = 2 * np.pi * np.fft.fftfreq(N, dx)
    Kk = kernel(k)
    Gamma = 1.0 / D111
    A = Gamma / (2 * np.pi / DHOP) ** 2

    def f(phi):
        ph = np.fft.fft(phi) * dx
        Eel = 0.5 * np.sum(Kk * np.abs(ph) ** 2) / Lb
        gel = np.real(np.fft.ifft(Kk * ph)) * N / Lb
        s = 2 * np.pi * phi / DHOP
        return (Eel + np.sum(A * (1 - np.cos(s))) * dx,
                (gel + A * (2 * np.pi / DHOP) * np.sin(s)) * dx)
    return f, x, Gamma, A


def arctan_pair(x, Lb, w):
    xa = x - Lb / 2
    return np.where(np.abs(x) < np.abs(xa), DHOP * (0.5 + np.arctan(x / w) / np.pi),
                    DHOP * (0.5 - np.arctan(xa / w) / np.pi))


def free_profile(N=2 ** 16, dx=0.01, w0=0.45, kernel=cosserat_kernel):
    """Kink-antikink pair on a periodic grid, profile free; returns x, the slip, convergence flag,
    Gamma, the misfit amplitude and the relaxed energy."""
    f, x, Gamma, A = pair_energy(N, dx, kernel)
    res = minimize(f, arctan_pair(x, N * dx, w0), jac=True, method='L-BFGS-B',
                   options={'maxiter': 80000, 'gtol': 1e-14, 'ftol': 1e-18, 'maxcor': 50})
    return x, res.x, res.success, Gamma, A, res.fun


def harmonics(x, phi, dx, A, Gamma):
    """First-harmonic exponents W in the misfit-energy-density reading of the chapter and in the
    slip-density reading, and the misfit energy read as an arctangent width (all over d)."""
    Lb = len(x) * dx
    near = np.abs(x) < Lb / 4
    k1 = 2 * np.pi / DHOP
    g = A * (1 - np.cos(2 * np.pi * phi / DHOP))
    g0 = np.sum(g[near]) * dx
    g1 = abs(np.sum(g[near] * np.exp(-1j * k1 * x[near])) * dx)
    r1 = abs(np.sum(np.gradient(phi, dx)[near] * np.exp(-1j * k1 * x[near])) * dx)
    return (-np.log(g1 / g0) / (2 * np.pi), -np.log(r1 / DHOP) / (2 * np.pi),
            g0 / (Gamma * DHOP ** 2 / (2 * np.pi)) / DHOP)


def part3():
    print('\nPart 3. Free minimiser of the continuum functional (N^2 = 1/pi, gamma = mu ell^2, Frenkel)')
    w_family = 0.83664                                   # arctangent-family equilibrium (chapter)
    for N, dx in ((2 ** 15, 0.02), (2 ** 16, 0.01)):
        x, phi, ok, Gamma, A, E = free_profile(N, dx)
        Wm, Ws, wE = harmonics(x, phi, dx, A, Gamma)
        print('  grid %6d x %.3f: W (misfit-density reading) = %.5f   W (slip-density reading) = %.5f'
              '   misfit energy as arctangent width = %.5f   converged %s' % (N, dx, Wm, Ws, wE, ok))
    # energy against the best arctangent pair on the same grid
    f, xg, _, _ = pair_energy(N, dx, cosserat_kernel)
    from scipy.optimize import minimize_scalar
    best = minimize_scalar(lambda w: f(arctan_pair(xg, N * dx, w))[0], bounds=(0.3, 0.7), method='bounded',
                           options={'xatol': 1e-8})
    print('  free minimiser lies %.2e per kink below the best arctangent pair (w/d = %.4f),'
          ' %.2f per cent of the misfit energy per kink'
          % ((best.fun - E) / 2, best.x / DHOP, 100 * (best.fun - E) / 2 / (wE * Gamma * DHOP ** 2 / (2 * np.pi))))
    # singularity: local decay rate of |rho_hat(k)| between k = 4 and 24 (units 1/ell)
    Lb = N * dx
    near = np.abs(x) < Lb / 4
    rho_hat = np.fft.fft(np.where(near, np.gradient(phi, dx), 0.0)) * dx
    kk = 2 * np.pi * np.fft.fftfreq(N, dx)
    sel = (kk > 4) & (kk < 24)
    slope, icpt = np.polyfit(kk[sel], np.log(np.abs(rho_hat[sel])), 1)
    print('  slip density: |rho_hat(k)| = R exp(-k zeta) with zeta/d = %.4f and R/d = %.4f'
          % (-slope / DHOP, np.exp(icpt) / DHOP))
    rho = np.gradient(phi, dx) / DHOP
    i0 = np.argmin(np.abs(x))
    j = i0
    while rho[j] > rho[i0] / 2:
        j += 1
    w = w_family * DHOP
    print('  shape: peak rho(0) d = %.4f (family %.4f); half-width at half maximum %.4f d (family %.4f d)'
          % (rho[i0] * DHOP, DHOP / (np.pi * w), x[j] / DHOP, w_family))
    for xx in (20, 40):
        i = np.argmin(np.abs(x - xx * DHOP))
        print('  far tail x^2 rho/d at x = %2d d: %.4f (family %.4f; Cauchy core of width d111/2: %.4f)'
              % (xx, x[i] ** 2 * rho[i] / DHOP, w / np.pi / DHOP, (D111 / 2) / np.pi / DHOP))
    # control: the same construction with the Cauchy kernel returns the arctangent of width d111/2
    x, phi, ok, Gamma, A, E = free_profile(2 ** 16, 0.01, kernel=cauchy_kernel)
    Wm, Ws, wE = harmonics(x, phi, 0.01, A, Gamma)
    print('  Cauchy control: W (misfit-density) = %.5f, W (slip-density) = %.5f, 1/sqrt2 = %.5f'
          % (Wm, Ws, 1 / np.sqrt(2)))


def misfit_of_translated_core(D, gam, X):
    """Misfit row sum of a fixed core translated rigidly by X rows, through its band-limited
    interpolant (the Nyquist component kept real)."""
    N = len(D)
    th = 2 * np.pi * np.fft.fftfreq(N)
    ph = np.exp(-1j * th * X)
    ph[N // 2] = np.cos(np.pi * X)
    return np.sum(gam(np.real(np.fft.ifft(np.fft.fft(D) * ph))))


def part4(N=4096):
    print('\nPart 4. The continuum functional sampled on its own rows (spacing d), three discretisations')
    th = 2 * np.pi * np.fft.fftfreq(N)
    gam, dgam = frenkel(1.0 / D111)                     # misfit per row in units of d (Delta/d)

    def scheme(name, K):
        if name == 'kernel cut at the zone boundary':
            return K(th / DHOP)
        if name == 'kernel at the lattice wavenumber':
            return K((2 / DHOP) * np.abs(np.sin(th / 2)))
        A = np.zeros_like(th)                           # slip interpolated linearly between rows
        for m in range(-300, 301):
            t = th + 2 * np.pi * m
            A += K(t / DHOP) * np.sinc(t / (2 * np.pi)) ** 4
        return A

    W = lambda a: -np.log(abs(a)) / (2 * np.pi)
    for kname, K in (('Cauchy', cauchy_kernel), ('Cosserat', cosserat_kernel)):
        for name in ('kernel cut at the zone boundary', 'kernel at the lattice wavenumber',
                     'slip linear between rows'):
            A = scheme(name, K)
            out = {}
            for c in ('row', 'bond'):
                res, D = relax(A, gam, dgam, N, c, arctan_base(N, c, 0.8))
                out[c] = (res.fun, D)
            Emean = 0.25 * (np.sum(gam(out['row'][1])) + np.sum(gam(out['bond'][1])))
            full = 0.5 * (out['row'][0] - out['bond'][0]) / (4 * Emean)
            fixed = []
            for c in ('row', 'bond'):
                e0 = misfit_of_translated_core(out[c][1], gam, 0.0)
                e1 = misfit_of_translated_core(out[c][1], gam, 0.5)
                fixed.append(abs(e0 - e1) / (2 * (e0 + e1)))
            print('  %-8s %-34s full barrier W = %.4f; misfit harmonic of the fixed core, row-centred'
                  ' W = %.4f, between rows W = %.4f' % (kname, name, W(full), W(fixed[0]), W(fixed[1])))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--check2d', action='store_true', help='relax the two-dimensional lattice too')
    ap.add_argument('--skip-free', action='store_true', help='skip part 3')
    a = ap.parse_args()
    part1()
    part2(check2d=a.check2d)
    if not a.skip_free:
        part3()
    part4()
