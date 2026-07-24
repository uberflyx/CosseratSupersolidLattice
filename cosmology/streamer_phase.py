#!/usr/bin/env python3
"""CMB acoustic-peak phase coefficient for a free-streamer of arbitrary speed.

Monograph: superfluid sector, "The phase coefficient at arbitrary streaming
speed" (eq:streamer_hamiltonian .. eq:streamer_phase_result, fig:streamer_phase).

Generalises the Bashinsky-Seljak / Baumann-Green-Meyers-Wallisch first-order
computation to streaming speed v. The acoustic-metric Hamiltonian of a
linear-dispersion quantum is w = v q (1 + Phi_+), exact at linear order for
any v, so the free-streaming kernels generalise by k*tau -> v*k*tau alone:
  D2(y) = -3 zeta j2(s y) + 3 s Int_0^y Phi_+(y') [ (2/5) j1(s(y-y'))
          - (3/5) j3(s(y-y')) ] dy',    s = v/c_gamma.
The kernel's mean vanishes exactly ((2/5)*1 - (3/5)*(2/3) = 0), which is why
the coefficient collapses at large s: the streamer's drive averages out over
its own crossing time. Anchor: f(c) = 0.602 vs literature 0.600.
Headline: f(3sqrt2 c) = 0.062, suppression 9.7, dl = 0.03 at dNeff = 0.029.
"""
import numpy as np
from scipy.special import spherical_jn
from scipy.signal import fftconvolve

C_GAMMA = 1.0/np.sqrt(3.0)     # plasma sound speed / c
A_NU = (8.0/7.0)*(11.0/4.0)**(4.0/3.0)


def B_at(s, Y, dy=0.0015):
    """Out-of-phase amplitude B/(zeta*eps) with integration cutoff Y."""
    N = int(Y/dy)
    y = np.linspace(1e-6, Y, N); dy = y[1] - y[0]
    Phi0 = 4.0*(np.sin(y) - y*np.cos(y))/y**3          # zeroth-order potential
    g = 0.4*spherical_jn(1, s*(y - y[0])) - 0.6*spherical_jn(3, s*(y - y[0]))
    conv = fftconvolve(Phi0, g)[:N]*dy                 # streaming convolution
    Phim = 4.0*spherical_jn(2, s*y)/y**2 - 4.0*s*conv/y**2
    d1 = np.gradient(Phim, dy, edge_order=2)
    d2 = np.gradient(d1, dy, edge_order=2)
    S = d2 + 2.0*d1/y + 3.0*Phim                       # source S~[Phi_-]
    # Phi_+^(1) via the radiation-era Green's function, trig-expanded cumulants
    A1 = np.cumsum(S*y*np.cos(y))*dy;   A2 = np.cumsum(S*y*np.sin(y))*dy
    A3 = np.cumsum(S*y*y*np.cos(y))*dy; A4 = np.cumsum(S*y*y*np.sin(y))*dy
    c, sn = np.cos(y), np.sin(y)
    Phi1 = (c*(A3 - y*A1 - A2 - y*A4) + sn*(A4 - y*A2 + A1 + y*A3))/y**3
    Bc = np.cumsum(Phi1*np.cos(y))*dy
    m = y > (Y - 20*np.pi)                              # period-averaged tail
    return float(np.mean(Bc[m]))


def f_of_v(v):
    """Phase coefficient f(v): Richardson-extrapolated in the 1/Y cutoff error
    (the 1/Y scaling is verified over three cutoff doublings)."""
    s = v/C_GAMMA
    return 2.0*B_at(s, 6400.0) - B_at(s, 3200.0)


if __name__ == "__main__":
    f_nu = f_of_v(1.0)
    print(f"anchor  f(c)     = {f_nu:.4f}   (literature 0.600)")
    f2 = f_of_v(18.0**0.5)
    print(f"result  f(3sqrt2 c)= {f2:.4f}   suppression {f_nu/f2:.2f}")
    dneff = (4/7)*(11/4)**(4/3)*(1/18)**1.5   # v2 = 3 sqrt2 c (C11 = 3 mu, f_n = 1/6)
    eps2 = (dneff/A_NU)/(1.0 + (3.046 + dneff)/A_NU)
    th2 = f2*eps2
    print(f"eps_2 = {eps2:.5f}  theta_2 = {th2:.3e} rad  dl = {th2/np.pi*330:.3f}")
    print(f"phase-active fraction = {f2/f_nu:.3f}  ->  (dNeff, Nfluid) = "
          f"({f2/f_nu*dneff:.4f}, {(1-f2/f_nu)*dneff:.4f})")
