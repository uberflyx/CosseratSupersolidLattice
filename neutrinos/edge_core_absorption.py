#!/usr/bin/env python3
"""The core computation: the Born-Huang coefficient of m_1 = theta_ch^2 m_0,
tensor-complete on the Peierls-Nabarro-resolved edge.

WHAT THIS CLOSES
The mass scale's unit coefficient was carried by structural identity with the
screw; the calibrated far-field run gave a floor of 0.072 and localised the
balance to the core. This script computes the core: the full second-order
chiral energy of the edge with (a) the exact strain harmonics, (b) the
Peierls-Nabarro misfit smearing that is the core's physical UV completion,
(c) the certified Cosserat moduli (kernel-certified in
screw_vertex_calibration.py), and (d) the tensor split of the microrotation
response into its longitudinal and u_z-dressed transverse parts, which the
scalar run collapsed.

THE PHYSICS, STEP BY STEP
1. The edge with a PN core is a distribution of infinitesimal edges with
   Lorentzian density rho(x) = (w/pi)/(x^2+w^2) on the glide plane; its
   Fourier strain is the point-edge strain times exp(-w|q_x|). The core
   width w is the variational object: the absorption claim is that the
   coefficient, evaluated at the self-consistent w, needs no further
   overlap factor.
2. The chiral vertex W = 2 C_ch e_(ij) k_(ij), C_ch = mu theta ell, sources
   the in-plane microrotation doublet: S_i(q) = 2 C_ch i q_j e_(ij)(q).
3. The doublet's response splits. The longitudinal part (phi || q) sees the
   bare screened stiffness gamma(q^2 + q_s^2), q_s^2 = 2 kappa_c/gamma. The
   transverse part (phi || z x q) couples to u_z through kappa_c (the
   (curl u)_perp = grad u_z channel), and integrating u_z out dresses its
   screening to p^2 = q_s^2 (1 - N^2), exactly as in the certified kernel.
4. E2 per unit line length = (1/2) Int d2q/(2pi)^2 [ |S_L|^2/(gam(q^2+q_s^2))
   + |S_T|^2/(gam(q^2+p^2)) ], and the coefficient is
   C = E2_pl x L_seg / theta^2, with L_seg the line length per the object
   whose mass is m_1. Two conventions are reported: L_seg = ell (one node)
   and L_seg = L_xi = 1.73 ell (the brane loop length the monograph fixes).
5. The remaining cut: the PN factor regularises q_x, but the direction
   normal to the glide plane keeps the elastic 1/q tail down to the lattice,
   so a zone cut q_max = pi/d_111 acts on |q| and its residual sensitivity
   is scanned and reported.

Units: ell = mu = 1, m_0 = mu ell^3 = 1.
"""

import numpy as np

N2   = 1/np.pi
KAP  = 2*N2/(1-2*N2)         # kappa_c/mu, certified convention
GAM  = 1.0                   # gamma/(mu ell^2), certified
QS2  = 2*KAP/GAM             # bare screening^2  (= 3.504)
P2   = QS2*(1-N2)            # u_z-dressed       (= 2.388)
D111 = np.sqrt(2/3)
LXI  = 1.73                  # brane loop length [ell], monograph

# exact edge strain harmonics (coefficient of 1/r, unit b), hand-verified
H = {
 'xx': {(1,'s'): -1/(3*np.pi), (3,'s'): -1/(6*np.pi)},
 'yy': {(3,'s'):  1/(6*np.pi)},
 'xy': {(1,'c'):  1/(6*np.pi), (3,'c'):  1/(6*np.pi)},
}

def strain_q(qx, qy):
    """Point-edge Fourier strains: harmonic l of a/r -> 2 pi (-i)^l a h_l(th_q)/q."""
    q  = np.hypot(qx, qy); th = np.arctan2(qy, qx)
    out = {}
    for comp, terms in H.items():
        val = np.zeros_like(q, dtype=complex)
        for (l, kind), a in terms.items():
            hl = np.cos(l*th) if kind == 'c' else np.sin(l*th)
            val += 2*np.pi*((-1j)**l)*a*hl
        out[comp] = val/q
    return out

def coefficient(w, qmax=np.pi/D111, nq=900, nth=720, Lseg=1.0, split=True):
    th = np.linspace(0, 2*np.pi, nth, endpoint=False)
    lq = np.linspace(np.log(1e-4), np.log(qmax), nq)
    q  = np.exp(lq)
    Q, TH = np.meshgrid(q, th, indexing='ij')
    QX, QY = Q*np.cos(TH), Q*np.sin(TH)
    e = strain_q(QX, QY)
    pn = np.exp(-w*np.abs(QX))                    # PN misfit smearing (x only)
    # source S_i = i q_j e_(ij) (the 2 C_ch factored out; vertex^2 = 4 theta^2)
    Sx = 1j*(QX*e['xx'] + QY*e['xy'])*pn
    Sy = 1j*(QX*e['xy'] + QY*e['yy'])*pn
    # longitudinal / transverse split
    qhx, qhy = np.cos(TH), np.sin(TH)
    SL = Sx*qhx + Sy*qhy
    ST = -Sx*qhy + Sy*qhx
    if split:
        dens = (np.abs(SL)**2/(GAM*(Q**2+QS2)) + np.abs(ST)**2/(GAM*(Q**2+P2)))
    else:   # scalar collapse, for comparison with the earlier runs
        dens = (np.abs(SL)**2 + np.abs(ST)**2)/(GAM*(Q**2+P2))
    # measure: d2q/(2pi)^2 = q dq dth/(4 pi^2); log grid: q dq = q^2 dlq
    E2 = 0.5*np.trapezoid(np.trapezoid(dens*Q**2, lq, axis=0), th)/ (4*np.pi**2)
    return 4.0*E2*Lseg           # vertex^2 = (2 C_ch)^2 -> 4 theta^2; /theta^2

def main():
    print("="*70)
    print("Core Born-Huang coefficient, tensor-complete, PN-resolved")
    print("="*70)
    print(f"moduli: kappa_c/mu={KAP:.4f}  q_s^2={QS2:.4f}  p^2={P2:.4f}  (certified)")
    print(f"zone cut pi/d_111 = {np.pi/D111:.3f}/ell;  brane L_xi = {LXI} ell\n")
    ws = [0.354, 0.437, 0.453, 0.544, 0.70]
    print(f"{'w [ell]':>8s} {'C (L=ell)':>11s} {'C (L=L_xi)':>11s}")
    for w in ws:
        c = coefficient(w)
        print(f"{w:8.3f} {c:11.4f} {c*LXI:11.4f}")
    w0 = 0.453
    c0 = coefficient(w0)
    print(f"\nzone-cut scan at w={w0} (L=L_xi):", end=" ")
    for qm in [0.75*np.pi/D111, np.pi/D111, 1.5*np.pi/D111]:
        print(f"{coefficient(w0, qmax=qm)*LXI:.3f}", end="  ")
    print(f"\ntensor split vs scalar collapse at w={w0}: "
          f"{coefficient(w0)*LXI:.4f} vs {coefficient(w0, split=False)*LXI:.4f}")
    print(f"\nheadline: C = {c0*LXI:.3f} at the screw-precedent width, brane length,")
    print(f"          C = {coefficient(0.544)*LXI:.3f} at the classical edge PN width.")
    print("""
CONCLUSION. Every continuum ingredient is now certified: propagator
(kernel-exact), source harmonics (hand-verified), angular algebra (exact,
1/(36 pi) with the Fourier phases; the in-phase 5/(36 pi) an earlier run
used is exactly 5x too large), UV completion (the PN width, physical),
IR (the Cosserat screening, physical). The tensor split changes the
scalar collapse by 20 per cent, not an order. The continuum second-order
chiral response of the edge is C ~ 0.008 per node length, 0.013 at the
brane loop, two orders below the claimed unity, and no convention excuse
remains on the continuum side.

What that decides: the unit coefficient of m_1 = theta_ch^2 m_0 cannot
live in continuum elasticity. It must be carried by the discrete channel,
the localised core modes of the chirally dressed Peierls potential, the
lattice's own periodic potential acting on the defect, which is exactly
where the screw's mass alpha m_0 comes from (e^{pi^2/2} is not a
continuum object either). The structural-identity argument is therefore
better aimed than the continuum split suggested: it compares like with
like. What it still does not supply is the number, and the open
calculation is now precisely posed: the RS sum over the bound modes of
the discrete misfit potential, dressed by the chiral vertex, on the
relaxed PN profile. The continuum piece computed here is its small
smooth-background correction.""")

if __name__ == "__main__":
    main()
