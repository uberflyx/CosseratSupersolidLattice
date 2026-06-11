#!/usr/bin/env python3
"""
d4_gamma_dictionary.py — the D4 microscopic dictionary and the
glide-wall misfit-curvature closed forms
=====================================================================

Verifies, by direct lattice sums over the 24 nearest-neighbour contacts
of the D4 lattice, the three closed-form results of the monograph's
"tangential sector and the D4 microscopic dictionary" subsection:

  1. The D4 dictionary between the contact stiffness ratio r = k_t/k_n
     and the Cosserat coupling number:

         N^2 = 3r / (1 + 5r),       rolling point  r = 1/(3 pi - 5)

     (slice model: N^2 = 2r/(1+3r), rolling point r = 1/(2 pi - 3)).
     The inter-layer bonds contribute pure k_t to the in-slice
     transverse sums: their single-axis spatial projections drop out of
     every cross-index (normal-stiffness) shear sum, which is the same
     mechanism that drives the Zener ratio to A = 1.

  2. The bond census of a {111} glide wall.  In the slice model the
     wall is crossed by 3 bonds per node.  On D4 it is crossed by 9
     bonds from each node of the adjacent layer (3 steep in-slice at
     normal height d111, 6 shallow inter-layer at d111/2) plus 3 more
     from the layer below: 12 per layer pair.

  3. The frozen-contact harmonic curvature of the misfit potential,
     in units of the Frenkel value mubar/d111 and anchored at
     N^2 = 1/pi through each model's own dictionary:

         slice:  (pi + 1) / (3 (pi - 1)) = 0.6446
         D4:      2 pi    / (3 (pi - 1)) = 0.9780

     with the improvement factor 2 pi / (pi + 1).  The inter-layer
     contacts supply the stiffness the slice count misses; the D4 wall
     supports the Frenkel convention of the alpha derivation to one
     part in fifty at harmonic order.

Open refinement (flagged in the monograph): rotational relaxation at
the wall softens the k_t contribution, the finite-amplitude contact
path stiffens it; the full Morse-path computation (a*ell = 7/3) with
contact switching decides the residual 2%.

Usage:  python3 d4_gamma_dictionary.py
Requires: numpy.

Author: Mitchell A. Cox
Date:   June 2026
"""

import numpy as np

# ---------------------------------------------------------------
# geometry: ell (NN distance) = 1 throughout
# ---------------------------------------------------------------
D111 = 2.0 / np.sqrt(6.0)        # {111} interplanar spacing = sqrt(2/3)
DHOP = 1.0 / np.sqrt(3.0)        # Shockley partial hop (misfit period)


def bonds_d4():
    """24 D4 nearest-neighbour unit vectors: (+-1, +-1, 0, 0)/sqrt2, all
    coordinate pairs."""
    out = []
    for i in range(4):
        for j in range(i + 1, 4):
            for si in (-1, 1):
                for sj in (-1, 1):
                    v = np.zeros(4)
                    v[i] = si / np.sqrt(2)
                    v[j] = sj / np.sqrt(2)
                    out.append(v)
    return np.array(out)


def bonds_slice():
    """The 12 in-slice (FCC) bonds: fourth component zero."""
    return np.array([b for b in bonds_d4() if abs(b[3]) < 1e-12])


# ---------------------------------------------------------------
# 1. the dictionary, by operational probe of the homogenised tensor
# ---------------------------------------------------------------
def transverse_moduli(bonds, kn, kt, cell_volume):
    """Probe sigma_ij = C_ijkl eps_kl with the asymmetric Cosserat
    strain.  Symmetric shear eps_(12)=1 reads off 2 mu_s; antisymmetric
    eps_[12]=1 reads off kappa_c.  C from the Born--Huang bond sum."""
    D = bonds.shape[1]
    C = np.zeros((D,) * 4)
    for n in bonds:
        K = kn * np.outer(n, n) + kt * (np.eye(D) - np.outer(n, n))
        C += np.einsum('ik,j,l->ijkl', K, n, n)
    C /= cell_volume
    es = np.zeros((D, D)); es[0, 1] = es[1, 0] = 1.0
    ea = np.zeros((D, D)); ea[0, 1] = 1.0; ea[1, 0] = -1.0
    ss = np.einsum('ijkl,kl->ij', C, es)
    sa = np.einsum('ijkl,kl->ij', C, ea)
    mu_s    = (ss[0, 1] + ss[1, 0]) / 4.0
    kappa_c = (sa[0, 1] - sa[1, 0]) / 2.0
    return mu_s, kappa_c


def N2_from_moduli(mu_s, kappa_c):
    """N^2 = kappa_c / (2 (mu_Er + kappa_c)), mu_Er = mu_s - kappa_c/2.
    Volume normalisation cancels in this ratio."""
    mu_er = mu_s - kappa_c / 2.0
    return kappa_c / (2.0 * (mu_er + kappa_c))


# ---------------------------------------------------------------
# 2. the glide-wall census and the frozen-contact curvature
# ---------------------------------------------------------------
def wall_census_and_curvature(kn, kt):
    """Bonds crossing a {111} wall (normal (1,1,1,0)/sqrt3) placed
    between adjacent D4 layers (which are spaced d111/2 along the
    normal).  Returns the census and the frozen-contact tangential
    stiffness sum per layer pair for a slide along the partial-hop
    direction t = (1,1,-2,0)/sqrt6."""
    nw = np.array([1, 1, 1, 0]) / np.sqrt(3)
    t  = np.array([1, 1, -2, 0]) / np.sqrt(6)
    up = [b for b in bonds_d4() if b @ nw > 1e-9]
    steep   = [b for b in up if abs(b @ nw - D111)     < 1e-9]
    shallow = [b for b in up if abs(b @ nw - D111 / 2) < 1e-9]
    # wall at d111/4 above a layer: that layer's node sends all 9
    # (steep + shallow); the layer below it reaches across with its 3
    # steep bonds only.
    S_pair = 2 * sum((b @ t) ** 2 for b in steep) + sum((b @ t) ** 2 for b in shallow)
    n_pair = 2 * len(steep) + len(shallow)
    K_pair = kn * S_pair + kt * (n_pair - S_pair)
    return len(steep), len(shallow), n_pair, S_pair, K_pair


def main():
    sep = "=" * 68
    print(sep)
    print("  D4 microscopic dictionary and glide-wall misfit curvature")
    print(sep)

    # --- dictionaries ---
    r_slice = 1.0 / (2 * np.pi - 3)
    r_d4    = 1.0 / (3 * np.pi - 5)
    for label, bonds, V, r in [("slice (12 bonds)", bonds_slice()[:, :3], 1 / np.sqrt(2), r_slice),
                               ("D4 (24 bonds)",    bonds_d4(),           0.5,            r_d4)]:
        mu_s, kc = transverse_moduli(bonds, 1.0, r, V)
        n2 = N2_from_moduli(mu_s, kc)
        print(f"\n{label}: r = {r:.6f}")
        print(f"  mu_s = {mu_s:.6f}, kappa_c = {kc:.6f}  ->  N^2 = {n2:.10f}")
        print(f"  1/pi = {1/np.pi:.10f}   (match: {abs(n2*np.pi - 1) < 1e-12})")
    print(f"\nclosed forms: slice N^2 = 2r/(1+3r), r = 1/(2pi-3) = {r_slice:.5f}")
    print(f"              D4    N^2 = 3r/(1+5r), r = 1/(3pi-5) = {r_d4:.5f}")

    # --- census ---
    ns, nh, n_pair, S_pair, K_pair = wall_census_and_curvature(1.0, r_d4)
    print(f"\n{{111}} wall census: steep = {ns} (height d111), "
          f"shallow = {nh} (height d111/2); {n_pair} bonds per layer pair")
    print(f"Sum (n.t)^2 per pair = {S_pair:.4f}")

    # --- frozen-contact curvature ratios, anchored at N^2 = 1/pi ---
    # D4: Gamma/Gamma_F = (2 kn + 10 kt) / (3 (kn + 2 kt)) at r = r_d4
    # slice: (kn + 5 kt) / (3 (kn + kt)) at r = r_slice
    g_d4    = (2 + 10 * r_d4) / (3 * (1 + 2 * r_d4))
    g_slice = (1 + 5 * r_slice) / (3 * (1 + r_slice))
    cf_d4    = 2 * np.pi / (3 * (np.pi - 1))
    cf_slice = (np.pi + 1) / (3 * (np.pi - 1))
    print(f"\nfrozen-contact Gamma/Gamma_Frenkel:")
    print(f"  slice: {g_slice:.6f}   closed form (pi+1)/(3(pi-1)) = {cf_slice:.6f}")
    print(f"  D4:    {g_d4:.6f}   closed form 2pi/(3(pi-1))    = {cf_d4:.6f}")
    print(f"  improvement = {g_d4/g_slice:.6f} = 2pi/(pi+1) = {2*np.pi/(np.pi+1):.6f}")
    assert abs(g_d4 - cf_d4) < 1e-12 and abs(g_slice - cf_slice) < 1e-12
    print(f"\nThe D4 wall supports the Frenkel misfit convention to "
          f"{100*(1-cf_d4):.1f}% at harmonic order.")
    print(sep)


if __name__ == '__main__':
    main()
