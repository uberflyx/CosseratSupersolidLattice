#!/usr/bin/env python3
"""The mass-scale coefficient in closed form, and the width it demands.

Established this session (sympy-verified, see transcript of derivation):
1. hbar-cancellation: linear coupling to the core libration gives the
   classical static response g^2/(2K); no quantisation convention enters.
2. Exact misfit integral: Int eps_F^2 dx = b^2 w/(2 pi d^2), with
   eps_F = (b/2 pi d) sin(2 pi phi/b) the Frenkel interplanar strain on the
   PN slip profile.
3. The coefficient of theta_ch^2 m_0 in m_1, with the microrotation gradient
   at the Cosserat length l_c = l/2 (exact) and the two bounding node sheets:

       coeff(w) = (6 w / ell) (1 - 2N^2)/(pi N^2)

   One line. No cutoff, no hbar, no fitted constant.

Consequence: both measured splittings ride w^2, so the data fix the edge
core width to 0.35 per cent:

       w_edge(data) = 0.4585 +/- 0.0016 ell
       w_star(unity) = pi N^2 ell / (6(1-2N^2)) = 0.45866 ell  (closed form)

Candidate widths and their verdicts are printed below. The screw-precedent
width misses at -2.5 sigma on each splitting; the brane width L_b is
excluded outright. NOT adopted: any choice among them. The open input is
the edge's own Cosserat Peierls-Nabarro equilibrium width, derivable by
the same machinery that gave the screw w/d = pi/4 and alpha to ten digits.

Exposure list (recorded so nothing moves quietly): the sheet count
lambda_L = 2 could be argued 1 (relative mode of the bounding pair), and
the local stiffness K = kappa_c ell^3 comes from (9-3)/12 contact counting;
each has a physical argument, each moves the required width linearly.
"""
import numpy as np

N2 = 1/np.pi
d111 = np.sqrt(2/3)
def coeff(w): return 6*w*(1-2*N2)/(np.pi*N2)

D21P, D21E, D21S = 74.711, 75.37, (0.94, 1.00)   # meV^2 (hi, lo sigma)
D31P, D31E, D31S = 2523.4, 2511.0, (21.0, 20.0)

def pull(p, e, s):
    return (p-e)/(s[0] if p > e else s[1])

def main():
    wstar = np.pi*N2/(6*(1-2*N2))
    print(f"coeff(w) = 6 w (1-2N^2)/(pi N^2 ell);  w_star(unity) = {wstar:.5f} ell")
    print(f"{'w/ell':>8s} {'provenance':>24s} {'coeff':>8s} {'Dm21':>8s} {'pull':>6s} {'Dm31':>8s} {'pull':>6s}")
    for w, lab in [(np.pi/(4*np.sqrt(3)), "screw pi d_p/4"),
                   (wstar,                "unity, closed form"),
                   (0.589,                "brane L_b")]:
        c = coeff(w); s = c*c
        print(f"{w:8.4f} {lab:>24s} {c:8.4f} {D21P*s:8.2f} {pull(D21P*s,D21E,D21S):+6.1f}"
              f" {D31P*s:8.1f} {pull(D31P*s,D31E,D31S):+6.1f}")
    # data-preferred width: splittings ~ w^2, joint fit
    ws = np.linspace(0.44, 0.48, 40001)
    chi = [(pull(D21P*coeff(w)**2, D21E, D21S)**2 +
            pull(D31P*coeff(w)**2, D31E, D31S)**2) for w in ws]
    chi = np.array(chi); i = chi.argmin()
    band = ws[chi <= chi[i]+1]
    print(f"\ndata-preferred width: w = {ws[i]:.4f} (+{band.max()-ws[i]:.4f}/-{ws[i]-band.min():.4f}) ell")
    print("the open calculation: the edge's Cosserat PN equilibrium width.")

if __name__ == "__main__":
    main()
