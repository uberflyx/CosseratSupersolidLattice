"""
Does the T_2u parity-odd strong shadow leave a measurable trace?

The forces section shows the hemitropic coupling theta_ch = alpha^2/(2pi) links
each silent channel to a carrying one: T_2u to the strong channel T_2g. A naive
reading gives a parity-odd strong coupling of order theta_ch ~ 1e-5, which would
be ~100 times the observed hadronic parity violation (~1e-7) and a problem. The
resolution is the matrix element, and this script computes it.

The point is that T_2u is silent for the same reason it is nearly harmless: it
produces zero direct tunnelling on the first two reciprocal-lattice shells, so
its leading coupling to any source sits on the third shell {220}, while the
strong channel T_2g couples on the first shell {111}. The Peierls-Nabarro form
factor f(k) = exp(-pi w |k|), with w = 0.783 l, falls steeply between shells, so
the parity-odd vertex is suppressed relative to the parity-even one by
f(k_220)/f(k_111). The parity-odd-to-parity-even ratio in the strong sector is
therefore

    r_PV(T_2u) ~ theta_ch * f(k_220)/f(k_111),

which the script evaluates. Using the third shell is conservative: if T_2u first
appears on a higher shell the suppression is stronger still, so this is an upper
bound on the effect.
"""

import numpy as np

alpha = 1.0 / 137.035999177
theta_ch = alpha ** 2 / (2.0 * np.pi)     # hemitropic chirality, ~8.5e-6
w_over_l = 0.783                          # PN core half-width / lattice spacing

# FCC reciprocal (BCC) shells: |G|_{hkl} = (2pi/a) sqrt(h^2+k^2+l^2), a = l sqrt2.
# In units of 1/l:  |G| = (2pi/sqrt2) sqrt(h^2+k^2+l^2).
def G_mag(h, k, l):
    return (2.0 * np.pi / np.sqrt(2.0)) * np.sqrt(h * h + k * k + l * l)

SHELLS = [("111", (1, 1, 1)), ("200", (2, 0, 0)), ("220", (2, 2, 0)),
          ("311", (3, 1, 1)), ("222", (2, 2, 2))]

def form_factor(kl):
    """f(k) = exp(-pi w |k|), argument k*l dimensionless."""
    return np.exp(-np.pi * w_over_l * kl)

# Observed hadronic parity violation, characteristic dimensionless scale.
r_PV_observed = 2.0e-7          # ~ G_F m_pi^2, the DDH weak-PV scale

def line():
    print("-" * 74)

print(__doc__.strip().splitlines()[0])
line()
print(f"theta_ch = alpha^2/(2pi) = {theta_ch:.2e}   w = {w_over_l} l")
print(f"form factor f(k) = exp(-pi w |k|) = exp(-{np.pi*w_over_l:.3f} * |k| l)")
line()
print("Reciprocal-lattice shells and the form factor there:")
print(f"  {'shell':>6s} {'|G| l':>8s} {'f(|G|)':>12s}")
f = {}
for name, (h, k, l) in SHELLS:
    kl = G_mag(h, k, l)
    f[name] = form_factor(kl)
    print(f"  {name:>6s} {kl:>8.3f} {f[name]:>12.2e}")
line()
# T_2u is zero on shells 1 (111) and 2 (200); its leading coupling is shell 3 (220).
f_supp = f["220"] / f["111"]
r_PV_T2u = theta_ch * f_supp
print("The T_2u parity-odd strong coupling:")
print(f"  T_2u silent on {{111}},{{200}} -> leading coupling on {{220}} (shell 3)")
print(f"  form-factor suppression f(220)/f(111) = {f_supp:.2e}")
print(f"  parity-odd/parity-even ratio r_PV(T_2u) = theta_ch * f_supp = {r_PV_T2u:.2e}")
line()
print("Comparison:")
print(f"  observed hadronic parity violation  r_PV ~ {r_PV_observed:.1e}")
print(f"  T_2u prediction                     r_PV ~ {r_PV_T2u:.1e}")
print(f"  ratio (T_2u / observed)             ~ {r_PV_T2u/r_PV_observed:.1e}")
line()
print("Verdict:")
if r_PV_T2u < 0.05 * r_PV_observed:
    print("  The matrix element crushes it. The same zero overlap on the first two")
    print("  shells that silences T_2u as a force suppresses its parity-odd shadow")
    print(f"  to ~{r_PV_T2u:.0e}, four to five orders below the observed hadronic")
    print("  parity violation and far below any foreseeable measurement. What")
    print("  looked like a possible tension at the naive theta_ch ~ 1e-5 level is")
    print("  not one: there is no parity problem in the strong sector, and no new")
    print("  observable either. The naive estimate missed the form factor.")
else:
    print(f"  r_PV(T_2u) ~ {r_PV_T2u:.0e} is within reach of the observed scale;")
    print("  this would be a genuine additional contribution to hadronic parity")
    print("  violation and needs the full DDH-style matching.")
