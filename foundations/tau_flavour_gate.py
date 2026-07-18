"""Tau elastic stability: the flavour gate, quantified.

The generations are ring (Koide Z_3) states over the three compact stacking
registries; flavour change is Delta n != 0 mod 3. Wound elastic quanta
(k_4 = +/-1, bare shear rung E_1 = 3 m_0 c^2 / sqrt(2) = 148.546 MeV) carry
exactly that charge, so elastic de-excitation ell -> e + X(k_4) is
Z_3-allowed. Kinematics then splits the leptons:

  electron: lightest, stable trivially.
  muon:     m_mu - m_e = 105.147 < E_1 = 148.546 MeV. CLOSED. The muon is
            exactly elastically stable by kinematics, and mu -> e gamma is
            forbidden by compact crystal-momentum conservation (photon is
            k_4 = 0). Margin 43.4 MeV; curiously E_1/(m_mu - m_e) = sqrt(2)
            to 0.10 percent (equivalently m_mu - m_e vs 3 m_0/2), an
            observation logged but not used.
  tau:      m_tau - m_e = 1776.3 MeV >> E_1. OPEN in naive harmonic
            elasticity, and dangerously so: a Holstein-type analysis shows a
            wound phonon modulating the three registry energies with phases
            1, w, w^2 connects ring Bloch states n -> n +/- 1 at elastic
            strength WITHOUT the tunnelling overlap. Ungated width scale
            ~ alpha m_tau ~ 13 MeV; observed total tau width 2.27 meV:
            excluded by ten orders of magnitude.

THE GATE (the requirement this makes sharp). The tau's narrowness requires
that every Delta n != 0 coupling to bulk radiation carry the chiral
(hemitropic) insertion theta_ch = alpha^2/(2 pi) per amplitude, i.e. that
no free Delta-n-carrying quantum couples at elastic strength. The framework
supplies this structure independently three times: (i) it derives the same
gate for neutrinos, Gamma(nu_j -> nu_i gamma) ~ theta_ch^4 m^5 ~ G_F^2 m^5;
(ii) it places the wound sector at the hadronic stiffness gap whose cold
asymptotic states are composites (confinement read as temporal binding), so
one-quantum wound states are not asymptotic and Delta n transfer requires
composite/defect-pair emission, which is the weak channel (the observed tau
decays: tau -> pi nu, ell nu nu); (iii) it makes G_F^2 proportional to
theta_ch^4 for the whole weak sector. With the gate, the elastic-flavour
width scale is G_F^2 m_tau^5/(192 pi^3) = 0.405 meV, inside the observed
weak width: no conflict, and the tau becomes a ten-orders-of-magnitude
consistency lock between the generation sector and the weak sector's
chirality structure. REMAINING DERIVATION (assigned, not done): show from
the temporal-binding/confinement machinery that the registry-diagonal
deformation coupling to wound phonons is absent for asymptotic states, or
equivalently derive the co-sliding covariance argument that non-chiral
couplings cannot change registry.
"""
import numpy as np

m0, me = 0.51099895069*137.035999177, 0.51099895069        # MeV
mmu, mtau = 105.6583755, 1776.86                            # MeV, PDG
E1 = 3*m0/np.sqrt(2)
alpha = 1/137.035999177
GF = 1.1663788e-5                                           # GeV^-2

if __name__ == "__main__":
    for name, m in [("electron", me), ("muon", mmu), ("tau", mtau)]:
        thr = m - me
        state = "CLOSED (stable)" if thr < E1 else "OPEN (gate required)"
        print(f"{name:8s}: m - m_e = {thr:9.3f} MeV vs E_1 = {E1:.3f} -> {state}")
    print(f"muon ratio E_1/(m_mu-m_e) = {E1/(mmu-me):.5f} (sqrt2 to "
          f"{E1/(mmu-me)/np.sqrt(2)-1:+.2%})")
    print(f"ungated width ~ alpha m_tau = {alpha*mtau:.1f} MeV; observed "
          f"{6.582119569e-16/2.903e-13*1e3:.2f} meV; exclusion ~1e10")
    print(f"gated width ~ G_F^2 m_tau^5/(192 pi^3) = "
          f"{GF**2*(mtau/1e3)**5/(192*np.pi**3)*1e12:.3f} meV: consistent")
