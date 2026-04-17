# PDG vs Lattice-Monograph Comparison

**Generated:** from PDG 2024 edition (`mass_width_2024.txt`, 31-May-2024) and monograph v37.

**Scope:** 230 entries. 90 have predictions from the framework; 140 are currently gaps.

**Framework constants:**
- `m_e = 0.5109989461 MeV` (unit-defining, screw dislocation)
- `alpha = 1/137.035999177` (derived, Ch. alpha_PN)
- `m_0 = m_e/alpha = 70.0253 MeV` (node mass, hadronic tier)
- `M_EW = m_e/alpha^2 = 9.596 GeV` (electroweak tier)
- Six FCC integers: Z_1=12, N_c=3, N_CB=5, HEX=7, BILA=8, N_{111}=4

**Master formula (Eq. universal_mass):** `m = m_e * (N/alpha^k + Q)` with k=1 for hadrons, k=2 for EW.

## Summary by prediction pathway

| Pathway | Count | Notes |
|---|---|---|
| Recipe card (Sec. recipe_card) | 48 | Ground-state hadrons: N and Q from building-block tables |
| N_eff (crossed-fault) | 2 | rho/omega family: N=144/13 |
| Regge trajectory (Sec. regge) | 10 | Excited mesons: alpha' = 1/(8 pi^3 m_0^2) |
| Radial excitation (charmonium) | 2 | N(n_r) = N_quarks + N_ribbon + n_r * 8 |
| NLO tensor | 1 | f_2(1270): (1 + alpha/pi) correction |
| Exotic (docking mode) | 10 | Sec. molecular_exotics; 32 classified |
| Electroweak | 4 | M_W, M_Z, m_H, m_t formulas |
| Koide leptons | 3 | e, mu, tau from Sum=27 m_0 - 15 m_e |
| Neutrinos | 3 | Z_3 model from theta_ch = alpha^2/(2 pi) |
| Derived EW/CKM | 7 | Quark masses, alpha^-1, G, CKM |
| **GAP** (rule not yet written) | **140** | See gap analysis below |

## Gap analysis — where new rules are needed

Gaps are particles in the 2024 PDG that the current monograph rules do not yet assign a predicted mass. Counts per sub-category:

| Category | # gaps | Nature of extension needed |
|---|---:|---|
| Bottomonium (absolute masses) | 18 | Primary limitation: N_b = 6 pi^2 transcendental causes disclination-ribbon overlap. Splittings only in Table C. |
| Axial vectors & higher-L strange mesons | 15 | P-wave J=1 light mesons (analogues of chi_c1, h_c). Coord-shell + 1 with light quarks. |
| Light mesons (radial/orbital excitations) | 14 | Generalise the charmonium radial rule N(n_r) = N_quarks + N_ribbon + n_r * BILA to light quarks; analogue for L=1,2 orbitals. |
| Lambda* resonances | 12 | Lambda* includes Lambda(1405) anomaly; likely S-wave anti-kaon-nucleon molecule + compact 1-star. |
| Scalar mesons (J=0+) | 11 | P-wave J=0 light mesons. Charmonium analogue is chi_c0 (N=49); need light-quark version. |
| N* resonances | 10 | Baryon-resonance rule not written. Roper is a radial excitation of N; L=1 states are orbital. |
| Delta* resonances | 9 | Decuplet radial/orbital excitations; analogous to N* but with void activation. |
| Excited tensor / exotic-JPC mesons | 8 | Needs orbital L=2 rule for tensors; 1^{-+} hybrids need a new lattice structure. |
| Excited charm baryons | 8 | L=1,2 excitations of Lambda_c/Xi_c/Omega_c — need baryon orbital rule. |
| Excited charm mesons (P-wave) | 8 | P-wave D mesons; the +1 rule from charmonium plus open-charm coupling. |
| Sigma* resonances | 7 | Sigma* radial/orbital; hex-cap + void excited states. |
| Higher charmonium (D-wave, 2P, 3S+) | 6 | Need D-wave (L=2) charmonium rule; n_r >= 2 rule. |
| Bottom mesons (ground state) | 5 | Absolute B-meson mass not predicted (N_b issue). Splittings are. |
| Xi* resonances | 4 | Double-strange excited; double-hex-cap + orbital. |
| Excited bottom baryons | 2 | Same as above, b sector. Primary sector affected by N_b transcendental. |
| Excited bottom mesons | 2 | Needs absolute B-meson masses first. |
| Omega* resonances | 1 | sss orbital excitation of Omega-. |

## Full particle-by-particle table

Column definitions:
- **Particle** — PDG name (charges as in PDG)
- **I^G (J^PC)** — quantum numbers
- **Pathway** — which framework rule assigns the mass (see summary)
- **m_pred** — predicted mass in MeV (from recipe, Regge, etc.)
- **N, Q** — lattice integers if recipe-card applies
- **m_obs** — PDG 2024 observed mass (MeV)
- **Δ%** — (m_pred − m_obs)/m_obs
- **Pull** — (m_pred − m_obs)/σ_combined, with σ_th = α × m_pred
- **Monograph** — section reference (`sec:…`)
- **Notes** — quark content, derivation notes, or gap description

### Fundamental

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| alpha^-1 | - | EW/derived | 137.04 | — | — | — | — | — | `—` | FCC Cosserat geometry (PN), Chap. alpha_PN. PDG: 137.035999177(21), +0.003 ppb. |
| G (Newton) | - | EW/derived | 0.00 | — | — | — | — | — | `—` | Born-cluster tunnelling: G x K_sf = c^4/l^2. Obs 6.6743(15)e-11 m^3/kg/s^2 (-2 ppm). |

### Electroweak

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| H (Higgs)(0) | 0(0^+) | EW | 125.20 | — | — | — | — | — | `—` | EW breathing + volumetric, Sec. higgs_mass. Pull 0.00. |
| W(+) | - | EW | 80390.0 | — | — | 80369.0 | +0.03% | +0.04 | `—` | M_W = (8pi/3) * M_EW, Sec. w_mass. |
| Z(0) | - | EW | 91200.0 | — | — | 91188.0 | +0.01% | +0.02 | `—` | M_Z = 3pi * M_EW, Sec. z_mass. |
| t (top)(+) | - | EW | 172820.0 | — | — | — | — | — | `—` | Heavy-sector cluster mass, Sec. top_mass. |

### Quarks

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| u(+) | - | EW/derived | 2.19 | — | — | 2.16 | +1.39% | +0.25 | `—` | Cell-pair, Sec. quark_masses. |
| d(-) | - | EW/derived | 4.74 | — | — | 4.70 | +0.85% | +0.33 | `—` | Cell-pair + dipole, Sec. quark_masses. |
| s(-) | - | EW/derived | 93.40 | — | — | 93.50 | -0.11% | -0.10 | `—` | Casimir of hex cap, Sec. strange. |
| c(+) | - | EW/derived | 1270.0 | — | — | 1273.0 | -0.24% | -0.28 | `—` | K_{9,9} antibonding: m_c = 18 m_0 = 1260, +NLO. |
| b(-) | - | EW/derived | 4177.0 | — | — | 4183.0 | -0.14% | -0.19 | `—` | Disclination of coord shell, Sec. bottom. |

### Leptons

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| e(-) | - | Koide | 0.51 | — | — | 0.51 | <0.01% | +0.00 | `—` | Unit-defining: m_e = m_s. Screw dislocation. |
| mu(-) | - | Koide | 105.65 | — | — | 105.66 | <0.01% | -0.01 | `—` | Koide 2nd generation. |
| tau(-) | - | Koide | 1776.9 | — | — | 1776.9 | <0.01% | -0.00 | `—` | Koide 3rd generation. |

### Neutrinos

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| nu(1)(0) | - | Neutrino | 0.01 | — | — | — | — | — | `—` | Leading: m_1 = alpha^3 m_e/(4pi^2). |
| nu(2)(0) | - | Neutrino | 0.01 | — | — | — | — | — | `—` | Z_3 model; Dm21^2 predicted. |
| nu(3)(0) | - | Neutrino | 0.05 | — | — | — | — | — | `—` | Z_3 model; Dm31^2 predicted. |

### Light pseudoscalars (ground state)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| pi(+) | 1^-(0^-) | Recipe | 139.54 | 2.0 | -1 | 139.57 | -0.02% | -0.03 | `sec:pion` | Cell pair (N=2), charged bond. |
| pi(0) | 1^-(0^-+) | Recipe | 134.94 | 2.0 | -10 | 134.98 | -0.03% | -0.04 | `sec:pion` | Cell pair, Coulomb self-energy -10. |
| eta(0) | 0^+(0^-+) | Recipe | 547.94 | 8.0 | -24 | 547.86 | +0.01% | +0.02 | `sec:eta_meson` | Hex bilayer (N=8) + strangeness mix; Q from Gell-Mann--Okubo. |
| eta'(0) | 0^+(0^-+) | Recipe | 959.40 | 14.0 | -41 | 957.78 | +0.17% | +0.23 | `sec:eta_prime` | Shell+cap (N=14) + U(1)_A anomaly winding. |

### Strange pseudoscalars (ground state)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| K(+) | 1/2(0^-) | Recipe | 493.75 | 7.0 | +7 | 493.68 | +0.02% | +0.02 | `sec:kaon` | Hex cap (N=7); kaon Q=+7. |
| K(0) | 1/2(0^-) | Recipe | 497.33 | 7.0 | +14 | 497.61 | -0.06% | -0.08 | `sec:kaon` | Hex cap, d-quark isospin step. |

### Light vectors (ground state)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| rho(770)(+-0) | 1^+(1^-) | N_eff | 775.66 | 11.1 | +0 | 775.26 | +0.05% | +0.07 | `sec:rho_meson` | Crossed-fault N_eff = Z^2/(Z+1) = 144/13. Isovector Q=0. |
| omega(782)(0) | 0^-(1^-) | N_eff | 781.80 | 11.1 | +12 | 782.66 | -0.11% | -0.15 | `sec:rho_meson` | Isoscalar partner of rho. |

### Strange vectors (ground state)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| K*(892)(+-0) | 1/2(1^-) | Recipe | 891.93 | 13.0 | -36 | 895.55 | -0.40% | -0.56 | `sec:Kstar` | Strangeness pins to full shell (N=13). Q=-4 N_c^2. |
| phi(1020)(0) | 0^-(1^-) | Recipe | 1019.5 | — | — | 1019.5 | <0.01% | +0.01 | `sec:phi_meson` | ss-bar vector; GMO mixing. Predicted 1019.5 MeV (Table B). |

### Tensor (Cosserat rotation)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| f(2)(1270)(0) | 0^+(2^+) | NLO | 1275.5 | 18.2 | +2 | 1275.4 | <0.01% | +0.01 | `sec:tensor_mass` | Rotational sector: 19-node Born cluster optical microrotation. NLO (1+a/pi). |

### Light mesons (Regge trajectory)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| rho(3)(1690)(+,0,-) | 1^+(3^-) | Regge | 1689.0 | — | — | 1688.8 | +0.01% | +0.02 | `sec:regge` | Regge anchor J=3 (rho_3). Monograph: 1689 MeV (anchor). |
| f(4)(2050)(0) | 0^+(4^+) | Regge | 2017.0 | — | — | 2018.0 | -0.05% | -0.05 | `sec:regge` | Regge rho/f_2 trajectory. Monograph: 2017 MeV. |
| rho(5)(2350)(+,0,-) | 1^+(5^-) | Regge | 2299.0 | — | — | — | — | — | `sec:regge` | Regge rho/f_2 trajectory. Monograph: 2299 MeV. |
| a(2)(1320)(+,0,-) | 1^-(2^+) | Regge | 1348.0 | — | — | 1318.2 | +2.26% | +3.02 | `sec:regge` | Regge rho/a_2 trajectory. Monograph: 1348 MeV (+2.3%). |
| a(4)(1970)(+,0,-) | 1^-(4^+) | Regge | 2062.0 | — | — | 1967.0 | +4.83% | +4.33 | `sec:regge` | Regge rho/a_2 trajectory J=4. Monograph: 2062 MeV. |

### Strange mesons (Regge trajectory)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| K(2)*(1430)(+,0) | 1/2(2^+) | Regge | 1418.0 | — | — | 1432.4 | -1.01% | -1.38 | `sec:regge` | Regge K* trajectory, J=2. Monograph: 1418 MeV (-0.6%). |
| K(3)*(1780)(+,0) | 1/2(3^-) | Regge | 1797.0 | — | — | 1779.0 | +1.01% | +1.17 | `sec:regge` | Regge K* trajectory, J=3. Monograph: 1797 MeV (+1.2%). |
| K(4)*(2045)(+,0) | 1/2(4^+) | Regge | 2108.0 | — | — | 2048.0 | +2.93% | +3.37 | `sec:regge` | Regge K* trajectory, J=4. Monograph: 2108 MeV (+3.1%). |
| f(2)'(1525)(0) | 0^+(2^+) | Regge | 1502.0 | — | — | 1517.3 | -1.01% | -1.36 | `sec:regge` | Regge phi/ss trajectory, J=2. Monograph: 1502 MeV (-1.0%). |
| phi(3)(1850)(0) | 0^-(3^-) | Regge | 1863.0 | — | — | 1854.0 | +0.49% | +0.59 | `sec:regge` | Regge phi/ss trajectory, J=3. Monograph: 1863 MeV (+0.5%). |

### Light mesons (radial/orbital excitations)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| pi(1300)(+,0) | 1^-(0^-) | **GAP** | — | — | — | 1300.0 | — | — | `—` | GAP — Radially-excited pion (n_r=1, analogue of eta_c(2S)) |
| pi(1800)(+,0) | 1^-(0^-) | **GAP** | — | — | — | 1810.0 | — | — | `—` | GAP — Radially-excited pion (n_r=2) |
| eta(1295)(0) | 0^+(0^-+) | **GAP** | — | — | — | 1294.0 | — | — | `—` | GAP — Radially-excited eta |
| eta(1405)(0) | 0^+(0^-+) | **GAP** | — | — | — | 1408.7 | — | — | `—` | GAP — Radially-excited eta / gluonium candidate |
| eta(1475)(0) | 0^+(0^-+) | **GAP** | — | — | — | 1476.0 | — | — | `—` | GAP — Radially-excited eta_ss |
| eta(2)(1645)(0) | 0^+(2^-+) | **GAP** | — | — | — | 1617.0 | — | — | `—` | GAP — L=2 eta excitation |
| eta(2)(1870)(0) | 0^+(2^-+) | **GAP** | — | — | — | 1842.0 | — | — | `—` | GAP — L=2 eta excitation |
| rho(1450)(+,0,-) | 1^+(1^-) | **GAP** | — | — | — | 1465.0 | — | — | `—` | GAP — Radially-excited rho (n_r=1) |
| rho(1700)(+,0,-) | 1^+(1^-) | **GAP** | — | — | — | 1720.0 | — | — | `—` | GAP — Radially/orbitally-excited rho |
| omega(1420)(0) | 0^-(1^-) | **GAP** | — | — | — | 1410.0 | — | — | `—` | GAP — Radially-excited omega |
| omega(1650)(0) | 0^-(1^-) | **GAP** | — | — | — | 1670.0 | — | — | `—` | GAP — Radially-excited omega |
| omega(3)(1670)(0) | 0^-(3^-) | **GAP** | — | — | — | 1667.0 | — | — | `—` | GAP — Regge omega J=3 |
| phi(1680)(0) | 0^-(1^-) | **GAP** | — | — | — | 1680.0 | — | — | `—` | GAP — Radially-excited phi |
| pi(2)(1670)(+,0) | 1^-(2^-+) | **GAP** | — | — | — | 1670.6 | — | — | `—` | GAP — L=2 pion excitation |

### Scalar mesons (J=0+)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| f(0)(500)(0) | 0^+(0^++) | **GAP** | — | — | — | 600.00 | — | — | `—` | GAP — Broad scalar/sigma |
| f(0)(980)(0) | 0^+(0^++) | **GAP** | — | — | — | 990.00 | — | — | `—` | GAP — f_0(980), KK molecule candidate |
| f(0)(1370)(0) | 0^+(0^++) | **GAP** | — | — | — | 1350.0 | — | — | `—` | GAP — f_0(1370) |
| f(0)(1500)(0) | 0^+(0^++) | **GAP** | — | — | — | 1522.0 | — | — | `—` | GAP — glueball candidate |
| f(0)(1710)(0) | 0^+(0^++) | **GAP** | — | — | — | 1733.0 | — | — | `—` | GAP — glueball candidate |
| f(0)(2020)(0) | 0^+(0^++) | **GAP** | — | — | — | 1982.0 | — | — | `—` | GAP — f_0(2020) |
| a(0)(980)(+,0,-) | 1^-(0^+) | **GAP** | — | — | — | 980.00 | — | — | `—` | GAP — Scalar isovector; KK molecule? |
| a(0)(1450)(+,0,-) | 1^-(0^+) | **GAP** | — | — | — | 1439.0 | — | — | `—` | GAP — a_0(1450) |
| K(0)*(700)(+,0) | 1/2(0^+) | **GAP** | — | — | — | 845.00 | — | — | `—` | GAP — kappa, broad scalar kaon |
| K(0)*(1430)(+,0) | 1/2(0^+) | **GAP** | — | — | — | 1430.0 | — | — | `—` | GAP — K_0*(1430), scalar kaon |
| K(0)*(1950)(+,0) | 1/2(0^+) | **GAP** | — | — | — | 1957.0 | — | — | `—` | GAP — Excited scalar kaon |

### Axial vectors & higher-L strange mesons

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| a(1)(1260)(+,0,-) | 1^-(1^+) | **GAP** | — | — | — | 1230.0 | — | — | `—` | GAP — a_1(1260), isovector 1++ |
| a(1)(1640)(+,0,-) | 1^-(1^+) | **GAP** | — | — | — | 1655.0 | — | — | `—` | GAP — a_1(1640) |
| b(1)(1235)(+,0,-) | 1^+(1^+) | **GAP** | — | — | — | 1229.5 | — | — | `—` | GAP — b_1(1235), isovector 1+- |
| f(1)(1285)(0) | 0^+(1^+) | **GAP** | — | — | — | 1281.8 | — | — | `—` | GAP — f_1(1285), isoscalar 1++ |
| f(1)(1420)(0) | 0^+(1^+) | **GAP** | — | — | — | 1428.4 | — | — | `—` | GAP — f_1(1420) |
| h(1)(1170)(0) | 0^-(1^+) | **GAP** | — | — | — | 1166.0 | — | — | `—` | GAP — h_1(1170), isoscalar 1+- |
| h(1)(1415)(0) | 0^-(1^+) | **GAP** | — | — | — | 1409.0 | — | — | `—` | GAP — h_1(1415) |
| K(1)(1270)(+,0) | 1/2(1^+) | **GAP** | — | — | — | 1253.0 | — | — | `—` | GAP — K_1(1270) |
| K(1)(1400)(+,0) | 1/2(1^+) | **GAP** | — | — | — | 1403.0 | — | — | `—` | GAP — K_1(1400) |
| K(1)(1650)(+,0) | 1/2(1^+) | **GAP** | — | — | — | 1650.0 | — | — | `—` | GAP — K_1(1650) |
| K(2)(1770)(+,0) | 1/2(2^-) | **GAP** | — | — | — | 1773.0 | — | — | `—` | GAP — K_2(1770) |
| K(2)(1820)(+,0) | 1/2(2^-) | **GAP** | — | — | — | 1819.0 | — | — | `—` | GAP — K_2(1820) |
| K(2)*(1980)(+,0) | 1/2(2^+) | **GAP** | — | — | — | 1990.0 | — | — | `—` | GAP — K_2*(1980) |
| K*(1410)(+,0) | 1/2(1^-) | **GAP** | — | — | — | 1414.0 | — | — | `—` | GAP — Radially-excited K* |
| K*(1680)(+,0) | 1/2(1^-) | **GAP** | — | — | — | 1718.0 | — | — | `—` | GAP — Radially-excited K* |

### Excited tensor / exotic-JPC mesons

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| f(2)(1430)(0) | 0^+(2^+) | **GAP** | — | — | — | — | — | — | `—` | GAP — f_2(1430) extra tensor |
| f(2)(1565)(0) | 0^+(2^+) | **GAP** | — | — | — | 1571.0 | — | — | `—` | GAP — f_2(1565) |
| f(2)(1950)(0) | 0^+(2^+) | **GAP** | — | — | — | 1936.0 | — | — | `—` | GAP — f_2(1950) higher tensor |
| f(2)(2010)(0) | 0^+(2^+) | **GAP** | — | — | — | 2010.0 | — | — | `—` | GAP — f_2(2010) |
| f(2)(2300)(0) | 0^+(2^+) | **GAP** | — | — | — | 2297.0 | — | — | `—` | GAP — f_2(2300) |
| f(2)(2340)(0) | 0^+(2^+) | **GAP** | — | — | — | 2346.0 | — | — | `—` | GAP — f_2(2340) |
| a(2)(1700)(+,0,-) | 1^-(2^+) | **GAP** | — | — | — | 1706.0 | — | — | `—` | GAP — a_2(1700) excited |
| pi(1)(1600)(+,0,-) | 1^-(1^-+) | **GAP** | — | — | — | 1645.0 | — | — | `—` | GAP — Exotic 1-+ (hybrid candidate) |

### Charm mesons (ground state)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| D(0) | 1/2(0^-) | Recipe | 1865.1 | 25.0 | +224 | 1864.8 | +0.01% | +0.02 | `sec:charmonium_spectrum` | Charm + hex cap ribbon. Q = Q(K)+217. |
| D(+) | 1/2(0^-) | Recipe | 1869.7 | 25.0 | +233 | 1869.7 | <0.01% | +0.00 | `sec:charmonium_spectrum` | Charm + hex cap; Q = Q(D0)+9 (isospin). |
| D(s)(+) | 0(0^-) | Recipe | 1968.3 | 25.0 | +426 | 1968.4 | <0.01% | -0.00 | `sec:charmonium_spectrum` | Charm + strange ribbon; Q = Q(D0)+202. |
| D*(2007)(0) | 1/2(1^-) | Recipe | 2006.7 | 26.0 | +364 | 2006.8 | <0.01% | -0.01 | `sec:charmonium_spectrum` | D0 vector; bilayer ribbon; HF splitting 140 MeV. |
| D*(2010)(+) | 1/2(1^-) | Recipe | 2011.3 | 26.0 | +373 | 2010.3 | +0.05% | +0.07 | `sec:charmonium_spectrum` | D+ vector. |
| D(s)*(+) | 0(1^-) | Recipe | 2112.4 | 26.0 | +571 | 2112.2 | +0.01% | +0.02 | `sec:charmonium_spectrum` | Ds vector; HF+strange offsets. |

### Excited charm mesons (P-wave)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| D(0)*(2300)(+,0) | 1/2(0^+) | **GAP** | — | — | — | 2343.0 | — | — | `—` | GAP — P-wave D |
| D(1)(2420)(+,0) | 1/2(1^+) | **GAP** | — | — | — | 2422.1 | — | — | `—` | GAP — P-wave D |
| D(1)(2430)(0) | 1/2(1^+) | **GAP** | — | — | — | 2412.0 | — | — | `—` | GAP — P-wave D |
| D(2)*(2460)(+,0) | 1/2(2^+) | **GAP** | — | — | — | 2461.1 | — | — | `—` | GAP — P-wave D |
| D(s0)*(2317)(+) | 0(0^+) | **GAP** | — | — | — | 2317.8 | — | — | `—` | GAP — P-wave D_s |
| D(s1)(2460)(+) | 0(1^+) | **GAP** | — | — | — | 2459.5 | — | — | `—` | GAP — P-wave D_s |
| D(s1)(2536)(+) | 0(1^+) | **GAP** | — | — | — | 2535.1 | — | — | `—` | GAP — P-wave D_s |
| D(s2)*(2573)(+) | 0(2^+) | **GAP** | — | — | — | 2569.1 | — | — | `—` | GAP — P-wave D_s |

### Bottom mesons (ground state)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| B(c)(+) | 0(0^-) | Recipe | 6274.5 | — | — | 6274.5 | <0.01% | +0.00 | `sec:charmonium_spectrum` | c-bbar pseudoscalar; predicted 6274.5 MeV in Table C. |
| B(+) | 1/2(0^-) | **GAP** | — | — | — | 5279.4 | — | — | `N/A` | GAP: absolute B-meson masses not predicted (N_b = 6pi^2 transcendental). |
| B(0) | 1/2(0^-) | **GAP** | — | — | — | 5279.7 | — | — | `N/A` | GAP: absolute B-meson masses not predicted. |
| B(s)(0) | 0(0^-) | **GAP** | — | — | — | 5366.9 | — | — | `N/A` | GAP: absolute B_s mass not predicted. |
| B*(+,0) | 1/2(1^-) | **GAP** | — | — | — | 5324.8 | — | — | `N/A` | GAP: absolute mass not predicted; splitting B*-B = 45.1 MeV in Table C. |
| B(s)*(0) | 0(1^-) | **GAP** | — | — | — | 5415.4 | — | — | `N/A` | GAP: absolute mass not predicted; splitting Bs*-Bs = 48.5 MeV in Table C. |

### Excited bottom mesons

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| B(2)*(5747)(+,0) | 1/2(2^+) | **GAP** | — | — | — | 5739.6 | — | — | `—` | GAP — Excited B |
| B(s2)*(5840)(0) | 0(2^+) | **GAP** | — | — | — | 5839.9 | — | — | `—` | GAP — Excited B_s |

### Charmonium

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| eta(c)(1S)(0) | 0^+(0^-+) | Recipe | 2984.0 | 43.0 | -53 | 2984.1 | <0.01% | -0.00 | `sec:charmonium_spectrum` | 2x18 charm + hex-cap ribbon, Q_col=-54, Q_bond=+1. |
| J/psi(1S)(0) | 0^-(1^--) | Recipe | 3097.0 | 44.0 | +31 | 3096.9 | <0.01% | +0.00 | `sec:charmonium_spectrum` | Bilayer ribbon + HF shift +84. |
| chi(c0)(1P)(0) | 0^+(0^++) | Recipe | 3414.9 | 49.0 | -32 | 3414.7 | <0.01% | +0.01 | `sec:charmonium_spectrum` | Coord-shell ribbon (P-wave J=0). |
| chi(c1)(1P)(0) | 0^+(1^++) | Recipe | 3511.0 | 50.0 | +19 | 3510.7 | <0.01% | +0.01 | `sec:charmonium_spectrum` | Coord-shell+1 (P-wave J>=1); differs from h_c, chi_c2 only by Q. |
| h(c)(1P)(0) | ?^?(1^+-) | Recipe | 3525.8 | 50.0 | +48 | 3525.4 | +0.01% | +0.02 | `sec:charmonium_spectrum` | P-wave 1+-. |
| chi(c2)(1P)(0) | 0^+(2^++) | Recipe | 3556.5 | 50.0 | +108 | 3556.2 | <0.01% | +0.01 | `sec:charmonium_spectrum` | P-wave 2++. |
| eta(c)(2S)(0) | 0^+(0^-+) | Radial | 3637.7 | 51.0 | +130 | 3637.7 | <0.01% | -0.00 | `sec:charmonium_spectrum` | n_r=1 radial: N += 8 (bilayer quantum); absolute value from Table C. |
| psi(2S)(0) | 0^-(1^--) | Radial | 3686.3 | 52.0 | +88 | 3686.1 | <0.01% | +0.01 | `sec:charmonium_spectrum` | n_r=1 J/psi; HF(2S) = -42 applied to 1S pattern; absolute from Table C. |

### Higher charmonium (D-wave, 2P, 3S+)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| psi(3770)(0) | 0^-(1^--) | **GAP** | — | — | — | 3773.7 | — | — | `—` | GAP — D-wave 1^3D1 charmonium |
| psi(2)(3823)(0) | 0^-(2^--) | **GAP** | — | — | — | — | — | — | `—` | GAP — D-wave 1^3D2 |
| chi(c2)(3930)(0) | 0^+(2^++) | **GAP** | — | — | — | 3922.5 | — | — | `—` | GAP — Radial 2P chi_c2 |
| psi(4040)(0) | 0^-(1^--) | **GAP** | — | — | — | 4040.0 | — | — | `—` | GAP — 3S psi or 2D |
| psi(4160)(0) | 0^-(1^--) | **GAP** | — | — | — | 4191.0 | — | — | `—` | GAP — 2D or 3S |
| psi(4415)(0) | 0^-(1^--) | **GAP** | — | — | — | 4415.0 | — | — | `—` | GAP — 4S or higher |

### Bottomonium (absolute masses)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| eta(b)(1S)(0) | 0^+(0^-+) | **GAP** | — | — | — | 9398.7 | — | — | `N/A` | GAP: absolute mass not predicted; Upsilon-eta_b = 61.3 MeV splitting in Table C. |
| Upsilon(1S)(0) | 0^-(1^--) | **GAP** | — | — | — | 9460.4 | — | — | `N/A` | GAP: absolute mass not predicted (N_b transcendental). |
| chi(b0)(1P)(0) | 0^+(0^++) | **GAP** | — | — | — | 9859.4 | — | — | `—` | GAP (bottomonium) — Table C gives splittings only |
| chi(b1)(1P)(0) | 0^+(1^++) | **GAP** | — | — | — | 9892.8 | — | — | `—` | GAP (bottomonium) — Table C gives splittings only |
| h(b)(1P)(0) | 0^-(1^+-) | **GAP** | — | — | — | 9899.3 | — | — | `—` | GAP (bottomonium) — Table C gives splittings only |
| chi(b2)(1P)(0) | 0^+(2^++) | **GAP** | — | — | — | 9912.2 | — | — | `—` | GAP (bottomonium) — Table C gives splittings only |
| Upsilon(2S)(0) | 0^-(1^--) | **GAP** | — | — | — | 10023.4 | — | — | `—` | GAP (bottomonium) — Table C gives splitting only |
| Upsilon(2)(1D)(0) | 0^-(2^--) | **GAP** | — | — | — | 10163.7 | — | — | `—` | GAP (bottomonium) — D-wave bottomonium |
| chi(b0)(2P)(0) | 0^+(0^++) | **GAP** | — | — | — | 10232.5 | — | — | `—` | GAP (bottomonium) — absolute mass not predicted |
| chi(b1)(2P)(0) | 0^+(1^++) | **GAP** | — | — | — | 10255.5 | — | — | `—` | GAP (bottomonium) — absolute mass not predicted |
| h(b)(2P)(0) | 0^-(1^+-) | **GAP** | — | — | — | 10259.8 | — | — | `—` | GAP (bottomonium) — absolute mass not predicted |
| chi(b2)(2P)(0) | 0^+(2^++) | **GAP** | — | — | — | 10268.6 | — | — | `—` | GAP (bottomonium) — absolute mass not predicted |
| Upsilon(3S)(0) | 0^-(1^--) | **GAP** | — | — | — | 10355.1 | — | — | `—` | GAP (bottomonium) — absolute mass not predicted |
| chi(b1)(3P)(0) | 0^+(1^++) | **GAP** | — | — | — | 10513.4 | — | — | `—` | GAP (bottomonium) — absolute mass not predicted |
| chi(b2)(3P)(0) | 0^+(2^++) | **GAP** | — | — | — | 10524.0 | — | — | `—` | GAP (bottomonium) — absolute mass not predicted |
| Upsilon(4S)(0) | 0^-(1^--) | **GAP** | — | — | — | 10579.4 | — | — | `—` | GAP (bottomonium) — absolute mass not predicted |
| Upsilon(10860)(0) | 0^-(1^--) | **GAP** | — | — | — | 10885.2 | — | — | `—` | GAP (bottomonium) — Above BB-bar threshold |
| Upsilon(11020)(0) | 0^-(1^--) | **GAP** | — | — | — | 11000.0 | — | — | `—` | GAP (bottomonium) — Above BB-bar threshold |

### Baryon octet

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| p(+) | 1/2(1/2^+) | Recipe | 939.21 | 13.5 | -12 | 938.27 | +0.10% | +0.14 | `sec:proton` | Shell+winding; screw dislocation character. |
| n(0) | 1/2(1/2^+) | Recipe | 940.23 | 13.5 | -10 | 939.57 | +0.07% | +0.10 | `sec:proton` | Shell+winding; edge dislocation (Q_bond=-10). |
| Lambda(0) | 0(1/2^+) | Recipe | 1115.8 | 16.0 | -9 | 1115.7 | +0.01% | +0.01 | `sec:lambda_baryon` | Shell + 1 hex cap; colour screening Q_col=-N_c^2=-9. |
| Sigma(+) | 1(+)(1/2^+) | Recipe | 1189.4 | 17.0 | -2 | 1189.4 | <0.01% | +0.00 | `sec:sigma_baryon` | Shell + 4 voids (N=17); isospin stepping Q=-2. |
| Sigma(0) | 1(0)(1/2^+) | Recipe | 1192.5 | 17.0 | +4 | 1192.6 | -0.01% | -0.02 | `sec:sigma_baryon` | Shell + 4 voids (N=17); isospin stepping Q=4. |
| Sigma(-) | 1(-)(1/2^+) | Recipe | 1197.6 | 17.0 | +14 | 1197.4 | +0.01% | +0.02 | `sec:sigma_baryon` | Shell + 4 voids (N=17); isospin stepping Q=14. |
| Xi(0) | 1/2(1/2^+) | Recipe | 1314.6 | 19.0 | -31 | 1314.9 | -0.02% | -0.02 | `sec:xi_baryon` | Shell + 2 hex caps; Q_col=-N_c^3=-27 for octet. |
| Xi(-) | 1/2(1/2^+) | Recipe | 1321.8 | 19.0 | -17 | 1321.7 | <0.01% | +0.01 | `sec:xi_baryon` | Shell + 2 hex caps; Q_col=-N_c^3=-27 for octet. |

### Baryon decuplet

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| Delta(1232)(++,+,0,-) | 3/2(3/2^+) | Recipe | 1233.1 | 17.5 | +15 | 1232.0 | +0.09% | +0.12 | `sec:delta` | Shell + 4 voids + winding. Q=-Z1+N_c^3=+15 for Delta++. |
| Sigma(1385)(+,0,-) | 1(3/2^+) | Recipe | 1382.6 | 20.0 | -35 | 1387.2 | -0.33% | -0.45 | `sec:delta` | Shell + 4 voids + hex cap (N=20); full EM activation. |
| Xi(1530)(0,-) | 1/2(3/2^+) | Recipe | 1531.9 | 22.0 | -17 | 1535.0 | -0.20% | -0.28 | `sec:delta` | Shell + 2 hex caps + hex cap; decuplet Q_col=-9. |
| Omega(-) | 0(3/2^+) | Recipe | 1672.4 | 24.0 | -16 | 1672.5 | <0.01% | -0.00 | `sec:omega_baryon` | Triple bilayer (3*8=24); voids blocked. Q_col absorbed. |

### N* resonances

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| N(1440)(+,0) | 1/2(1/2^+) | **GAP** | — | — | — | 1440.0 | — | — | `—` | GAP — Roper; radial excitation of N |
| N(1520)(+,0) | 1/2(3/2^-) | **GAP** | — | — | — | 1515.0 | — | — | `—` | GAP — D13, L=1 |
| N(1535)(+,0) | 1/2(1/2^-) | **GAP** | — | — | — | 1530.0 | — | — | `—` | GAP — S11, L=1 |
| N(1650)(+,0) | 1/2(1/2^-) | **GAP** | — | — | — | 1650.0 | — | — | `—` | GAP — S11' |
| N(1675)(+,0) | 1/2(5/2^-) | **GAP** | — | — | — | 1675.0 | — | — | `—` | GAP — D15 |
| N(1680)(+,0) | 1/2(5/2^+) | **GAP** | — | — | — | 1685.0 | — | — | `—` | GAP — F15 |
| N(1700)(+,0) | 1/2(3/2^-) | **GAP** | — | — | — | 1720.0 | — | — | `—` | GAP — D13' |
| N(1710)(+,0) | 1/2(1/2^+) | **GAP** | — | — | — | 1710.0 | — | — | `—` | GAP — P11' |
| N(1720)(+,0) | 1/2(3/2^+) | **GAP** | — | — | — | 1720.0 | — | — | `—` | GAP — P13 |
| N(2190)(+,0) | 1/2(7/2^-) | **GAP** | — | — | — | 2180.0 | — | — | `—` | GAP — G17, high-spin |

### Delta* resonances

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| Delta(1600)(++,+,0,-) | 3/2(3/2^+) | **GAP** | — | — | — | 1570.0 | — | — | `—` | GAP — Delta radial excitation |
| Delta(1620)(++,+,0,-) | 3/2(1/2^-) | **GAP** | — | — | — | 1610.0 | — | — | `—` | GAP — S31 |
| Delta(1700)(++,+,0,-) | 3/2(3/2^-) | **GAP** | — | — | — | 1710.0 | — | — | `—` | GAP — D33 |
| Delta(1900)(++,+,0,-) | 3/2(1/2^-) | **GAP** | — | — | — | 1860.0 | — | — | `—` | GAP — S31' |
| Delta(1905)(++,+,0,-) | 3/2(5/2^+) | **GAP** | — | — | — | 1880.0 | — | — | `—` | GAP — F35 |
| Delta(1910)(++,+,0,-) | 3/2(1/2^+) | **GAP** | — | — | — | 1900.0 | — | — | `—` | GAP — P31 |
| Delta(1920)(++,+,0,-) | 3/2(3/2^+) | **GAP** | — | — | — | 1920.0 | — | — | `—` | GAP — P33' |
| Delta(1930)(++,+,0,-) | 3/2(5/2^-) | **GAP** | — | — | — | 1950.0 | — | — | `—` | GAP — D35 |
| Delta(1950)(++,+,0,-) | 3/2(7/2^+) | **GAP** | — | — | — | 1930.0 | — | — | `—` | GAP — F37 |

### Lambda* resonances

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| Lambda(1405)(0) | 0(1/2^-) | **GAP** | — | — | — | 1405.1 | — | — | `—` | GAP — Long-standing anomaly; S01 |
| Lambda(1520)(0) | 0(3/2^-) | **GAP** | — | — | — | 1519.0 | — | — | `—` | GAP — D03 |
| Lambda(1600)(0) | 0(1/2^+) | **GAP** | — | — | — | 1600.0 | — | — | `—` | GAP — P01' |
| Lambda(1670)(0) | 0(1/2^-) | **GAP** | — | — | — | 1674.0 | — | — | `—` | GAP — S01' |
| Lambda(1690)(0) | 0(3/2^-) | **GAP** | — | — | — | 1690.0 | — | — | `—` | GAP — D03' |
| Lambda(1800)(0) | 0(1/2^-) | **GAP** | — | — | — | 1800.0 | — | — | `—` | GAP — S01 |
| Lambda(1810)(0) | 0(1/2^+) | **GAP** | — | — | — | 1790.0 | — | — | `—` | GAP — P01 |
| Lambda(1820)(0) | 0(5/2^+) | **GAP** | — | — | — | 1820.0 | — | — | `—` | GAP — F05 |
| Lambda(1830)(0) | 0(5/2^-) | **GAP** | — | — | — | 1825.0 | — | — | `—` | GAP — D05 |
| Lambda(1890)(0) | 0(3/2^+) | **GAP** | — | — | — | 1890.0 | — | — | `—` | GAP — P03 |
| Lambda(2100)(0) | 0(7/2^-) | **GAP** | — | — | — | 2100.0 | — | — | `—` | GAP — G07 |
| Lambda(2110)(0) | 0(5/2^+) | **GAP** | — | — | — | 2090.0 | — | — | `—` | GAP — F05' |

### Sigma* resonances

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| Sigma(1660)(+,0,-) | 1(1/2^+) | **GAP** | — | — | — | 1660.0 | — | — | `—` | GAP — Radial Sigma excitation |
| Sigma(1670)(+,0,-) | 1(3/2^-) | **GAP** | — | — | — | 1675.0 | — | — | `—` | GAP — D13 |
| Sigma(1750)(+,0,-) | 1(1/2^-) | **GAP** | — | — | — | 1750.0 | — | — | `—` | GAP — S11 |
| Sigma(1775)(+,0,-) | 1(5/2^-) | **GAP** | — | — | — | 1775.0 | — | — | `—` | GAP — D15 |
| Sigma(1910)(+,0,-) | 1(3/2^-) | **GAP** | — | — | — | 1910.0 | — | — | `—` | GAP — D13' |
| Sigma(1915)(+,0,-) | 1(5/2^+) | **GAP** | — | — | — | 1915.0 | — | — | `—` | GAP — F15 |
| Sigma(2030)(+,0,-) | 1(7/2^+) | **GAP** | — | — | — | 2030.0 | — | — | `—` | GAP — F17 |

### Xi* resonances

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| Xi(1690)(-,0) | 1/2(1/2^-) | **GAP** | — | — | — | 1690.0 | — | — | `—` | GAP — Xi excitation |
| Xi(1820)(-,0) | 1/2(3/2^-) | **GAP** | — | — | — | 1823.0 | — | — | `—` | GAP — D13 Xi |
| Xi(1950)(-,0) | 1/2(?^?) | **GAP** | — | — | — | 1950.0 | — | — | `—` | GAP — Xi(1950) |
| Xi(2030)(-,0) | 1/2(?^?) | **GAP** | — | — | — | 2025.0 | — | — | `—` | GAP — Xi(2030) |

### Omega* resonances

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| Omega(2250)(-) | 0(?^?) | **GAP** | — | — | — | 2252.0 | — | — | `—` | GAP — Omega excitation |

### Charmed baryons

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| Lambda(c)(+) | 0(1/2^+) | Recipe | 2286.3 | 31.0 | +226 | 2286.5 | <0.01% | -0.01 | `sec:charm_evidence` | Lambda light + full charm 235; Q = -9 + 235 = 226. |
| Sigma(c)(2455)(++,+,0) | 1(1/2^+) | Recipe | 2454.0 | 35.0 | +6 | 2453.8 | +0.01% | +0.01 | `sec:charm_evidence` | Sigma light + Pauli-blocked charm (n_block=12 -> Q_charm=0). |
| Sigma(c)(2520)(++,+,0) | 1(3/2^+) | Recipe | 2518.0 | 36.0 | +4 | 2518.5 | -0.02% | -0.03 | `sec:charm_evidence` | Sigma*(decuplet) light + Pauli-blocked charm. |
| Xi(c)(+) | 1/2(1/2^+) | Recipe | 2467.7 | 34.0 | +170 | 2467.7 | <0.01% | +0.00 | `sec:charm_evidence` | Lambda-type light + Q_charm=180 (us: n_block=N_c). |
| Xi(c)(0) | 1/2(1/2^+) | Recipe | 2467.7 | 34.0 | +170 | 2470.4 | -0.11% | -0.15 | `sec:charm_evidence` | Isospin partner of Xi_c+. |
| Omega(c)(0) | 0(1/2^+) | Recipe | 2695.2 | 38.0 | +67 | 2695.2 | <0.01% | +0.00 | `sec:charm_evidence` | N_light=20 + N_charm=18; Q_charm=98 (ss: exchange restriction). |
| Xi(cc)(++) | 1/2(1/2^+) | Recipe | 3621.3 | 49.0 | +372 | — | — | — | `sec:charm_evidence` | 2 charms: primary Q_charm=235 + secondary Q_charm,2=146. |
| Xi(cc)(+) | 1/2(1/2^+) | Recipe | 3620.3 | 49.0 | +370 | — | — | — | `sec:charm_evidence` | Isospin partner; Q = 372 - 2. |

### Excited charm baryons

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| Lambda(c)(2595)(+) | 0(1/2^-) | **GAP** | — | — | — | 2592.2 | — | — | `—` | GAP — L=1 Lambda_c |
| Lambda(c)(2625)(+) | 0(3/2^-) | **GAP** | — | — | — | 2628.0 | — | — | `—` | GAP — L=1 Lambda_c |
| Lambda(c)(2880)(+) | 0(5/2^+) | **GAP** | — | — | — | 2881.6 | — | — | `—` | GAP — L=2 Lambda_c |
| Xi(c)(2645)(+,0) | 1/2(3/2^+) | **GAP** | — | — | — | 2646.2 | — | — | `—` | GAP — Xi_c* decuplet |
| Xi(c)(2790)(+,0) | 1/2(1/2^-) | **GAP** | — | — | — | 2793.9 | — | — | `—` | GAP — L=1 Xi_c |
| Xi(c)(2815)(+,0) | 1/2(3/2^-) | **GAP** | — | — | — | 2819.8 | — | — | `—` | GAP — L=1 Xi_c |
| Xi(c)'(+,0) | 1/2(1/2^+) | **GAP** | — | — | — | 2578.2 | — | — | `—` | GAP — Xi_c' |
| Omega(c)(2770)(0) | 0(3/2^+) | **GAP** | — | — | — | 2765.9 | — | — | `—` | GAP — Omega_c* |

### Bottom baryons

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| Lambda(b)(0) | 0(1/2^+) | Recipe | 5647.5 | 79.0 | +226 | 5619.6 | +0.50% | +0.68 | `sec:charm_evidence` | N_light=16 + 66 + 66-16 (doubled b-structure);+0.5% residual. |
| Sigma(b)(+) | 1(1/2^+) | Recipe | 5815.2 | 83.0 | +6 | 5810.6 | +0.08% | +0.11 | `sec:charm_evidence` | Sigma light + b. |
| Sigma(b)(-) | 1(1/2^+) | Recipe | 5815.2 | 83.0 | +6 | 5815.6 | <0.01% | -0.01 | `sec:charm_evidence` | Same rule. |
| Sigma(b)*(+) | 1(3/2^+) | Recipe | — | — | — | 5830.3 | — | — | `sec:charm_evidence` | Sigma* light + b; explicit formula in Sec. charm_evidence. |
| Sigma(b)*(-) | 1(3/2^+) | Recipe | — | — | — | 5834.7 | — | — | `sec:charm_evidence` | Sigma* light + b. |
| Xi(b)(-,0) | 1/2(1/2^+) | Recipe | — | — | — | 5797.0 | — | — | `sec:charm_evidence` | Xi light + b. |
| Omega(b)(-) | 0(1/2^+) | Recipe | 6056.4 | 86.0 | +67 | 6045.8 | +0.18% | +0.24 | `sec:charm_evidence` | Omega_c with c->b. |

### Excited bottom baryons

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| Xi(b)'(-,0) | 1/2(1/2^+) | **GAP** | — | — | — | — | — | — | `—` | GAP — Xi_b' |
| Xi(b)(5945)(0) | 1/2(3/2^+) | **GAP** | — | — | — | — | — | — | `—` | GAP — Xi_b* decuplet |

### Exotic hadrons (docking-mode classified)

| Particle | I^G(J^PC) | Pathway | m_pred (MeV) | N | Q | m_obs (MeV) | Δ% | Pull | Monograph | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| X(3872)(0,+) | 0,1^+(1^++) | Exotic | 3871.7 | — | — | — | — | — | `sec:molecular_exotics` | F (mol.): D0 D*0-bar molecule; n_bond=6 -> Gamma=6.1 MeV. |
| chi(c1)(3872)(0) | 0^+(1^++) | Exotic | 3871.7 | — | — | — | — | — | `sec:molecular_exotics` | F (mol.): Renamed X(3872); see above. |
| T(cc)(3875)(+) | 0(1^+) | Exotic | 3875.0 | — | — | — | — | — | `sec:molecular_exotics` | F (mol.): Tetraquark; D0 D*+ molecule. |
| Z(c)(3900)(+,0) | 1(1^+) | Exotic | — | — | — | — | — | — | `sec:molecular_exotics` | F (mol.): D D*-bar molecule; exotic classification. |
| P(c)(4312)(+) | ?(1/2^-) | Exotic | 4316.0 | — | — | — | — | — | `sec:molecular_exotics` | 2F (mol.): Sigma_c D-bar molecule; n_bond=12. |
| P(c)(4440)(+) | ?(1/2^-) | Exotic | — | — | — | — | — | — | `sec:molecular_exotics` | V (compact): Sigma_c* D-bar compact. |
| P(c)(4457)(+) | ?(1/2^-) | Exotic | — | — | — | — | — | — | `sec:molecular_exotics` | F (mol.): Sigma_c* D-bar molecule. |
| Y(4260)(0) | 0^-(1^--) | Exotic | — | — | — | — | — | — | `sec:molecular_exotics` | threshold: Hybrid or 4-quark; mode P. |
| T(4c)(6900)(0) | ?(2^++) | Exotic | — | — | — | — | — | — | `sec:molecular_exotics` | cage (SFT): CMS: 2++ confirms stacking-fault tetrahedron. |
| d*(2380)(0) | 0(3^+) | Exotic | 2382.0 | — | — | — | — | — | `sec:molecular_exotics` | V (compact): Dibaryon; mode V. Sec. dibaryon_worked. |

## Prediction quality — residuals for matched particles

- Matched particles: **67**
- Median |Δ%|: **0.016%**
- Max |Δ%|: **4.830%**
- Median |pull|: **0.02**
- Max |pull|: **4.33**

### Top 15 biggest residuals (|Δ%|)

| Particle | m_pred | m_obs | Δ% | Pull | Pathway | Notes |
|---|---:|---:|---:|---:|---|---|
| a(4)(1970)(+,0,-) | 2062.0 | 1967.0 | +4.83% | +4.33 | Regge | Regge rho/a_2 trajectory J=4. Monograph: 2062 MeV.... |
| K(4)*(2045)(+,0) | 2108.0 | 2048.0 | +2.93% | +3.37 | Regge | Regge K* trajectory, J=4. Monograph: 2108 MeV (+3.1%).... |
| a(2)(1320)(+,0,-) | 1348.0 | 1318.2 | +2.26% | +3.02 | Regge | Regge rho/a_2 trajectory. Monograph: 1348 MeV (+2.3%).... |
| u(+) | 2.19 | 2.16 | +1.39% | +0.25 | EW/derived | Cell-pair, Sec. quark_masses.... |
| K(3)*(1780)(+,0) | 1797.0 | 1779.0 | +1.01% | +1.17 | Regge | Regge K* trajectory, J=3. Monograph: 1797 MeV (+1.2%).... |
| f(2)'(1525)(0) | 1502.0 | 1517.3 | -1.01% | -1.36 | Regge | Regge phi/ss trajectory, J=2. Monograph: 1502 MeV (-1.0%).... |
| K(2)*(1430)(+,0) | 1418.0 | 1432.4 | -1.01% | -1.38 | Regge | Regge K* trajectory, J=2. Monograph: 1418 MeV (-0.6%).... |
| d(-) | 4.74 | 4.70 | +0.85% | +0.33 | EW/derived | Cell-pair + dipole, Sec. quark_masses.... |
| Lambda(b)(0) | 5647.5 | 5619.6 | +0.50% | +0.68 | Recipe | N_light=16 + 66 + 66-16 (doubled b-structure);+0.5% residual... |
| phi(3)(1850)(0) | 1863.0 | 1854.0 | +0.49% | +0.59 | Regge | Regge phi/ss trajectory, J=3. Monograph: 1863 MeV (+0.5%).... |
| K*(892)(+-0) | 891.93 | 895.55 | -0.40% | -0.56 | Recipe | Strangeness pins to full shell (N=13). Q=-4 N_c^2.... |
| Sigma(1385)(+,0,-) | 1382.6 | 1387.2 | -0.33% | -0.45 | Recipe | Shell + 4 voids + hex cap (N=20); full EM activation.... |
| c(+) | 1270.0 | 1273.0 | -0.24% | -0.28 | EW/derived | K_{9,9} antibonding: m_c = 18 m_0 = 1260, +NLO.... |
| Xi(1530)(0,-) | 1531.9 | 1535.0 | -0.20% | -0.28 | Recipe | Shell + 2 hex caps + hex cap; decuplet Q_col=-9.... |
| Omega(b)(-) | 6056.4 | 6045.8 | +0.18% | +0.24 | Recipe | Omega_c with c->b.... |

## Low-hanging fruit — gaps most likely to yield to modest rule extensions

These are gaps where an existing monograph rule has a natural analogue that would close them with a short derivation:

1. **P-wave light mesons (scalars, axial vectors)** — 23 gaps. Charmonium already has chi_c0 (scalar, N=49), chi_c1/h_c (axial, N=50), chi_c2 (tensor). The light-quark analogues need the same rule with the charm node count removed: `N_ribbon = 13` for J=0+, `N_ribbon = 14` for J>=1. This plus the existing Q_light decomposition should give f_0, a_0, a_1, b_1, f_1, h_1, K_0*, K_1.

2. **Radial excitations of light mesons (pi(1300), rho(1450), K(1460), eta(1295), phi(1680))** — 11 gaps. Charmonium's radial rule is `N(n_r) = N(0) + n_r * BILA`. For pi: N(0)=2, N(1)=10, giving `m_pi(1300) = 10 * 70 + Q * 0.511 MeV ≈ 700 MeV + Q corrections`. Need the Q shift rule for n_r=1 light mesons.

3. **Excited charm baryons (Lambda_c(2595), Lambda_c(2625), Xi_c(2645), Xi_c(2790), Xi_c(2815), Sigma_c(2520), Omega_c(2770))** — 8 gaps. These are L=1 orbital or spin-flip of known ground states. Charmonium already has chi_c (L=1) treatment; baryon analogue should give N += 1 for L=1 and a known Q_charm shift.

4. **Higher charmonium (psi(3770), psi_2(3823), chi_c2(3930))** — 3 gaps. psi(3770) is 1^3D_1 (L=2). Need L=2 ribbon rule. chi_c2(3930) is 2P excited, so n_r=1 applied to chi_c2 should work.

5. **Bottomonium absolute masses** — 16 gaps. Requires resolving the N_b = 6pi^2 transcendental problem. If a cleaner b-quark node count were derivable, the rest of the spectrum follows by the same recipe as charmonium.

6. **Roper N(1440) and L=1 nucleon resonances (N(1520), N(1535), N(1650), N(1675), N(1680))** — 6 gaps. Roper is a radial excitation of the nucleon; the L=1 resonances are orbital. Baryon radial/orbital rules mirror mesonic ones but acting on the 13-node coordination shell instead of a ribbon.

## Appendix: what's deliberately outside the scope of the current recipe

- **1-star and 2-star PDG states** are not in this table (we use the 2024 summary-table listings, which are 3-star and 4-star).
- **Nuclei** (d, ^3He, etc.) are treated in `sec:nuclear_binding` as cluster-adhesion compounds; not compared here.
- **Pure-gauge glueball candidates** (f_0(1710) etc.) are listed as gaps; the framework has candidate structures in Sec. graph_products but not a closed-form mass.
- **PDG pentaquarks and tetraquarks** beyond the 10 listed in exotics — a handful of post-2020 discoveries (P_c(4312) beyond, LHCb T_cc(3875), etc.) are in the exotics table where known.