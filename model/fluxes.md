# Baseline flux models

## Sign and unit policy

Flux directions below follow the verbal and balance-equation descriptions in the published 2018 model. Before code is written, Phase 01 must establish one explicit convention for positive transport, current (`I`) versus molar flux (`J`), membrane area/density scaling, and every conversion by Faraday’s constant. The paper is not fully consistent on these points.

## Transmembrane transporters, pump, and channels

| Symbol | Mechanism and direction in the baseline balances | Published constitutive location | Published provenance | Calibration status |
|---|---|---|---|---|
| `J_Nkcc1` | Imports 1 Na+, 1 K+, and 2 Cl− from interstitium | Appendix 1 Eq. (17), Table 3 | Two-state simplification from Palk et al. (2010), based on Benjamin and Johnson (1997) | `α_Nkcc1` fitted in 2018 paper |
| `J_NaK` / `I_NaK` | Extrudes 3 Na+ and imports 2 K+ per ATP | Appendix 2 Eq. (18), Table 4; current contribution in Eq. (7) | Palk et al. (2010), based on Smith and Crampin (2004) | `α_NaK` fitted in 2018 paper; current conversion not explicit |
| `J_CaKC` / `I_CaKC` | Basolateral K+ efflux | Appendix 3 Eqs. (19)–(21) | Hill open probability and ohmic/Nernst relation based on Takahata et al. (2003), with `K_CaKC` from Palk et al. (2010) | `G_CaKC` fitted in 2018 paper |
| `J_CaCC` / `I_CaCC` | Apical Cl− efflux to lumen | Appendix 4 Eqs. (22)–(24) | Paper-selected simplification similar to Takahata et al. (2003); compared with Arreola et al. (2002) | `G_CaCC` fitted in 2018 paper |
| `J_Ae2` | Imports Cl− while exporting HCO3−, bidirectionally | Appendix 5 Eq. (25), Table 5 | Falkenberg and Jakobsson (2010) Michaelis-Menten difference form | `G_Ae2` fitted in 2018 paper |
| `J_Nhe1` | Imports Na+ while exporting H+ | Appendix 6 Eq. (26), Table 6 | Falkenberg and Jakobsson (2010) | `G_Nhe1` fitted in 2018 paper |
| `J_Ae4` | Imports 1 Cl− while exporting 2 HCO3− and one monovalent cation, partitioned between Na+ and K+ by intracellular mole fractions | Appendix 7 Fig. 10, Eqs. (27)–(31), Table 7 | Stoichiometry from Peña-Münzenmayer et al. (2016); Markov construction method from Dupont et al. (2016) | All Table 7 values marked model-determined; units are inconsistent across text/table/balances |

## CO2 and buffer terms

| Symbol | Published relation | Location and provenance |
|---|---|---|
| `J_CO2_a` | Proportional to `[CO2]_l - [CO2]_i` | Appendix 8; permeability from Sharp et al. (2015) |
| `J_CO2_b` | Proportional to `[CO2]_e - [CO2]_i` | Appendix 8; permeability from Sharp et al. (2015) |
| `J_CO2` | Conservation convention is proportional to `[CO2]_l + [CO2]_e - 2[CO2]_i`; the paper prints its negative | Appendix 8; Phase 10 selects the sum of the two individually printed influxes and retains the aggregate printed sign only as a sensitivity branch |
| `J_Buffer` | Forward CO2 hydration minus reverse H+/HCO3− association | Appendix 8, Table 8; Sharp et al. (2015) |

## Tight-junction transport

Appendix 10 defines linear current-voltage relations for paracellular K+ and Na+ transport using `V_t = V_a - V_b` and Nernst potentials computed from lumen/interstitium concentration ratios. These terms supply lumen Na+ and K+ and enter both membrane-potential balances.

The appendix labels the expressions `J_K^t` and `J_Na^t`, while the main model uses `I_K^t` and `I_Na^t`. Table 2 uses uppercase `G` for maximum conductance while Appendix 10 uses lowercase `g`. These are not silently equated here.

## Water and convective fluxes

| Symbol | Function | Published location | Status |
|---|---|---|---|
| `q_a` | Osmotic water flux from cell to lumen, including `Ψ_l` and `x_i/ω_i` | §2.5 Eq. (10) | Direction and volume contribution resolved by conservation in Phase 10 |
| `q_b` | Osmotic water flux from interstitium to cell | §2.5 Eq. (11) | Direction and volume contribution resolved by conservation in Phase 10 |
| `q_t` | Paracellular osmotic water flux from interstitium/lumen gradient | §2.5 Eq. (12) | Published |
| `q_tot` | `q_a + q_t`; total lumen inflow and, under constant `ω_l`, ductal outflow | §2.5 Eq. (13) | Primary fluid-secretion observable |
| Convective ion loss | `q_tot [X]_l` for luminal Na+, K+, and Cl− | §2.6 Eqs. (14)–(16) | Published; Eq. (16) drops the `l` subscript in print |

## Required Phase 01 checks

1. Resolve all `I`/`J` conversions and dimensions from published sources.
2. Use the conservation-resolved `dω_i/dt=q_b-q_a` and verify that steady results remain invariant to the isolated printed row sign.
3. Obtain a numerical, citable `[Ca2+]_i(t)` input or document a new independent input model.
4. Reconcile the Ae4 table units with Eq. (31) and the amount balances.
5. Confirm the intended CO2 net-flux sign.
