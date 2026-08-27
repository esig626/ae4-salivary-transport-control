# Baseline parameters

## Policy

Values below are transcribed from the published 2018 paper’s parameter tables. “Paper-derived” means the value was calibrated or solved within that paper, not independently measured. Such values are nominal reconstruction targets, not validated biological constants.

## Global, geometry, water, and conductance parameters

Published location: Appendix 11, Table 2.

| Symbol | Description | Value | Units | Provenance in paper |
|---|---|---:|---|---|
| `G_CaCC` | Maximum CaCC conductance | 71.3 | nS | Paper-derived |
| `G_CaKC` | Maximum CaKC conductance | 30.4 | nS | Paper-derived |
| `G_Na^t` | Maximum paracellular Na+ conductance | 12.46 | nS | Paper-derived |
| `G_K^t` | Maximum paracellular K+ conductance | 0.9 | nS | Paper-derived |
| `ω_i0` | Cell volume | 1.3 | pL | Palk et al. (2010) |
| `ω_l/ω_i0` | Lumen-to-cell volume ratio | 0.02 | dimensionless | Palk et al. (2010) |
| `P_a` | Apical water permeability | 4.32 × 10^-12 | L^2 mol^-1 s^-1 | Palk et al. (2010) |
| `P_b` | Basolateral water permeability | 5.15 × 10^-11 | L^2 mol^-1 s^-1 | Palk et al. (2010) |
| `P_t` | Tight-junction water permeability | 2.6 × 10^-13 | L^2 mol^-1 s^-1 | Palk et al. (2010) |
| `R` | Universal gas constant | 8.3144621 | J mol^-1 K^-1 | Moldover et al. (1988) |
| `T` | Temperature | 310 | K | Table 2; no row-level citation marker |
| `F` | Faraday constant | 96485.3365 | C mol^-1 | DeBiévre et al. (1994) |

## Nkcc1

Published location: Appendix 1, Table 3.

| Symbol | Value | Units | Provenance |
|---|---:|---|---|
| `α_Nkcc1` | 2.15 | amol/μm^3 | Paper-derived |
| `a_1` | 157.5 | s^-1 | Palk et al. (2010) |
| `a_2` | 2.0096 × 10^7 | mM^-4 s^-1 | Palk et al. (2010) |
| `a_3` | 1.0306 | s^-1 | Palk et al. (2010) |
| `a_4` | 1.3852 × 10^6 | mM^-4 s^-1 | Palk et al. (2010) |

## NaK-ATPase

Published location: Appendix 2, Table 4.

| Symbol | Value | Units | Provenance |
|---|---:|---|---|
| `α_NaK` | 4.84 | amol/μm^3 | Paper-derived |
| `r` | 1.305 × 10^6 | mM^-3 s^-1 | Palk et al. (2010) |
| `α_1` | 0.641 | mM^-1 | Palk et al. (2010) |

## Channel gating constants

Published locations: Appendices 3–4.

| Symbol | Value | Units | Provenance |
|---|---:|---|---|
| `η_2` | 2.54 | dimensionless | CaKC model attributed to Takahata et al. (2003) |
| `K_CaKC` | 0.182 | μM | Palk et al. (2010) |
| `η_1` | 1.46 | dimensionless | CaCC simplification in 2018 paper |
| `K_CaCC` | 0.26 | μM | CaCC simplification in 2018 paper |
| `z_K` | +1 | dimensionless | Ionic valence |
| `z_Cl` | -1 | dimensionless | Ionic valence |
| `z_Na` | +1 | dimensionless | Ionic valence |

## Ae2 and Nhe1

Published locations: Appendix 5 Table 5 and Appendix 6 Table 6.

| Mechanism | Symbol | Value | Units | Provenance |
|---|---|---:|---|---|
| Ae2 | `G_Ae2` | 0.01807 | fmol/s | Paper-derived |
| Ae2 | `K_Cl` | 5.6 | mM | Falkenberg and Jakobsson (2010) |
| Ae2 | `K_B` | 10^4 | mM | Falkenberg and Jakobsson (2010) |
| Nhe1 | `G_Nhe1` | 0.0305 | fmol/s | Paper-derived |
| Nhe1 | `K_H` | 4.5 × 10^-4 | mM | Falkenberg and Jakobsson (2010) |
| Nhe1 | `K_Na` | 15 | mM | Falkenberg and Jakobsson (2010) |

## Ae4

Published location: Appendix 7, Table 7.

| Symbol | Value | Units printed in table | Provenance |
|---|---:|---|---|
| `G_Ae4` | 0.66 | amol/μm^3 | Paper-derived |
| `k_+` | 1.92 × 10^-2 | mM^-4 s^-1 | Paper-derived |
| `k_-` | 1.3 × 10^-5 | mM^-4 s^-1 | Paper-derived |

The prose after Eq. (31) instead describes `G_Ae4` as having units of fmol and `J_Ae4` as concentration/time, while the main balances require amount/time. This family is not implementation-ready.

## CO2 transport and bicarbonate buffer

Published location: Appendix 8, Table 8; values attributed to Sharp et al. (2015).

| Symbol | Value | Units printed in table |
|---|---:|---|
| `k_1` | 11 | s^-1 |
| `k_-1` | 2.6 × 10^4 | s^-1 |
| `P_CO2` | 1.97 × 10^-13 | s^-1 |

The dimensional adequacy of `k_-1` for a term proportional to `[H+][HCO3−]` must be checked against the upstream source.

## Resting bath values and calibration targets

Published location: Table 1. These rows mix fixed inputs, experimental targets, and paper-derived values.

| Quantity | Nominal value | Source status in Table 1 |
|---|---:|---|
| `[Ca2+]_i` | 58 nM | Foskett and Melvin (1989) |
| `[Cl−]_i` | experimental 50.1 mM; model 50 mM | Peña-Münzenmayer et al. (2015) |
| `[K+]_i` | experimental/model 120 mM | Pedersen and Petersen (1973) |
| `[Na+]_i` | experimental 20 mM; model 25 mM | Grinstein and Foskett (1990) |
| `[HCO3−]_i` | experimental 12 mM; model 12.1 mM | Crampin et al. (2006) |
| `pH_i` | experimental 6.8; model 6.91 | Peña-Münzenmayer et al. (2015) |
| `[CO2]_i` | 6.6 mM | Paper-derived |
| `x_i/ω_i0` | 82.8 mM | Paper-derived |
| `V_a` | experimental -50.2 mV; model -50.24 mV | Lau and Case (1988) |
| `V_b` | experimental -61.8 mV; model -62.8 mV | Lau and Case (1988) |
| `[Cl−]_l` | 124.3 mM | Palk et al. (2010) |
| `[K+]_l` | 5.6 mM | Mangos et al. (1973) |
| `[Na+]_l` | 118.7 mM | Palk et al. (2010) |
| `[HCO3−]_l` | 1.5 × 10^-4 mM | Paper-derived |
| `pH_l` | 6.81 | Jayaraman et al. (2001); model column incorrectly appends “mM” |
| `Ψ_l` | 48.8 mM | Paper-derived from steady-state volume/flow |
| `[Cl−]_e` | 102.6 mM | Mangos et al. (1973) |
| `[K+]_e` | 5.3 mM | Mangos et al. (1973) |
| `[Na+]_e` | 140.2 mM | Mangos et al. (1973) |
| `[HCO3−]_e` | 42.9 mM | Paper-derived |
| `pH_e` | 7.4 | Mangos et al. (1973) |
| `[CO2]_e` | 1.9 mM | Paper-derived |
| `ω_i` | 1.3 pL | Palk et al. (2010) |

## Missing or non-numeric inputs

- `C_m` appears in Eqs. (7)–(8) but is not listed in Table 2; QSS removes it from the reported simulations.
- The prescribed `[Ca2+]_i(t)` is shown only as a qualitative figure.
- `ω_l` is available only through its ratio to `ω_i0`.
- Solver tolerances, DAE initialization controls, and any parameter-fitting objective/weights are not reported.
