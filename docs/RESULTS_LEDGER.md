# Results ledger

Only validated results belong in this file.

For each result record

- result identifier
- concise statement
- status as reproduced or new
- analysis path
- code commit
- input data or parameter set
- numerical tolerances
- figure or table pointer
- robustness status
- caveats

## Phase 10 validated results

| ID | Statement | Status | Analysis / code | Input and tolerance | Pointer | Robustness and caveats |
| --- | --- | --- | --- | --- | --- | --- |
| R10-01 | At a regular steady state, `D_theta u*=-F_u^-1 F_theta` and `D_theta H=h_theta-h_u F_u^-1 F_theta`; a locally constant rank-`r` observation has fiber dimension `2-r` | Independently derived for this project; standard theorem | `analysis/10_identifiability_discrimination/theory.md`; [code/result commit](https://github.com/esig626/ae4-salivary-transport-control/commit/6248e1eb43a21cd89f39f9c3f79f57e9c013c0f6) | `C1` maps, nonsingular `F_u`, locally constant rank | Theory, Proposition 1 | A regular scalar leaves a curve; critical rank alone is inconclusive; no global claim |
| R10-02 | The printed AE2/AE4 activity-forcing columns have rank two when both turnover factors are nonzero; a fixed complete state cannot be shared by distinct activities under fixed conditions | New model-specific conditional proposition | `analysis/10_identifiability_discrimination/theory.md`; commit `6248e1eb43a21cd89f39f9c3f79f57e9c013c0f6` | Common flux units, multiplicative activities, fixed nuisances; Cl/HCO3 minor `-phi2 phi4` | `panel_rank_summary.csv`; exact unit tests | Full-state result only; it does not certify any partial physiological panel |
| R10-03 | With `A=C-2N=J2+J4` and `B=J2+2J4`, the exchanger fluxes are exactly `J2=2A-B`, `J4=B-A`; nonnegative forward fluxes give `A<=B<=2A` | New checked reduction of the published balances | `reduction.md`, `theory.md`, `src/ae4_identifiability.py`; commit `6248e1eb43a21cd89f39f9c3f79f57e9c013c0f6` | Common signed amount-flux convention; inverse tested to absolute tolerance `1e-14` | `stoichiometric_geometry.csv`; `identifiability_stoichiometric_geometry.png` | Identifies cycle fluxes from internal aggregates, not activities from measured physiology |
| R10-04 | The printed Table 1 row is not a simultaneous steady state of the printed water, apical-current, and acid/base equations | Published-baseline reproduction audit | `reduction.md`, `src/ae4_identifiability.py`; commit `6248e1eb43a21cd89f39f9c3f79f57e9c013c0f6` | Direct substitution: `q_b/q_a=59.606`, printed-flow ratio `117.6`, apical residual `-67.839 pA`, equilibrium CO2 `3.51857 mM` versus `6.6 mM` | `baseline_printed_equation_audit.json`; regression tests | Diagnostic of the citable specification, not a recalibration or reproduced full model |
| R10-05 | The upstream Palk corrigendum converts the NKCC1 turnover at the printed state from approximately `-14.51 s^-1` under literal 2018-table units to `+0.598 s^-1` | Reproduced source-backed correction | `model/parameters.md`, `src/ae4_identifiability.py`; commit `6248e1eb43a21cd89f39f9c3f79f57e9c013c0f6` | Palk M-based coefficients converted to mM; regression interval `(0.59,0.61) s^-1` | `baseline_printed_equation_audit.json` | Does not resolve the fitted density scale, NaK convention, Ae2 value, or full baseline |
| R10-06 | A consistent complete residual-row sign change leaves implicit sensitivity invariant | New derivation and reproduced software check | `theory.md`, `src/identifiability_discrimination/`; commit `6248e1eb43a21cd89f39f9c3f79f57e9c013c0f6` | Diagonal signs `+/-1`; acceptance `1e-14` | `derivative_verification.json` | Does not make a physical flux-direction redefinition invariant |
| R10-07 | The retained affine fixture's implicit sensitivity matches independently differenced equilibria | Reproduced computational check | `src/identifiability_discrimination/`; commit `6248e1eb43a21cd89f39f9c3f79f57e9c013c0f6` | Declared diagnostic `A,S`; central step `1e-6`; acceptance `2e-10` absolute | `derivative_verification.json` | Software validation only; proxy outputs and panel ranks are non-physiological |

All 20 Phase 10 unit tests passed before commit `6248e1eb43a21cd89f39f9c3f79f57e9c013c0f6`. No validated result identifies a named physiological `Q`-plus-ion panel or claims global partial-observation identifiability.
