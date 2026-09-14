# Task 37 WT-only NBC validation

`WT_NBC_VALIDATION_PASS_PARTITION_OUTSIDE_CONTEXT`

One WT CCh + isoproterenol trajectory completed to 600 s and passed every specified physical, conservation and secretion gate. The positive chloride-loading partition was **79.9248% NKCC1, 20.0752% AE4 and 0% AE2** over 60–600 s. NKCC1 is above the 65–75% experimental context band; this is a validation discrepancy, not a failed secretion or physiology gate. No parameter was tuned.

Work is confined to `codex/task-37-wt-nbc-validation`, prepared head `d925785512a33b241ec24843660ff3ef88458549`. `AGENTS.md` was read before `prompts/37_wt_nbc_validation.md`. All 82 prepared-head files were verified against their GitHub blob hashes. Publication uses the connected GitHub integration; this task does not merge or modify `main`.

The accepted Task 31 R09 checkpoint is `R09_WT_density_5`, with core-state SHA-256 `d77c1907e3e7124061721d606148a733f858cce2ef591c4e9f16e5cb7fc2469d`. Its 12 core coordinates were copied exactly, with the inherited resting regulatory coordinate equal to zero. NHE1 carrier amount remains `2.339370005697548e-05 fmol`; NBC capacity remains `0.11570913197464398 fmol/s`. AE4 expression is exactly `1.0`. The historical `AE4NA05` text in the R09 root/routing identifier is not an AE4-expression perturbation.

Phase A passed all nine selected existing tests: six focused NBC tests and three directly relevant Cha NHE1 regressions. One literal correction changed a stale import in `tests/test_task36_minimal_nbc.py` from the absent `task31_nhe1_mechanistic` module to the existing `task31_nhe1_repair` module. The focused tests were rerun once and passed. The historical expression-prediction test was excluded under the WT-only instruction. Pytest was unavailable in the runtime, so the exact selected test bodies ran through the standard-library unittest runner. The initial environment/import failures and successful rerun are retained in [focused_test.log](focused_test.log). No production model source was edited.

The corrected inward NBC cycle has conserved source `(Na, K, Cl, TIC, TA) = (+1, 0, 0, +2, +2)`, transports −1 charge equivalent into the cell, and contributes positive conventional cell-to-bath current `I_NBC = F J_NBC`. The full affinity includes `+V_b/(RT/F)`; forward/reverse current and flux signs and electrochemical reversal passed. Isolated TIC/TA identities and charge/current cancellation closed at round-off. Independent whole-cell source assembly differed from the wrapper RHS by at most `2.7755575615628914e-17 fmol/s`.

At REST, NBC activation, flux and current are exactly zero and the NHE1 multiplier is exactly one. At the central stimulus they are activation one and NHE1 multiplier 2.3. Direct REST nesting used no new root solve:

| Comparison | Observed difference | Required tolerance |
|---|---:|---:|
| All amount RHS rows | 0 fmol/s | ≤ 1e-12 fmol/s |
| Both volume RHS rows | 0 pL/s | ≤ 1e-14 pL/s |
| Apical, basolateral and transepithelial potentials | 0 V | ≤ 1e-12 V |
| Resting outflow | 0 pL/s | ≤ 1e-12 pL/s |

The complete RHS was bitwise identical. Both models returned resting outflow `0.0010757073853493032 pL/s`. Full checks and vectors are in [rest_nesting.json](../../results/37_wt_nbc_validation/rest_nesting.json).

The single production Radau integration used `rtol=1e-7`, amount and regulatory `atol=1e-10`, volume `atol=1e-12 pL`, and maximum step `2 s`. The fixed stimulus was CCh `0.3 µM`, isoproterenol `5 µM`, Ca `0.25 µM`, beta occupancy `1.0`, and the inherited AE4 beta/PKA and NKCC1 regulation. The inherited protocol is resting at exactly `t=0` and active for `0<t≤600`; its production integration starts at `1e-6 s` with the unchanged saved state. Low-level solver stepping allowed gates to be checked before each next step without restarting the integration.

There were 333 accepted Radau steps, 2336 RHS evaluations, and 934 checked states including every accepted endpoint and every integer-second dense-output point. States, observables, tracked amounts and volumes passed finite/positive checks at all monitored points. These are numerical trajectory checks, not a formal enclosure between samples.

| Intracellular readout | Monitored minimum | Monitored maximum | Gate |
|---|---:|---:|---|
| Na (mM) | 11.6361257 | 21.8257348 | < 40 |
| K (mM) | 106.522832 | 116.763209 | 50–200 |
| Cl (mM) | 54.0795196 | 60.3049669 | 30–80 |
| pH | 6.91000022 | 7.03482684 | 6.6–7.3 |
| TIC (mM) | 5.67174623 | 8.70517426 | > 0 |
| HCO3 (mM) | 4.90937455 | 7.70215442 | > 0 and < 100 |
| Cell volume (pL) | 1.39831725 | 1.4796298 | > 0 and < 3 |

All eleven inherited conservation diagnostics stayed within their declared production tolerances. The largest tolerance-normalized residual was 0.00598632255, from cell speciation alkalinity. The two NBC source diagnostic fields are physical sources (`−J_NBC` charge and `+2 J_NBC` carbon), and were checked against those values rather than incorrectly treated as zero residuals.

| Conservation diagnostic | Maximum absolute residual | Production tolerance |
|---|---:|---:|
| Apical current (A) | 2.6495629e-25 | 1e-20 |
| Basolateral current (A) | 2.30899713e-23 | 1e-20 |
| Cell charge rate (fmol/s) | 2.40363285e-13 | 1e-10 |
| Lumen charge closure (fmol/s) | 2.74866935e-15 | 1e-10 |
| Carbon accounting (fmol/s) | 2.77555756e-17 | 1e-10 |
| Water accounting (pL/s) | 4.33680869e-19 | 1e-12 |

| Secretion gate | WT result | Required threshold | Result |
|---|---:|---:|---|
| Mean outflow, 60–600 s | 0.00162624074525298 pL/s | > 0.00118327812388423 | PASS |
| Cumulative outflow, 0–600 s | 0.986995055625111 pL | > 0.70996687433054 | PASS |
| Outflow at 600 s | 0.00156842215184375 pL/s | ≥ 0.0010757073853493 | PASS |

Mean stimulated outflow was 1.511787 times resting outflow; cumulative outflow was 1.529219 times the resting 600-second reference.

Flux integrals use trapezoidal quadrature on the inherited one-second output grid. Positive chloride loading includes only `max(J_NKCC1,Cl,0)`, `max(J_AE4,Cl,0)` and `max(J_AE2,Cl,0)` before integration. AE2 remained a signed chloride exporter and contributed zero to the positive pool. CaCC export is cell-to-lumen; paracellular chloride return is lumen-to-bath.

| Integral over 60–600 s | Value | Unit |
|---|---:|---|
| positive_nkcc1_cl_loading | 155.723018967 | fmol |
| positive_ae4_cl_loading | 39.1139293575 | fmol |
| positive_ae2_cl_loading | 0 | fmol |
| signed_ae2_cl_flux | -0.464994757386 | fmol |
| nbc_cycles | 36.9765137919 | fmol |
| nbc_bicarbonate_equivalent_influx | 73.9530275839 | fmol |
| nhe1_flux | 3.90465701384 | fmol |
| cacc_cl_export | 205.923012665 | fmol |
| paracellular_cl_return | 76.3142356764 | fmol |
| water_outflow | 0.878170002437 | pL |

The positive basolateral chloride-loading pool was 194.836948324 fmol. AE4 supplied 20.0752% versus 1.8903% at REST, a 10.6202-fold increase in share. AE4 is therefore a substantial positive stimulated loader in this WT trajectory. No predeclared numerical AE4-share gate was imposed. NKCC1 exceeded the upper context boundary by 4.9248 percentage points. No inference about an AE4 knockout or reduced-expression phenotype follows from this WT-only run.

All 11 requested time points (0, 60, …, 600 s) and all 23 time/readout columns are in [wt_timeseries.csv](../../results/37_wt_nbc_validation/wt_timeseries.csv). Integrals are in [wt_integrated_fluxes.csv](../../results/37_wt_nbc_validation/wt_integrated_fluxes.csv), shares in [wt_chloride_partition.json](../../results/37_wt_nbc_validation/wt_chloride_partition.json), and complete gate evidence in [wt_verification.json](../../results/37_wt_nbc_validation/wt_verification.json).

Budget: **one intended and completed WT integration, zero numerical retries, zero stationary solves, zero optimisation calls, zero sweeps, zero genotype tests/integrations, one scientific agent/worker, and one thread in each BLAS pool**. Numerical execution used 3.740857 s of the 600 s ceiling. [budget.json](../../results/37_wt_nbc_validation/budget.json) records the counts; [source_verification.json](../../results/37_wt_nbc_validation/source_verification.json) verifies the unchanged model source. The executable audit and run scripts are retained alongside this report. No further scientific execution was performed. Task 37 ends here.
