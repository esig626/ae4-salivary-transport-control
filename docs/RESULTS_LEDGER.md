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

## Phase 11 validated results

| ID | Statement | Status | Analysis / code | Input and tolerance | Pointer | Robustness and caveats |
| --- | --- | --- | --- | --- | --- | --- |
| R11-01 | Every one of the 67 files under `archive/legacy-2017/` is accounted for by content, role, dependency, chronology, and integrity evidence | Reproduced forensic inventory | `analysis/11_forensic_reconstruction/archive_audit.md`; [evidence commit](https://github.com/esig626/ae4-salivary-transport-control/commit/80921227fad5823e29ff12d522aeb72fbf602c2e) | Exact 67-path set; 65 byte-exact local SHA-256 checks; two connector exceptions tied to authoritative remote Git blobs/sizes | `archive_hashes.csv`; file-level ledger | The archived article payload and normalized log require remote Git-object evidence; local transport bytes are explicitly not substituted for them |
| R11-02 | `Par.mat` contains only one `1 x 1` struct with 45 scalar-double fields; all 45 exactly match `Par.m`, and Original calibration/transform formulas reproduce every stored field with maximum relative difference `2.51e-15` | Reproduced historical parameter lineage | `code_lineage.md`, `tools/dump_par_mat.py`, `tools/reproduce_historical_matlab.py`; [evidence commit](https://github.com/esig626/ae4-salivary-transport-control/commit/80921227fad5823e29ff12d522aeb72fbf602c2e) | Exact field/value/shape/type comparison; calibration-transform acceptance below `2.6e-15` relative | `par_mat_inventory.json`; `reproduction_summary.json` | The MAT header date is 2018-10-10; it stores no trajectories or hidden solver state and does not prove publication-generating identity |
| R11-03 | The Original historical calibration closes its 12 displayed residuals to maximum absolute `1.39e-17`; the later seven-state WT initial row has maximum raw RHS `5.24e-12` under its own conventions | Reproduced historical constructed baseline | `reproduction.md`, forensic evaluator and tests; [evidence commit](https://github.com/esig626/ae4-salivary-transport-control/commit/80921227fad5823e29ff12d522aeb72fbf602c2e) | Direct archived-expression arithmetic; test bounds `1.4e-17` and `5.3e-12` | `reproduction_summary.json`; `reproduction_residuals.csv` | Closure is achieved by sequentially fitting conductances/activities to the same row; HCO3, gates, water, acid/base, and other conventions differ from the published specification |
| R11-04 | Under the archived calcium-step protocol, translated seven-state endpoint flow ratios are `0.999052` for AE2 KO and `0.998556` for AE4 KO; the AE4 first post-step sample is `1.101877`. The five-state AE4 endpoint is `1.059019` and still moving | Reproduced translated-equation nonmatch to the central phenotype | `reproduction.md`, forensic evaluator and tests; [evidence commit](https://github.com/esig626/ae4-salivary-transport-control/commit/80921227fad5823e29ff12d522aeb72fbf602c2e) | SciPy BDF, archived `1e-6` tolerances and RHS scales; flow-ratio regression tolerance `2e-8` | `reproduction_summary.json`; `reproduction_trajectories.csv` | MATLAB/Octave was unavailable; this is translated execution, not native `ode15s`. Neither variant reproduces the published approximately 24% AE4-KO flow reduction, and the calcium/time protocols differ |
| R11-05 | Common model-family lineage is confirmed, but the exact archived snapshot is graded `UNLIKELY` to be the publication-generating 2018 implementation; the clean-reconstruction and physiological-identifiability gates fail, so the project remains `STOP` | New evidence-based forensic classification | `final_assessment.md`, `identity_evidence.md`, `paper_code_concordance.md`, `blocker_resolution.md`; [evidence commit](https://github.com/esig626/ae4-salivary-transport-control/commit/80921227fad5823e29ff12d522aeb72fbf602c2e) | Evidence matrix over equation/parameter identity, output reproduction, figures, chronology, and dependencies; prompt threshold `VERY LIKELY` not met | Final assessment decision table | `UNLIKELY` applies to exact snapshot identity, not common ancestry. No clean implementation, physiological sensitivity map, or manuscript was produced |

All 33 repository tests passed on the evidence state committed as `80921227fad5823e29ff12d522aeb72fbf602c2e`. The historical trajectory checks explicitly distinguish translated SciPy-BDF evidence from native MATLAB execution.

## Phase 12 validated results

| ID | Statement | Status | Analysis / code | Input and tolerance | Pointer | Robustness and caveats |
| --- | --- | --- | --- | --- | --- | --- |
| R12-01 | At the frozen historical WT row, exact closure requires the AE4 source to equal the Na-only `1:1:2` source vector. Of 768 Round-2 rows, only 10 close at `1e-8`, all effectively on the pure-Na boundary (`f_K <= 2.5718e-7`); none is source-admissible. | New exact obstruction plus reproduced enumeration | `analysis/12_ae4_mechanism_reconstruction/structural_explanation.md`, `model_comparison.md`; [scientific commit](https://github.com/esig626/ae4-salivary-transport-control/commit/1a154ad2e65833294aa240c08a9370208a630319) | Fixed seven-state WT row; amount-balance tolerance `1e-8`; 17 primary candidates and 768-row structural grid | `round2_envelope.csv`; `model_comparison.csv`; `summary.json` | Frozen-row theorem only. Moved roots are not excluded. C1's `0.144386%` decrement and the unsupported `~0.1476%` envelope are inadmissible controls, not maxima over an admissible AE4 family. |
| R12-02 | KO-blind continuation finds rest-compatible `2:1:3` C6, C7, and C8-SAT bands, but all 30 source-default, rest-screened positive-PKA rows admitted by the new volume sensitivity predict endpoint and integral `AE4-KO/WT > 1`. The minimum production ratios are `1.0137777584` and `1.0209390021`, so the maximum positive reduction on this grid is `0%`. | New tested-domain exclusion; independently reproduced | `analysis/12_ae4_mechanism_reconstruction/model_comparison.md`, `adversarial_audit.md`; commit `1a154ad2e65833294aa240c08a9370208a630319` | Source-default weights `1.5:1.6`, `K_eq=1`, `sigma=1`, bath HCO3 `21 mM`; finite `1e-8`-`1e4` capacity continuation with `0.025`-log-decade local refinement; WT Cl/pH within `+/-2` reported SEM; post-hoc `+/-10%` volume sensitivity; AE4-specific folds `1.25` and `3`; Radau/BDF `rtol=1e-7`, `atol=1e-8` | `rest_capacity_envelope.csv`; `pka_activation_envelope.csv`; `solver_crosscheck.json`; `trajectories.csv` | Independent BDF reproduces all 30 rows: maximum BDF/Radau differences `1.36e-9` endpoint and `3.37e-9` integral; conservation `<=4.98e-14`; thermodynamic violation `0`. The volume band is a new sensitivity, not experimental uncertainty; unsupported axes and disconnected roots are not globally excluded; physical time remains uncertified. |
| R12-03 | The AE4-null vector field is candidate/PKA-independent; its screened local branch predicts `Cl_i=48.7822 mM` and `pH_i=7.3785`, respectively `+7.68` and `+24.4` reported SEM from direct AE4-null targets. All 15 candidate-specific surviving `2:1:3` AE2-null roots also fail the direct AE2-null chloride target, with z-scores `-4.736` to `-4.401`. | New fixed-chassis genotype-rest falsification | `analysis/12_ae4_mechanism_reconstruction/model_comparison.md`, `adversarial_audit.md`; commit `1a154ad2e65833294aa240c08a9370208a630319` | Direct genotype targets; `+/-2` reported-SEM acceptance; deterministic 20-start AE4-null screen; AE4-null raw residual `5.16e-14`; 15 AE2-null roots close to at most `5.52e-13`; no knockout target used in calibration | `genotype_resting_equilibria.csv`; `genotype_resting_equilibria.json`; `summary.json` | The AE4-null result is a common screened local branch, not a proof of global equilibrium uniqueness. AE2-null roots are candidate-dependent. Near-unity acute AE2 secretion does not override the chronic resting-Cl contradiction. |
| R12-04 | On the fixed steady chassis, a flow ratio of `0.65` requires a `35.851%` apical chloride-flux loss. The exact cation/current identity localizes the missing support to pump plus K efflux. The Na/K-pump direction has `abs(cos(theta))=0.98058` but leaves a `19.612%` residual; pump-only and NHE1-only corrections each fail, while pump plus an independent Na-loading direction spans the frozen residual. | New exact pathway localization | `analysis/12_ae4_mechanism_reconstruction/structural_explanation.md`; commit `1a154ad2e65833294aa240c08a9370208a630319` | Fixed steady water and amount-balance equations; pathway signatures projected onto the frozen mixed-cation correction | `pathway_localization.csv`; structural equations 13-17 | Analytical localization, not a phenotype-tuned simulation. It does not establish a changed pump, K channel, or NHE1 protein, provide physiological flux scaling, or prove a two-pathway correction sufficient. |
| R12-05 | No declared C1-C8 implementation passes all source, WT/rest, genotype-rest, thermodynamic, AE2, conservation, complexity, and held-out AE4-flow gates. The controlling classification is `AE4 ALONE INSUFFICIENT; MISSING PATHWAY LOCALIZED`; identifiability is not reopened. | Independently reproduced and adversarially audited final classification | `analysis/12_ae4_mechanism_reconstruction/final_answer.md`, `adversarial_audit.md`; commit `1a154ad2e65833294aa240c08a9370208a630319` | Held-out AE4 ratio `0.65 +/- 0.047`, interpreted as `[0.603,0.697]` with endpoint and integral both required; new `+/-5%` qualitative AE2-dynamic sensitivity; genotype rest within `+/-2` reported SEM; zero complete hard-gate successes | `summary.json`; `model_comparison.json`; final-answer controlling decision | Exhaustive only over the declared fixed chassis and tested finite domains, not all AE4 biology. Absolute flow/time remain uncertified, and no physiological identifiability or manuscript claim follows. |

All 126 repository tests passed on the scientific state committed as `1a154ad2e65833294aa240c08a9370208a630319`. The Phase 12 tests cover the independent chassis, candidate registry and signatures, thermodynamics, moving/genotype resting gates, knockout-data isolation, source-bounded stimulation, solver agreement, generated-artifact schemas, and archive immutability. Passing supports the stated tested-domain exclusion; it does not prove global root uniqueness.

## Task 13B pre-reveal reconstruction results (current, incomplete)

These rows record source and WT-gate results only.  They do not constitute a
substantive Task-13B classification and do not authorize held-out reveal.

| ID | Statement | Status | Analysis / code | Input and tolerance | Pointer | Robustness and caveats |
| --- | --- | --- | --- | --- | --- | --- |
| R13B-P01 | The old G5 pre-reveal ensemble is `REJECTED_SUPERSEDED` under two independent hard failures: `ABSOLUTE_GLAND_GEOMETRY_SCALE_FAIL` and `NATIVE_CONDUCTANCE_SOURCE_MAP_FAIL` | Current controlling pre-reveal decision | `analysis/13B_modern_full_model/dynamic_validation.md`, `generation_log.md` | Immutable old gate SHA-256 `19eb2300...`; manifest `6a820b35...`; selected counts `294,706,403` and `260,704,116` | `g5_supersession.json`; immutable old G5 gate/manifest | The old files retain their historical internal eligibility field and are not rewritten; the separate supersession artifact makes current eligibility false |
| R13B-P02 | The native production conductance map uses literal whole-cell maxima `31.4 nS` CaCC and `14 nS` total K with the Palk gate `C_eff^1.46/(0.26^1.46+C_eff^1.46)` and zero K backgrounds | Source-map correction; replacement hierarchy complete through the preliminary absolute screen | `evidence_ledger.md`, `cation_topology.md`, `dimensional_ledger.md` | `C_eff` in `uM`; gates at `0.058/0.10/0.12/0.15 uM` are `0.10062056/0.19860330/0.24437206/0.30936970` | `native_source_map.json`; `native_source_panel.csv` | Rat-parotid-lineage transfer; `0.10/0.12/0.15 uM` are assumed WT spatial-average sensitivities, not direct matched-protocol calcium measurements and not genotype fits |
| R13B-P03 | A generous one-SMG geometry permits at most `51,469,231` modeled `1.30 pL` acinar cells; `9--10 uL/min` therefore requires `0.00291436--0.00323818 pL/s` per cell | Independent WT absolute-scale gate | `evidence_ledger.md`, `dynamic_validation.md`, `dimensional_ledger.md` | Kondo `66.91%` acinar fraction; `100 mg` gland; unit-density transfer | `g5_supersession.json` | Deliberately generous upper bound, not a point estimate; exact gland mass, collection efficiency, and local-unit map remain uncertain; paired-gland interpretation is sensitivity-only |
| R13B-P04 | Uniform pump surface density plus published rat-parotid areas gives an apical transfer `0.075075`, with `0.056878/0.094586` sensitivity; K `0.20/0.30/0.40` and AE4 Na `0.05/0.10/0.20` are assumption panels | Source/assumption separation | `evidence_ledger.md`, `cation_topology.md`, `wt_calibration.md` | Conservative low/high area combinations; total K fixed at `14 nS` | `native_source_map.json`; `native_source_panel.csv` | Cross-gland transfer and uniform-density assumption remain explicit; “even distribution” is not interpreted as 50% of total capacity |
| R13B-P05 | The frozen native regulatory contract contains dynamic beta/cAMP/PKA/AE4 members with R1 `10/30/90 s`, R2 `(3,10)/(10,30)/(30,90) s`, and symmetric R3 isomorphs | Mandatory 2021-evidence interface retained; kinetics unidentified | `camp_pka_model.md`, `dynamic_validation.md`, `dimensional_ledger.md` | Capacity sensitivities `1.25/1.70`; exact R2/R3 nesting where symmetric | `src/modern_full_model/native_dynamic_contract.py`; `tests/test_modern_native_dynamic_contract.py` | The 2021 study establishes pathway dependence and construct logic but no acute kinetic trace; R0 is a negative control and cannot alone satisfy the dynamic-regulation requirement |
| R13B-P06 | A beta/volume-sensitive anion-exit diagnostic `gV*beta*max(Vi/Vrest-1,0)` never opens on the tested IPR-only trajectory: volume changes `1.3000 -> 1.273730 pL` with zero swelling | Rejected quantitative rescue; required topology remains unresolved | `evidence_ledger.md`, `dynamic_validation.md`; `src/modern_full_model/vbeta_diagnostic.py` | Exact-rest WT input; R1 `tau=30 s`; `Ca=0.058 uM`; any nonnegative `gV` | `vbeta_native_exact_rest_diagnostic.csv`; `vbeta_native_exact_rest_summary.json` | Catalan et al. license a distinct beta/volume-sensitive exit topology but not an absolute `gV`; even an exogenous swelling diagnostic cannot license a production amplitude |
| R13B-P07 | The complete canonical H1 checkpoint has 45 native/pump-sensitivity panels, 135 attempts, 35 numerical roots, and zero WT-passing roots; every root fails independent WT Cl and K gates | Completed production-root generation result; not a final Task-13B classification | `generation_log.md`; `src/modern_full_model/native_source_panel.py` | 15 native-pump panels plus 30 conditional pump-capacity `0.5/2` panels; WT ionic gates fixed before scanning | `native_source_root_attempts.csv`; `native_source_wt_summary.json` | Some roots also fail Na and/or pH; the later absolute tier is recorded separately in R13B-P09 |
| R13B-P08 | HW and HWQ each screen 225 diagnostic-only panels/675 attempts, produce 175 numerical roots, and produce zero WT-passing roots; combined H1+HW+HWQ is 495 panels/1,485 attempts/385 roots/zero passes | Completed negative hydraulic diagnostic | `generation_log.md`, `wt_calibration.md`, `dynamic_validation.md`, `independent_reproduction.md` | `H=2.5/3.5/5/6.75/8`; HW scales `Pa/Pb/Pt`, HWQ also scales outflow; scale-one forms are exact H1 aliases | `native_source_panel.csv`; `native_source_wt_summary.json` | At each H, HWQ matches HW displayed 35-root intracellular ranges (Na `8.107--29.441`, K `70.840--92.515`, Cl `23.589--29.243 mM`, pH `6.7522--6.9440` over the panel); numerical diagnostic, not an analytic identity; H>1 cannot rescue production |
| R13B-P09 | The complete three-start absolute-source tier has 630 panels/1,890 attempts/354 numerical roots; `N_ABS_NKCC` contributes 229 roots and ten preliminary WT passes, all at scale 4, whereas `AN_ABS_AE4_NHE` contributes 125 roots and zero passes | Completed preliminary production screen; confirmation pending | `generation_log.md`, `wt_calibration.md`, `dynamic_validation.md`; `src/modern_full_model/native_source_panel.py` | N counts by scales `0.5/2/4/5/6/7/8`: `40/35/33/31/30/30/30`; AN: `45/25/20/10/10/10/5`; 1,847 optimizer-success, 1,890 admissible, 981 numerical attempts; 276 panels rootless | `native_source_panel.csv`; `native_source_root_attempts.csv`; `native_source_wt_roots.csv`; `native_source_wt_summary.json` | All ten passes cover all five topologies at pump 1/`AE4NA05` and pump 2/`AE4NA20`; Cl `51.0172--53.0228 mM` (minimum margin `0.07723 mM`), pH `6.83860--7.04827` (minimum margin `0.001728`), Na `14.3306--18.9814`, K `105.935--118.116 mM`; residual `1.01e-12--2.41e-11`, rank 11/nullity 0, no boundary hits. These rows are pre-seven-start/pre-33-start and not dynamics- or reveal-eligible. |
| R13B-P10 | The combined native hierarchy through the absolute screen has 1,125 panels/3,375 attempts/739 numerical roots and ten preliminary WT passes; promotion is paused at a code-audit engineering gate | Current execution-integrity disposition; no scientific retuning | `generation_log.md`, `independent_reproduction.md`; `src/modern_full_model/run_native_source_panel.py`; `tests/test_modern_native_source_panel.py` | Panel-atomic rerun replacement; both quartiles retained in three-/seven-start designs; duplicate ID and eligibility fail closed; exact scale-one alias coverage; complete upstream panel-ID checks; 33-start intake rejects duplicates, failed/ineligible rows, and unconfirmed rows; no first-passer pruning | `native_source_panel.csv`; `native_source_wt_summary.json`; `native_source_scale_one_identity.csv` | The audit changes persistence/selection safeguards only, not equations, source-class scales, WT gates, parameter bounds, or calibration. Seven-start confirmation and 33-start geometry remain incomplete; `eligible_dynamic_root_ids=[]` and held-out reveal remains prohibited. |

| R13B-F01 | All ten N-scale-4 survivors pass seven-start confirmation and 33-start geometry with one retained branch per panel; independent equation transcription agrees for all 739 numerical roots and all ten WT gates | Final native-root geometry result | `wt_calibration.md`, `independent_reproduction.md` | 33 starts per survivor; full-rank 11-coordinate Jacobians; exact table hashes | `native_source_wt_summary.json`; `independent_native_summary.json` | Root robustness does not imply a valid stimulated WT model |
| R13B-F02 | The native WT contract has 800 production and 320 nearby cases; every numerical, solver, conservation, sustainment, co-stimulation, and regulatory gate passes, but zero of ten roots passes the absolute one-SMG scale | Final pre-reveal WT dynamic result | `dynamic_validation.md`, `heldout_validation.md` | Required scale `6362.31--8210.72`; one-SMG ceiling `3088.15386` | `native_dynamic_contract_gate.json`; `native_dynamic_contract_group_gate.csv` | This shared failure blocks genotype evaluation and holdout reveal |
| R13B-F03 | The complete diagnostic high-capacity tier produces zero additional WT-rest-passing roots | Targeted-repair exclusion | `generation_log.md`, `wt_calibration.md` | Predeclared N and AN stress scales `12/16/24/32`; no phenotype access | `native_source_wt_summary.json`; `native_source_panel.csv` | Diagnostic tiers cannot rescue production and were not promoted |
| R13B-F04 | Final classification is `FULL MODEL NOT IDENTIFIABLE FROM EXISTING DATA — ONE DECISION-CRITICAL MEASUREMENT SPECIFIED`; heldout reveal count remains zero | Controlling Task-13B outcome | `final_answer.md`, `heldout_validation.md`, `adversarial_audit.md` | Decision threshold `0.00291436--0.00323818 pL/s` per modeled 1.30-pL cell | `final_classification.json`; `reveal_log.json` | The result is scoped to the tested source classes; it does not argue against AE4 biology or the historical article |

## Task 50A formal audit qualification

**R50A-01:** The generic beta conditioned construction has formal WT, REST,
AE4 null CCh only and AE2 knockout nesting. Its mandated
`lambda=0.89488127156712` does not have the exact complement of the actual frozen
Task 41 `b=0.10511872843288446`. Their factors differ by
`-4.46e-15*a*(1-e)` at beta one under exact decimal arithmetic. The displayed
report value `b=0.10511872843288` does satisfy the requested identity, but is a
shortened coefficient.

Status: new provenance qualification, not a numerical prediction. Scope:
FORMAL DEDUCTION concerning a TARGET-CALIBRATED CONSTRUCTION. Evidence:
`analysis/50_minimal_beta_conditioned_effective_coupling/EQUIVALENCE_PROOF.md`
and `output/coefficient_audit.json`; start commit
`8f5fefd9fc0b563a6bd963895211b4acf76972b0`; frozen Task 41 publication
`86f135e075aa91a70a752b650d7bb894353fee04`. Seven frozen execution hashes pass.
There is no numerical tolerance in the exact literal comparison. No model
evaluation or production calculation was performed and no biological
significance is inferred. Task 50 stops at its prescribed exact identity gate;
the master ledger R50A records the qualification and preserves prior results.

## Task 50B reduced construction implementation

**R50B-01:** Under the user's explicit numerical treatment of the coefficient
rounding difference, ten tests pass for the fixed Task 50 law. There are 21
bit identical parent comparisons and six beta one fixed state comparisons to
frozen Task 41. Largest factor/RHS differences are `4.5103e-15` and `8.0908e-15`
respectively, with unit specific tolerances predeclared in D50-07. Current,
charge, chloride source and unrelated fixed state checks pass. This is software
verification of a TARGET-CALIBRATED CONSTRUCTION, not biological validation or a
global trajectory error bound. Evidence: Task 50 `IMPLEMENTATION_VERIFICATION.md`
and `output/verification_attempt_02.json`. Scientific source hashes are in that
record. No production trajectory, stationary solve or fit was run. The earlier
governance hash setup failure is preserved separately with zero evaluations.

## Task 50C qualified reuse of the frozen phenotype

**R50C-01:** Twenty eight artefacts match original Task 41 publication
`86f135e075aa91a70a752b650d7bb894353fee04` byte for byte. Reused cumulative
deficits are `23.163365263893244%` at 5% AE4 and `30.26115064458511%` at null,
with endpoint flow deficits `20.108927720908575%` and `27.038018389118477%`.
Status: inherited TARGET-CALIBRATED CONSTRUCTION, under the authorised numerical
treatment of coefficient rounding and the verified 50B comparisons at
`537ce3204d11a1a6cc51b1dd39032d12cb42df5e`. Saved observations establish finite
window persistence only. Evidence: Task 50 `FROZEN_RESULT_REUSE.md` and its two
machine readable output records. No new model evaluation or integration was
run; no global trajectory error bound, independent validation, chronic rest
repair, IPR only uptake repair or microscopic mechanism is established.

## Task 50D final reduced construction

**R50D-01:** Task 50 is complete with the user authorised numerical precision
qualification, preserving exactly `lambda=0.89488127156712`. Exact parent nesting
and the ten passing software tests support one fixed beta conditioned scalar.
Twenty eight frozen artefacts supply the inherited 23.1634% and 30.2612%
cumulative deficits, with finite window persistence only. Classification:
TARGET-CALIBRATED CONSTRUCTION. No new production trajectory, stationary solve
or fit was run; chronic resting adaptation, IPR only uptake and microscopic
identity remain unresolved. The original exact coefficient discrepancy remains
true and qualified reuse is not a global trajectory error bound.

Evidence: Task 50 `MINIMAL_EFFECTIVE_COUPLING_REPORT.md` and `output/` receipts.
Implementation was verified at `537ce3204d11a1a6cc51b1dd39032d12cb42df5e`;
reuse was verified at `88997ef6888498e6ea6bdd9bd4acd7e5ec108dbc`. Tolerances
remain those predeclared in D50-07. Master ledger R50D is the controlling result.
Stop at 50D; Tasks 51 and 52 remain unstarted.
