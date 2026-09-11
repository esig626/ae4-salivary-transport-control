# Task 15: AE4 sign and mixed-cation slip audit

## Mission

Task 14B produced a robust but biologically wrong-direction result: reducing AE4 expression to 5% increased total 0 to 600 s secretion in all 30 retained root/calcium cases.

Before proposing any missing biology, compensation, retuning, or alternative transporter mechanism, determine whether the modern model uses the correct transport signs and whether its reversible shared Na/K AE4 construction generates a large mixed-cation slip mode under physiological conditions.

The controlling completed Task 14B commit is:

`4d4403f279fc36aec940c9e4759486d02f07e3f6`

Work on branch:

`codex/task-15-ae4-sign-slip-audit`

Read `AGENTS.md` first, then this prompt, then the Task 13B AE4 module/evidence audit, Task 14 final report, and Task 14B report and saved trajectories.

Issue: #9.

Do not edit `archive/` and do not merge to `main`.

## Scientific status

The AE4 phenotype is already revealed. Nothing in this task is blind.

The observed experiment reports lower secretion after AE4 loss. Task 14B instead found an 8.23% to 21.66% increase at 5% AE4 expression. That wrong-direction result is the reason for this audit, but it must not be used to choose signs, tune parameters, select roots, or select calcium values.

Task 15 is diagnostic only.

Do not repair the model in this task.

If a sign error is found, document it and stop before changing the scientific model. If no sign error is found but the shared-carrier architecture produces large opposing Na and K cycles with small net chloride loading, document that result and leave any no-slip or alternative-cycle model to a later task.

## Hard freeze

Base all scientific evaluation on exactly commit `4d4403f279fc36aec940c9e4759486d02f07e3f6`.

Do not recalibrate, optimise, retune, refit, change a root, change a calcium value, change a regulatory model, change any transporter capacity, alter a source sign, change a stoichiometry, remove a branch, suppress slip, or introduce genotype-specific compensation.

Retain all ten roots and calcium values:

`0.10, 0.25, 0.50 uM`.

Retain regulatory member:

`R1_G125_P21_PROSE_TREFERENCE`.

Retain the 5% AE4 states and trajectories from Task 14B and the matching WT trajectories from Task 14.

Do not attempt exact AE4 zero in Task 15.

## Stage A: provenance and sign truth table

Before any new diagnostic calculation, freeze a Task 15 audit manifest containing:

- controlling commit and tree;
- hashes of every source file used in the audit;
- hashes of Task 14 WT trajectories and Task 14B 5% AE4 trajectories;
- the ten retained root IDs;
- the three calcium values;
- the AE4 expression values being compared, WT 1.0 and reduced 0.05;
- a declaration that no model parameter or scientific source file will be changed.

Construct an explicit sign truth table for every transport/water component relevant to secretion.

At minimum audit:

1. AE4 Na branch;
2. AE4 K branch;
3. net AE4 chloride source;
4. AE4 bicarbonate/TIC/alkalinity source;
5. NKCC1;
6. NHE1;
7. AE2;
8. apical CaCC chloride flux;
9. apical and basolateral K flux;
10. Na/K pump on each membrane;
11. paracellular Na, K, Cl and HCO3 fluxes;
12. basolateral water flux;
13. apical water flux;
14. paracellular water flux;
15. lumen outflow.

For each quantity report:

- positive mathematical direction;
- biological direction represented by that positive sign;
- intracellular source sign;
- luminal source sign where applicable;
- exact code path where the sign is created;
- exact code path where it is inserted into the whole-cell RHS;
- whether the direction agrees with the retained literature evidence or is only a modeling convention.

Do not infer correctness merely because conservation tests pass. A consistently reversed electroneutral cycle can conserve charge perfectly while being biologically backwards.

## Stage B: independent AE4 reaction-direction derivation

Independently derive the source vector from the declared positive AE4 reaction

`Cl_o + C_i + 2 HCO3_i <-> Cl_i + C_o + 2 HCO3_o`, with `C = Na or K`.

Verify algebraically that a positive branch flux implies:

- intracellular Cl increases;
- intracellular cation decreases;
- intracellular HCO3/TIC/alkalinity decreases;
- net transported charge is zero.

Then verify that the implementation in `src/modern_full_model/transporters.py` and its insertion in `src/modern_full_model/model.py` reproduce exactly that source vector with no additional sign reversal or double scaling.

Explicitly audit expression scaling. Distinguish raw full-capacity branch values retained inside AE4 diagnostics from the actual expression-scaled source applied to the whole-cell model. Confirm that WT uses scale 1.0 and Task 14B uses scale 0.05 exactly once.

## Stage C: reversal tests independent of the fitted roots

Add focused unit tests using constructed concentration states, not phenotype-selected states.

Required AE4 tests:

1. global equilibrium gives zero Na, K and net Cl flux;
2. a constructed state strongly favouring the declared positive Na cycle gives positive `J_Na`, positive intracellular Cl source, negative intracellular Na source and negative intracellular HCO3 source;
3. reversing the same thermodynamic gradient reverses all those source signs;
4. the equivalent pair of tests for the K branch;
5. branch affinity sign, branch flux sign and entropy production are thermodynamically consistent;
6. expression scaling by 0.05 preserves sign and scales the actual whole-cell AE4 source by exactly 0.05;
7. the whole-cell RHS receives the AE4 source with no second sign inversion.

Add analogous minimal direction tests for NKCC1, NHE1, AE2, CaCC, Na/K pump, K channels and water/outflow where existing tests do not already prove the source direction.

Do not change production equations to make these tests pass. If a test exposes a genuine sign error, preserve the failing evidence, classify it, and stop before repair.

## Stage D: full frozen-ensemble AE4 branch audit

Using the already saved WT and 5% AE4 trajectories, evaluate the frozen model at every retained time point. Do not rerun the ODEs unless an integrity failure makes the saved trajectory unusable; if so, stop and report the integrity failure rather than silently regenerating data.

For every root, calcium value, genotype state and time point record at minimum:

- `affinity_na`;
- `affinity_k`;
- actual expression-scaled `J_na`;
- actual expression-scaled `J_k`;
- `J_total = J_na + J_k`;
- actual intracellular Na source from AE4;
- actual intracellular K source from AE4;
- actual intracellular Cl source from AE4;
- actual intracellular HCO3/TIC/alkalinity source from AE4;
- AE4 entropy production;
- AE4 capacity/regulatory multiplier;
- expression scale;
- intracellular Na, K, Cl, HCO3, pH and volume;
- bath Na, K, Cl and HCO3 used in the affinity calculation.

Report values at rest, onset/right-limit, 60, 180, 300 and 600 s in a compact table, while retaining the full time-grid data machine-readably.

## Stage E: quantify productive transport versus cation slip

For each time point define

`branch_turnover = abs(J_na) + abs(J_k)`.

If `branch_turnover > 0`, define

`productive_fraction = abs(J_total) / branch_turnover`

and

`cancellation_fraction = 1 - productive_fraction`.

When `J_na * J_k < 0`, also define

`slip_pair_magnitude = min(abs(J_na), abs(J_k))`

and report the direction of the implied Na/K exchange.

Interpretation discipline:

- `J_total > 0` means net AE4 chloride loading under the declared sign convention;
- `J_total < 0` means net AE4 chloride unloading;
- large `cancellation_fraction` means opposing branch cycles dominate the total branch turnover;
- do not call this cation slip experimentally established physiology; it is a property of the current model architecture;
- do not claim that slip causes the wrong secretion phenotype unless a later counterfactual mechanism test demonstrates causality.

For all 30 WT cases summarize:

- fraction of time points with Na and K branch fluxes of opposite sign;
- fraction of time points with net AE4 Cl loading versus unloading;
- median and range of `cancellation_fraction`;
- median and range of `productive_fraction`;
- integrated positive/negative branch turnover over 0 to 600 s;
- integrated net AE4 Cl source over 0 to 600 s.

Repeat the same summaries for the 30 5% AE4 cases using actual expression-scaled fluxes.

## Stage F: compare against the retained experimental constraints

Use only claims already captured in the repository evidence ledger and primary-source audit.

At minimum distinguish:

- 2016 evidence that AE4 transports Na and K and behaves reversibly under manipulated ionic gradients;
- 2025 evidence that Na-supported and K-supported function can respond differently to mutations;
- evidence for macroscopic electroneutrality;
- absence of a direct physiological mixed-Na/K experiment establishing large simultaneous opposing Na and K branch cycles in native salivary acinar cells.

Do not upgrade pure-cation assays into evidence for the model's mixed-bath slip cycle.

Compare with the 2018 model only at the level of sign and coupling structure. Determine whether the older implementation permitted independent opposite-sign Na and K AE4 cycles or instead used one net AE4 exchange direction with cation partitioning.

## Stage G: limited downstream sanity check

This is not yet the full causal mechanism diagnosis, but report the already available WT versus 5% changes that make the sign result interpretable:

- resting and endpoint intracellular Cl;
- resting and endpoint Na;
- resting and endpoint K;
- resting and endpoint pH;
- resting and endpoint cell volume;
- net 0 to 600 s secretion ratio already established in Task 14B.

Where feasible without a new ODE run, also report at rest and endpoint the frozen NKCC1, AE2, NHE1 and apical Cl flux signs for WT and 5% AE4.

Do not attribute causality in this task.

## Required classifications

Choose exactly one primary classification:

1. `CLERICAL OR ALGEBRAIC AE4 SIGN ERROR FOUND`
2. `WHOLE-CELL AE4 SOURCE MAPPING OR EXPRESSION-SCALING ERROR FOUND`
3. `SIGNS CORRECT; AE4 NET LOADS CL BUT LARGE OPPOSING NA/K SLIP IS PRESENT`
4. `SIGNS CORRECT; AE4 RUNS NET REVERSE IN THE FROZEN PHYSIOLOGICAL STATES`
5. `SIGNS CORRECT; LARGE MIXED-CATION SLIP IS NOT PRESENT`
6. `SIGN OR SLIP AUDIT INCONCLUSIVE`

If classification 1 or 2 applies, do not repair the model in Task 15. Identify the smallest exact defect and the files/functions involved, then stop.

If classification 3 applies, explicitly state that the current shared-carrier architecture permits a large model-generated Na/K exchange cycle that is not established by the available mixed-bath experimental evidence. Do not yet replace it.

If classification 4 applies, state clearly that the implemented thermodynamics place AE4 on the opposite net transport side under the frozen physiological states.

If classification 5 applies, the wrong secretion sign must be sought elsewhere in a later task.

## Required artifacts

Create at minimum:

- `results/15_ae4_sign_slip_audit/audit_manifest.json`
- `results/15_ae4_sign_slip_audit/sign_truth_table.csv`
- `results/15_ae4_sign_slip_audit/ae4_branch_timecourse.csv`
- `results/15_ae4_sign_slip_audit/slip_summary.csv`
- `results/15_ae4_sign_slip_audit/downstream_state_summary.csv`
- `results/15_ae4_sign_slip_audit/test_results.txt`
- `analysis/15_ae4_sign_slip_audit/sign_derivation.md`
- `analysis/15_ae4_sign_slip_audit/evidence_comparison.md`
- `analysis/15_ae4_sign_slip_audit/final_answer.md`

The final answer must state plainly:

- whether any sign is actually wrong;
- whether AE4 is net loading or unloading Cl in the frozen WT states;
- whether Na and K branches oppose one another;
- how large the cancellation/slip is;
- whether that behaviour is experimentally established, merely permitted, or contradicted;
- what the next task should test, without performing the repair here.

## Final actions

Run focused tests and any existing inherited tests needed to establish that Task 15 added diagnostics only and did not change the scientific model.

Commit and push all Task 15 work to `codex/task-15-ae4-sign-slip-audit`.

Do not edit `archive/`.

Do not merge to `main`.

Do not draft manuscript text.
