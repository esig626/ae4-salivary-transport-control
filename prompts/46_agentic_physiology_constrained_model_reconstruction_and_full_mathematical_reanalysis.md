# Task 46: physiology-constrained model reconstruction and full mathematical reanalysis

## Objective

Build a scientifically defensible replacement for the current reconstructed AE4 whole-cell model architecture that reproduces the key physiological AE4 phenotype for the right mechanistic reasons, rather than by arbitrary fitting, and then subject the accepted reconstruction to the same level of full-system mathematical analysis completed in Tasks 43–44.

This is an agentic scientific reconstruction task. It is not a cosmetic retuning exercise.

The accepted model must retain the conservation-explicit architecture: carbon, alkalinity, finite cell and lumen volumes, finite lumen composition, algebraic pH, electrical closure, water transport and AE4 regulation. Do not revert to the Na/K/Cl-only reduction in `Marty.tex`.

## Starting state

Read `AGENTS.md` first.

Then inspect and use, in order:

1. `analysis/task-43-parameter-provenance` at remote head `547f113d9123ab1976774639a69483faf40cef34`.
2. `analysis/task-44-full-system-mathematics` at remote head `e5fa9840bf96be3147c6118daee4942468c78f8b`.
3. Tasks 37–42 and their frozen outputs.
4. `analysis/43_parameter_provenance/report.tex` and `AE4_parameter_provenance.pdf`.
5. `analysis/44_full_system_mathematics/report.tex` and `AE4_full_system_mathematics.pdf`.
6. The current production model and manuscript only as read-only scientific context.

Create and work only on a new branch:

`analysis/task-46-physiology-constrained-model-reconstruction`

Do not modify `main`, manuscript text, Tasks 1–45, existing frozen outputs or the two completed analysis branches.

## Scientific premise established by Task 44

Treat the following as findings to explain, not defects to hide:

- stationary AE4-null secretion deficit of approximately 0.94617%;
- local NKCC chloride compensation of approximately 92.901%;
- finite-time replacement of missing AE4 chloride of approximately 88.507%;
- integrated NKCC increase of approximately 23.163%;
- WT/null slow modes of approximately 345.82/459.92 s;
- equal-routing cancellation does not imply zero total AE4 sensitivity;
- a physiological NBC-absent equilibrium can exist if AE4 loading is substantially reduced;
- Task 41 requires approximately 89.51% AE4-dependent stimulated CaCC recruitment in its prescribed inverse family to reach the selected 600 s phenotype threshold;
- that Task 41 crossing is not a global lower bound and is not independent validation.

The reconstruction must explain why the previous architecture overcompensates AE4 loss and must remove or constrain that overcompensation through biologically defensible mechanisms.

## Core principle

Do not ask “what parameter values make the knockout equal 35%?”

Ask instead:

1. Which model assumptions produce the excessive compensation?
2. Which of those assumptions are weakly sourced, phenomenological, inherited from another preparation, or contradicted/tensioned by salivary experiments?
3. What is the smallest literature- and physiology-consistent architectural repair?
4. Does that repair reproduce the principal phenotype while also satisfying independent observables that were not used to choose it?

If no defensible model can meet the constraints, report that result rather than manufacture one.

## Phase A: diagnose the current model before modifying it

1. Reproduce the Task 44 reference calculations exactly from the saved source and verify the numerical checkpoint.
2. Trace the full derivative/network pathway responsible for large NKCC compensation after AE4 loss.
3. Use the Task 43 provenance inventory and Task 44 sensitivity crosswalk to rank the assumptions responsible for compensation.
4. Separate:
   - published kinetic coefficients;
   - published model laws inherited from other systems;
   - salivary-specific measured quantities;
   - closure-derived quantities;
   - effective capacities/abundances;
   - regulatory multipliers;
   - numerical regularisation;
   - target-selected parameters.
5. Produce a concise machine-readable causal diagnosis linking each high-impact assumption to its effect on AE4-null compensation, chloride storage, sodium, pH, volume and flow.

Do not alter the model until this diagnosis is committed and pushed.

## Phase B: predeclare calibration and holdout evidence

Before testing new architectures, write a frozen `EVIDENCE_SPLIT.md` and machine-readable equivalent.

Use a conservative split that prevents circular validation.

### Calibration / construction constraints

The model may use the following to select a reconstruction:

- normal WT resting physiology and conservation;
- normal WT stimulated physiology;
- experimentally established qualitative direction that AE4 loss causes a substantial sustained salivary secretion deficit;
- the central AE4-null secretion magnitude only as an explicit calibration target if needed, labelled as such;
- literature-supported transporter stoichiometry, thermodynamic direction, localisation and known regulation;
- biochemical/electrophysiological constraints needed to define physically admissible kinetics.

### Holdout / validation evidence

Do not use these to choose or tune the accepted model unless a source review proves they are inseparable from the construction constraint; if that happens, document the change before fitting:

- temporal development of the AE4-null phenotype;
- the 5% AE4-expression trajectory;
- isolated NKCC1 assay behaviour in AE4/AE2 knockout preparations;
- intracellular sodium phenotype;
- intracellular chloride phenotype;
- resting and stimulated pH observations not already required for WT calibration;
- initial CaCC/TMEM16A behaviour;
- Task 42 Catalán mechanism-class outcomes;
- any independent flux/current or cell-volume observation not explicitly used for construction.

If the literature gives uncertainty ranges, preserve them. Do not turn point estimates into fake confidence intervals.

Commit and push the evidence split before model search.

## Phase C: literature- and mechanism-constrained candidate architectures

Investigate candidate repairs in a hierarchy from least to most structurally invasive.

At minimum test, where scientifically meaningful:

1. corrected/limited NKCC1 capacity or regulation using salivary evidence rather than unconstrained compensatory scaling;
2. corrected Na/K pump capacity or distribution if Task 43/44 shows it controls compensation;
3. corrected NHE1 and NBC effective capacities/activation where the current values are weakly constrained;
4. explicit saturation or thermodynamic limits missing from inherited laws;
5. physiologically justified coupling between AE4 state and one or more apical/basolateral processes;
6. slow regulatory or trafficking states only if demanded by independent timing evidence and not already explained by existing modes;
7. gland-level or ductal coupling only if a single-acinar-cell model cannot reproduce the organ-level phenotype without violating independent acinar observations.

Do not add arbitrary mechanism classes merely because they are convenient. Every added state, flux or coupling must have:

- a biological interpretation;
- literature or experimental motivation;
- thermodynamic/charge consistency;
- parameter provenance;
- an explicit statement of whether it is measured, inherited, derived, assumed or calibrated.

Prefer the smallest model that passes the evidence hierarchy.

## Phase D: selection rule for the accepted reconstruction

The accepted model must satisfy all of the following:

1. Exact conservation/charge/electrical closure checks.
2. Physiological WT rest.
3. Physiological WT 600 s stimulation.
4. A substantial AE4-null secretion deficit compatible with the construction target.
5. No unphysiological chloride depletion, sodium explosion, pH drift, volume failure or charge failure during the calibration horizon.
6. Compensation behaviour that is compatible with the independent NKCC evidence or, where protocol mismatch prevents direct comparison, does not rely on clearly extreme compensation without explicit qualification.
7. Better agreement with the holdout observations than the Task 40/44 architecture.
8. No hidden dependence on Task 41’s target-selected CaCC construction unless an independently justified apical coupling survives holdout testing.

Use a lexicographic scientific selection rule, not a single arbitrary weighted loss:

- first reject physically/physiologically invalid models;
- then reject models violating held-out qualitative facts;
- then prefer fewer new assumptions/parameters;
- only then compare quantitative phenotype fit.

If multiple models remain, retain the set and report non-identifiability rather than choosing one by aesthetic preference.

## Phase E: validation panel

For every surviving candidate and the final accepted model, run at least:

- WT rest;
- WT 600 s stimulus;
- 5% AE4 expression;
- AE4 null;
- relevant AE2 comparison if available;
- isolated or mimicked NKCC protocol if the experimental protocol can be represented honestly;
- the seven Task 42 source vectors under the accepted kinetic architecture, where meaningful;
- longer-time continuation beyond 600 s sufficient to detect delayed pH/volume/chloride failure;
- perturbations of the key uncertain parameters identified by Task 43/44.

Preserve failed cases. Do not rescue a candidate after seeing a holdout failure without returning to a newly versioned construction stage and documenting that the holdout has been consumed.

## Phase F: parameter provenance and robustness for the accepted model

Create an updated parameter inventory for all changed/new parameters.

For each changed parameter record:

- old value/law;
- new value/law;
- source/provenance;
- mathematical role;
- calibration target, if any;
- admissible range if genuinely supported;
- sensitivity of principal conclusions;
- whether the result survives reasonable variation.

Perform local sensitivity at minimum and global/structured robustness where justified by supported ranges. Do not invent ranges solely to enable a global sensitivity plot.

## Phase G: repeat the full Task 44 mathematical analysis on the accepted model

Perform the same standard of analysis used in Task 44, adapted to the new architecture.

This includes:

1. Nondimensionalise the full retained system without freezing carbon, alkalinity, volume, lumen composition or voltage.
2. Verify dimensional/dimensionless trajectory equivalence.
3. Re-establish carbonate closure and electrical closure results for the accepted model.
4. Derive the exact intracellular and whole-system conservation identities.
5. Re-derive the alkalinity storage/budget result and any capacity obstruction under the new architecture.
6. Re-derive the equal-routing or general AE4 source-vector projections; state explicitly which cancellations survive and which do not.
7. Solve constant-input equilibria on the exact neutral charge manifold.
8. Compute local Jacobians, eigenvalues, slow modes and nonnormal/transient qualifications.
9. Compute implicit-function equilibrium sensitivities and validate them against nearby independently solved roots.
10. Define chloride compensation correctly with stoichiometric factors, distinguishing local derivative compensation, finite-time replacement and relative integrated increases.
11. Reanalyse all seven Task 42 source classes under the accepted architecture if scientifically coherent.
12. If an AE4-dependent regulatory coupling remains, analyse it structurally and determine what is exact, local numerical, finite-time numerical and unresolved.
13. Re-run the clean numerical replay and verification suite, recording exact matching counts and all gate checks.

Do not claim global uniqueness/stability/positivity unless proved.

## Required scientific comparison

Produce a direct old-vs-new table comparing at least:

- WT physiology;
- WT chloride-loading partition;
- AE4-null stationary deficit;
- AE4-null 600 s cumulative deficit;
- 5% AE4 result;
- local NKCC chloride compensation;
- finite-time replacement of missing AE4 chloride;
- integrated NKCC increase;
- intracellular sodium/chloride response;
- pH range;
- volume range;
- slowest local modes;
- timing of secretion divergence;
- agreement/tension with isolated NKCC data;
- CaCC behaviour;
- number of new assumptions/parameters.

The final report must explain exactly what changed the phenotype and why.

## Deliverables

Under `analysis/46_physiology_constrained_model_reconstruction/` produce at minimum:

- `PROTOCOL.md`
- `EVIDENCE_SPLIT.md`
- machine-readable evidence split
- causal diagnosis of Task 44 overcompensation
- candidate-model ledger
- rejected-candidate ledger with reasons
- accepted model source/configuration
- updated parameter provenance inventory
- all calibration and holdout results
- robustness/sensitivity outputs
- exact mathematical derivations in TeX
- equilibrium/eigenvalue/sensitivity outputs
- old-vs-new comparison tables
- verification scripts and machine-readable verification record
- `report.tex`
- compiled standalone PDF report
- `FINAL_STATUS.md`

The standalone report must be suitable as the technical basis for later rewriting of the main BMB manuscript, but do not edit the manuscript in this task.

## Checkpoint discipline

At every safe checkpoint:

1. fetch remote branch state;
2. pull with `--ff-only` where direct Git is available;
3. independently verify the completed stage;
4. commit;
5. push immediately.

Never force-push.

Required remote checkpoints include at least:

1. reproduced Task 44 baseline + causal diagnosis;
2. frozen evidence split;
3. candidate architecture panel defined;
4. each major candidate-family evaluation;
5. accepted model frozen before mathematical reanalysis;
6. completed calibration/holdout validation;
7. full mathematical reanalysis;
8. parameter robustness/provenance update;
9. compiled report and final verification.

Do not wait until the end to publish work.

## Stop conditions

Do not stop for routine clarification. Work autonomously and make conservative choices.

Stop without declaring success if:

- no biologically defensible candidate can reproduce the calibration phenotype without violating holdout evidence or physiology;
- the phenotype can only be recovered by extreme unsupported couplings comparable to or worse than Task 41;
- available evidence is insufficient to distinguish multiple materially different models.

In those cases, publish the negative result and the smallest discriminating experiments required.

## Final response

Report:

- final branch SHA;
- accepted model architecture and exactly what changed;
- which data were calibration vs holdout;
- WT, 5% and null headline results;
- old vs new compensation metrics;
- holdout successes and failures;
- changed/new parameter provenance;
- mathematical conclusions classified as proved/exact, local numerical, finite-time numerical or unresolved;
- robustness limitations;
- whether one accepted model or an equivalence class survived;
- paths to `report.tex` and compiled PDF;
- verification/replay counts;
- confirmation that `main`, manuscript, Tasks 1–45, production source outside the Task 46 reconstruction and prior frozen outputs were untouched;
- confirmation that the branch was not merged.
