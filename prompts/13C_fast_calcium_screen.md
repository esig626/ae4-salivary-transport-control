# Task 13C fast calcium screen

## Mission

Task 13B already established the WT result at stimulated calcium `0.10 uM`. Do not rerun that case merely to reproduce information already frozen in the repository.

This task asks the smallest useful remaining question:

> With the entire Task 13B WT model family frozen, does increasing stimulated calcium to `0.25 uM` or `0.50 uM` close the absolute WT secretion deficit?

The calculation must be designed to minimise runtime while preserving the scientific conclusion.

Work on branch `codex/task-13c-calcium-input-sweep`.

Read `AGENTS.md` first. Then read the final Task 13B analysis and freeze artifacts, especially:

* `analysis/13B_modern_full_model/final_answer.md`
* `analysis/13B_modern_full_model/dynamic_validation.md`
* `results/13B_modern_full_model/native_dynamic_contract_manifest.json`
* `results/13B_modern_full_model/native_dynamic_contract_group_gate.csv`
* `results/13B_modern_full_model/native_source_wt_roots.csv`

Do not edit anything under `archive/`.

Do not merge to `main`.

Do not draft a manuscript.

## Scientific freeze

All Task 13B model equations, parameters, resting roots, geometry, regulatory parameters, transporter parameters, hydraulic parameters, solver acceptance gates, WT observation envelope, and holdout firewall remain frozen.

The only scientific input that may change is stimulated calcium.

Use exactly these calcium values:

* existing frozen reference: `0.10 uM`, reuse saved Task 13B results and do not recompute it
* new screen point: `0.25 uM`
* new screen point: `0.50 uM`

Do not add intermediate calcium values and do not perform threshold refinement. If `0.25 uM` fails and `0.50 uM` passes, report only that the transition lies somewhere between those values. Do not launch another grid.

The primary protocol is WT `CCH_IPR` only. Do not run `CCH_ONLY`, `IPR_ONLY`, AE2 null, AE4 null, or any other genotype or stimulus arm in this task.

The strict AE4 null phenotype remains sealed. Do not open, parse, inspect, decode, simulate against, or otherwise access the held out target values.

## Runtime objective

The historical Task 13B dynamic contract was intentionally exhaustive. This task is not.

The primary screen contains exactly:

`10 frozen WT roots x 2 new calcium values x 1 regulatory member x 1 protocol arm x 1 screening solver = 20 trajectories`

Use only the preferred Task 13B regulatory member:

`R1_G125_P21_PROSE_TREFERENCE`

Do not run the 20 member regulatory ensemble during the primary screen.

### Parallel execution

The ten roots are independent and should be run in parallel.

Reuse the existing `ProcessPoolExecutor` pattern or an equivalent process based implementation. Expose a `--workers` argument. Use a sensible default based on available CPU cores, capped at ten workers because there are ten roots. Do not use nested process pools.

Document the recommended invocation with BLAS thread counts fixed to one per worker, for example:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m src.modern_full_model.run_calcium_fast_screen --workers 10
```

Do not hard code ten if fewer CPU cores are available.

### Screening solver

Use the existing declared `LOOSE_RADAU` solver for the 20 trajectory primary screen. Do not invent a new tolerance set unless the existing loose solver is demonstrably unusable.

The screening calculation is exploratory. It must still retain positivity, conservation, flow, sustainment, physiology, and WT scale diagnostics, but it need not pay production solver cost before the screen shows that a case matters.

Reuse the existing physical 600 second protocol and existing output logic. Do not shorten the biological protocol to save time.

### Production confirmation

Escalate selectively.

1. If `0.25 uM` produces a loose screen WT absolute flow pass for any root, or comes sufficiently close that solver accuracy could change the classification, rerun all ten roots at `0.25 uM` with `PRODUCTION_RADAU`.
2. If `0.25 uM` clearly fails by a large margin, do not production rerun all ten roots at `0.25 uM`.
3. Evaluate `0.50 uM` for all ten roots in the loose screen.
4. If `0.50 uM` passes for any root, or lies near the WT scale boundary, rerun all ten roots at `0.50 uM` with `PRODUCTION_RADAU`.
5. If `0.50 uM` clearly fails for every root, production confirm only the best flowing `0.50 uM` root to verify the negative endpoint.
6. Use `PRODUCTION_BDF` only as a final crosscheck on the representative scientifically decisive case or cases. Do not run BDF across the whole screen.

Define and record the rule used for "sufficiently close" before inspecting the screen classification. A reasonable default is within ten percent of the one SMG flow scale boundary, but use an existing repository convention if one is already defined and better justified.

The saved Task 13B Radau versus BDF agreement may be cited as justification for selective rather than exhaustive solver crosschecking.

## Do not rerun old work

Do not:

* recalibrate or re solve resting roots
* rerun the Task 13B `0.10 uM` production ensemble
* rerun nearby state robustness
* rerun the full regulatory ensemble during screening
* rerun both stiff solvers for every case
* perform a new model reconstruction
* add a new transporter mechanism
* change the calcium gate
* change the 600 second protocol
* change the WT 9 to 10 uL/min observation envelope
* change the one SMG geometry ceiling

Validate the inherited freeze using saved hashes and artifacts instead of regenerating them.

## Runtime profiling

Record runtime information because runtime itself is now an engineering concern.

For each new trajectory record, where available:

* wall clock seconds
* solver method
* `nfev`
* `njev`
* `nlu`
* number of output samples
* success or failure

Before launching the full 20 trajectory screen, run one representative `0.25 uM` loose Radau trajectory and record its runtime and solver counters. If this single trajectory remains unexpectedly expensive, profile it before launching the rest.

Do not immediately rewrite validated model equations for speed.

If further optimisation is genuinely required, identify the bottleneck first. In particular inspect the cost of repeated acid base speciation and numerical Jacobian construction. Any optimisation that changes the numerical implementation of the validated model must be checked against the existing implementation on representative trajectories before being used for a scientific conclusion.

Do not alter the physical equations merely to gain speed.

## Required outputs

Prefer a small new runner such as:

`src/modern_full_model/run_calcium_fast_screen.py`

that calls the existing Task 13B implementation rather than duplicating the full model.

Create:

### Analysis

* `analysis/13C_calcium_fast_screen/evidence_and_freeze.md`
* `analysis/13C_calcium_fast_screen/runtime.md`
* `analysis/13C_calcium_fast_screen/final_answer.md`

### Machine readable results

* `results/13C_calcium_fast_screen/calcium_screen.csv`
* `results/13C_calcium_fast_screen/calcium_confirmatory.csv`
* `results/13C_calcium_fast_screen/frozen_manifest.json`
* `results/13C_calcium_fast_screen/artifact_hashes.csv`

The screen CSV must contain one row for every new root and calcium combination and must include:

* root ID
* frozen parameter hash
* calcium value
* calcium gate activation
* solver label
* runtime and solver counters where available
* flow at 60, 120, ..., 600 seconds
* 0 to 600 second mean flow
* cumulative flow
* minimum and maximum post stimulation flow
* sustainment ratio
* shared WT observation scale interval
* one SMG geometry pass or fail
* endpoint intracellular Na, K, Cl, pH and volume
* endpoint membrane potentials if already available through the existing diagnostics
* numerical and conservation gate status

Reuse the saved Task 13B `0.10 uM` results in the final comparison table but do not place them in the new trajectory count as newly simulated data.

## Tests

Add focused tests proving:

1. the ten resting roots are byte identical or hash identical to Task 13B
2. all non calcium model parameters are unchanged
3. only `0.25 uM` and `0.50 uM` are newly simulated
4. the `0.10 uM` row is loaded from frozen Task 13B artifacts rather than recomputed
5. the primary screen uses only WT, `CCH_IPR`, `R1_G125_P21_PROSE_TREFERENCE`, and `LOOSE_RADAU`
6. the holdout code path is never invoked
7. parallel and serial execution give identical scientific outputs within existing numerical tolerances
8. production confirmation is triggered only by the predeclared decision rule
9. final artifacts are deterministic apart from wall clock runtime fields

## Interpretation

Answer these questions directly:

1. Relative to the existing `0.10 uM` result, how much does WT per cell secretion increase at `0.25 uM` and `0.50 uM`?
2. Does `0.25 uM` close the absolute one SMG WT flow deficit?
3. Does `0.50 uM` close it?
4. Does increased calcium cause another physiological or numerical gate to fail before the flow deficit closes?
5. Are conclusions qualitatively consistent across all ten roots?
6. How much wall clock runtime was saved relative to the originally proposed nine point Task 13C screen?
7. Is a matched SMG calcium measurement still decision critical after these two tests?

Do not make an AE4 phenotype claim. Task 13C ends with WT calcium only.

## Final classification

End with exactly one of:

`WT ABSOLUTE FLOW GATE PASSES AT 0.25 UM`

`WT ABSOLUTE FLOW GATE FAILS AT 0.25 UM BUT PASSES AT 0.50 UM`

`CALCIUM 0.50 UM REMAINS INSUFFICIENT TO CLOSE WT ABSOLUTE FLOW DEFICIT`

`CALCIUM INCREASE EXPOSES PHYSIOLOGICAL OR NUMERICAL FAILURE BEFORE FLOW CLOSURE`

If roots disagree, state that explicitly before giving the classification and do not hide the heterogeneity behind an average.

## Repository actions

Run the focused new tests and only the inherited regression tests needed to prove that the frozen Task 13B implementation remains intact.

Commit all new Task 13C fast screen code, analysis, and results to `codex/task-13c-calcium-input-sweep`.

Push the branch.

Do not merge to `main`.