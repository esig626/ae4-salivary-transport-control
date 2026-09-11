# Task 14B: five percent AE4 secretion check

## Mission

Task 14 established that all ten frozen AE4 continuation branches contain valid connected steady states at 5% AE4 activity, but its exact-zero rule prevented any stimulated AE4 trajectory from being run when continuation failed at zero.

That exact-zero requirement is not the scientific question in this task.

The question is:

> With the model otherwise completely frozen, how much does total stimulated secretion change when AE4 activity is reduced to the already validated 5% state?

The controlling completed Task 14 commit is:

`a2fde97e30d79601e7455d9f1d77d8d67609cb49`

Work on branch:

`codex/task-14b-five-percent-ae4-check`

Read `AGENTS.md`, the Task 14 final report, and the Task 14 continuation artifacts before running anything.

Do not edit `archive/` and do not merge to `main`.

## Scientific interpretation fixed in advance

Five percent residual AE4 activity is being used as a near-complete functional-loss model state. Do not call it literal genetic zero.

The primary result is descriptive, not a newly fitted acceptance test. Report the magnitude of the reduction and its robustness across the frozen model ensemble.

The experimental approximately 35% AE4-null reduction may be shown as biological context, but:

- do not fit to 35%;
- do not choose a root or calcium value because it is closer to 35%;
- do not change any model parameter;
- do not introduce a new observation scale;
- do not demand exact agreement with 35%;
- do not use an exact minute-by-minute experimental time signature as an acceptance criterion.

The model is not required to reproduce a plateau followed by a particular decline, an exact onset minute, or any other detailed shape for which the underlying biological mechanism is unresolved. Preserve time traces as diagnostic information only.

The central qualitative question is whether near-complete AE4 loss produces a substantial loss of secretion while the already established AE2 result remains essentially neutral.

## Hard freeze

Do not recalibrate, optimise, retune, repair or replace anything.

Keep exactly as inherited from Task 14:

- all ten WT roots;
- all equations and parameters;
- the R1 regulatory member `R1_G125_P21_PROSE_TREFERENCE`;
- NKCC1 stimulation and normalisation;
- calcium inputs `0.10`, `0.25`, `0.50 uM`;
- CCh + IPR protocol and 600 s duration;
- solver tolerances;
- the selected connected 5% AE4 resting state for each root.

Do not solve a new 5% root. Use the already retained selected state at expression `0.05` from the Task 14 continuation artifacts.

Do not try to reach AE4 expression zero in this task.

## Required simulations

For every one of the ten frozen roots and each calcium input

`0.10`, `0.25`, `0.50 uM`,

run the validated 5% AE4 resting state under the frozen stimulated protocol with production Radau.

This is exactly 30 new AE4 trajectories.

Use the matching frozen Task 14 WT production Radau trajectory as the denominator. Do not recompute or replace the WT denominator unless a file-integrity check proves it is unusable; if that happens, stop and report the problem rather than silently substituting a new denominator.

For each pair compute at minimum

`R_total = total_secretion_AE4_5pct / total_secretion_WT`

and

`D_total = 1 - R_total`.

Report the ratio and percentage reduction.

Also retain numerical, positivity and conservation diagnostics for the 5% trajectories.

Minute flows, cumulative curves, ionic states and pH may be recorded for scientific context, but they are not phenotype acceptance gates in Task 14B.

## Required summaries

Report:

1. the total secretion ratio and percentage reduction for every root/calcium pair;
2. range and median by calcium value;
3. overall range across all 30 pairs;
4. whether every 5% trajectory remains numerically valid;
5. whether the effect is consistent in direction across all roots;
6. whether AE4 loss is clearly much larger than the Task 14 AE2 effect of approximately 0.03% to 0.16%;
7. the experimental approximately 35% reduction only as contextual scale, not an exact target;
8. any major disagreement between the two AE4 routing families.

Do not discard any root.

## Existing runner

A small runner has already been added on this branch:

`src/modern_full_model/task14b_five_percent.py`

Audit it before use. Fix it only if required for correctness. It is intended to:

- read the retained Task 14 5% states;
- use the frozen WT production totals;
- run only the thirty missing 5% AE4 trajectories;
- write `results/14B_five_percent_ae4_check/five_percent_ae4_flow.csv`;
- write `results/14B_five_percent_ae4_check/summary.json`;
- write `analysis/14B_five_percent_ae4_check/final_answer.md`.

## Interpretation

Do not turn exact-zero continuation failure into the headline conclusion if the 5% functional-loss state already gives the biological result of interest.

Conversely, do not call the model successful merely because 5% is numerically valid. The magnitude of the predicted secretion change must be reported plainly.

The scientifically useful comparison is the general phenotype:

- AE4 near-complete loss: substantial secretion effect or not;
- AE2 loss: essentially no secretion effect.

No detailed time-course mimicry is required.

## Final actions

Run focused tests sufficient to verify that the runner does not alter the frozen model and that all thirty rows correspond to the declared roots/calcium values.

Commit and push all Task 14B results to the current branch.

Do not merge to `main`.
