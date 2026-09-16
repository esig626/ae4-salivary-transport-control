# Frozen Task 50 effective-coupling benchmark

## Status

This document is binding project memory for every post-Task-50 mechanistic reconstruction.

The Task 50 model is the **accepted working phenotype benchmark**. It is deliberately classified as a **TARGET-CALIBRATED CONSTRUCTION**: it demonstrates the magnitude and stimulus/genotype specificity that the reconstructed network must be capable of producing, but it does not identify the molecular mechanism that produces that coupling.

The exact preserved implementation is archived at:

- branch: `archive/task-50-working-effective-coupling-benchmark`
- pinned commit: `b16c30094b95f61a73d8f1977cd58e79c7bb50f6`
- scientific publication boundary: `7cb032e8bb9372c910d0e5400712ea79beb23e8f`
- implementation: `src/modern_full_model/task50_effective_coupling.py`
- implementation blob SHA: `e3a79489d20530d7e977a3f6c37f08c222060444`
- controlling report: `analysis/50_minimal_beta_conditioned_effective_coupling/MINIMAL_EFFECTIVE_COUPLING_REPORT.md`
- controlling ledger result: R50D in `docs/MANDATORY_RESEARCH_LEDGER.md`

The archive branch is a preservation branch. Future scientific work must not rewrite it or repoint it.

## Frozen law

The effective apical chloride conductance factor is

`m50 = 1 - lambda * a_Ca * beta * (1 - e_AE4)`

with

`lambda = 0.89488127156712`.

Here:

- `a_Ca` is the inherited normalised calcium/channel activation;
- `beta` is the inherited effective beta/IPR protocol input;
- `e_AE4` is the inherited AE4 expression coordinate.

The multiplier changes only the existing apical chloride conductance before the inherited electrical closure. It does not directly multiply water or secretion, add a signalling state, change AE4 transport stoichiometry, or alter other transporter parameters.

## Exact nesting properties

The frozen implementation establishes:

- WT (`e_AE4 = 1`): exact parent model;
- REST (`a_Ca = 0`): exact parent model;
- AE4-null CCh-only (`beta = 0`): exact parent model;
- AE2 knockout with AE4 retained (`e_AE4 = 1`): exact parent model.

The beta-containing mutant comparison to the original Task 41 coefficient retains the authorised numerical precision qualification recorded in R50A-R50D. This qualification does not alter the phenotype benchmark below.

## Frozen phenotype benchmark

Under the inherited combined CCh+IPR protocol (`0.3 uM` CCh, `5 uM` IPR, beta occupancy `1`, stimulated calcium `0.25 uM`, 600 s):

| Case | Cumulative secretion 0-600 s | Cumulative deficit vs WT | Flow deficit at 600 s |
| --- | ---: | ---: | ---: |
| WT | `0.992524544081835 pL` | `0%` | `0%` |
| AE4 5% | `0.7626224586023684 pL` | `23.163365263893244%` | `20.108927720908575%` |
| AE4 null | `0.6921751966127494 pL` | `30.26115064458511%` | `27.038018389118477%` |

Across the saved cumulative observations from 60 to 600 s:

- AE4 5% deficit remains approximately `23.16-28.82%`;
- AE4-null deficit remains approximately `30.26-37.05%`.

Thus the effect is not confined to one cherry-picked time point.

The frozen phenotype record is:

`analysis/50_minimal_beta_conditioned_effective_coupling/output/inherited_phenotype.json`

and the frozen common rest hash is:

`a0b96ef16e5967de4dfe2bdabb882f2a28d478e211f920bc19cd735a5a143acb`.

## Interpretation boundary

The benchmark establishes the following engineering/scientific requirement:

> The reconstructed transport network needs an additional beta-conditioned effect of sufficient strength to convert the otherwise strongly compensated AE4 deletion into a substantial stimulated secretion deficit of the experimentally observed order.

It does **not** establish that AE4 molecularly recruits TMEM16A, that `lambda` is a biological constant, or that the Task 50 algebra is the true molecular mechanism.

The unresolved molecular identity is the reason to test a source-backed mechanistic replacement.

## Requirement for any successor mechanistic model

A future mechanistic reconstruction must be evaluated against this frozen benchmark while keeping the benchmark itself untouched.

The mechanistic model must **not inherit or secretly retain the Task 50 multiplier**. Instead, independently source-backed equations/parameters must generate the missing effective coupling through their own dynamics.

At minimum, a successful successor must demonstrate:

1. a substantial AE4-null combined-stimulus secretion reduction of the same order as the Task 50 benchmark and the approximately 35% experimental phenotype;
2. persistence over a broad stimulation window rather than a single selected instant;
3. preservation of WT secretion physiology and AE2 specificity to the extent supported by the primary data;
4. no arbitrary genotype-specific secretion multiplier or post-reveal retuning;
5. explicit comparison against the Task 50 benchmark values above;
6. clear reporting of whether it also improves the still-unresolved IPR-only AE4-null response and chronic resting Cl/pH phenotype.

The Task 50 values are a **benchmark/reference result**, not an independent dataset to be fitted exactly. The successor model succeeds scientifically if independently constrained physiology reproduces a comparably substantial and robust AE4-loss phenotype without the Task 50 effective multiplier.

## Preservation rule

Do not delete, rewrite, repoint, or reinterpret the archive branch as an independently validated molecular model. If a later mechanistic model fails, Task 50 remains the accepted proof-of-sufficiency benchmark showing quantitatively how much missing beta-conditioned network coupling is required by the current reconstructed architecture.
