# Task 27: inverse AE4-null phenotype repair and mechanistic diagnosis

Work on branch `codex/task-27-ae4-inverse-phenotype-repair`.

Scientific baseline is the current `main` lineage through Task 26, starting from commit `b9d4e769747e636ad79d44a23c64c729f31b506d`.

Read `AGENTS.md` first.

## Critical execution and publication restriction

This task is computation and analysis first.

Do **not** use the GitHub blob API. Do not manually create Git blobs or trees. Do not encode files as Base64. Do not create or upload `.npz`, `.tar`, `.gz`, `.zip`, bundles, trajectory archives, compressed results packages, or other binary artefacts. Do not upload trajectories individually.

Full trajectory arrays may exist locally if required for computation, but leave them local and untracked. The GitHub deliverable, if publication is attempted later, is text only.

If ordinary Git authentication or `git push` is unavailable, do not fall back to any alternate publication mechanism. Finish the computation, preserve the text outputs locally, and report the publication failure separately.

## Scientific objective

Interrogate the present routed-AE4 model rather than guessing another biological fix.

Determine the **smallest intervention to the existing model** that makes complete AE4 loss:

1. possess an admissible resting equilibrium under the inherited REST gates; and
2. produce the experimentally observed approximately 35% reduction in secretion.

Then determine exactly which altered mechanism causes that reduction.

This is an inverse diagnostic of model architecture. A successful parameter intervention is not automatically a biological explanation.

## Fixed phenotype target

For each retained root at the central stimulation level,

`r_i(theta) = Q_AE4-null(theta) / Q_WT`.

Target `r = 0.65` with diagnostic acceptable band

`0.60 <= r_i <= 0.70`.

The band is fixed before optimisation and must not be adjusted after seeing results.

## Important Task 26 baseline result

Use `analysis/ae4_exact_null/ae4_exact_null_report.md` and `results/ae4_exact_null/` as the prior diagnostic.

Task 26 found no admissible exact-null phenotype. Five null roots numerically closed but exceeded inherited pH/volume bounds, five remained unresolved. The closed formal subset increased secretion. Exact-null states showed severe TIC/HCO3 accumulation, alkalinisation and swelling, and NKCC1 replaced about 138-140% of lost AE4 chloride loading.

Do not reinterpret those formal out-of-domain secretion values as accepted physiological roots.

## Genotypes

Use only:

- WT AE4 expression = `1.0`
- exact-null AE4 expression = `0.0`

Do not rerun the 5% experiment as an optimisation arm.

AE4 itself is not an optimisation lever.

## Candidate intervention mechanisms

Freeze this set before optimisation. Each lever is a multiplicative scale on the existing implementation only:

1. NKCC1 capacity/activity
2. apical CaCC / chloride conductance
3. total K-channel conductance
4. Na/K pump capacity
5. AE2 capacity
6. NHE1 capacity
7. basolateral CO2 permeability
8. apical CO2 permeability
9. paracellular conductance
10. water permeability

Map each conceptual lever to the exact code parameter path and save the mapping in `contract.json`.

Do not introduce new equations or transport mechanisms.

Search bounds for every scale:

`0.05 <= theta_j <= 5.0`.

Optimise in log parameter space.

Water permeability and CaCC are output-proximal controls. Keep them in the diagnostic screen, but identify them explicitly as such in interpretation.

## Hard REST requirement

No phenotype objective may be evaluated from an unresolved or physiologically inadmissible null REST state.

For every candidate parameter point and every root, re-equilibrate exact-null REST before stimulation using inherited numerical and conservation gates.

A candidate is globally phenotype-feasible only if all ten roots have admissible null REST states and all ten central-calcium secretion ratios are in the target band.

Do not loosen the inherited pH, volume, positivity, conservation, current, charge or boundary gates to make a solution feasible.

Record failed and unresolved candidates. Never substitute the last solver iterate for a resting state.

## Stage 1: local mechanism screen

For each candidate mechanism, use symmetric small log perturbations around the current model to estimate local control over:

- REST solvability/admissibility;
- secretion at Ca = 0.25 µM;
- intracellular Cl, TIC, HCO3, pH and volume;
- NKCC1 chloride loading.

Where exact-null REST is unavailable at the unperturbed point, report the derivative as unavailable rather than fabricating it. Also record whether the perturbation moves the failed REST state toward or away from the inherited domain.

## Stage 2: single-mechanism inverse search

For every candidate mechanism individually, determine whether varying only that mechanism can produce ten admissible exact-null REST states and put all ten Ca = 0.25 µM null/WT secretion ratios inside 0.60-0.70.

Use a deterministic bounded one-dimensional optimiser and cache every evaluated scale. Do not rerun identical points.

Among feasible single-mechanism solutions, minimise `|log(theta)|`.

Validate every feasible solution at Ca = 0.10, 0.25 and 0.50 µM.

## Stage 3: sparse pair search

If no satisfactory single mechanism exists, rank mechanisms by REST rescue, secretion control and single-mechanism performance, then search all pairs among the six strongest candidates.

A pair is feasible only if all ten null REST states are admissible and all ten central-calcium ratios lie within 0.60-0.70.

Among feasible pairs minimise

`|log(theta_1)| + |log(theta_2)|`.

Do not prefer a numerically closer 0.650 fit once all ten roots are already inside the target band.

Validate feasible pair solutions at all three calcium levels.

## Stage 4: sparse multidimensional repair

Only if singles and pairs fail, search across the full candidate set.

Use a fixed-seed or deterministic optimisation method and retain every evaluated point.

Optimisation is lexicographic:

1. obtain admissible null REST states for all ten roots;
2. minimise phenotype-band violations at Ca = 0.25 µM;
3. among globally phenotype-feasible solutions, minimise the number of altered mechanisms;
4. among equally sparse solutions, minimise `sum_j |log(theta_j)|`.

Freeze the definition of an altered mechanism and all numerical tolerances before optimisation.

## Stage 5: evidence-constrained repair

Starting from the best sparse unconstrained repair, solve the same inverse problem with additional diagnostic constraints based on experimental observations:

- intracellular Cl should decrease substantially: diagnostic target `Cl_null / Cl_WT ≈ 0.73`, acceptable range `0.65-0.81`;
- resting pH should remain approximately unchanged: `|pH_null - pH_WT| <= 0.10`;
- NKCC1 activity should remain approximately unchanged: `0.85 <= NKCC1_null / NKCC1_WT <= 1.15`.

These are fixed diagnostic tolerances, not claimed experimental confidence intervals.

Retain all inherited REST admissibility requirements.

Determine whether an evidence-compatible phenotype repair exists within the fixed parameter bounds.

## Stage 6: causal decomposition

For the best unconstrained repair and, if it exists, the best evidence-constrained repair:

1. fully re-equilibrate REST;
2. run Ca = 0.10, 0.25 and 0.50 µM;
3. record the complete state/flux ledger;
4. individually revert each altered mechanism to baseline while keeping the others altered;
5. re-equilibrate and rerun stimulation.

This is the leave-one-mechanism-out causal ablation.

If the winning repair contains five or fewer active mechanisms, additionally calculate exact Shapley contributions to the secretion deficit across the active mechanisms.

Quantify how many percentage points of the secretion reduction are attributable to each mechanism and to interactions.

## Mechanistic ledger

For WT, untreated exact null where available, best unconstrained repair and best evidence-constrained repair, record at REST and during stimulation:

- intracellular Na, K, Cl, TIC, HCO3 and pH
- cell volume
- apical and basolateral voltage
- AE4 flux
- NKCC1 cycle flux and chloride loading
- NHE1 flux
- AE2 flux
- apical and basolateral CO2 flux
- CaCC chloride efflux
- K-channel fluxes
- Na/K pump flux
- paracellular flux
- water flow
- cumulative secretion

Explain the causal chain:

`intervention -> resting state -> stimulated transport -> osmotic driving -> fluid secretion`.

## Questions the report must answer

A. What interventions rescue admissible exact-null REST?

B. Which model mechanisms have the largest control over secretion and over the exact-null homeostasis failure?

C. Can any single mechanism produce the required approximately 35% secretion deficit while restoring all ten resting states?

D. If not, what is the sparsest successful combination?

E. Which altered mechanism actually causes most of the secretion reduction after causal ablation?

F. Can the model reproduce the secretion phenotype while simultaneously giving reduced intracellular Cl, nearly unchanged pH and approximately unchanged NKCC1 activity?

G. Does the minimal unconstrained repair contradict known experimental evidence?

H. Does the evidence-constrained inverse problem have a satisfactory solution?

I. Based strictly on these diagnostics, is the current failure more consistent with a missing parameter effect or structural/model-architecture misspecification?

## Required local text outputs

Create `results/ae4_inverse_phenotype_repair/` with at least:

- `local_sensitivity.csv`
- `single_mechanism_search.csv`
- `pair_search.csv`
- `multi_mechanism_search.csv`
- `evidence_constrained_search.csv`
- `search_evaluations.csv`
- `best_solutions.json`
- `causal_ablation.csv`
- `shapley.csv` if applicable
- `resting_states.csv`
- `trajectory_summary.csv`
- `flux_ledger.csv`
- `verification.json`
- `contract.json`
- `summary.json`
- scientific Markdown report

No binary trajectory files are part of the deliverable.

When computation and validation are complete, report the main scientific conclusion and paths to the text outputs. Do not use a publication failure as a reason to rerun simulations or create binary/archive fallback uploads.
