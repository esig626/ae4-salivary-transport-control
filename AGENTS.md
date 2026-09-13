# Active phase: Task 33

Work only on `codex/task-33-nkcc70-chloride-partition`.
Read and execute `prompts/33_nkcc70_chloride_partition.md`.

This task calibrates the WT chloride transport partition, not the absolute WT intracellular chloride concentration.

## Scientific target

For WT R09 at the central stimulated condition (Ca = 0.25 uM):

- NKCC1 must carry 70% +/- 2% of positive basolateral chloride uptake;
- the remaining approximately 30% must come from positive AE4 + AE2 chloride uptake;
- CaCC must carry the corresponding apical chloride efflux required by stationary whole-cell balance;
- total positive chloride throughput must remain within 10% of the uncalibrated reference so the target cannot be met by shutting transport down.

The Task 31 Cha NHE1 mechanism and frozen NHE1 density are retained. Do not introduce another NHE1 model.

## Allowed calibration magnitudes

Only these may change:

1. NKCC1 capacity;
2. AE4 activity/capacity magnitude;
3. AE2 capacity;
4. CaCC conductance.

Use one deterministic bounded local constrained calibration. No parameter grids or combinatorial search.

R09 is calibration. Freeze the accepted parameter set before any genotype evaluation. R10 is out-of-sample validation only.

Do not fit WT Cl_i to 50.1 mM. Do not fit knockout Cl, pH, or secretion. Those are validation observations.

## Hard restrictions

- No Cartesian grid, pairwise/factorial enumeration, random, global, Bayesian, evolutionary or high-dimensional search.
- No automatic escalation to other transporters or mechanisms.
- Do not alter transporter stoichiometry or direction laws.
- Do not alter NHE1, pump, K, CO2, paracellular or water laws.
- Do not use the AE4-null phenotype during calibration.
- Absolute numerical limits: 12 stationary solver calls, 20,000 residual evaluations, 6 local calibration iterations/evaluations, 4 stimulated integrations, one worker, one BLAS thread and 30 minutes numerical execution.
- Do not run dynamics unless an admissible matched R09 WT/exact-null REST pair exists.
- Do not merge or modify `main`.

## GitHub publication in ChatGPT Work

Use the connected GitHub integration for all remote writes.

- Do not run `gh auth`, request a PAT, configure SSH or depend on shell `git push` credentials.
- Shell Git is allowed only for local status/diff inspection.
- Publish completed source edits, focused tests and compact UTF-8 outputs directly to `codex/task-33-nkcc70-chloride-partition` through the connected GitHub integration.
- Do not use manual Git objects, Base64 upload loops, archives or binary trajectory publication as a workaround.
- A publication problem must never trigger scientific recomputation.
