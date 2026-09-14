# Active phase: Task 38 AE4 perturbation validation

Work only on `codex/task-38-ae4-perturbation-validation`.
Read and execute `prompts/38_ae4_perturbation_validation.md`.

This task begins from the successful frozen Task 37 WT NBC model. It is a prediction task only. The model is not to be recalibrated, rescued, expanded or searched.

## Frozen reference

- Task 37 WT reference commit: `2ba0e4e0c2a9e940aa96ed86c3c72456afd0e3aa`.
- Use exactly the same accepted Task 31 R09 WT resting initial state as Task 37.
- Keep the Task 31 Cha et al. mechanistic NHE1 carrier amount unchanged.
- Keep routed AE4 with `1 Cl : 1 monovalent cation : 2 HCO3` working stoichiometry.
- Keep the corrected electrogenic NBC exactly unchanged: `1 Na : 2 HCO3`, stimulus-recruited, with the same fixed capacity and voltage-dependent affinity/current closure.
- Keep the source-fixed stimulated NHE1 factor exactly unchanged.
- Keep all Task 37 AE4 beta/PKA, NKCC1, AE2, pump, channels, CO2, paracellular, water, bath, geometry and solver settings unchanged.
- Do not solve genotype-specific REST states.

## Exactly two intended runs

Run only:

1. AE4 expression = `0.05`
2. AE4 expression = `0.00`

For both, use the same WT initial state and identical Task 37 standard CCh + isoproterenol stimulus for 600 s.

The only scientific change between runs is `genotype.ae4_expression`.

## No search

- no optimisation
- no stationary solve
- no parameter sweep
- no alternative NBC capacity or stoichiometry
- no alternative NHE1 gain/model
- no alternate AE4 regulation
- no Ca/CCh/IPR/beta changes
- no AE2 or NKCC1 perturbations
- no intermediate AE4 expression values
- no R10
- no phenotype fitting
- no rescue attempt after failure
- no automatic fallback mechanism
- no pairwise, factorial, random, global, evolutionary, Bayesian, Latin-hypercube or Shapley analysis

## Gates

Apply the same physical and conservation gates used in Task 37.

If AE4=0.05 fails a scientific gate, stop before AE4=0.00 with `AE4_5PCT_PERTURBATION_FAILED`.

If AE4=0.05 passes but AE4=0.00 fails, stop with `AE4_NULL_PERTURBATION_FAILED`.

If both pass, report `AE4_PERTURBATION_VALIDATION_COMPLETE`.

The experimental AE4-null secretion deficit is held-out validation only. Do not tune anything to match it.

## Hard compute limit

- exactly two intended perturbation integrations
- at most one purely numerical solver retry per perturbation
- zero optimisation calls
- zero stationary solves
- zero additional scientific parameter evaluations
- one worker
- one BLAS thread
- maximum 10 minutes numerical execution

## GitHub publication

Use the connected GitHub integration for remote writes.
Do not use shell Git authentication, PATs, SSH configuration or `gh auth`.
Do not merge to `main`.
Publish compact UTF-8 results only.

The task ends after the two perturbation reports. Do not start another mechanism or calibration.
