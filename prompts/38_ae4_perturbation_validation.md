# Task 38: AE4 5% and 0% perturbation validation

## Scope

This task starts from the successful Task 37 WT NBC validation at commit `2ba0e4e0c2a9e940aa96ed86c3c72456afd0e3aa`.

The model is frozen.

Run exactly two additional stimulated trajectories from the same accepted Task 31 R09 WT resting state used in Task 37:

- AE4 expression = 0.05
- AE4 expression = 0.00

Use the identical Task 37 stimulus, NBC mechanism, NHE1 stimulation, AE4 beta/PKA regulation, NKCC1, AE2, pump, channels, CO2, paracellular transport, water, geometry, bath and numerical solver settings.

Do not solve genotype-specific resting states.
Do not alter any parameter.
Do not rerun a WT trajectory unless needed only to load or verify the frozen Task 37 reference outputs. The Task 37 WT trajectory is the reference.

## Absolute prohibitions

- no optimisation
- no stationary solve
- no parameter sweep
- no alternative NBC capacity
- no alternative NBC stoichiometry
- no alternative NHE1 multiplier
- no Ca sweep
- no beta/IPR sweep
- no AE2 perturbation
- no NKCC1 perturbation
- no intermediate AE4 expression values
- no R10
- no fitting to the experimental AE4-null secretion deficit
- no rescue attempt if a perturbation fails
- no automatic addition of another mechanism
- no pairwise, factorial, random, Latin-hypercube, global, evolutionary, Bayesian or Shapley analysis

## Initial state and stimulus

Use exactly the same saved Task 31 R09 WT resting core state and resting regulatory coordinate as Task 37.

The only changed quantity at t=0 is `genotype.ae4_expression`.

Use the identical standard stimulated protocol:

- CCh = 0.3 uM
- isoproterenol = 5 uM
- stimulated Ca = 0.25 uM
- beta occupancy = 1.0
- duration = 600 s
- production Radau settings identical to Task 37

The corrected electrogenic 1 Na : 2 HCO3 NBC and source-fixed stimulated NHE1 factor remain unchanged.

## Runs

Run exactly:

1. AE4 = 0.05
2. AE4 = 0.00

Maximum two intended integrations total. One retry for a purely numerical solver failure is permitted per trajectory, but only if no scientific gate has failed. A retry must use the already-declared independent stiff fallback solver or equivalent pre-existing production fallback, with no parameter changes.

## Physical and conservation gates

Apply the same Task 37 monitored-state and conservation gates to both trajectories.

At minimum require:

- all tracked amounts, concentrations and volumes finite and positive
- Na_i < 40 mM
- 50 <= K_i <= 200 mM
- 30 <= Cl_i <= 80 mM
- 6.6 <= pH_i <= 7.3
- HCO3_i > 0 and < 100 mM
- V_i > 0 and < 3 pL
- inherited current, charge, carbon, buffer and water residual tolerances pass

If either perturbation violates a physical or conservation gate, report the failure and STOP. Do not alter anything.

## Required outputs

For each perturbation report:

- cumulative water outflow over 0-600 s
- ratio of cumulative secretion to frozen Task 37 WT
- mean q_out over 60-600 s
- q_out at 600 s
- q_out at 0, 60, 120, 180, 240, 300, 360, 420, 480, 540 and 600 s
- Na_i, K_i, Cl_i, pH_i, HCO3_i, TIC_i and cell volume at those same times
- integrated positive NKCC1 chloride loading over 60-600 s
- integrated positive AE4 chloride loading over 60-600 s
- signed AE2 chloride flux over 60-600 s
- integrated NBC cycles and bicarbonate-equivalent influx
- integrated NHE1 flux
- integrated CaCC chloride export
- integrated paracellular chloride return
- chloride-loading shares in the positive pool
- WT minus perturbation cumulative secretion versus time at 60 s intervals

Also report:

- whether AE4 = 0.05 and AE4 = 0.00 are materially different from each other
- whether the secretion deficit is small early and becomes larger after approximately 2-3 minutes
- whether loss of AE4 lowers intracellular chloride
- whether pH remains reasonably controlled
- whether NKCC1 changes substantially despite no parameter change
- whether AE2 becomes a compensatory positive chloride loader or remains negligible/reverse

## Interpretation rules

The experimental approximately 35% AE4-null secretion deficit is held-out validation only.
Do not tune to it.

Do not call a mismatch a reason to change parameters inside this task.

The purpose is only to measure the model prediction after freezing the successful Task 37 WT architecture.

## Compute budget

- exactly two intended perturbation integrations
- at most one purely numerical retry per perturbation
- zero optimisation calls
- zero stationary solves
- zero parameter evaluations beyond the two fixed AE4 expressions
- one worker
- one BLAS thread
- maximum 10 minutes numerical execution

## Stop conditions

If 5% fails a physical/conservation gate, stop before 0% and report:

`AE4_5PCT_PERTURBATION_FAILED`

If 5% passes but 0% fails, report:

`AE4_NULL_PERTURBATION_FAILED`

If both pass, report:

`AE4_PERTURBATION_VALIDATION_COMPLETE`

and include the two secretion ratios and all required diagnostics.

Do not start another task or mechanism.

## GitHub publication

Use the connected GitHub integration for remote writes.
Do not use shell Git authentication, PATs, SSH configuration or `gh auth`.
Do not merge to `main`.
Publish compact UTF-8 result tables/reports and any minimal execution script needed for this task only.
