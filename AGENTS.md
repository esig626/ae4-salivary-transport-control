# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The active task is Task 14. The executable prompt is

`prompts/14_scale_free_genotype_holdout_validation.md`.

Task 13B produced a modern conservation explicit WT model family with ten retained native WT roots. Task 13C then tested whether increasing only stimulated calcium could close the inherited absolute one SMG WT flow discrepancy. It could not. The controlling completed Task 13C commit is

`2f54e7c4f87b6da746d3a8427c7bdca40a77c3de`.

Task 14 asks the central genotype question with the model frozen:

> Does exact AE4 deletion predict the held out WT relative secretion phenotype and its time structure, while exact AE2 deletion produces a much smaller secretion effect?

This is a blind prediction and validation task. It is not a calibration task and not a reconstruction task.

## Scientific interpretation fixed for Task 14

The Task 13B and Task 13C absolute whole gland secretion discrepancy remains recorded, but the inherited 9 to 10 uL/min one SMG mapping is now a nonblocking observation scale diagnostic.

Do not delete, hide or rewrite that discrepancy.

Do not use it to veto genotype prediction.

A common multiplicative cellular to gland scale cancels from null to WT secretion ratios. The primary Task 14 quantities are therefore paired WT relative secretion measures, not absolute whole gland flow.

Do not rescale WT and null separately. Do not fit an observation scale after phenotype reveal.

## Hard freeze

The Task 14 branch begins exactly from the completed Task 13C commit above.

Keep frozen exactly as inherited:

* all ten retained WT root payloads;
* all whole cell equations;
* all transporter capacities and source multipliers;
* NKCC1 family, source scale and kinetics;
* AE4 carrier model, rates, stoichiometry and Na K routing;
* AE2 and NHE1;
* pump and K channel capacities and membrane fractions;
* CaCC and calcium gate;
* paracellular pathways;
* acid base and carbon chemistry;
* hydraulic parameters, lumen parameters and outflow law;
* geometry and bath composition;
* beta input and regulatory gains and kinetics;
* resting OTHER osmoles;
* solver and conservation tolerances;
* exact genotype deletion semantics.

No model parameter may be optimised or changed in Task 14.

No genotype specific compensation is permitted.

## Calcium panel

The exact matched SMG calcium amplitude remains uncertain. Carry forward exactly the three already evaluated calcium values

`0.10, 0.25, 0.50 uM`.

Use all three as a fixed sensitivity panel.

Do not add a new calcium value and do not choose the best calcium after holdout reveal.

The preferred dynamic regulatory member is

`R1_G125_P21_PROSE_TREFERENCE`.

The primary protocol remains WT or null `CCH_IPR`, 0.3 uM CCh plus 5 uM IPR, over 600 s with the inherited beta input.

## Genotype resting states

Do not evaluate deletion dynamics from the WT resting state.

For every retained WT root, continue AE4 expression and AE2 expression independently from one to exact zero using the existing target free continuation framework.

Fit zero model parameters during continuation.

Keep the connected branch by the predeclared numerical distance rule only.

Require exact deleted transporter flux at zero expression and retain alternate mathematical roots for audit.

## Blind holdout firewall

Before the blind prediction checkpoint is committed and pushed, do not inspect the exact held out AE4 secretion magnitude or time course.

Do not open, parse, digitise, search history for, or otherwise recover concealed phenotype values.

Do not use exact target values in code, tests, thresholds, root selection or calcium selection.

Before simulation create the machine readable prediction contract required by the prompt.

Generate and save all blind predictions, tests and hashes. Then commit and push the blind checkpoint.

Only after that pushed commit exists may the held out phenotype be opened for one time evaluation.

After reveal, no model change is permitted.

## Primary prediction quantities

For every valid root and calcium value report matched WT relative quantities for AE4 null and AE2 null.

At minimum report

* total 0 to 600 s secretion ratio;
* minute specific flow ratios;
* minute specific cumulative secretion ratios;
* early 0 to 180 s integrated ratio;
* sustained 180 to 600 s integrated ratio;
* the change from early to sustained ratio;
* AE4 effect, AE2 effect and their contrast;
* resting and endpoint intracellular Cl, pH, Na, K and volume differences.

Absolute cellular and mapped gland flows remain secondary diagnostics only.

All ten roots and all three calcium values remain in the report regardless of agreement.

No best root or best calcium selection is allowed after reveal.

## Numerical standard

Use production Radau for the primary blind dynamic prediction.

Use BDF confirmation according to the target independent representative rule in the Task 14 prompt.

Preserve positivity, exact deletion, current closure, charge consistency, carbon and water conservation and the inherited numerical tolerances.

Parallel execution across roots is encouraged, with numerical library thread counts limited to avoid oversubscription.

Any continuation or solver anomaly must be resolved or recorded before holdout reveal without phenotype information.

## Post reveal discipline

Compare the frozen blind predictions with the held out data at the resolution actually supported by the experiment.

Experimental uncertainty is evidence, not a universal biological law.

Where uncertainty intervals are reported, use them. Where they are absent, report residuals without inventing intervals.

Do not fit a scale, time shift, gain, offset or smoothing parameter.

Do not change a parameter, root, calcium input, regulation rule, topology or mechanism after reveal.

If the prediction fails, diagnose the failure and stop. Any repair belongs in a later independently designed task.

## General discipline

Do not edit `archive/`.

Do not merge to `main`.

Do not draft a manuscript.

Do not perform model reduction, GSPT, bifurcation analysis or a new identifiability study.

Do not state that the 2018 paper was wrong merely because the absolute gland scale remains unresolved.

Do not state that AE4 is mechanistically validated unless the blind genotype comparison supports that conclusion.
