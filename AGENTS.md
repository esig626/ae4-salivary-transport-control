# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The current branch is not another parameter sweep of the historical AE4 model. Its purpose is to determine what mechanistic structure is required for a modern, internally coherent salivary acinar-cell model to reproduce the experimentally observed AE4-null phenotype.

The executable task for this branch is `prompts/13_state_resolved_ae4_cycle_multiagent.md`.

The central question is:

> Can the AE4 knockout phenotype be reconstructed from a source-constrained, state-resolved AE4 transport cycle embedded in a physiologically defensible whole-cell transport network, and if not, what is the smallest additional network component required?

## Required correction to the interpretation of Task 12

Do **not** state or imply that “AE4 itself is insufficient.” That is not supported by the biology and is not what Task 12 established.

Task 12 established only the following narrower result:

> None of the tested coarse-grained AE4 flux formulations reproduced the AE4 knockout phenotype on the fixed seven-state historical chassis while satisfying the declared source and resting-state constraints.

This leaves open at least two scientifically distinct explanations:

1. the AE4 representation remains too coarse and omits an important internal transport-cycle feature; or
2. the surrounding cell model omits the network structure through which loss of AE4 produces the observed chloride, pH, volume, and secretion phenotype.

The present task must distinguish these possibilities rather than assuming either one.

## Experimental evidence has priority over historical model details

The biological knockout phenotype is independent evidence. Historical code, the 2018 mathematical article, and later reduced models are model records, not ground truth.

When evidence conflicts, use the following hierarchy.

1. Conservation laws, electroneutrality, thermodynamic feasibility, dimensional consistency, and sign conventions.
2. Primary experimental measurements of AE4 transport and salivary phenotypes.
3. Structural and mutagenesis evidence that constrains AE4 ion coordination and transport-cycle architecture.
4. Reproducible historical implementation behavior.
5. Published mathematical equations and parameter tables.
6. New modeling decisions, which must be labelled and justified.

The old article must not be discredited merely because its generating code is unavailable. Treat it as a legitimate historical model under the assumptions and information available at the time. The new work is a revisit and extension, not an attempted repudiation.

## Primary AE4 evidence that must inform the new model

At minimum, audit and use the following primary sources.

- Peña-Münzenmayer et al. 2015, JBC, DOI `10.1074/jbc.M114.612895` for the AE4-null salivary phenotype and resting ionic phenotypes.
- Peña-Münzenmayer et al. 2016, JGP, DOI `10.1085/jgp.201611571` for direct Na and K transport, electroneutrality, cation dose responses, Hill coefficients, and the alternative stoichiometries discussed by the authors.
- The AE4 PKA-regulation paper cited and audited in Task 12, DOI `10.1152/ajpgi.00145.2021`, for beta/cAMP/PKA activation constraints.
- Catalán et al. 2025, AJP Cell Physiology, DOI `10.1152/ajpcell.00346.2024` for the cation-coordination site, Na/K-specific mutant behavior, and the proposed alternative state-level transport-cycle hypotheses.
- The 2018 apical Na/K-pump and K-channel salivary paper already present in the project literature for experimentally supported cation-homeostasis topology.

Do not treat any proposed molecular cycle in these papers as established merely because it is drawn or discussed. Distinguish direct measurements, thermodynamic deductions, molecular-dynamics hypotheses, and speculation.

## State-resolved AE4 requirement

Task 13 must move beyond a single scalar flux law whenever the evidence supports doing so.

A valid state-resolved AE4 model should, where possible, distinguish:

- inward-facing and outward-facing transporter conformations;
- separate Na and K binding or coordination behavior rather than replacing both by a single pooled cation concentration;
- anion binding and release steps;
- reversible transport and transporter reversal;
- thermodynamic cycle constraints;
- cation-specific effects suggested by the T448/T756/D709/T713 coordination evidence;
- experimentally supported cAMP/PKA modulation;
- alternative HCO3 versus carbonate hypotheses only when the required chemical states are explicitly represented.

Do not invent microscopic detail merely to make the model elaborate. Prefer the smallest Markov or alternating-access scheme capable of representing the experimentally supported distinctions.

## Held-out validation discipline

The approximately 30–35 percent reduction in stimulated salivary secretion after AE4 deletion must not be directly fitted in the primary reconstruction.

Use staged validation.

### Stage A

Use WT physiology, transporter-level experiments, thermodynamic constraints, and source-supported model parameters to define and calibrate candidate AE4 mechanisms. Keep all AE4-null salivary outcomes held out.

### Stage B, only if Stage A fails

The AE4-null chloride, pH, and volume phenotypes may then be used to localize and reconstruct a missing whole-cell homeostasis module. The **magnitude and time course of the AE4-null saliva deficit remain held out**.

This is important. A model that is tuned directly to a 0.65 secretion ratio has not explained that ratio.

## Whole-cell reconstruction discipline

If state-resolved AE4 alone cannot rescue the phenotype on the current chassis, do not stop. Reconstruct the surrounding module systematically from independent evidence.

Candidate modules include, but are not limited to:

- apical and basolateral Na/K ATPase localization and activity;
- apical and basolateral K conductances and K recycling;
- NHE1 and acid-base coupling;
- CO2/HCO3 buffering and carbonic chemistry;
- NKCC1, while respecting evidence that AE4 deletion does not simply increase NKCC1 activity;
- membrane-potential closure;
- osmotic water transport and luminal outflow.

Add source-supported structure before adding free parameters. Record exactly which new module removes which residual or phenotype discrepancy.

## Multi-agent standard

Task 13 must be run as a genuine multi-agent investigation. Independent agents must cover at least:

1. experimental evidence and data extraction;
2. structural AE4 transport-cycle inference;
3. thermodynamic and Markov-model derivation;
4. whole-cell epithelial transport reconstruction;
5. numerical inference and continuation;
6. adversarial scientific audit;
7. independent numerical reproduction;
8. final synthesis.

Agents must record disagreements rather than averaging them away.

## Do not stop at the first failed family

The task must proceed hierarchically.

1. Build and test state-resolved AE4 cycles.
2. If they fail on the current chassis, mathematically localize the missing source/balance direction in the AE4-null system.
3. Add the smallest independently supported whole-cell module capable of supplying that direction.
4. Recalibrate only against non-secretion evidence.
5. Retest the held-out AE4-null secretion magnitude and time course.
6. Iterate until either a source-supported reconstruction succeeds or the remaining missing information is localized to a specific experimentally unconstrained flux, state, or regulatory process.

A generic `STOP` is not an acceptable final answer.

## Mathematical standard

Use exact balance relations, stoichiometric source vectors, thermodynamic cycle affinities, structural rank, monotonicity, residual projection, and sparse/minimal module reconstruction before broad parameter search.

When fitting is unavoidable, distinguish parameter estimation from model selection. Penalize unnecessary new degrees of freedom. Prefer nested model comparison and out-of-sample prediction over unconstrained optimization.

If a successful model is found, determine what feature is structurally necessary and what is merely one parameterization. Only then reopen identifiability, model-class discrimination, experimental design, or hypothesis testing.

## Numerical standard

Every accepted result must be reproducible with at least two independent numerical routes where feasible. For stiff dynamics use independent solvers or tolerances. Verify conservation and positivity. Report failed runs. Do not discard inconvenient branches silently.

All code must live outside `archive/`. Never edit the archive.

## Claims discipline

Separate throughout:

1. direct experimental fact;
2. source-supported mechanistic constraint;
3. model-derived result;
4. historical implementation evidence;
5. new modeling assumption.

Do not say that the old paper was wrong. Do not say that AE4 is insufficient. State precisely which model family, chassis, or mechanistic assumption failed.

## Completion standard

Task 13 is complete only when it reaches one of the following substantive outcomes with an evidence trail.

1. `PHENOTYPE RECONSTRUCTED — STATE-RESOLVED AE4 SUFFICIENT ON VALIDATED CHASSIS`
2. `PHENOTYPE RECONSTRUCTED — AE4 PLUS MINIMAL CHASSIS COUPLING REQUIRED`
3. `MULTIPLE SOURCE-CONSISTENT MECHANISMS RECONSTRUCT PHENOTYPE`
4. `PHENOTYPE UNRESOLVED — SPECIFIC MISSING EXPERIMENT OR FLUX IDENTIFIED`

For outcomes 1–3, the secretion deficit must be a held-out prediction. For outcome 4, the missing information must be localized narrowly enough that a specific experiment, parameter, or process would decide the question.

## Manuscript discipline

Do not draft a manuscript during Task 13. First solve the reconstruction problem, freeze the accepted model and validation results, and perform an adversarial audit.