# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The current branch asks a narrower and more consequential question than the previous reconstruction tasks.

The biological AE4 knockout phenotype is independent evidence. The published 2018 model and the archived MATLAB implementations are both imperfect records of the intended model. Do not assume that either record is correct.

The immediate objective is to determine whether the missing approximately 30 percent AE4 knockout reduction can be explained by the **structure of the AE4 transport model itself**, while changing as little else in the salivary secretion model as possible.

The executable task for this branch is `prompts/12_ae4_mechanism_reconstruction_multiagent.md`.

## Established inputs

Phase 00, Task 10, and Task 11 are complete in this branch history.

Task 11 established the following.

- The archived code has direct lineage with the published project but is unlikely to be the exact publication generating implementation.
- The archived seven state model closes its own calibrated WT baseline.
- Its AE4 knockout changes endpoint flow by only about minus 0.144 percent rather than the experimentally observed approximately 30 percent reduction.
- The archived AE4 law differs materially from the published AE4 mechanism, including its cation treatment.
- The exact publication generating implementation remains unavailable.

Treat those findings as evidence, not as a reason to stop the present task.

## Evidence hierarchy

Do not treat the 2018 article as ground truth. Use the following hierarchy when records conflict.

1. Conservation laws, electroneutrality, dimensional consistency, and thermodynamic feasibility.
2. Primary experimental evidence for AE4 transport and salivary knockout phenotypes.
3. Reproducible behaviour of historical implementations.
4. Published equations, parameter tables, and prose.
5. New modelling choices, which must be explicitly labelled and justified.

The published paper may contain errors. The archived code may also contain errors or represent a different development branch.

## Non negotiable provenance rules

1. Never edit, rename, delete, reformat, or overwrite files under `archive/`.
2. Historical files may be read, parsed, visually inspected, hashed, and executed unchanged for forensic purposes.
3. Any new implementation must live outside `archive/` and must not be a silent copy of unpublished code.
4. Every new equation or parameter convention must be labelled as `published`, `primary experimental evidence`, `historical implementation evidence`, or `new modelling decision`.
5. Update `docs/PROVENANCE.md` when historical scientific content is consulted.
6. Do not attribute sole authorship to historical collaborative work.

## Reconstruction discipline

The knockout phenotype must **not** be used as a calibration target in the primary model comparison.

For each AE4 mechanism candidate

- keep the non AE4 chassis fixed in the primary comparison
- fit or normalize only AE4 specific quantities using WT baseline closure and independent AE4 transport evidence
- then predict the AE4 knockout phenotype as held out validation
- report failure as failure

A secondary robustness stage may vary non AE4 quantities within source supported uncertainty ranges, but it must remain separate from the primary comparison and must not hide parameter compensation.

Do not tune arbitrary parameters until the desired knockout ratio appears.

## Scientific question

Determine which structural properties of AE4 are required, if any, for an otherwise fixed salivary secretion model to reproduce the experimentally observed AE4 knockout phenotype while retaining the negligible AE2 knockout phenotype and other source supported physiological constraints.

Possible structural differences include, but are not limited to

- Na only versus Na and K coupling
- explicit Na and K branches rather than a pooled cation term
- cation selectivity
- cooperative cation dependence
- transported cations versus cation dependent allosteric gating
- alternative electroneutral stoichiometries supported or discussed by the primary AE4 experiments
- reversible thermodynamically consistent transport laws
- transporter reversal or saturation regimes

Do not assume any one of these is correct before testing it.

## Multi agent standard

The executable prompt uses independent agents for experimental evidence, model lineage, AE4 mechanism derivation, numerical reconstruction, and adversarial validation. The lead agent must reconcile their conclusions and must not suppress disagreements.

## Completion standard

Do not stop at the first candidate that fails or approximately matches the knockout result. Systematically exhaust the source justified candidate family defined in the task.

The task is complete only when one of the following evidence backed answers is reached.

1. A minimal AE4 mechanism class reproduces the held out phenotype and the required structural ingredient is identified.
2. No source justified AE4 mechanism in the tested family can reproduce the phenotype on the fixed chassis, with a documented exclusion result and the next missing model component localized.
3. Multiple mechanisms remain observationally equivalent, in which case the task must identify the experiment or observable required to distinguish them.

A vague `STOP` is not an acceptable answer for this branch.

## Mathematical standard

Use analytical reasoning where it helps isolate mechanism. Mass balance, stoichiometry, thermodynamic driving force, monotonicity, limiting cases, rank, and exact algebra should be used before broad numerical search.

If a successful mechanism family is found, then and only then revisit identifiability, mechanism discrimination, experimental design, or hypothesis testing.

## Manuscript discipline

Do not draft a manuscript during Task 12. First obtain the mechanism level answer, validate it independently, and record it in the results ledger.