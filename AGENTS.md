# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The current branch is a forensic reconstruction of the salivary AE4 project. The immediate objective is to determine whether the historical files in `archive/legacy-2017/` implement, directly precede, or otherwise explain the published 2018 AE4 model, and whether they resolve the published inconsistencies that caused Task 10 to classify the identifiability paper as `STOP`.

The executable task for this branch is `prompts/11_full_forensic_reconstruction.md`. It supersedes the earlier quick-paper workflow for this branch only.

Do not assume either that the archived MATLAB files are the 2018 implementation or that they are not. Establish provenance and model identity from equations, parameters, dependencies, figures, and numerical reproduction.

## Non negotiable provenance rules

1. Treat historical unpublished material as immutable evidence.
2. Never edit, rename, delete, reformat, or overwrite files under `archive/`.
3. Historical files may be read, parsed, visually inspected, hashed, and executed unchanged for forensic purposes.
4. If historical execution requires writable files, use a temporary or clearly labelled forensic working copy outside `archive/`.
5. Do not copy unpublished prose, figures, derivations, or source code into a new manuscript or independent scientific implementation.
6. Published sources remain the primary scientific sources. Historical material may resolve implementation provenance or ambiguous conventions, but must be labelled as historical implementation evidence rather than published evidence.
7. If a clean implementation is created, derive it independently and record whether each convention is `published`, `historical implementation evidence`, or a `new modelling decision`.
8. Update `docs/PROVENANCE.md` whenever historical scientific content is consulted.
9. Do not attribute sole authorship to historical collaborative work. Any new project must stand on newly derived results.

## Established inputs

Phase 00 is complete at commit `2f71f23098037037836853fd96814a0e286cc0af` and is part of this branch history.

Task 10 is also complete on the parent branch. Its `STOP` classification is provisional for the present forensic task because it used the published specification rather than the historical implementation as the executable source. Preserve its analytical results and blocker diagnostics, but recompute model-specific conclusions if a validated implementation is recovered.

## Forensic standard

Inspect all available evidence before concluding provenance or model identity. Account for every archive file. Do not rely on filenames or directory names alone.

For historical MATLAB variants, establish the dependency graph and likely lineage using

- state dimensions and algebraic eliminations
- equation structure
- numerical parameter matches
- unit and scaling conventions
- solver and calcium protocols
- output figures
- references in historical TeX
- reproducibility of published baseline and knockout outputs

For every mismatch between publication and code, distinguish

- algebraic equivalence
- likely publication typo
- undocumented implementation convention
- calibration choice
- substantive model difference
- unresolved inconsistency

Do not repair a model silently.

## Model implementation rules

- Preserve mass balance, charge conventions, stoichiometry, and units explicitly.
- Keep parameter definitions separate from solver code.
- Keep observables separate from state equations.
- Add regression tests for any claimed baseline reproduction.
- Verify historical-code residuals under the historical code's own conventions before judging them against the publication.
- Do not retune historical parameters merely to force agreement with the paper.
- If a clean reconstructed model is created, verify it independently against both the historical implementation and published outputs.
- Do not claim global identifiability from local rank alone.
- Do not claim mechanism discrimination when physiological observation sets overlap over the stated domains.

## Mathematical standard

The project is not a sensitivity-analysis or bifurcation-catalogue exercise. If a trustworthy physiological steady-state map is recovered, resume the identifiability programme using analytical rank, level-set, elimination, stoichiometric, monotonicity, global-geometry, perturbation, and testing arguments where justified.

Do not stop at routine local rank if stronger structure is present. Conversely, do not force a general theorem or hypothesis-testing section when the recovered model does not support one.

## Claims discipline

Separate three evidence levels throughout

1. published scientific evidence
2. historical implementation or provenance evidence
3. new inference or modelling decisions

Do not promote implementation details to biological facts. Do not hide contradictions. A failed reproduction is a result.

## Manuscript discipline

Do not draft manuscript prose during Task 11. First establish what the historical materials are, whether the model can be reproduced, whether the previous blockers are resolved, and whether the project should be classified `QUICK PAPER`, `PUSH HARDER`, or `STOP`.
