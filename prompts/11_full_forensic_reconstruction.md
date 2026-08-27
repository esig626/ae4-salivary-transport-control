# Codex task. Full forensic reconstruction of the salivary model

Work in repository `esig626/ae4-salivary-transport-control` on branch `codex/full-forensic-reconstruction`.

Read `AGENTS.md`, `README.md`, `model/SOURCE_MAP.md`, every file in `model/`, `analysis/00_inventory/`, every file in `analysis/10_identifiability_discrimination/`, `docs/PROVENANCE.md`, `docs/DECISIONS.md`, and `docs/RESULTS_LEDGER.md` before doing anything else.

This branch inherits the completed Phase 00 audit and the completed identifiability attempt. The previous classification was `STOP` because the published 2018 specification could not be made into a dimensionally coherent, baseline-consistent executable model. That classification is now provisional because historical MATLAB files may be the actual implementation used for, or immediately ancestral to, the 2018 article.

The user has specifically flagged

- `archive/legacy-2017/Ae4_Dynamics_Project/Original/Parameters.m`
- `archive/legacy-2017/Ae4_Dynamics_Project/Original/Saliva_Ae4.m`

as possibly belonging to the 2018 implementation. Do not assume that they do. Establish this from evidence.

## Core objective

Perform an exhaustive forensic reconstruction of the project using **all available evidence in the repository**. Determine what the historical files actually are, whether any of them generated or directly encode the published 2018 AE4 model, which implementation conventions resolve the published inconsistencies, and whether a reproducible physiological model can be recovered strongly enough to reopen the identifiability project.

This task is intentionally broader than Phase 00. Phase 00 classified the archive without opening historical scientific content. This task must inspect it.

## Non-negotiable handling rules

- Do not edit, rename, delete, reformat, or overwrite anything under `archive/`.
- Historical files may be read, parsed, executed unchanged, compared, hashed, and inspected as forensic evidence.
- If execution requires a writable copy, copy files to an untracked temporary directory or a new clearly labelled forensic workspace. Never modify the archived originals.
- Do not treat unpublished code, manuscript text, or figures as published scientific evidence. Label conclusions based on historical material as implementation evidence or provenance evidence.
- Do not copy old unpublished prose, derivations, or source code into a new manuscript or into the new independent implementation.
- If a clean implementation is eventually created, derive it independently and document exactly which conventions come from the published paper and which are recovered from historical implementation evidence.
- Do not draft manuscript prose in this task.

## 1. Inspect every archive item

The archive inventory contains 67 files. Account for every one.

For each file under `archive/legacy-2017/`, inspect as much content as the format allows and record

- file path and type
- scientific role
- dependencies on other archive files
- whether it appears to predate, implement, analyse, postdate, or merely build the 2018 model
- direct links to equations, parameters, figures, simulations, or analysis in the 2018 publication
- direct links to the later unpublished Golubitsky manuscript
- confidence level and evidence for the classification

Do not skip build logs, auxiliary files, bibliography files, figures, MAT files, or converted PDFs merely because Phase 00 classified them as artefacts. They may contain useful chronology, figure names, variable names, dependencies, or parameter provenance.

Create `analysis/11_forensic_reconstruction/archive_audit.md` with a complete file-level ledger.

## 2. Establish code lineage rather than guessing it

Inspect all historical MATLAB files, including at least

- `Original/Parameters.m`
- `Original/Saliva_Ae4.m`
- `Par.m`
- `Par.mat`
- `Salivary.m`
- `Salivary2.m`
- `Salivary_ex.m`

Reconstruct their dependency graph and likely development lineage.

Determine whether `Original/` is genuinely an earlier/full implementation and whether the outer files are later reductions or reorganisations. Base the conclusion on structural evidence such as

- state dimension
- equation structure
- parameter names and numerical values
- algebraic eliminations
- calcium protocol
- solver choice and time scaling
- AE2 and AE4 parameterisation
- membrane-potential treatment
- water-flow conventions
- output figures
- comments and variable names
- references from the TeX manuscript

Do not rely on directory names alone.

Create `analysis/11_forensic_reconstruction/code_lineage.md`.

## 3. Build an equation-by-equation paper/code concordance

Compare the historical implementation against the published 2018 article equation by equation and parameter by parameter.

For every published state equation, algebraic constraint, constitutive flux, water flux, transporter model, current equation, parameter table, baseline state, and reported observable, identify

- the corresponding historical code expression if one exists
- whether it is algebraically identical
- whether signs differ
- whether units or scale factors differ
- whether numerical constants differ
- whether the code uses a calibrated value not transparently reported in the paper
- whether the code contains a correction to an apparent publication typo
- whether the code contains a different model rather than merely a corrected implementation

Pay particular attention to

- `d omega_i / dt` sign
- definitions and directions of `q_a`, `q_b`, `q_t`, and `q_tot`
- water permeability rescaling or calibration
- NKCC1 units and the Palk corrigendum
- NaK scaling
- current versus particle-flux conversion
- membrane areas and volumes
- AE4 fourth-order kinetics and density/activity units
- AE2 and NHE1 units
- CO2 transport sign and buffer scaling
- luminal chloride balance
- electroneutrality and eliminated H+, CO2, and luminal chloride variables
- quasi-steady membrane potentials
- calcium input
- solver time scaling
- knockout implementation

Create `analysis/11_forensic_reconstruction/paper_code_concordance.md` with a table whose status categories include at least `exact match`, `equivalent after algebra`, `publication typo likely`, `implementation convention`, `substantive model difference`, and `unresolved`.

## 4. Inspect the historical TeX, bibliography, build records, and figures

Read the unpublished TeX files and bibliography only as provenance evidence.

Determine

- which code outputs are referenced by which figures
- whether figure captions, variable names, parameter ranges, or numerical values identify the MATLAB version used
- whether the unpublished mathematical-analysis manuscript explicitly describes a reduction of the published 2018 model or a different implementation
- whether any figures in the archive duplicate or closely correspond to figures in the 2018 paper
- whether build logs or auxiliary files reveal chronology or dependencies

Inspect figure files visually or numerically where feasible. Compare axes, ranges, curve shapes, knockout labels, and baseline values with the published 2018 figures. Do not infer identity from visual similarity alone.

Create `analysis/11_forensic_reconstruction/figure_and_manuscript_provenance.md`.

## 5. Inspect `Par.mat` and all hidden numerical state

Load `Par.mat` using MATLAB, Octave, Python `scipy.io.loadmat`, or another non-destructive reader.

Record every stored variable, value, shape, and type. Compare them against `Par.m`, `Parameters.m`, the 2018 tables, and the later MATLAB functions.

Search the repository for any other binary or generated parameter/state files that contain implementation values not obvious from the text sources.

Create machine-readable dumps in `results/11_forensic_reconstruction/`.

## 6. Reproduce the historical implementation unchanged

If the required runtime is available, execute the historical MATLAB code unchanged from a temporary working copy. If MATLAB is unavailable, use Octave if compatible. If neither is available, recreate only the execution harness needed to evaluate the archived equations without changing them, and state exactly what could not be executed literally.

First determine the correct dependency combination. Do not arbitrarily mix `Original/` files with later files unless the dependency graph supports it.

For every plausible runnable version

- compute the baseline residuals
- verify whether the prescribed baseline is actually a steady state under the code's own conventions
- reproduce steady or transient state values
- reproduce WT, AE2 knockout, and AE4 knockout results if the historical code supports them
- compare generated outputs numerically with the 2018 tables and figures
- record solver settings and any hidden time scaling

Do not retune parameters merely to obtain a match. Run the historical implementation first as found.

Save scripts used only for forensic execution under `analysis/11_forensic_reconstruction/tools/` or `src/forensics/`, clearly labelled as forensic tooling rather than the new scientific implementation.

Create `analysis/11_forensic_reconstruction/reproduction.md` and machine-readable residual/output files under `results/11_forensic_reconstruction/`.

## 7. Revisit every blocker from the previous STOP decision

The prior analysis reported, among other issues

- approximately `59.6 x` water-balance mismatch
- approximately `67.8 pA` electrical residual
- acid/base inconsistency
- AE4 unit ambiguity
- missing calcium input
- unresolved flux/current conversions

For each previous blocker, classify it after the forensic audit as one of

- `resolved by historical implementation convention`
- `publication typo resolved by code and conservation`
- `implementation-specific calibration not stated clearly in paper`
- `historical code also inconsistent`
- `still unresolved`

Recompute every numerical discrepancy using the implementation's actual unit and scaling conventions. Explain exactly why the previous audit obtained a mismatch if the mismatch disappears.

Create `analysis/11_forensic_reconstruction/blocker_resolution.md`.

## 8. Decide whether the historical code is actually the 2018 implementation

Give a graded conclusion rather than a binary guess.

Use categories

- `CONFIRMED`
- `VERY LIKELY`
- `PLAUSIBLE BUT UNPROVEN`
- `UNLIKELY`
- `NOT THE 2018 IMPLEMENTATION`

Support the classification with an evidence matrix covering equation identity, parameter identity, output reproduction, figure correspondence, chronology, and dependency structure.

A directory name or the user's memory is not sufficient evidence. Conversely, do not reject the code merely because it contains algebraic reductions or calibration steps if those are mathematically equivalent to the published model.

## 9. Conditional clean reconstruction

Only if the forensic evidence reaches at least `VERY LIKELY` and the historical implementation reproduces the central published baseline/phenotype sufficiently well, reopen the executable-model work.

Then

1. write an independent clean implementation under a new namespace such as `src/reconstructed_2018/`;
2. derive it from the published equations plus explicitly documented historical conventions needed to disambiguate the paper;
3. do not port source code line by line;
4. add unit, balance, steady-state, WT, AE2-KO, and AE4-KO regression tests;
5. verify the clean implementation against the historical implementation and published outputs separately.

Every recovered convention must be labelled either `published`, `historical implementation evidence`, or `new modelling decision`.

Do not proceed if reproduction requires arbitrary tuning or undocumented choices that materially affect the phenotype.

## 10. Conditional reopening of identifiability and discrimination

If and only if a validated physiological steady-state map is recovered, revisit the previous `STOP` classification.

Using the reconstructed executable model

- compute `F_u` and `F_theta`
- verify implicit derivatives against independently re-solved finite differences
- certify or refute candidate panels such as `(Q, Na_i)`, `(Q, K_i)`, `(Q, Cl_i)`, and other biologically plausible outputs
- compute actual physiological observation geometry over a justified activity domain
- test whether one additional experimental condition rotates the sensitivity directions enough to restore identifiability
- inspect global equivalence, folds, self intersections, singular sets, and mechanism-class overlap
- reopen hypothesis testing only if the physiological class images support a nontrivial discrimination problem

Reuse the analytical results from Task 10 only after checking that their assumptions match the recovered model.

If the recovered model produces a genuinely stronger result, classify `PUSH HARDER`. If it supports only a compact but defensible result, classify `QUICK PAPER`. If the historical code cannot establish a trustworthy physiological map, retain `STOP`.

## 11. Literature and provenance cross-check

Use web/literature search where necessary to verify

- Palk et al. model/corrigendum details
- upstream transporter equations
- whether a corrected or supplementary implementation of the 2018 model exists publicly
- whether later salivary modelling papers reused or clarified these conventions

Keep published-source evidence separate from historical-code evidence.

Update `docs/PROVENANCE.md` with every historical scientific file consulted and its role.

## Required outputs

Create at minimum

- `analysis/11_forensic_reconstruction/archive_audit.md`
- `analysis/11_forensic_reconstruction/code_lineage.md`
- `analysis/11_forensic_reconstruction/paper_code_concordance.md`
- `analysis/11_forensic_reconstruction/figure_and_manuscript_provenance.md`
- `analysis/11_forensic_reconstruction/reproduction.md`
- `analysis/11_forensic_reconstruction/blocker_resolution.md`
- `analysis/11_forensic_reconstruction/final_assessment.md`
- machine-readable outputs under `results/11_forensic_reconstruction/`
- forensic helper scripts if needed
- tests for every reproducible numerical claim
- updates to `docs/PROVENANCE.md`, `docs/DECISIONS.md`, and `docs/RESULTS_LEDGER.md` only where justified

If conditional reconstruction succeeds, also create the independent clean implementation and its regression tests.

## Final assessment requirements

`final_assessment.md` must answer directly

1. What are these historical MATLAB files?
2. Which files are most likely to have generated the 2018 results?
3. What exact evidence supports that conclusion?
4. Which published inconsistencies are resolved by the implementation?
5. Which inconsistencies remain real?
6. Does the implementation reproduce the published baseline and AE2/AE4 phenotype?
7. Is there now a trustworthy physiological parameter-to-observation map?
8. Does the previous `STOP` classification change to `QUICK PAPER`, `PUSH HARDER`, or remain `STOP`?
9. What is the strongest mathematical result now available?
10. What are the next three steps, if any?

## Acceptance criteria

This task is complete only when

- every archive file has been accounted for
- every historical scientific source has actually been inspected rather than merely classified by filename
- all MATLAB variants have a documented dependency and lineage analysis
- `Par.mat` has been inspected
- the code/paper comparison is equation-by-equation and parameter-by-parameter
- the previous numerical blockers have been recomputed using the historical conventions
- at least one serious attempt has been made to execute/reproduce the historical implementation
- the conclusion about whether these files belong to the 2018 model is evidence based
- `archive/` hashes are unchanged
- no manuscript has been drafted

Commit all work to `codex/full-forensic-reconstruction`. Do not merge to `main`.
