# Phase 00 prompt. Inventory and source map

Work in repository `esig626/ae4-salivary-transport-control` on a new branch named `codex/phase-00-inventory`.

Read `AGENTS.md`, `README.md`, `docs/RESEARCH_PLAN.md`, and `analysis/00_inventory/README.md` before doing anything else.

Do not edit any file under `archive/`.

## Task

Complete Phase 00 only.

1. Inventory the contents of `archive/legacy-2017/` and classify each item as published source, unpublished manuscript material, historical code, figure or build artefact.
2. Use the published 2018 AE4 paper under `literature/core/2018_AE4_model.pdf` as the primary source for reconstructing the baseline model.
3. Create `model/SOURCE_MAP.md` mapping every state equation, flux model, parameter table, and observable required for the baseline reconstruction to a published location.
4. Create skeleton files `model/states.md`, `model/parameters.md`, `model/fluxes.md`, `model/equations.md`, and `model/observables.md` using only information supported by published sources.
5. Create `analysis/00_inventory/inventory.md` describing the historical material without copying unpublished prose or derivations.
6. Create `analysis/00_inventory/open_questions.md` listing unresolved model definitions, parameter provenance problems, and literature gaps that block independent reconstruction.
7. Update `docs/PROVENANCE.md` with any archive files consulted and the reason they were consulted.

## Prohibitions

- Do not port historical MATLAB code into `src/`.
- Do not copy text, derivations, tables, or figures from the unpublished manuscript.
- Do not begin Phase 01 implementation.
- Do not draft manuscript prose.
- Do not infer missing equations or parameter values without labelling them as unresolved.

## Acceptance criteria

Phase 00 is complete only if a new researcher could reconstruct the baseline model from the published source map without relying on the unpublished manuscript or historical code, or if every obstacle preventing that reconstruction is explicitly documented.
