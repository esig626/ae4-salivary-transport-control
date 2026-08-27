# Provenance record

## Historical material

The repository was initialised with a directory labelled `M. Golubitsky` containing an older AE4 analysis project. It includes an unpublished manuscript and associated build products, figures, MATLAB files, and a copy of the published AE4 model paper.

On 2026-08-27 the directory was moved without content modification to `archive/legacy-2017/` to separate historical collaborative material from the new work. The original state remains available in Git history.

A second repository pointer to the published 2018 AE4 paper is stored at `literature/core/2018_AE4_model.pdf` because that published paper is an appropriate scientific source for independent baseline reconstruction.

## Use policy

Historical unpublished material may be inspected to understand what was attempted and to prevent accidental duplication of old reasoning. It is not a source for new sole authored text or code.

The new implementation must be derived from published model equations and documented published evidence.

## Consultation log

### 2026-08-27 — Phase 00 inventory and source map

- Consulted `archive/README.md` to confirm the archive’s immutable historical status and its broad content categories.
- Consulted Git tree metadata (paths, file formats, and sizes) for `archive/legacy-2017/` to produce a complete file-level inventory. This metadata inspection did not open any archived file payload.
- Did **not** consult the contents of `Ae4_Basis.tex`, `Marty.tex`, `mybib.bib`, any historical MATLAB/MAT file, any archived figure, or any historical build output.
- Did **not** use the duplicate archived article. The scientific source was the separate published pointer `literature/core/2018_AE4_model.pdf` (cross-checked against the same open published article at PMCID PMC5792321).
- Purpose: classify historical holdings, prevent accidental reuse, and establish that the new `model/` specification depends only on published sources.
- Independently produced: state inventory, balance-equation map, flux/source map, parameter provenance table, observable definitions, and a list of published-record ambiguities. No Phase 01 code or manuscript prose was created.

## Future entries

Whenever archived material materially informs a research decision, append the date, file consulted, purpose, and what was independently rederived.
