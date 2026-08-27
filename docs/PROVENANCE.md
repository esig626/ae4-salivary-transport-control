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

### 2026-08-27 — Phase 10 identifiability and discrimination audit

- Consulted the open published article by Vera-Sigüenza et al., *Bulletin of
  Mathematical Biology* 80 (2018), DOI
  `10.1007/s11538-017-0370-6`, PMCID `PMC5792321`, as the primary equation and
  parameter source. The public PMC HTML was used because the repository PDF
  pointer was not locally retrievable in this environment.
- Consulted the Phase 00 `model/` specification, inventory, open questions,
  provenance record, and decision log as established inputs.
- Consulted the published upstream Palk et al. salivary model (2010), DOI
  `10.1016/j.jtbi.2010.06.027`, and its 2013 corrigendum to resolve the volume
  sign and NKCC1 coefficient units. Consulted Sharp et al. (2015), DOI
  `10.1016/j.jtbi.2015.06.050`, for the individually defined inward CO2
  fluxes, and Peña-Münzenmayer et al. (2016), DOI
  `10.1085/jgp.201611571`, for the experimental AE4 transport context.
- Did **not** open or inspect the payload of any file under `archive/`, and did
  not use unpublished historical prose, equations, code, data, or figures.
- Independently derived the steady balance identities, regular-rank results,
  two-condition criterion, and diagnostic code from the published equations.
  The implementation is deliberately narrower than a full simulator because
  the printed model does not specify a reproducible numerical baseline.
- The earlier affine diagnostic already present on the working branch was
  retained as an explicitly uncalibrated software fixture. It was not used as
  physiological evidence, and its activity domain was narrowed to normalized
  knockout-to-WT values.
- For the targeted novelty audit, consulted only public scholarly records and
  primary publications listed with DOI links in
  `analysis/10_identifiability_discrimination/novelty_audit.md`. No text or
  code was copied from those publications.

## Future entries

Whenever archived material materially informs a research decision, append the date, file consulted, purpose, and what was independently rederived.
