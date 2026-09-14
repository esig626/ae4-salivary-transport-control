# AE4 manuscript

Working title:

**AE4 control of salivary secretion nearly a decade later: acid-base balance, transporter stoichiometry and missing network regulation**

This folder contains the current source-audited manuscript built around the final reconstruction and mechanism tests. It is organised by scientific argument rather than by task chronology.

## Current PDF

[`AE4_manuscript_current.pdf`](AE4_manuscript_current.pdf) is the current compiled manuscript for reading and author review. It is generated from the LaTeX source in this directory and the frozen figure values in `data/figure_source_values.csv`; it is not an independent source of numerical results.

## Reproducible build

From the repository root:

```bash
bash manuscript/build.sh
```

The build first regenerates all six manuscript figures from frozen values, then runs LaTeX, BibTeX and the final LaTeX passes. The GitHub Actions workflow `.github/workflows/build-manuscript.yml` performs the same complete build on manuscript changes and exposes the compiled PDF as a workflow artifact.

## Files

- `main.tex`: title, abstract, Introduction and submission statements
- `model.tex`: complete model definition and numerical protocol
- `results.tex`: results organised around the figures
- `discussion.tex`: interpretation, scope, limitations and experimental predictions
- `references.bib`: consolidated BibTeX database, including the corrected 2018 archive bibliography
- `figure_plan.md`: figure order and the argument carried by each figure
- `claims_and_evidence.md`: manuscript claim and evidence ledger
- `corrected_archive_audit.md`: audit of the corrected 2018 source archive against the current article
- `data/figure_source_values.csv`: frozen values used by the plotting script
- `scripts/build_figures.py`: deterministic figure builder; it does not rerun the model
- `figures/`: generated PDF figures for LaTeX and PNG previews
- `AE4_manuscript_current.pdf`: complete compiled manuscript for review

## Submission status

The scientific draft and corrected-archive audit are complete enough for author review, but several submission items remain deliberately unresolved: final author list and order, affiliations, funding, contribution statements, competing-interest confirmation, the final placement of the target-selected Task 41 construction, Springer/BMB production formatting, final figure styling, and a final line-by-line numerical and citation audit.

The live checklist is maintained in [`../docs/MANUSCRIPT_STATUS.md`](../docs/MANUSCRIPT_STATUS.md).
