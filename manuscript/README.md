# Bulletin of Mathematical Biology manuscript draft

Working title:

**AE4 control of salivary secretion nearly a decade later: acid-base balance, transporter stoichiometry and missing network regulation**

This folder contains a complete first manuscript draft built around the final Tasks 37-42 results. It is deliberately arranged by scientific story rather than by task chronology.

## Build

From the repository root:

```bash
python manuscript/scripts/build_figures.py
cd manuscript
latexmk -pdf -bibtex main.tex
```

The manuscript uses BibTeX, `cleveref`, `equation`, and `align`. Equations are numbered only when they are referenced in the text.

## Files

- `main.tex`: title, abstract, Introduction and submission statements
- `model.tex`: complete model definition and numerical protocol
- `results.tex`: results organised around the figures
- `discussion.tex`: interpretation, limitations and experimental predictions
- `references.bib`: draft BibTeX database
- `figure_plan.md`: proposed figure order and the argument carried by each figure
- `claims_and_evidence.md`: manuscript claim ledger
- `data/figure_source_values.csv`: frozen values used by the plotting script
- `scripts/build_figures.py`: deterministic figure builder; it does not rerun the model
- `figures/`: PDF figures for LaTeX and PNG previews

## Before submission

The draft is scientifically complete enough for author editing, but it is not submission-ready. The author list, affiliations, funding, contributions and competing-interest statements need agreement. The Springer/BMB production class should replace the portable `article` class once the author group begins submission formatting. The manuscript should also receive a full line-by-line source audit and a final figure-design pass after the scientific text is agreed.
