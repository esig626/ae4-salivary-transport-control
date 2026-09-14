# Manuscript: AE4 control nearly a decade later

This folder contains the first complete article draft built around the final Task 37, 40, 41 and 42 results.

## Working title

**AE4 Control of Salivary Secretion Nearly a Decade Later: Acid–Base Constraints, Transporter Compensation and Mechanism Discrimination**

The draft is written in British English and follows a figure-led structure. The scientific distinction between source-constrained prediction and target-selected construction is preserved throughout.

## Files

- `main.tex`: complete first manuscript draft, including the reconstructed model equations.
- `references.bib`: BibTeX database.
- `FIGURE_PLAN.md`: the intended story for each main and supplementary figure.
- `REVISION_CHECKLIST.md`: unresolved scientific, authorship and submission items.
- `make_figures.py`: deterministic script used to create the current figures.
- `figures/`: generated PDF figures and PNG previews. These are recreated by `make_figures.py` from the frozen results and are not treated as primary source files.
- `manuscript.pdf`: locally compiled review copy. It is supplied separately for author review rather than tracked as a source file.

## Build

From this directory:

```bash
python make_figures.py
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

The current draft and all five figures were compiled and visually checked locally. The source uses a deliberately plain `article` class so it compiles without journal-specific dependencies while the scientific text and figures are being edited. Before submission, convert `main.tex` to the current Springer Nature `sn-jnl` template and flatten figures into the submission directory. The source is already kept as one main TeX document, uses BibTeX, `cleveref`, `equation` and `align`, and numbers only equations that are cited in the text.

## Scientific status

This is a genuine first draft, not a submission-ready final article. In particular:

- the final author list and order must be agreed;
- the precise experimental comparator and uncertainty should be checked against the original data;
- the full parameter table should be moved from the repository audit into a publication supplement;
- figure styling and panel selection should be revised after the first author edit;
- Task 41 remains a constructive, target-selected comparison and must not be presented as independent validation;
- Task 42 does not falsify Catalán et al. (2025), because their experiments do not determine a complete salivary whole-cell kinetic law.
