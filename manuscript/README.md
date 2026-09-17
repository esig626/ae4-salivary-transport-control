# AE4 manuscript

Revised article: **AE4 control of salivary secretion nearly a decade later: chloride availability, conservation and the limits of transporter compensation**.

The source baseline for this manuscript-only revision is scientific main `764e32a648ee3b69265364788df5a59f76cb2843`, with existing manuscript decisions retained. No scientific task is reopened and no model is evaluated to build the article.

## Build and deliverables

Run `bash manuscript/build.sh` from the repository root. It reads recorded arrays and CSVs, regenerates eight time series figures, compiles LaTeX/BibTeX and checks unresolved references. The reading copy is `AE4_manuscript_current.pdf`. Editable sources are `main.tex`, `model.tex`, `results.tex`, `discussion.tex` and `references.bib`. `AE4_manuscript_standalone.tex` is a flattened source; it uses the same figures and bibliography.

Every active figure has EPS, PDF and PNG versions in `figures/`. `data/time_series/` contains exact plotting values, input SHA256 hashes and a file-only arithmetic audit. There are no active bar plots or schematic replacement figures. Historical endpoint-only outcomes remain tables.

## Evidence scope

The verified central paired trajectories show 20.04–20.10% cumulative secretion deficits and independently closed chloride budgets. Historical figures use saved outputs, including the accepted `candidate_08` benchmark rather than intermediate root-level trials. The original model comparison uses its corrected archive record. The article distinguishes measured findings, modelling assumptions, exact identities, numerical predictions, calibrated constructions and failed or unidentified models.

The additional supply/projection/single-agonist outcomes are reported continuation summaries preserved in `revision_task52/DECISION_LOG.md`. Their original raw arrays are not in this snapshot; they are not reverified or plotted. See `reported_continuation_summary.md` and the article's data availability statement. These raw files, the complete author list, affiliations, contributions, funding and disclosures require author confirmation/deposition before submission.

The generated PDF is a complete revised reading manuscript, not a claim that those submission requirements have been met.
