# AE4: complete mathematical analysis and model diagnosis

**Completed technical report: 154 pages.**

[Read the full report](AE4_mathematical_analysis_report.pdf)

The report consolidates the completed Codex Task 43 parameter audit and Task 44 full system mathematical analysis. It adds a detailed account of the numerical methods, a careful interpretation of model inadequacy, complete numerical appendices and an unchanged archive of the scientific evidence. It does not construct a new model, alter the original results or claim new experimental validation.

## Reading guide

- Part I: purpose, evidence classes, complete model scope, interpretation of the failed phenotype and expanded computational methods.
- Part II: complete final Task 44 exposition, including all included mathematical fragments, conditional proofs, constitutive laws, conservation identities, source projections, equilibrium and mode analysis, sensitivities, inverse construction and literature comparison.
- Part III: complete final Task 43 parameter provenance report, including all inventory rows and its sensitivity crosswalk.
- Part IV: all 26 equilibrium cases, 338 dimensional state entries, 26 complete Jacobians (3,146 entries), 286 eigenvalues, modal projections, expression derivatives, 208 parameter and metric sensitivity rows, the complete 17 point inverse scan, both inverse sensitivity observables and 137 detailed provenance records.
- Part V: article level conclusions, unsupported claims to avoid, source limitations, original verification receipts, numerical controls and the complete electronic evidence catalogue.

There are 137 inventory records with 109 marked active. These are not 137 independently fitted parameters. The audit distinguishes physical constants, effective settings, initial states, protocol labels, legacy fields and other types of record.

## Files

- `report.tex`: main TeX entry point.
- `interpretation.tex`, `methods.tex`, `conclusions.tex`: new explanatory text.
- `assemble.py`: builds the combined exposition and expanded tables from saved sources only.
- `generated/`: complete expanded TeX, tables, bibliography and coverage record.
- `evidence/task43/`, `evidence/task44/`: 126 unchanged source files, including the original reports, all published numerical payloads, mathematical fragments, scripts and verification receipts.
- `SOURCE_MANIFEST.json`: original commit, path, byte count and SHA256 for every archived source file.
- `COVERAGE.md`: detailed coverage counts and where full precision electronic records are retained.
- `BUILD_VERIFICATION.json`: report compilation, citation, coverage and source integrity checks.
- `VISUAL_REVIEW.md`: rendered page review and comparison of the local and GitHub builds.
- `ASSEMBLY_CORRECTIONS.patch`: audit record of formatting and source schema corrections already applied to `assemble.py`; no separate application is needed to build the final report.

The PDF prints all state vectors, Jacobians and eigenvalues. Full precision modal vectors, every stored nearby trial, diagnostic residual, source location and verification exclusion remain in the accompanying electronic evidence. This preserves detail without presenting rounded PDF tables as the numerical authority.

## Exact scientific inputs

- Task 43: `547f113d9123ab1976774639a69483faf40cef34`.
- Task 44: `e5fa9840bf96be3147c6118daee4942468c78f8b`.
- Original independently replayed numerical checkpoint: `f783785a863440df459e9c5530beba4c7eb59110`.

The original Task 44 clean replay and its supplementary checks are documented, not rerun by this report assembly. The report has its own independent build and source integrity checks. No production model, original output, source analysis branch or manuscript text was changed.

The original BibTeX databases are preserved under `evidence/`. The combined report uses separate source namespaces rather than silently reconcile differing metadata. Identical underlying papers can consequently appear with both original records. The manuscript bibliography can be reconciled separately when this technical report is integrated into an article.

## Build

Requirements: Python 3, pdfLaTeX, BibTeX, latexmk, the LaTeX packages listed in `report.tex`, and Poppler utilities. No scientific Python packages or model reruns are required for this report build.

From a full repository clone containing the two pinned source commits:

```sh
bash report/build.sh
```

From an archive or shallow clone, use the already included exact evidence snapshots. From the repository root:

```sh
bash report/build.sh \
  --source43 "$PWD/report/evidence/task43" \
  --source44 "$PWD/report/evidence/task44"
```

The build checks all preserved source hashes, coverage counts, references, citations and overflow warnings. It does not fit parameters or execute the model.

The published PDF is 999,015 bytes. SHA256:

`15fc598194fdf92c3df3307a7988e7cae64f5582e8acc800ee19b84e3931c68c`

This is a technical basis for a serious article. It distinguishes what is proved, what is numerically demonstrated and what remains biologically unresolved; it does not label an unvalidated model as correct.
