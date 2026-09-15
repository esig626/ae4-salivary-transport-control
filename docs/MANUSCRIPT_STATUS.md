# Manuscript status

## Working title

**AE4 control of salivary secretion nearly a decade later: acid-base balance, transporter stoichiometry and missing network regulation**

## Current state

The repository contains a complete source-audited manuscript under `manuscript/`, together with a reproducible build and a tracked compiled reading copy at `manuscript/AE4_manuscript_current.pdf`.

The article contains:

- an Introduction centred on the 2018 paper and the intervening experimental work;
- the complete conservation-explicit model;
- the NBC and pH reconstruction;
- Task 40 equal-routing results;
- the Task 41 target-selected CaCC construction;
- the predeclared Task 42 Catalán mechanism-class test;
- a Discussion organised around what survived, what failed and what remains unresolved;
- an explicit scope and limitations section;
- six reproducible figures built only from frozen source values;
- a consolidated BibTeX bibliography;
- a claim-and-evidence ledger.

The corrected 2018 article archive has been audited against the current manuscript. That audit restored historical constraints that were under-emphasised in the first draft, including the original approximately 24% modelled AE4-null deficit, the roughly 2 min versus roughly 10 min timing discrepancy, minimal NKCC1 compensation in the 2018 model, and the isolated 2015 NKCC1/NHE observations. The audit also sharpened the distinction between primary acinar secretion and downstream ductal modification. Details are recorded in `manuscript/corrected_archive_audit.md`.

The old and new bibliographies have been consolidated to 56 unique entries with duplicate legacy keys collapsed and malformed legacy metadata corrected where verified.

## Scientific story

1. The 2018 model identified AE4 as the correct biological actor and cation dependence as an important property.
2. Rewriting acid-base chemistry in conserved TIC and TA coordinates shows that an additional inward bicarbonate source is required if AE4 is to carry a substantial stimulated chloride load while pH remains physiological.
3. With reconstructed NHE1, NKCC1 and NBC-like transport included, equal AE4 cation routing produces only a small secretion deficit because NKCC1 compensates strongly in the model.
4. The exact equal-routing balance removes the explicit AE4 term from sustained apical chloride export.
5. The Catalán 2025 source classes do not recover the phenotype under the inherited whole-cell kinetic law, and three carbonate classes fail stimulated wild-type pH.
6. A target-selected AE4-dependent CaCC coupling recovers phenotype magnitude but not the full ionic or temporal signature.
7. The strong Task 40 NKCC1 compensation is a model prediction in experimental tension with the isolated 2015 assay, not an experimentally established knockout response.

## Build and reproducibility

The canonical local build is:

```bash
bash manuscript/build.sh
```

`.github/workflows/build-manuscript.yml` independently regenerates the six figures, builds the complete LaTeX/BibTeX manuscript and uploads the resulting PDF as a workflow artifact whenever manuscript source changes.

No model source, frozen trajectory, simulation checkpoint or figure source value was altered by the corrected-archive audit or the submission-facing repository cleanup.

## Remaining submission work

- confirm the complete author list, order and affiliations;
- confirm funding, author contributions and competing-interest statements;
- decide whether the target-selected Task 41 result remains fully in the main text or moves partly to supplementary material;
- perform the final figure-design pass;
- migrate the portable `article` draft to the current Springer/BMB submission class;
- perform the final line-by-line numerical, source and citation audit of the submission version;
- write the cover letter and final author-facing claim summary.

These items require scientific or author decisions and are not silently filled in by repository housekeeping.
