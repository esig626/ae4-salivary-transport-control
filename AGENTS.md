# Active phase: manuscript development

`main` is the canonical scientific lineage and now contains the completed analyses through Task 42.

The active writing branch is:

`manuscript/bmb-ae4-decade-reassessment`

## Primary task

Develop the *Bulletin of Mathematical Biology* manuscript under `manuscript/` without changing frozen scientific outputs.

Read in this order:

1. `README.md`
2. `manuscript/README.md`
3. `manuscript/figure_plan.md`
4. `manuscript/claims_and_evidence.md`
5. `analysis/37_wt_nbc_validation/final_answer.md`
6. `analysis/40_ae4_equal_cation_routing/final_answer.md`
7. `analysis/41_ae4_loss_algebraic_design/final_answer.md`
8. `analysis/42_catalan_2025_ae4_mechanism_classes/final_answer.md`

## Writing constraints

- Use British English.
- Preserve the distinction between source-derived experimental findings, model results and target-selected constructions.
- Do not claim that Catalán et al. 2025 established the proposed whole-cycle stoichiometries experimentally.
- Do not describe Task 41 as independent validation.
- Do not identify the NBC-like structural pathway with a specific salivary isoform unless new evidence is added.
- Keep the full reconstructed model in the manuscript. Use `equation` and `align`; number only equations that are referenced and use `cleveref` for cross-references.
- Maintain BibTeX references rather than a handwritten bibliography.
- Manuscript figures must be generated from frozen repository values by `manuscript/scripts/build_figures.py`. Do not rerun the scientific model merely to redraw a figure.

## Scientific state that must remain unchanged

- Task 40 equal-routing null deficit: 3.8575%.
- Task 41 target-selected null deficit: 30.2612%.
- Task 42 pre-reveal checkpoint: `4cb75308bd81d06f064bb80b11831985fa4bb219`.
- Task 42 conclusion: none of the seven source classes recovers the held-out phenotype under inherited Task 40 kinetics.
- Task 41 and Task 42 failures, limitations and thermodynamic qualifications must remain visible.

## Repository constraints

- Do not rewrite or delete `archive/`, numbered analyses or frozen results.
- Keep LaTeX build products out of version control.
- New scientific calculations require a new numbered task and branch. Manuscript editing alone must not alter model code or result files.
- Do not merge a manuscript change until the LaTeX source compiles, citations and references resolve, and the figures render without clipping.
