# AE4 Salivary Transport Control

This is the working repository for the **AE4 salivary transport project**. It contains the current reconstruction, validation and mechanism analysis built from the published AE4 secretion model, while preserving the earlier AE4 project under `archive/` as historical and provenance material.

The published AE4 model remains the scientific starting point: [A Mathematical Model Supports a Key Role for Ae4 (Slc4a9) in Salivary Gland Secretion](https://link.springer.com/article/10.1007/s11538-017-0370-6).

## Central question

Why does loss of AE4 produce a substantial secretory phenotype while loss of AE2 produces little or no effect, and does that distinction survive a scale invariant analysis of transporter control?

The project separates three possibilities.

1. The apparent AE4 dominance is a genuine structural property of the transport network.
2. It is partly or largely an artefact of parameter scaling in the older sensitivity calculation.
3. The phenotype depends specifically on AE4 monovalent cation coupling and can therefore discriminate among alternative microscopic transport mechanisms.

## Repository map

- `archive/` immutable historical material and provenance notes
- `literature/` published sources, evidence tables, and literature audit
- `model/` canonical equations, parameters, observables, and source mapping
- `src/` newly written reusable model and analysis code
- `analysis/` numbered research analyses in execution order
- `tests/` regression, conservation, and reproducibility tests
- `data/` curated or derived data that may legally be version controlled
- `results/` machine readable outputs used for figures and tables
- `figures/` publication figures and their source notes
- `manuscript/` the new manuscript only
- `docs/` research plan, provenance, decisions, and result ledgers
- `prompts/` reproducible Codex task prompts

## Scientific rule

No result from the historical material is automatically a result of this project. New claims must be independently derived from published sources and newly written code, then recorded with reproducibility and provenance.

## Current status

Task 14, the scale free genotype holdout validation, is complete and merged to `main`. The scientific model state remained frozen at commit `2f54e7c4f87b6da746d3a8427c7bdca40a77c3de` during the blind experiment.

All ten AE4 continuation branches reached 5% expression but the declared exact zero continuation failed numerically in every branch, so no valid exact AE4 null secretion trajectory was produced. All ten AE2 continuations succeeded, and AE2 deletion predicted only a 0.03% to 0.16% secretion increase. The blind checkpoint was committed before phenotype reveal at `607ac21279d821b95d44ba63f35df27078503cf6`.

After reveal, the retained experiment reported 35 ± 4.7% lower total AE4 null secretion. Because no exact zero model trajectory existed, Task 14 did not establish either agreement or disagreement for the AE4 phenotype. The final classification was **`GENOTYPE CONTINUATION OR NUMERICS PREVENT A VALID HOLDOUT TEST`**. See `analysis/14_scale_free_genotype_holdout/` and `results/14_scale_free_genotype_holdout/` for the complete audit trail.

The next planned diagnostic is a post reveal, frozen model calculation using the already valid 5% AE4 states. It is not a second blind holdout and must not be used to retune the model.
