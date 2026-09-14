# AE4 Salivary Transport Control

This is the working repository for the **AE4 salivary transport project**. It contains the current reconstruction, validation and mechanism analysis built from the published AE4 secretion model, while preserving the earlier AE4 project under `archive/` as historical and provenance material.

The published AE4 model remains the scientific starting point: [A Mathematical Model Supports a Key Role for Ae4 (Slc4a9) in Salivary Gland Secretion](https://link.springer.com/article/10.1007/s11538-017-0370-6).

## Central question

Why does loss of AE4 produce a substantial secretory phenotype while loss of AE2 produces little or no effect, and can a mechanistically constrained transport model reproduce that distinction without fitting the held out genotype phenotype?

The project separates three possibilities.

1. The apparent AE4 dominance is a genuine structural property of the transport network.
2. It is partly or largely an artefact of parameter scaling or transporter allocation in the older model.
3. The phenotype depends on transport structure or whole cell coupling that is absent from the current model.

## Repository map

* `archive/` immutable historical material and provenance notes
* `literature/` published sources, evidence tables and literature audit
* `model/` canonical equations, parameters, observables and source mapping
* `src/` newly written reusable model and analysis code
* `analysis/` numbered research analyses in execution order
* `tests/` regression, conservation and reproducibility tests
* `data/` curated or derived data that may legally be version controlled
* `results/` machine readable outputs used for figures and tables
* `figures/` publication figures and their source notes
* `manuscript/` the new manuscript only
* `docs/` research plan, provenance, decisions and result ledgers
* `prompts/` reproducible Codex task prompts

## Scientific rule

No result from the historical material is automatically a result of this project. New claims must be independently derived from published sources and newly written code, then recorded with reproducibility and provenance.

The approximately 35% experimental AE4 loss secretion reduction is treated as held out phenotype information. It is not a calibration target.

## Current production state

The modern production lineage is consolidated on `main` through **Task 40**.

Tasks 30 to 38 repaired the modern whole cell model and consolidated the active lineage. Task 39 replaced the NKCC1 concentration response core with the source fixed Palk and Benjamin steady law. Task 40 then tested a literal equal Na/K routing rule for AE4 while preserving the inherited total AE4 cycle and every other scientific mechanism.

Task 40 used one WT resting solve and froze that state for WT, 5% AE4 and exact AE4 null 600 s trajectories. All focused tests, numerical gates, physiology gates and conservation checks passed without retries, parameter tuning or genotype specific resting solves.

The final Task 40 results are:

| Quantity | AE4 = 0.05 | AE4 = 0.00 |
| --- | ---: | ---: |
| Na_i difference from WT at 600 s | -0.138 mM | -0.093 mM |
| NKCC1 compensation | +20.82% | +23.16% |
| Final cumulative secretion ratio | 0.965824 | 0.961425 |
| Final cumulative secretion deficit | 3.4176% | 3.8575% |

Na_i rises transiently after AE4 loss but is below WT by 600 s. Intracellular chloride remains lower by about 1.89 mM and 2.20 mM. NBC and NHE1 activity decrease substantially, integrated Na/K pump activity increases slightly, and CaCC export and paracellular return decrease. AE2 supplies less than 0.27% of positive chloride loading.

The secretion deficit persists but decreases late in the simulation. The exact null deficit remains far below the held out experimental reduction of about 35%.

## Current interpretation

The present model is now a stable mechanistic reference chassis, but it does **not** reproduce the magnitude of the held out AE4 loss secretion phenotype.

This negative result has survived multiple changes that were plausible candidates for hidden rescue or failure, including transporter allocation, resting re equilibration, NHE1 repair, electrogenic NBC repair, source fixed NKCC1 kinetics and equal AE4 cation routing. In the current model, compensation remains too effective, principally through NKCC1 and the associated ionic rebalancing.

The result should therefore be read as a localisation result, not as proof that AE4 biology is unimportant. The remaining discrepancy points to missing or misrepresented coupling outside the tested AE4 routing and transporter response structures.

The complete Task 40 audit trail is in `analysis/40_ae4_equal_cation_routing/` and `results/40_ae4_equal_cation_routing/`.
