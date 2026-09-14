# AE4 Salivary Transport Control

This repository contains the current reconstruction, validation and mechanism analysis of AE4 control in salivary acinar fluid secretion. The published 2018 model remains the scientific starting point, while the earlier implementation is preserved under `archive/` as historical and provenance material.

## Central question

Why does loss of AE4 produce a substantial secretory phenotype while loss of AE2 produces little or no effect, and can a mechanistically constrained whole-cell model reproduce that distinction without fitting the held-out genotype phenotype?

The project now separates three claims that were previously entangled:

1. AE4 is an important biological contributor to stimulated chloride loading.
2. The original cation-balance explanation is not robust once wild-type acid–base physiology and source-constrained transporter laws are imposed.
3. The organ phenotype requires additional network regulation or coupling beyond the tested AE4 source stoichiometries.

## Repository map

- `archive/`: immutable historical material and provenance notes
- `literature/`: published sources, evidence tables and literature audit
- `model/`: canonical equations, parameters, observables and source mapping
- `src/`: reusable model and analysis code
- `analysis/`: numbered research analyses in execution order
- `tests/`: regression, conservation and reproducibility tests
- `data/`: curated or derived data that may legally be version controlled
- `results/`: machine-readable outputs used for figures and tables
- `figures/`: project figures and source notes
- `manuscript/`: current article draft, figure plan, BibTeX database and reproducible figure source
- `docs/`: research plan, provenance, decisions and result ledgers
- `prompts/`: reproducible task prompts

## Scientific rule

No result from the historical material is automatically a result of the reconstructed project. New claims must be independently derived from published sources and newly written code, then recorded with reproducibility and provenance.

The approximately 35% experimental AE4-loss secretion reduction is treated as held-out phenotype information. It was not used to calibrate the modern wild-type model or to select among the predeclared Catalán 2025 source classes.

## Current scientific state

The production lineage is consolidated on `main` through Tasks 41 and 42.

### Wild-type reconstruction and Task 40

Source-constrained NHE1 and NKCC1 laws, a stimulus-recruited electrogenic 1 Na : 2 HCO3 pathway, explicit total inorganic carbon and alkalinity balances, complete current closure and water accounting produce a physiological wild-type trajectory. AE4 supplies approximately 20% of positive stimulated chloride loading while NKCC1 remains dominant.

Under equal Na/K routing, AE4 loss causes only a small secretion phenotype:

| Quantity | AE4 = 0.05 | AE4 = 0.00 |
| --- | ---: | ---: |
| NKCC1 compensation | +20.82% | +23.16% |
| Final cumulative secretion deficit | 3.4176% | 3.8575% |

Intracellular chloride falls, but sodium does not remain above wild type. The network compensates too effectively for the cation split alone to explain the experimental phenotype.

### Task 42: Catalán 2025 mechanism classes

Seven fixed source stoichiometries motivated by Catalán et al. (2025) were tested behind a pre-reveal checkpoint. All seven wild-type resting states passed. Eighteen of 21 production trajectories passed; the three carbonate-class wild-type trajectories crossed the pH lower bound and were retained as failures without rescue.

None of the source classes recovered the approximately 35% phenotype under the inherited Task 40 scalar AE4 kinetic law. The NaCl-in/KHCO3-out class increased secretion after deletion because NKCC1 compensation approached 95%. The largest valid deficit was 7.67%, from a class whose forward affinity was opposed at the frozen reference state.

The result does not falsify the Catalán experiments. It shows that stoichiometric source accounting alone is insufficient in the whole-cell network when class-specific AE4 kinetics and regulation remain unknown.

### Task 41: constructive comparison

A separate target-selected construction made stimulated CaCC recruitment depend on AE4 expression. It produced 23.16% lower cumulative secretion at 5% AE4 and 30.26% lower secretion in the null while passing the declared 600 s gates. The required dependency is strong, sodium remains below wild type, and the temporal phenotype is incorrect. Task 41 is therefore a constructive existence result, not biological validation.

## Manuscript

The first complete article draft is in `manuscript/` under the working title:

**AE4 Control of Salivary Secretion Nearly a Decade Later: Acid–Base Constraints, Transporter Compensation and Mechanism Discrimination**

The draft is organised around five main figures, includes the full reconstructed model equations, uses British English, BibTeX and `cleveref`, and has been compiled and visually checked as a 20-page review PDF. The final author list, parameter supplement and journal-template conversion remain to be completed after author revision.
