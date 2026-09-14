# AE4 Salivary Transport Control

This repository contains a source-audited reconstruction of the AE4 salivary secretion model, its finite mechanism tests, and the manuscript now being prepared for the *Bulletin of Mathematical Biology*. The earlier project is retained under `archive/` as historical and provenance material.

The scientific starting point is [A Mathematical Model Supports a Key Role for Ae4 (Slc4a9) in Salivary Gland Secretion](https://doi.org/10.1007/s11538-017-0370-6).

## Scientific question

Why does loss of AE4 produce a substantial secretory phenotype while loss of AE2 produces little or no effect, and which parts of the mechanism proposed in 2018 survive a conservation-explicit reconstruction using the intervening experimental literature?

## Current conclusions

The rebuilt model gives a coherent wild type and passes the declared charge, current, carbon, water, positivity and physiology checks over the 600 s production protocol.

The main findings are:

1. **Acid-base supply constrains AE4.** Once AE4 is required to carry a substantial stimulated chloride load, the original NHE1/CO2 architecture cannot replace the bicarbonate and alkalinity exported by AE4. A stimulus-recruited electrogenic NBC-like pathway supplies the missing degree of freedom. In the validated wild type, positive chloride loading over 60-600 s is 79.9248% NKCC1, 20.0752% AE4 and 0% AE2.
2. **Equal AE4 cation routing is insufficient.** With the total AE4 cycle fixed and its cation source split equally between sodium and potassium, complete AE4 loss reduces cumulative secretion by only 3.8575%. NKCC1 compensation rises by 23.16%.
3. **Catalán 2025 stoichiometry alone is insufficient.** Seven predeclared source classes were frozen before phenotype comparison. None recovers the approximately 35% AE4-null secretion phenotype under the inherited scalar AE4 kinetic law. Three carbonate classes fail the stimulated wild-type pH gate.
4. **A constructive regulatory coupling can recover the magnitude.** A target-selected AE4-dependent reduction in stimulated CaCC recruitment gives 23.16% loss at 5% AE4 and 30.26% in the null while passing the declared 600 s gates. This is a modelling hypothesis, not independent biological validation, and it does not reproduce the complete sodium or temporal phenotype.

The resulting interpretation is narrow but useful: the biological importance of AE4 survives the decade, but cation stoichiometry alone does not explain the knockout. Acid-base supply, NKCC1 compensation and sustained regulation of apical chloride exit remain the decisive unresolved mechanisms.

## Manuscript

The first complete article draft is under `manuscript/`.

Working title:

> **AE4 control of salivary secretion nearly a decade later: acid-base balance, transporter stoichiometry and missing network regulation**

The manuscript is organised around six figures and includes the full reconstructed model, BibTeX references, `cleveref` cross-references, a claim ledger and deterministic figure-building code.

Build it with:

```bash
bash manuscript/build.sh
```

The author list, affiliations, funding and contribution statements still require agreement before submission.

## Repository map

- `archive/`: immutable historical material and provenance notes
- `literature/`: published sources and evidence tables
- `model/`: canonical equations, observables and source mapping
- `src/`: reusable model and analysis code
- `analysis/`: numbered scientific analyses in execution order
- `tests/`: regression, conservation and reproducibility tests
- `results/`: machine-readable outputs used for figures and tables
- `figures/`: project figures outside the article draft
- `manuscript/`: article source, figures, data and build scripts
- `docs/`: research plans, decisions and result ledgers
- `prompts/`: reproducible scientific task prompts

## Scientific rules

- Historical outputs are not automatically results of the reconstructed project.
- The approximately 35% experimental AE4-loss secretion reduction is not a calibration target unless a task is explicitly labelled target-selected.
- Numerical failure, physiological failure and biological mismatch are reported separately.
- Failed mechanism classes are retained rather than removed after comparison.
- `archive/` is not rewritten.
