# AE4 Salivary Transport Control

This repository contains a reconstruction of the Solute Carrier Family 9 Anion Transporter 4, in a salivary secretion model, its finite mechanism tests, and the current manuscript being prepared for the *Bulletin of Mathematical Biology*. The earlier project is retained under `archive/` as historical and provenance material.

The scientific starting point is [A Mathematical Model Supports a Key Role for Ae4 (Slc4a9) in Salivary Gland Secretion](https://doi.org/10.1007/s11538-017-0370-6).



## Why does loss of AE4 produce a substantial secretory phenotype while loss of AE2 produces little or no effect, and which parts of the mechanism proposed in 2018 survive a conservation-explicit reconstruction using the intervening experimental literature?

## Current conclusions

The rebuilt model gives a coherent wild type and passes the declared charge, current, carbon, water, positivity and physiology checks over the 600 s production protocol.

The main findings are:

1. **Acid-base supply constrains AE4.** Once AE4 is required to carry a substantial stimulated chloride load, the original NHE1/CO2 architecture cannot replace the bicarbonate and alkalinity exported by AE4. A stimulus-recruited electrogenic NBC-like pathway supplies the missing degree of freedom. In the validated wild type, positive chloride loading over 60-600 s is 79.9248% NKCC1, 20.0752% AE4 and 0% AE2.
2. **Equal AE4 cation routing is insufficient.** With the total AE4 cycle fixed and its cation source split equally between sodium and potassium, complete AE4 loss reduces cumulative secretion by only 3.8575%. NKCC1 compensation rises by 23.16%.
3. **Catalán 2025 stoichiometry alone is insufficient.** Seven predeclared source classes were frozen before phenotype comparison. None recovers the approximately 35% AE4-null secretion phenotype under the inherited scalar AE4 kinetic law. Three carbonate classes fail the stimulated wild-type pH gate.
4. **A constructive regulatory coupling can recover the magnitude.** A target-selected AE4-dependent reduction in stimulated CaCC recruitment gives 23.16% loss at 5% AE4 and 30.26% in the null while passing the declared 600 s gates. This is a modelling hypothesis, not independent biological validation, and it does not reproduce the complete sodium or temporal phenotype.


## Manuscript

The current source-audited article is under `manuscript/`.

Working title:

> **AE4 control of salivary secretion nearly a decade later: acid-base balance, transporter stoichiometry and missing network regulation**

**Current compiled manuscript:** [`manuscript/AE4_manuscript_current.pdf`](manuscript/AE4_manuscript_current.pdf)

The tracked PDF is a convenience copy of the complete manuscript. The LaTeX source and frozen figure values remain authoritative. A GitHub Actions workflow independently rebuilds the figures and manuscript whenever manuscript source changes.

To build locally from the repository root:

```bash
bash manuscript/build.sh
```

The corrected 2018 source archive has been audited against the new article, and its bibliography has been consolidated into `manuscript/references.bib`. The remaining submission work is listed in [`docs/MANUSCRIPT_STATUS.md`](docs/MANUSCRIPT_STATUS.md). Author order, affiliations, funding, contributions and competing-interest statements still require agreement and are deliberately not invented here.

## Repository map

- `archive/`: immutable historical material and provenance notes
- `literature/`: published sources and evidence tables
- `model/`: canonical equations, observables and source mapping
- `src/`: reusable model and analysis code
- `analysis/`: numbered scientific analyses in execution order
- `tests/`: regression, conservation and reproducibility tests
- `results/`: machine-readable outputs used for figures and tables
- `figures/`: project figures outside the article
- `manuscript/`: article source, frozen figure data, build scripts and current compiled PDF
- `docs/`: research plans, decisions and result ledgers
- `prompts/`: reproducible scientific task prompts





## Revised AE4 manuscript, 18 September 2026

The revised article in `manuscript/` centres the recorded chloride availability result, compares it with the original model and later failed/compensating constructions, and discusses Catalán 2015 and 2025. All eight active figures are time series, with EPS copies under `manuscript/figures/`. Build using `bash manuscript/build.sh`; the stable reading copy is `manuscript/AE4_manuscript_current.pdf`. See the manuscript README and M20260918 ledger entry for source scope and explicitly summary-only continuation sensitivities. No frozen science was changed or rerun.
