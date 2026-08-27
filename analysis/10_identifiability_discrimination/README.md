# Phase 10 identifiability/discrimination audit

This directory executes `prompts/10_identifiability_discrimination_quickpaper.md` without modifying the historical archive or starting a manuscript.

The audit separates two levels of result:

1. analytical statements that hold for any regular steady-state model `F(u; G2, G4) = 0`; and
2. an explicitly uncalibrated, baseline-anchored stoichiometric diagnostic used to test code paths and visualize rank geometry.

The published 2018 model could not be independently reproduced from the citable specification. Therefore, no numerical result from the diagnostic is presented as a salivary-gland prediction. The final classification is recorded in `viability.md`.

Reproduce the executable checks from the repository root with:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python -m identifiability_discrimination.run_analysis
```

Generated data are in `results/10_identifiability_discrimination/`; figures are in `figures/10_identifiability_discrimination/`.
