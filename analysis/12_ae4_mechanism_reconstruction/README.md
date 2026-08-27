# Phase 12: AE4 mechanism reconstruction

## Outcome

**`AE4 ALONE INSUFFICIENT; MISSING PATHWAY LOCALIZED`** on the best available
internally coherent seven-state comparison chassis.

The result is a chassis-scoped exclusion, not a claim that AE4 is biologically
unimportant. Primary experiments require transported Na and K and show
beta/cAMP/PKA activation, but the historical chassis was calibrated around a
Na-only AE4 source. Source-justified mixed-cation, cooperative, allosteric,
reversible, saturating, alternative-stoichiometry, and stimulus-gated
candidates cannot satisfy the WT resting constraints and the held-out
knockout phenotypes together. Exact balance analysis localizes the first
missing model property to a two-dimensional cation-homeostasis closure module;
it does not validate a phenotype-tuned pump or NHE1 correction.

## Reading order

1. `experimental_evidence.md` — primary knockout, transport, and PKA evidence;
2. `chassis.md` — fixed non-AE4 equations, provenance, dimensions, and checks;
3. `ae4_candidate_family.md` — complete C1–C8 registry and rejected/unresolved
   modes;
4. `model_comparison.md` — knockout-blind calibration and iterative numerical
   comparisons;
5. `adversarial_audit.md` — independent hard-gate and solver audit;
6. `structural_explanation.md` — exact balance and thermodynamic explanation;
7. `final_answer.md` — controlling answers to the twelve Task-12 questions.

## Reproduce

From the repository root:

```bash
PYTHONPATH=src python -m ae4_mechanism_reconstruction.run_analysis
PYTHONPATH=src python -m unittest discover -s tests -p 'test*.py' -v
```

Machine-readable tables and metadata are in
`results/12_ae4_mechanism_reconstruction/`; generated diagnostics are in
`figures/12_ae4_mechanism_reconstruction/`. Absolute physical time and flow
remain uncertified, so the reported secretion outputs are dimensionless
within-protocol ratios. Nothing under `archive/` is read at runtime or changed
by this analysis.
