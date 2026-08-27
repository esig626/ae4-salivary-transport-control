# Phase 10: identifiability and discrimination

## Outcome

The requested quick-paper analysis terminates with **`STOP`** under the
currently published model specification. The branch contains checked
structural results and printed-equation diagnostics, but no full-model
baseline reproduction or physiological minimal-panel claim.

Read the files in this order:

1. `reduction.md` — independent symbolic reduction and reconstruction gate;
2. `theory.md` — regular local theory and exact conditional identities;
3. `observables.md` — measurement hierarchy and certification boundary;
4. `escalation_scan.md` — mandatory global, perturbation, design, testing, and
   singularity scan;
5. `novelty_audit.md` — focused comparison with primary literature; and
6. `viability.md` — objective paper decision and restart conditions.

## Reproduction

From the repository root:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test*.py' -v
python src/run_identifiability_audit.py
```

The generated JSON/CSV files are in
`results/10_identifiability_discrimination/`; the corresponding exact
flux-coordinate illustration is
`figures/identifiability_stoichiometric_geometry.png`.

The generated geometry is not a surrogate physiological model. It visualizes
only the certified linear aggregate map

\[
(J_2,J_4)\mapsto(J_2+J_4,J_2+2J_4).
\]

The pre-existing `identifiability_discrimination` package is retained only as
an affine derivative/rank software fixture. It can be regenerated with

```bash
PYTHONPATH=src python -m identifiability_discrimination.run_analysis
```

Its proxy outputs, panel ranks, and figures are explicitly uncalibrated and do
not enter the physiological or publication conclusions.
