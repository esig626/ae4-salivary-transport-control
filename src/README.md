# Source code

All code here must be newly written for this project from the published model specification in `model/`.

The implementation should eventually expose separate components for parameters and units, transporter and channel fluxes, state equations, steady state solving, continuation, observables, sensitivity and control coefficients, Jacobian decomposition, mechanism variants, uncertainty, and identifiability.

Do not port source code from the historical material line by line.

## Phase 10 audit code

`ae4_identifiability.py` contains independently written stoichiometric identities,
regular-rank utilities, and diagnostics obtained by direct substitution into the
published equations. It is deliberately **not** a full simulator: the published
baseline does not close under the printed water, current, acid/base, and unit
conventions.

`run_identifiability_audit.py` regenerates the checked Phase 10 JSON/CSV outputs
and the structural flux-geometry figure. Run it from the repository root with

```bash
python src/run_identifiability_audit.py
```

`identifiability_discrimination/` is a separate affine software fixture already
present on the working branch. It verifies derivative and rank code under
declared matrices; its `Q_star`, panel ranks, and geometry are not published
model outputs and are not used as evidence in the `STOP` decision.

Install the small analysis environment from `requirements-analysis.txt`.
