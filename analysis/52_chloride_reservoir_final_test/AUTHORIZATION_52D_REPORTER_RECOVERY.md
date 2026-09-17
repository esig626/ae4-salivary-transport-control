# Task52D reporting-only recovery authorisation

**User-authorised recovery scope:** continue Task52D only from the verified immutable 52C scientific freeze. The prior fresh 52D attempt failed before any verified case result was saved because the output recorder assumed identical onset/active diagnostic column sets. The failure evidence at commit `9c0fde42b5eb16746413801fae555c84cb32e201` remains preserved and must not be deleted or rewritten.

This file provides the explicit authority required by `NUMERICAL_FREEZE.md` / D52-15 for one bounded software correction and a fresh 52D attempt.

## What is authorised

1. Change **reporting/serialisation code only** so the recorder can preserve the union of onset and active diagnostic fields when some named diagnostics are absent at a sample.
2. Preserve every original named diagnostic. Use an explicit presence mask / missing-field representation; do not invent a physical value for a field that is absent at a sample.
3. The correction may touch only output-recording / schema-normalisation / stage-replay bookkeeping needed to permit this authorised rerun. It must not change scientific equations, model evaluation, diagnostic definitions, physical/conservation gates, solver settings, tolerances, quadrature, onset states, case definitions, transport laws or parameters.
4. Re-run Task52D from the immutable 52C inputs, exactly the three predeclared combined CCh+IPR cases:
   - case_01: matched supply, central projection, `G_aux=0 S`;
   - case_02: matched supply, central projection, `G_aux=2.32e-9 S`;
   - case_03: matched supply, central projection, `G_aux=4.49e-9 S`.
5. Use the exact 52C freeze SHA `9cebb8fcf85557643f6c7838ff57f9fb9603f5a7479cbe6fb63ba8c40b417ba2`, the exact frozen onset states, one scientific process / one BLAS thread, the frozen 600 s Radau settings, quadrature and all original gates.
6. Preserve the exact chloride-reservoir audit for each case:
   `D_J = delta(0) + 2D_N + D_A + D_E - delta(600)`.
7. After all three cases finish, create the normal 52D summaries/checkpoint/receipt, commit and push them immediately, independently verify the remote SHA/tree/output hashes, then STOP.

## What is not authorised

- no scientific-source or frozen-input change;
- no Palk/Benjamin NKCC;
- no Task41/Task50 coupling;
- no parameter change, fit, tuning, search or extra case;
- no altered physical/conservation tolerance;
- no changed diagnostic definition to make a case pass;
- no BDF unless the original frozen rule independently permits it for an explicit Radau solver-status failure;
- no 52E;
- no 52F;
- no Task53.

## Interpretation

The reporting correction is a software-compatibility repair only. The new trajectories are a fresh execution of the already frozen 52D scientific experiment, not a new hypothesis or a retuned experiment. The failed first attempt remains part of the audit trail.
