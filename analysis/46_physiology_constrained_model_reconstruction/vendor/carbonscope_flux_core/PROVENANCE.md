# Frozen local CarbonScope-derived flux-analysis core

This directory contains a self-contained minimal port of the generic FBA/FVA/VFFVA substrate required for the Task 46 AE4 constraint-based screening stage.

Source repository consulted read-only: `esig626/CarbonScope`

Pinned source commit:

`a11e25f176a2cc70de22c3f36a1c0e623248393a`

Relevant pinned source components used as the implementation reference were:

- `src/fluxemu/flux_analysis/highs.py`
- `src/fluxemu/flux_analysis/results.py`
- `src/fluxemu/model/schema.py`
- `src/fluxemu/model/validation.py`
- `src/fluxemu/exceptions.py`
- `THIRD_PARTY_NOTICES.md`

The AE4 copy is intentionally reduced to the functionality needed for the small transporter CBM: an immutable flux schema, validation, HiGHS FBA, cold-start reference FVA, reusable-worker VFFVA, result records and numerical checks. Isotope, EMU, MFA, hypothesis-testing, workflow, CLI and experiment code were not brought across.

This is therefore an adapted local port, not a byte-for-byte copy of the full CarbonScope modules. The scientific solver ideas and relevant implementation details are pinned to the source commit above. Any further modifications must occur only in the AE4 Task 46 branch and must be documented here or in Task 46 checkpoint records.

## Isolation rule

The Task 46 workflow must use this local core. It must not write to, branch, commit, push, open pull requests in, or otherwise modify the CarbonScope repository. CarbonScope may be consulted read-only only if source comparison is genuinely needed.

Runtime dependencies are NumPy, pandas and `highspy`; the pinned CarbonScope source specifies `highspy>=1.11,<1.13`.

The VFFVA worker design preserves the third-party attribution and MIT notice in `THIRD_PARTY_NOTICES.md`.
