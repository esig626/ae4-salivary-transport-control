# Frozen CarbonScope flux-analysis snapshot

This directory vendors only the generic FBA/FVA solver substrate required for the Task 46 AE4 constraint-based screening stage.

Source repository: `esig626/CarbonScope`

Pinned source commit:

`a11e25f176a2cc70de22c3f36a1c0e623248393a`

The Task 46 workflow must use this local snapshot and must not write to, branch, commit, push, open pull requests in, or otherwise modify the CarbonScope repository.

The intended copied source files are:

- `src/fluxemu/flux_analysis/highs.py`
- `src/fluxemu/flux_analysis/results.py`
- `src/fluxemu/model/schema.py`
- the flux-model validation subset from `src/fluxemu/model/validation.py`
- `src/fluxemu/exceptions.py`
- `THIRD_PARTY_NOTICES.md`

Only generic flux-model / HiGHS FBA/FVA functionality is needed. Do not copy CarbonScope isotope, EMU, MFA, hypothesis-testing, workflow, CLI, or experiment code.

Runtime dependencies for the vendored core are NumPy, pandas, and `highspy`; the pinned CarbonScope version specifies `highspy>=1.11,<1.13`.

The VFFVA implementation preserves the third-party attribution and MIT notice copied in `THIRD_PARTY_NOTICES.md`.
