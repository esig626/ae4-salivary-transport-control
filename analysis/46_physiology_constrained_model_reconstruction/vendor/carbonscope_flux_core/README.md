# Task 46 local CBM solver core

This directory contains a self-contained minimal FBA/FVA/VFFVA solver layer for the AE4 Task 46 constraint-based screening stage.

It is adapted from the generic flux-analysis layer of `esig626/CarbonScope` pinned at commit `a11e25f176a2cc70de22c3f36a1c0e623248393a`, but it is intentionally much smaller and contains no isotope, EMU, MFA, hypothesis-testing, CLI or workflow code.

## Isolation rule

The AE4 Task 46 workflow must use this local directory only. It must not clone CarbonScope as a working repository, must not create branches there, and must not commit, push, open pull requests, edit files or otherwise mutate `esig626/CarbonScope`.

If additional solver functionality is needed, implement it in the AE4 Task 46 branch and record the provenance. CarbonScope may be consulted read-only only if a source comparison is genuinely necessary.

## Runtime

Install the local requirements:

```bash
python -m pip install -r analysis/46_physiology_constrained_model_reconstruction/vendor/carbonscope_flux_core/requirements.txt
```

Make the vendor parent directory importable, then run:

```bash
python -m analysis.46_physiology_constrained_model_reconstruction.vendor.carbonscope_flux_core.smoke_test
```

If module-path execution is awkward because the numbered analysis directory is not a Python package, add the vendor directory's parent to `PYTHONPATH` and import `carbonscope_flux_core` directly.

## Intended use

Use this core for:

- FBA feasibility and objective solves;
- cold-start reference FVA;
- reusable-worker HiGHS VFFVA;
- comparison of reference and fast FVA before production use.

AE4-specific paired WT/KO constraints, flux coupling and sparse minimal-relaxation calculations belong in Task 46 analysis code above this generic layer.
