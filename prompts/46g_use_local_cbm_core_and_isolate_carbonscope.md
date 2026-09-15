# Task 46G: use the local AE4 CBM core; CarbonScope is read-only reference only

This addendum is mandatory and supersedes any earlier instruction that would require Task 46 to use `esig626/CarbonScope` as a working repository or runtime dependency.

## Hard repository-isolation rule

Do not modify `esig626/CarbonScope` in any way.

Do not:

- create or switch branches there;
- commit or push there;
- open pull requests there;
- edit or delete files there;
- run cleanup, reset, checkout, rebase or merge operations there;
- change its issues, workflows, tags or releases;
- use it as the Task 46 working tree.

The only permitted CarbonScope interaction is read-only source comparison when genuinely necessary.

## Local solver core

Task 46 now contains a self-contained local CarbonScope-derived FBA/FVA/VFFVA core at:

`analysis/46_physiology_constrained_model_reconstruction/vendor/carbonscope_flux_core/`

Its provenance is pinned to CarbonScope commit:

`a11e25f176a2cc70de22c3f36a1c0e623248393a`

Read:

- `README.md`
- `PROVENANCE.md`
- `THIRD_PARTY_NOTICES.md`

before using it.

Use this local core for the Task 46 CBM structural-screening stage. If additional generic LP/FVA functionality is required, implement it locally in the Task 46 branch, add tests, document the change and publish a checkpoint. Do not fix or extend CarbonScope itself as part of this task.

Before scientific CBM calculations:

1. install the local dependencies from its `requirements.txt` if necessary;
2. run the local smoke test;
3. construct one tiny independent analytical control LP and verify FBA plus reference-FVA/VFFVA agreement;
4. record that verification in the Task 46 checkpoint.

The AE4 transporter network, paired WT/KO formulation, flux coupling and sparse minimal-relaxation analyses belong in Task 46 code above this generic local solver layer.

This isolation rule has priority over Task 46F or any earlier wording suggesting direct operational reuse of the CarbonScope repository.
