# Post-freeze genotype interpretation

The WT-only candidate checkpoint was pushed and independently read back from the remote branch at `f5df88448b1a669c7c7b3986d83ac17360b0338d`. Frozen manifest SHA-256 is `44e98d3f9a67857a93cd88c6886f3bf0958e6a27ba52cc823ca1c5440bbcc84f`. Only afterward did the separate post-freeze evaluator read the genotype results. Every checkpoint hash remains unchanged.

The frozen set contains thirty inherited-baseline root/calcium combinations and zero modified combinations. Consequently, no new AE4 expression continuation, AE4 ODE, AE2 deletion continuation or AE2 ODE is licensed. No exact-zero AE4 attempt is made. No old Task 14B 5% state is used with changed capacities.

The baseline comparisons reuse all thirty Task 14B 5% AE4 production Radau trajectories and their exact matched Task 14 WT denominators. File hashes, parameter hashes, the exact saved 5% starting state, genotype expression, solver, time grid and integrated secretion arithmetic are checked. The validated baseline AE2 ratio table is an additional input outside the 187-file sign/slip inventory; its Git blob is independently matched to Task 14B and its SHA-256 is recorded in `genotype_results_frozen.json`. No baseline trajectory is recomputed.

| Routing family | Calcium (uM) | 5% AE4 secretion change (%) | AE2 loss secretion change (%) |
|---|---:|---:|---:|
| AE4NA05 | 0.10 | +14.0521 to +16.3180 | +0.12546 to +0.13919 |
| AE4NA05 | 0.25 | +17.9915 to +20.2714 | +0.14365 to +0.15598 |
| AE4NA05 | 0.50 | +19.1979 to +21.3305 | +0.14888 to +0.16008 |
| AE4NA20 | 0.10 | +8.2287 to +12.4933 | +0.03180 to +0.03868 |
| AE4NA20 | 0.25 | +14.5430 to +19.5148 | +0.03926 to +0.04616 |
| AE4NA20 | 0.50 | +16.6968 to +21.6630 | +0.04165 to +0.04821 |

The reported change is 100*(R-1); the required reduction fraction is D=1-R. Thus baseline D_AE4 is -0.216630 to -0.082287 and D_AE2 is -0.001601 to -0.000318. The AE4-minus-AE2 reduction contrast is -0.216148 to -0.081969. AE2 is comparatively neutral, although its small positive secretion change is not exactly zero. Baseline 5% AE4 increases secretion in every one of the thirty matched cases, across both routing families and all calcium values.

Each output CSV has ninety rows: thirty valid baseline comparisons and sixty explicit `NOT_EVALUATED_WT_INADMISSIBLE` rows. The latter have blank genotype ratios, not zeros and not borrowed baseline effects. The reference target and both licensed capacity scales remain visible even when the endpoint is unavailable. No frozen candidate is removed after genotype evaluation.

All three result tables were hashed and frozen before adding experimental context. The approximately 35% experimental AE4-null secretion reduction is now mentioned only as contextual direction and magnitude from the Task 16 directive. It is not an exact 5%-expression target, uncertainty band, acceptance criterion or ranking score. It played no role in candidate generation, WT admission or numerical stopping. No minute-by-minute null shape is scored or fitted.

There is no supported comparison of genotype effect across the modified shares, and no supported sign-change interval, because neither modified endpoint yielded an admissible WT state. Simple capacity allocation has supplied no admissible repair in this declared study. Increased cation traffic and cancellation on the saved continuation paths make the mixed-cation architecture a reasonable subject for a later mechanism task, but Task 16 does not establish that slip causes the wrong phenotype or that changing the architecture is mathematically necessary. No mechanism change or further share search is performed here.
