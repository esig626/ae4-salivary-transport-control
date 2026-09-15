# Task 46F addendum: use CarbonScope as the CBM/FVA computational substrate

This addendum is mandatory for the CBM-first Task 46 continuation.

## Pinned external codebase

Use repository `esig626/CarbonScope` as the preferred CBM/FBA/FVA engine at the exact pinned commit:

`a11e25f176a2cc70de22c3f36a1c0e623248393a`

Treat CarbonScope as read-only scientific infrastructure. Do not modify CarbonScope during Task 46.

Relevant native APIs include the canonical `FluxModel` representation and the HiGHS-native flux-analysis layer:

- `compile_flux_lp`
- `prepare_highs_flux_region`
- `run_highs_fba`
- `run_highs_fva_reference`
- `run_highs_vffva`
- `run_prepared_highs_vffva`
- `sample_highs_flux_states`
- `sample_prepared_flux_states`
- `validate_flux_states`

The Task 46 CBM should use CarbonScope directly where practical rather than reimplementing an LP/FVA engine.

## What to reuse

Reuse the solver and feasible-region machinery, not the metabolic/isotope biology.

Construct a small AE4 salivary-acinar `FluxModel` programmatically using only the transporter/network reactions required for the Task 46 structural question. Do not build a genome-scale metabolic model and do not invoke EMU/isotope machinery.

The CBM should retain at minimum the conserved balances needed for the structural screen, including Na, K, Cl, total inorganic carbon/alkalinity equivalents where representable as linear bookkeeping coordinates, and any explicit charge/pump/apical-export bookkeeping required by the chosen steady-state formulation.

Use explicit reaction bounds and shared WT/KO capacity constraints rather than hidden objectives.

## Scientific use of FVA

CarbonScope's own semantics are binding here: FVA endpoints are diagnostics only. Never combine independently optimised FVA coordinates into a fictitious feasible flux vector.

For any conclusion requiring a complete flux state, use either:

- an explicitly solved LP with the full joint constraints; or
- complete jointly feasible states sampled/validated with CarbonScope's feasible-region tools.

## Speed and verification

Because the AE4 CBM is expected to be small, speed is secondary to correctness. Use `run_highs_fva_reference` at least once as an oracle and compare its endpoints against `run_highs_vffva`/`run_prepared_highs_vffva`. After equivalence is verified, use the production VFFVA path for repeated FVA calculations.

Record:

- CarbonScope commit SHA;
- the exact `FluxModel` reaction/metabolite ordering;
- the S-matrix;
- flux bounds and provenance;
- objective or feasibility formulation;
- solver tolerances;
- reference-FVA versus VFFVA endpoint comparison;
- complete-state validation where used.

## Paired WT / AE4-null formulation

The key Task 46 questions should be posed with paired WT and KO feasible regions, with shared capacity assumptions stated explicitly. AE4-null must impose zero AE4 flux. NKCC1 genotype-specific capacity must not be allowed to increase unless a declared evidence source permits it.

Use the CBM to determine, before any new kinetic model is built:

1. WT feasible transport states;
2. AE4-null feasible transport states;
3. maximum achievable KO secretion/apical chloride export under the independent NKCC constraint;
4. whether high secretion requires NKCC flux to exceed its WT-supported range;
5. FVA ranges for NKCC1, AE4, pump, NHE1/NBC/AE2 and apical chloride exit;
6. targeted flux coupling/conditional range analysis relevant to AE4 loss;
7. sparse minimal-relaxation results identifying the smallest biological constraints that must change to recover the phenotype.

## Anti-combinatorial rule

Do not turn CarbonScope into a high-throughput architecture enumerator. Use one curated CBM network and a small number of scientifically motivated bound scenarios. Do not perform a Cartesian product of transporter bounds or mechanisms.

## Publication checkpoints

Publish the CarbonScope-backed CBM work incrementally on the existing Task 46 branch. At minimum publish separately:

1. pinned CarbonScope integration + AE4 `FluxModel` + S-matrix + reference/VFFVA equivalence;
2. WT feasible-space analysis;
3. AE4-null paired analysis under shared capacities;
4. NKCC-constrained KO analysis;
5. targeted coupling/minimal-relaxation analysis;
6. CBM decision report identifying which kinetic block(s) merit reconstruction.

At each checkpoint update `CURRENT_STATUS.md`, create an immutable checkpoint file, commit, push, and verify the remote SHA before continuing.

## Decision gate

Do not begin the new kinetic/ODE model search until the CarbonScope-backed CBM decision report is published remotely.

The purpose of the CBM is not to replace the dynamic model. It is to reduce the kinetic search space by identifying which structural degrees of freedom are actually necessary to explain the AE4 phenotype.