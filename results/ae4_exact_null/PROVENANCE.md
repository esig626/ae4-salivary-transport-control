# Task 26 exact-null provenance

This directory records the exact-AE4-null diagnostic supplied after Task 25. The uploaded archive was labelled generically as `Task.zip`, but its contents identify it as Task 26, not the Task 25 carbon/NHE1 sensitivity study.

The full uploaded archive contained 11 useful UTF-8 files and no binary trajectories. The scientific report, summary, verification record, matched secretion comparison table, and flux decomposition are committed here. The larger raw diagnostic tables remain checksum-addressable from the supplied archive and can be restored without ambiguity if needed.

## Original uploaded-file SHA256 values

- `Task.md` -> `aa84102695ed98f0255964cb4ceff1d229449767abfa20642720c4985c510f4f`
- `run_ae4_exact_null.py` -> `fd4ec66fa9424e77f92580eec57e900c16e5be8fa9229495728ff86bbc3a2f56`
- `test_ae4_exact_null.py` -> `5c765b71dc2cac974f7635b40bf18a8452b5ab74b02b17231b9774dbe7f79b3c`
- `resting_states.csv` -> `1a3c9d7bfc912e564e6ba37cc603a1aa245a6ecfb1759f1b2503b30aab50a5e0`
- `comparison.csv` -> `a4bacd537ca1e228ed540109720bcde03ee68e8e80ffc5e6adef0f5f0e0e43b7`
- `flux_decomposition.csv` -> `75b57271640e9d17bdea5f8b443dc7f159e9c503102206ca76864d90e96553ed`
- `trajectory_summary.csv` -> `db2b5c93f165cd58194a9bcd993d73103a1fb2ad703c2c68440312a03550e893`
- `warm_start_attempts.csv` -> `06ab0e6f2a2b78ff78c93683dca2769a9bc79042cd6bfeff69ff0d8884814366`
- `verification.json` -> `a56fba4e5746adc2a9820ff2ac823e133c7b85663c58d3bf669de698b6320952`
- `contract.json` -> `da88c48ca0d2bae0d4b2a2b60f4964e7293665903dd41db1a07f212d158e822b`
- `summary.json` -> `33df414fb09fbeb704577f725821d37ad7c69a2312528ffa3fb583d3ae715d45`

## Scientific status

Only five of ten exact-null roots reached numerically closed equilibria, and all five violate the inherited resting pH/volume domain. Five roots remained unresolved. Therefore no admissible exact-null secretion phenotype was established.

For the five closed formal contrasts, exact AE4 deletion increased secretion rather than lowering it, with null/WT ratios 1.0358-1.0430 at Ca 0.10 µM, 1.0823-1.1031 at 0.25 µM, and 1.1024-1.1273 at 0.50 µM. The closed null states were strongly alkalinised and swollen, with pH about 8.00-8.01, intracellular HCO3 increased by about 71-73 mM, and volume 7.0-8.6 pL. NKCC1 replaced about 138-140% of the lost resting AE4 chloride loading.

These are formal model diagnostics, not admissible physiological predictions.

## Current lineage

The routed-AE4 true-rest baseline and Task 25 carbon/NHE1 sensitivity have been merged into `main`. Older pooled-AE4 Task 21/22 pull requests were closed as superseded and were not merged. This Task 26 result therefore sits on the active routed-AE4 lineage and is intended as the starting diagnostic for the subsequent inverse phenotype-repair optimisation.
