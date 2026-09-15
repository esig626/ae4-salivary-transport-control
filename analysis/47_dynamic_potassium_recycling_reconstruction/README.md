# Task 47 dynamic potassium recycling

Read `DYNAMIC_POTASSIUM_RECYCLING_REPORT.md` for the article narrative, `EQUATION_AUDIT.md` for the exact equations and `BASELINE_DIAGNOSTICS.md` for detailed flux accounting. The answer is negative for a binding ceiling in the frozen model. No independently specified physiological correction was identified; biological channel and pump capacities remain uncertain.

The immutable prediction commit is `56ab6a00cb4b75c70e5857d2c58e7d46731d4d5f`. It precedes the 5% run and phenotype comparison. `PREDICTION_MANIFEST.json` and `REMOTE_PREDICTION_RECEIPT.json` record the freeze. The final validation commit is `56305f56c374e74a922e240016b1d4c349f00c96`. The containing report commit is the final Task 47 checkpoint.

## Check saved results

From the repository root, with its NumPy, SciPy and Matplotlib dependencies:

```
python analysis/47_dynamic_potassium_recycling_reconstruction/verify_results.py
python analysis/47_dynamic_potassium_recycling_reconstruction/validate_predictions.py
python analysis/47_dynamic_potassium_recycling_reconstruction/build_figures.py
python analysis/47_dynamic_potassium_recycling_reconstruction/write_report.py
```

These commands do not run scientific trajectories. The independent verifier does not import the production model. CSV archives under `output/wt`, `output/ae4_null` and `output/ae4_5pct` are lossless gzip files. `output/compression_manifest.json` contains their compressed and uncompressed hashes. Each trajectory has the exact onset, stimulated right limit and 600 integer second samples; complete state vectors are separate. SVG and PNG figures are under `figures/`.

## Dynamic reproduction

`run_diagnostics.py` is the frozen production runner. It refuses to overwrite an existing case. `reproduce_in_scratch.py` can rerun a chosen case into a newly created temporary output directory and compare it with the committed prediction, without touching any archived output. It changes only where the analysis runner saves files. All production source, parameter and initial state hash checks remain active. This optional reproduction is not executed as part of the completed Task 47 work.

The original WT metadata recovery is documented separately. Do not run `recover_wt_summary.py` against a completed case; it is preserved as the exact record of that one recovery. Task 43 and Task 44 source pins are in `PROTOCOL.md`; their active parameter payload and relevant provenance are copied with attribution, without merging their branches or repeating their calculations.
