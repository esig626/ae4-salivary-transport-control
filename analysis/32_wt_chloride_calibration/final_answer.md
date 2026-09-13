# Task 32: bounded WT chloride calibration

**WT-only calibration did not recover R09 pH and chloride simultaneously.** The initial AE4 scalar stage did not obtain a valid chloride sign bracket. It stopped at an unclosed positive-activity endpoint; this is a bounded numerical result, not a proof of model infeasibility.

Prepared head: `6027f68dfcf8113c13d726b91d0b38cd0349ef47`. Branch: `codex/task-32-wt-chloride-calibration`. One scientific agent; no other transporter law changed.

| R09 WT case | AE4 activity | pH_i | Cl_i (mM) | Na_i (mM) | Volume (pL) | Stationary closure |
|---|---:|---:|---:|---:|---:|---|
| Task31_cached_R09_WT | 1 | 6.910000 | 60.304967 | 11.636126 | 1.448426 | Pass |
| AE4_initial_1 | 0.1 | 6.969014 | 60.211276 | 11.527528 | 1.480329 | Pass |
| AE4_initial_2 | 0.01 | 7.081989 | 59.927872 | 11.378376 | 1.559098 | FAIL: candidate only |

WT acceptance intervals were pH 6.84-6.98 and chloride 48.6-51.6 mM. The 0.1 activity state passes full closure, charge, current, conservation, positivity, rank, capacity, ordinary sodium and volume gates, but exceeds the upper WT chloride limit by 8.6113 mM. Its pH lies inside the WT interval. Lowering capacity tenfold changed chloride by only -0.09369 mM from the cached WT state.

The 0.01 activity trial reached the unchanged Task 31 local ceiling of 2,400 actual residual evaluations. Its maximum scaled stationary residual is 0.00750789, versus the required 1e-7; the maximum amount residual is 0.00107400 fmol/s. Charge/current accounting and normal sodium/volume do not rescue missing stationary closure. The displayed pH and chloride for this trial are candidate diagnostics, not a phenotype or a scalar-bracket endpoint.

1. **Simultaneous WT recovery:** No. The first AE4 stage stopped; NHE1 correction and the final AE4 correction were not reached.
2. **Required scalars:** No calibrated pair was established. NHE1 remained at the Task 31 amount `2.339370005697548e-05` fmol (about 14,088 carriers). AE4 activity trials were 0.1 and 0.01. The scalar multiplies the existing AE4 carrier amount (0.1501264265 fmol at activity 1); WT genotype expression stayed exactly 1.0.
3. **AE4 chloride-loading fraction:** No post-calibration fraction exists. At the last closed, uncalibrated WT trial, AE4 supplies 1.334204% of positive basolateral chloride loading. This is an output, never a calibration target.
4. **Exact AE4-null closure with ordinary sodium and volume:** Not evaluated. There is no accepted WT calibration to freeze.
5. **Untuned null pH/chloride versus experiment:** Not evaluated. No null measurement or historical null result entered calibration.
6. **Secretion direction and magnitude:** Unproven; no stimulated integrations were permitted.
7. **R10 out-of-sample survival:** Not tested because no accepted frozen R09 WT/null pair exists.
8. **What remains unproven:** Whether any admissible simultaneous WT solution exists under the permitted two-scalar procedure, whether exact AE4-null REST closes, and whether AE4 loss predicts the experimental chloride, pH or secretion effects. This failed local endpoint does not establish absence of another stationary root or of a finite positive chloride bracket.

The last closed WT trial has signed AE4, NKCC1 and AE2 chloride fluxes 0.00348353194137, 0.257610844387, and -0.00138110401246 fmol/s, respectively. Total positive loading is 0.261094376328 fmol/s; net loading is 0.259713272316 fmol/s. AE4 / positive NKCC1 = 0.0135224585; AE4 / net basolateral = 0.0134129916. AE4 / positive AE2 is undefined because AE2 is a chloride counterflux, with zero positive loading. Signed AE2, cation, carbon, apical, paracellular and water fluxes are retained in `flux_ledger.csv`. The cached Task 31 row leaves paracellular chloride blank because that flux was not stored; it was not recomputed.

Budget: 2/18 stationary calls, 3,283/30,000 residual evaluations, 2/6 AE4 scalar evaluations, 0/4 NHE1 scalar evaluations and 0/2 integrations. Numerical execution took 0.969779 s, excluding development and algebra/control-flow tests. Global budgets remain unused in part because the predeclared unclosed-endpoint stop was reached. The factor-ten steps were directional scalar bracket discovery, with no factorial, grid or multidimensional parameter search. No sign interpolation was possible. Successful cached cases were never rerun.

`frozen_parameters.json` explicitly records **not frozen: calibration failed**, with no accepted parameters. No freeze hash is claimed. `cases.json` preserves both trials; `budget.json` preserves their counters. All 11 focused tests pass. All prepared model source files remain byte-for-byte unchanged.

Reproduce the focused checks with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m unittest discover -s tests -p test_task32_calibration.py -v`. Regenerate this report from saved checkpoints with `python scripts/report_task32_wt_chloride_calibration.py`; it never calls a stationary solver. The numerical calibration entry point is `PYTHONPATH=src python -m modern_full_model.task32_wt_chloride_calibration` with the same thread variables. An existing completed execution is returned read-only.
