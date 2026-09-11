# Task 14B execution audit

The controlling scientific source is `a2fde97e30d79601e7455d9f1d77d8d67609cb49`. The working branch is `codex/task-14b-five-percent-ae4-check`, starting at `645ee8b25c82db885dedb1884d69fb8f3e778d75`. `AGENTS.md`, the Task 14 final report, the selected continuation states and the Task 14B prompt were read before execution. Task 14B explicitly supersedes the older Task 14 exact zero and detailed time signature requirements.

The supplied runner already selected the correct saved 5% AE4 states, built the frozen model and used production Radau at the three declared calcium values. Its model construction and simulation call were retained. Corrections were limited to the prompt's required integrity checks, complete diagnostic retention and reporting:

* The Task 14 final artifact hash ledger is checked against its Git blob at the controlling commit. All 158 files, including the ledger, remain unchanged. The 22 frozen equation source hashes, native model manifest and all root parameter payloads are verified.
* Each of the ten initial states is the existing selected, connected, rank 10 state at nominal 5% AE4, with passing numerical and conservation checks. The inherited floating point expression is `0.04999999999999993`; it identifies the nominal `0.05` step within the existing selection tolerance. The core state is copied exactly without a new resting solve. The simulated AE4 expression is exactly `0.05` and AE2 expression remains `1.0`.
* Each of the 30 WT denominator totals is checked against the matching saved production Radau trajectory, its full time grid, solver specification and parameter hashes. Direct trapezoidal integration confirms the saved total. An integrity failure stops execution; there is no WT recomputation path.
* A manifest of hashes, saved initial states, the full request matrix, fixed solver and protocol is written before simulations. New trajectories explicitly use the matching saved WT time grid, which is the original frozen 602 point grid. Every complete trajectory and its diagnostics is saved. Failed cases would remain in the pair table with unavailable primary ratios, and exceptions would be recorded. An existing production checkpoint prevents silent reruns or replacement.
* Summaries cover every root and calcium, both routing families, overall and per calcium medians and ranges, and the matching saved AE2 effect. Negative percentage reductions are retained with their original sign. Time traces and experimental values are not numerical gates. The final prose interpretation was added from saved results after execution.

Exactly 30 new 5% AE4 production trajectories were run with four processes and one thread per numerical library. There were no new resting solves, WT trajectories, AE2 trajectories or additional solver runs. No attempt was made to reach zero AE4 expression. All 30 cases completed and passed the inherited positivity and conservation checks. The largest dimensionless conservation residual ratio was `0.002277289468111121`, below the inherited limit of 1. Maximum cell and lumen charge residuals were `6.394884621840902e-13` and `5.945244296867713e-14` fmol. All input checks were repeated after execution and matched the saved manifest.

The effect is an increase for every root and calcium. This is reported without retuning or selecting cases. Neither numerical validity nor the large AE4 to AE2 effect contrast is presented as evidence of reproducing the direction of the observed AE4 secretion phenotype. The approximately 35% experimental reduction is context only.

Execution command, run once from the repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=src python -m modern_full_model.task14b_five_percent
```

Focused verification command:

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=src:tests python -m unittest -v test_task14b_five_percent
```

Seven runner checks passed before execution; four artifact checks were then enabled, and all eleven passed. Tests use saved trajectories or mocked returns and do not run additional ODEs or resting solves. They verify intact model construction, failed input rejection, the expression and solver passed to the worker, failure retention, exact saved initial states, all 30 pair identities, primary integral arithmetic and unpruned ensemble summaries.

Only Task 14B orchestration, tests, results and analysis are changed. The completed Task 14 source and artifacts, `archive/` and `main` remain unchanged.
