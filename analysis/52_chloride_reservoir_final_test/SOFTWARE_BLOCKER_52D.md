# Task52D fresh execution — frozen output-writer blocker

**52D INCOMPLETE. SOFTWARE EXECUTION FAILURE; no verified numerical result.**

The recovery execution started on `analysis/task-52D-recovery-from-52C` at
exactly `36f7c3d0bb2353bf55ac0ed9e22582720097c829`. Independent remote branch
and tree verification, all 52 checkpoint artifact hashes, all 50 immutable
source/input hashes, the parameter-freeze SHA and ledger prefix passed before
numerical work. Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0 match the freeze.
The only launch adaptation replaced the old branch-name check with the user's
explicitly authorised branch; frozen files and scientific functions were not
changed. See D52-14 and `output/publication_52C.json`.

## Actual execution

The unmodified scientific runner was invoked once for stage 52D, with one
scientific process and one thread in each BLAS pool. It started `case_01`
(central projection, matched supply, combined CCh+IPR, G_aux=0 S).
It terminated with `AssertionError` in `run_frozen_cases.py`, at:

```python
assert all(tuple(r)==columns for s in samples.values() for r in s['rows'])
```

This assertion is outside the integration exception handler and precedes
trajectory and summary persistence. The process exit code is 1. No trajectory,
summary, physical-gate outcome, reservoir budget or terminal time was preserved.
Consequently no claim of 600-second completion can be made, and missing results
are not zero predictions. The scientific evaluation count cannot be recovered
from the persisted evidence and is not invented.

| Case | G_aux (S) | Verified outcome |
| --- | ---: | --- |
| case_01 | 0 | Started; output-schema assertion; no saved numerical result |
| case_02 | 2.32e-9 | Not launched because the process terminated in case_01 |
| case_03 | 4.49e-9 | Not launched because the process terminated in case_01 |

`output/execution_52D.json` remains unchanged with `complete=false` and an empty
completed-case list. `case_01/started.json` and the original execution traceback
are preserved. No `output/checkpoints/52D.json` is created: the requested
three-case dependency boundary has not been reached.

## Source-level cause and bounded correction proposal

Read-only inspection of the exact frozen files identifies a deterministic
schema mismatch. At time zero, the inactive NBC branch of
`frozen_task37/src/modern_full_model/nbc_minimal.py` retains the baseline current
dictionary. At positive time the active closure also reports `nbc_basolateral`
and `basolateral_old`. `diagnostics.row_for` emits each dictionary member as a
column, while the runner fixes its columns from the time-zero row and demands
the same tuple for every subsequent row. Thus active rows contain
`current_nbc_basolateral_A` and `current_basolateral_old_A` absent at onset.
This is a schema-set mismatch, so sorting column names or suppressing the
assertion would not safely preserve the output.

A bounded remedy would change only output serialisation to preserve all
original named diagnostics, with an explicit presence mask for fields absent
at a sample, and retain the exact model, diagnostics, gates, solver and
quadrature. It must not invent a measured diagnostic value for an absent field.
No such correction has been applied or evaluated. The first process did not
persist its numerical arrays, so producing complete case_01 outputs would
require a newly authorised fresh attempt, not postprocessing a saved trajectory.

## Stop and preservation

The immutable 52C failure policy permits one BDF fallback only for an explicit
Radau solver-status failure, not this reporting exception. The frozen runner
also rejects an attempted stage/case rerun. Combined with the user's exact-freeze
instruction, this blocks an unapproved correction or silent repeat.

All 50 frozen scientific/source/input hashes still match after the failure.
No 52A, 52B or 52C scientific work was repeated. No earlier lost 52D work was
recovered or reused. No scientific parameter, model, diagnostic definition,
solver tolerance or quadrature was changed. No additional case, BDF fallback,
fit, mechanism search, 52E or 52F was run. This failure is operational evidence,
not a biological/mechanistic classification or a completed 52D result.
