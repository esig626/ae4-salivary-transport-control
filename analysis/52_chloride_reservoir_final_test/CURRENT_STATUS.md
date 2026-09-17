# Task52 current status

**State: 52D COMPLETE — all three central cases reached 600 s; STOP after publication verification.**

The authorised reporting-only recovery started at
`f2b56f903d33a0ac5f5c292db92cead8b6157f1f`, using the unchanged scientific 52C
freeze at `36f7c3d0bb2353bf55ac0ed9e22582720097c829` on
`analysis/task-52D-recovery-from-52C`. All 50 frozen scientific/input hashes pass.
The original failure at `9c0fde42b5eb16746413801fae555c84cb32e201` is preserved.

The three central matched-supply CCh+IPR cases completed with Radau and passed
all frozen sampled physical, conservation and reservoir gates. Their cumulative
KO secretion deficits are 20.0394286%, 20.0757527% and 20.1015822% for
G_aux = 0, 2.32e-9 and 4.49e-9 S, respectively. No BDF fallback was needed.
The maximum independent reservoir error is 7.49e-8 fmol (limit 1e-5 fmol).

Complete compact outputs: `output/reporter_recovery_52D/`.
Stage summaries: `PRODUCTION_SUMMARY_52D.md`, `output/summary_52D.json` and CSV.
Scientific checkpoint receipt: `output/checkpoints/52D.json`.
Remote verification is recorded separately in `output/publication_52D.json`.

The new recorder preserves the union of original diagnostic names and explicit
presence masks; no scientific law, input, diagnostic definition, numerical
setting or acceptance gate changed. Three fresh paired production trajectories;
no tuning, fits, extra cases or mechanism search. The separate setup-path failure
occurred before model construction and is preserved in the audit trail.

52E and 52F remain unexecuted and unauthorised by this continuation. No further
scientific execution follows this 52D publication. The task has no final 52F
mechanistic classification at this boundary.
