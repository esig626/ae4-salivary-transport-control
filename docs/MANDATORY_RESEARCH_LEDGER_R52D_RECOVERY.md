# Mandatory research ledger continuation — Task52 52D recovery rerun

**Binding status:** append-only operational/scientific provenance continuation for Task52. Read with the master ledger and Task52 decision records.

**Classification:** RECOVERY / REPRODUCIBILITY RERUN AUTHORIZATION. This does not introduce a new mechanism, parameter, case, fit or phenotype-informed change.

The original Task52 52D runs were reported by the prior Codex workspace as completed to 600 s, but the completed outputs were not published to GitHub before that workspace became unavailable. Subsequent workspaces cannot access those unpublished files. The remote branch contains only the recorder compatibility boundary and no completed `52D` checkpoint.

The user has explicitly authorised one clean rerun of the frozen 52D production stage so the three predeclared cases can be regenerated and remotely published.

The recovery rerun must use verified 52C `36f7c3d0bb2353bf55ac0ed9e22582720097c829` unchanged and exactly the three 52D cases defined in `prompts/52_chloride_reservoir_final_test.md`:

1. matched supply, central projection, combined CCh+IPR, `G_aux=0`;
2. matched supply, central projection, combined CCh+IPR, `G_aux=2.32e-9 S`;
3. matched supply, central projection, combined CCh+IPR, `G_aux=4.49e-9 S`.

The published recorder-only compatibility adapter may be used. No scientific source/input, projection, parameter, solver setting, tolerance, case definition, observation window, gate or interpretation rule may change.

Do not run 52E or 52F as part of this recovery. After the three cases are regenerated, publish/verify `52D` remotely and stop.

If any case fails a frozen scientific/physical gate, record and publish the failure without rescue or retuning.

Controlling detailed instruction: `analysis/52_chloride_reservoir_final_test/52D_RECOVERY_RERUN_AUTHORIZATION.md`.
