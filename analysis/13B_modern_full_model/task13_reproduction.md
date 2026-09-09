# Task 13 numerical reproduction (G0)

## Scope and firewall

Task 13B began from remote commit `78be17254a167bd054448d089aafa247b1539e3f` on branch `codex/task-13b-modern-full-model`. The Task 13 state-resolved source, saved ledgers, and tests were materialized byte-for-byte before any Task 13B model file was created. Nothing under `archive/` was edited or used as executable ground truth. The strict AE4-null saliva magnitude and time-course targets were not read by the numerical reproduction code.

## Reproduction command

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test_state_resolved_ae4_*.py' -v
```

Run date: 2026-08-27 UTC. Result: **72 tests passed** in 11.5 s.

## Decisive reproduced controls

- Frozen evidence and artifact hashes pass, including the original pre-reveal transporter ledger and CTMC snapshot.
- Registered SR families preserve carrier closure, local detailed balance, exact branch electroneutrality, reversal signs, and nonnegative entropy production.
- The smallest shared-pool SR2 family and the SR5 coordination alternative reproduce the saved conditional transporter-level behavior.
- All five saved Stage-A transporter/gauge survivors still have no accepted WT capacity and no valid WT denominator on the inherited seven-state chassis.
- The fixed-chassis AE4-null field is candidate independent; its bounded physiological-domain root and the failed WT residual are independently reproduced.
- The M1 apical/basolateral pump and K-channel split remains exactly nested at zero apical fractions, recloses both membrane currents and lumen charge, but does not repair the inherited-chassis WT gate.
- Conserved carbon/buffer reaction signatures close mass and charge, while the released knockout Cl/pH projection leaves the documented free directions.
- The saved KO-only IVP is reproduced by Radau and BDF but remains an unphysical-time diagnostic, not a Task 13B secretion prediction.

## G0 disposition

`G0_REPRODUCED_NEGATIVE_CONTROL` — retained only as a regression target and source of transporter-level modules. The inherited seven-state whole-cell chassis is not promoted. Task 13B therefore proceeds to a conservation-correct full acid-base/carbon reconstruction before any new whole-cell calibration.

