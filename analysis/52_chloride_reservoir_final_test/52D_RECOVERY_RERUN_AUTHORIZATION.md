# Task52 52D recovery rerun authorisation

**Status:** binding recovery instruction for Task52 only.

The original Task52 52D scientific executions were reported by the prior Codex workspace as completed to 600 s, but their completed artefacts were not published to GitHub before that workspace became unavailable. The remote branch therefore contains only the recorder-compatibility boundary at `02823457251284d147fb143f6c8fee63b57e97a0`; it does not contain a complete `52D` checkpoint or the completed case payloads.

The user has explicitly authorised a recovery rerun of **52D only** because those unpublished artefacts are no longer recoverable from subsequent workspaces.

## What is authorised

Run exactly the three 52D cases frozen at verified 52C `36f7c3d0bb2353bf55ac0ed9e22582720097c829`:

1. matched supply, central projection, combined CCh+IPR, `G_aux=0`;
2. matched supply, central projection, combined CCh+IPR, `G_aux=2.32e-9 S`;
3. matched supply, central projection, combined CCh+IPR, `G_aux=4.49e-9 S`.

Use the exact frozen 52C:

- scientific source/input hashes;
- projected onset states;
- paired WT/KO model;
- matched non-AE4 supply rule;
- solver and thread settings;
- quadrature/observation windows;
- physical and conservation gates;
- diagnostic requirements;
- exact reservoir audit;
- recorder-only compatibility adapter already published at the 52D dependency boundary.

This is an **identical-input recovery/reproducibility rerun**, not a new mechanism, parameter choice, sensitivity or phenotype-informed scientific case.

## Absolute restrictions

- Do not change any 52C scientific input, source law, parameter, projection, solver tolerance, acceptance gate or case definition.
- Do not use Palk/Benjamin NKCC.
- Do not use Task41/Task50 effective coupling.
- Do not tune or inspect outputs to alter subsequent 52D cases.
- Do not add cases.
- Do not run 52E.
- Do not run 52F.
- Do not attempt to reconstruct or claim identity with the lost unpublished local files; treat them as unavailable historical workspace artefacts.

## Publication rule

Run and publish 52D only. After the three frozen cases complete:

1. write all required case outputs and compact trajectories;
2. verify the exact reservoir mass-balance audit and inherited physical/conservation gates;
3. create `output/checkpoints/52D.json` and the required publication receipt/status metadata;
4. commit and push the completed 52D payload to `analysis/task-52-chloride-reservoir-final-test`;
5. independently verify the remote SHA, tree, checkpoint/receipt and output hashes;
6. update `CURRENT_STATUS.md` to state 52D is remotely verified and 52E remains unrun;
7. STOP.

If any 52D case fails a predeclared scientific/physical gate, publish that failure exactly as observed and stop according to the frozen Task52 rules. Do not rescue it.
