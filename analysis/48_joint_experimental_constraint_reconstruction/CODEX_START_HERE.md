# Task 48 Codex handoff

Branch:

`analysis/task-48-joint-experimental-constraint-reconstruction`

Frozen scientific parent:

`0a6adda00f0bb57d527f341412ceffa348b67a76`

Read and execute, in order:

1. `AGENTS.md`
2. `prompts/48_joint_experimental_constraint_reconstruction.md`
3. `prompts/48a_parallel_agents_and_crash_safe_publication.md`
4. `analysis/48_joint_experimental_constraint_reconstruction/constraints_seed.json`
5. `analysis/48_joint_experimental_constraint_reconstruction/CURRENT_STATUS.md`

The task is intentionally hard-gated. Do not broaden it into another mechanism search.

Use one orchestrator and up to five parallel workers exactly as prescribed in the addendum. Workers may verify evidence, reconstruct protocols, build observation operators, audit parameter provenance and independently assess inherited-architecture feasibility in parallel. Only the orchestrator owns the canonical branch, integrates changes, runs authoritative inference, accepts scientific corrections, commits and pushes.

Order of attack:

1. Verify the complete experimental constraint ledger and construct the assay/protocol map.
2. Build genotype-specific resting states and explicit observation operators.
3. Quantify the residual vector of the inherited architecture under the correct protocols.
4. If the inherited AE4 law lacks the experimentally demonstrated beta-adrenergic/cAMP/PKA activation, add only the smallest AE4-specific activation mechanism and identify it from isolated exchanger data.
5. Permit at most one further equation-level repair, and only if a specific remaining residual and independent evidence uniquely require it.
6. Stop with a clean negative result if the joint data cannot be represented without unsupported mechanisms.

Use one shared parameter vector across genotypes and protocols. Do not fit the 35% saliva endpoint in isolation. Do not invent genotype multipliers, transporter caps, percentage tolerances or secretion penalties. Do not rerun Tasks 46 or 47 merely to reproduce known results.

## Mandatory crash protection

Publish in pieces. Do not wait until the end.

Every milestone 48A through 48F must be committed and pushed immediately. Each push must include an updated `CURRENT_STATUS.md` and a machine-readable receipt under `output/checkpoints/`. After pushing, verify the remote branch points to that exact SHA before dependent work continues.

Never force-push and never rewrite a published milestone. If the Codex workspace dies, resume from the latest verified remote checkpoint rather than repeating completed scientific work.
