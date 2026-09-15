# Task 48A addendum: parallel agents and crash-safe publication

This addendum is mandatory and supplements `prompts/48_joint_experimental_constraint_reconstruction.md`.

## Parallel-agent execution

Use one orchestrator and up to five parallel read/analysis workers. Parallelism is for independent evidence verification, protocol reconstruction, observation mapping, provenance audit and inherited-architecture feasibility. It is not permission to create competing models.

Only the orchestrator may own the canonical Task 48 branch history, integrate code, run authoritative inference, accept a model correction, commit, push or open/merge a pull request.

Freeze one common starting SHA before spawning workers. Every worker must use that same immutable scientific starting point. Workers must not recursively spawn more workers.

Recommended division:

1. **Primary-data worker**: verify every observation in `constraints_seed.json` against the primary paper/supplement. Record exact mean, uncertainty type, sample size, units, timing, agonists/inhibitors and genotype. Distinguish quantitative observations from qualitative or non-significant comparisons.
2. **Protocol worker**: map WT, AE4-KO and AE2-KO rest, CCh, IPR, CCh+IPR, bicarbonate-free NKCC assay, NHE assay, isolated exchanger assays, ±IPR activation and whole-gland secretion to explicit simulator interventions. Flag non-representable protocols.
3. **Observation-operator worker**: map each experimental measurement to the corresponding model observable. Do not silently equate SPQ fluorescence slope, intracellular chloride slope and a single transporter flux.
4. **Parameter/provenance worker**: classify candidate parameters as measured/fixed, literature-derived, calibrated/effective, model-determined/placeholder or unidentifiable. Produce the exact admissible inference set.
5. **Independent feasibility worker**: using inherited equations plus correct protocols/resting states, examine whether the full joint fingerprint is structurally feasible and identify equation-level contradictions. No new mechanism families.

Workers may use isolated worktrees or scratch branches, but their work is advisory. Never merge a worker branch wholesale. The orchestrator must inspect and integrate evidence-backed changes selectively.

## Serial scientific spine

The authoritative sequence remains serial even though preparation is parallel:

`verified constraints -> executable protocols/observables -> inherited-architecture joint feasibility -> one directly data-required correction if needed -> final joint reconstruction -> article-facing decision report`

Do not allow parallel workers to bypass this order.

## Crash protection: push in pieces

Do not wait until the end to push. This task must survive a dead Codex session or lost workspace.

After every milestone below, the orchestrator must:

1. fetch the remote Task 48 branch and confirm the expected parent SHA;
2. verify no unrelated files changed;
3. update `analysis/48_joint_experimental_constraint_reconstruction/CURRENT_STATUS.md` with completed work, key numerical findings, unresolved blockers, exact next step and authoritative input/output SHAs;
4. write a machine-readable receipt under `analysis/48_joint_experimental_constraint_reconstruction/output/checkpoints/` with milestone name, parent SHA, produced files, tests and key result;
5. commit the complete milestone;
6. push immediately to `analysis/task-48-joint-experimental-constraint-reconstruction`;
7. fetch the remote branch and verify it points to exactly that commit;
8. require a clean canonical worktree before dependent work begins.

Never force-push. Never rewrite a published milestone. Corrections to an earlier interpretation must be new commits that explicitly supersede it.

## Mandatory published milestones

### 48A: evidence/protocol freeze

Publish the verified experimental constraint ledger, protocol map, observation-map plan and parameter-admissibility table. No fitting and no model correction yet.

### 48B: executable protocol and observation layer

Implement the minimum code required to simulate all representable experimental protocols and observation operators. Add regression tests showing that interventions isolate the intended pathways. Push before any joint inference.

### 48C: inherited-architecture joint feasibility

Give established AE4-KO and AE2-KO animals genotype-appropriate resting states rather than instantaneous deletion from the WT state. Run the inherited architecture with one shared admissible parameterisation against the full constraint set. Publish the full residual table and a clear feasible/infeasible decision even if the result is poor.

### 48D: directly data-required correction, only if necessary

If 48C fails because a specific experimentally established dependency is absent or wrong, implement only the smallest such correction. The first pre-authorised example is β-adrenergic/cAMP/PKA activation of AE4 if audit confirms it is absent. Its parameters must be constrained by isolated exchanger ±IPR data, not saliva output. Push the correction and a pre-secretion prediction before final phenotype comparison.

If 48C is jointly feasible, skip 48D.

### 48E: final joint reconstruction

Run the authoritative joint inference/validation with the accepted architecture and one shared parameterisation. Publish all observables, residuals, uncertainty handling, genotype/protocol simulations and identifiability diagnostics. Do not retune after examining the final whole-gland result.

### 48F: article-facing decision report

Publish `analysis/48_joint_experimental_constraint_reconstruction/JOINT_EXPERIMENTAL_CONSTRAINT_REPORT.md` with measured constraints, observation assumptions, inferred parameters, any equation change, held-out predictions and residual contradictions. Push and stop.

## Important scientific constraints

- The full phenotype vector is joint. Do not fit only the approximately 35% secretion deficit.
- Established AE4-KO and AE2-KO animal data require genotype-specific resting states.
- Use reported uncertainties correctly. SE is not automatically a hard admissible interval.
- A reported non-significant difference is not exact equality and does not justify an invented ±5% band.
- One shared parameterisation must apply across genotypes and protocols, apart from the actual gene deletion and directly established protocol interventions.
- Do not allow genotype-specific NKCC, NHE, NBC, K-channel or pump multipliers without independent evidence.
- No Cartesian grids, broad random sweeps, mechanism enumeration, arbitrary caps or secretion penalties.
- AE2-KO is a mandatory negative control.
- The early-near-WT then later-divergent secretion time course is a constraint, not merely the 10-minute endpoint.
- If the full fingerprint remains incompatible after the permitted data-required correction path, report architectural incompatibility rather than adding another degree of freedom.

## Stop condition

Stop after Task 48. Do not automatically create Task 49. The next scientific step must follow from the published residual structure.
