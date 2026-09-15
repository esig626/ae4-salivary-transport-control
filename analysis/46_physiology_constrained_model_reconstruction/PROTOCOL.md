# Task 46 execution protocol

## Authority and scope

Execute the user's checkpoint directive supplied in `Pasted markdown.md` and all
three mandatory specifications at main commit
`850706dd6b64d68793d7bca89adc65679c73d4b0`:

- `prompts/46_agentic_physiology_constrained_model_reconstruction_and_full_mathematical_reanalysis.md`
- `prompts/46a_nkcc1_model_reconstruction_addendum.md`
- `prompts/46b_independent_model_discovery_and_synthesis_addendum.md`

`AGENTS.md` was read first. Its manuscript assignment is superseded by the
explicit user instruction to execute Task 46. Its applicable scientific
preservation, British English and bibliography rules remain in force.

Work only on `analysis/task-46-physiology-constrained-model-reconstruction`.
All new tracked work belongs under this directory. Do not modify production
source, main, manuscript, Tasks 1–45, frozen results or previous analysis
branches. Never merge or force-push. The current main commit is the branch base;
this does not merge either companion analysis.

## Fixed inputs

- Task 43: `547f113d9123ab1976774639a69483faf40cef34`.
- Task 44: `e5fa9840bf96be3147c6118daee4942468c78f8b`.

Read these exact commits, in that order, then Tasks 37–42 and their frozen
outputs, the Task 43 and 44 reports and PDFs, and production/manuscript context.
Use detached scratch source checkouts or `git show`; do not copy their result
trees into Task 46. Keep original source and result bytes unchanged. Runs must
write only new Task 46 outputs or disposable scratch outputs. Any redirection
wrapper must be explicit and must not alter scientific equations.

`source_pins.json` records source identities, required prompt blob SHAs,
initial branch refs and the protected base tree. Verify rather than trust a
pre-existing local copy.

## Publication and recovery

One main scientific objective or closely related calculation family per
milestone. Split an oversized milestone and publish the completed part first.
At the end of every milestone, before the next begins:

1. Fetch the current remote Task 46 branch (or verify absence for first push).
2. Reconcile only by fast-forward. Divergent history is a blocker; preserve it.
3. Independently verify the completed work and protected files.
4. Update `CURRENT_STATUS.md`; create immutable `CHECKPOINT_XX.md`.
5. Commit, then push immediately, without force.
6. Independently confirm the remote branch contains the commit and that the
   working tree is clean. Do not begin the next milestone until both pass.

Each status/checkpoint gives number/name, pre-milestone SHA, work completed,
verification performed, products, conclusions, unresolved issues, exact next
milestone and resume commands, and publication cleanliness. Since a commit
cannot contain its own SHA, its containing Git commit identifies the checkpoint;
the next checkpoint records that SHA. Final publication reporting must list all
checkpoint SHAs. Numerical outputs use immutable `output/checkpoint_XX_*` names.
Earlier checkpoint files and outputs are never overwritten.

On recovery, fetch and read the latest remotely published checkpoint, confirm
its commit, and execute only its declared next work. Do not restart completed
scientific milestones. An interrupted local stage is uncompleted until verified
and published.

## Ordered milestones

- 00: branch bootstrap and recovery record; no numerical claims.
- 01: exact relevant Task 44 baseline replay, including WT/5%/null, stationary
  and 600 s deficits, all three compensation definitions, states and slow modes.
- 02: causal diagnosis of excessive compensation; no equation changes.
- 03: freeze construction/calibration/holdout/context evidence in Markdown and
  machine-readable form, including provenance and uncertainties.
- 04: forensic audit of every current NKCC1-law term; no replacement.
- 05: independent primary-literature search and model-discovery ledger; no
  detailed knockout simulation.
- 06: freeze genuinely distinct candidate equations, first-principles
  derivations, provenance, parameters and permitted calibration before testing.
- 07 onward: one NKCC1 candidate family per published checkpoint, WT first,
  preserving failures. Apply the scientific rejection hierarchy.
- Freeze the accepted NKCC1 law (or record no accepted model) before any
  downstream architecture changes.
- Whole-cell hypotheses, then one predeclared architecture per checkpoint.
  Remotely freeze new equations before fitting or phenotype comparison.
- Freeze accepted whole-cell model before calibration/holdout validation and
  mathematical analysis. A verified implementation bug is documented; a
  scientific revision requires a new versioned candidate stage.
- Separate checkpoints for calibration, holdouts, provenance/robustness,
  exact structure/nondimensionalisation, equilibria/local dynamics,
  sensitivities/compensation, source classes/longer time, old/new synthesis,
  clean replay, report source, and final compiled/report verification.

Number these later stages sequentially when the candidate panel is known.
Do not skip publication because a candidate fails. A scientific negative
result is allowed; missing execution or missing evidence is not model failure.

## Scientific safeguards

Retain Na/K/Cl, total inorganic carbon, alkalinity, dynamic finite intracellular
and lumen volumes/compositions, algebraic pH and buffering, electrical closure,
water, AE4 regulation, exact stoichiometry and applicable reversal. The old
architecture is a failure baseline, not a template constraint. NKCC1 is
untrusted. An arbitrary cap or phenotype-selected equation is insufficient.

Independently search primary sources beyond the repository. Separate salivary
evidence, cross-tissue evidence, model-derived quantities and new hypotheses.
Derive at least one first-principles finite-carrier reversible saturating NKCC1
candidate if literature does not supply an adequate salivary model. Evaluate
at least three structurally distinct architectures unless evidence establishes
one uniquely appropriate model; parameter variants do not count.

Freeze the evidence split before model search. In milestone 03 explicitly
resolve the specifications' screening/final-holdout distinction: evidence
revealed for selection is no longer unseen final validation. Never silently
turn holdouts into fitting targets. Keep the isolated Peña NKCC1 assay as an
independent discriminator wherever possible, with honest protocol matching.

Reject mathematical/physical impossibility, thermodynamic invalidity, invalid
WT physiology, qualitative evidence failures and extreme unsupported
compensation/couplings in that order. Then prefer fewer unsupported assumptions
and fitted parameters; only then compare phenotype closeness. No weighted
all-observation loss and no unsupported Task 41-type rescue.

Keep exact/proved, local numerical, finite-time numerical and unresolved
conclusions distinct. Do not invent uncertainty ranges, experimental facts,
global proofs or a successful model. Preserve an equivalence class if
indistinguishable. Publish negative results and discriminating experiments if
no defensible candidate survives.
