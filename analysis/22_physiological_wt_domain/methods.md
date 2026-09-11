# Task 22 methods and validation

## Checkpoint order and inheritance

The branch is `codex/task-22-physiological-wt-domain`, continuing from Task 21
commit `db734f52b8d71d670fcf793beab922967212bdd6` with the user's Task 22
instructions at `ecafebd3e3fc4237aeac3b063eb03deb67165955`. The physiology
contract and its implementation were committed and pushed at
`a1ec05b706dfaf4c1574abe3b33ea39f0195d472`. Authenticated remote reads
verified that ref, commit and contract blob before the first optimization.
The receipt is `stage0_checkpoint.json`.

All inherited production equations, parameters outside the fifteen declared
capacity coordinates, and prior Task 21 evidence are unchanged. The contract
hashes the inherited source files and all seven witness payloads. The archive
tree is `ba31f58b023e28030515d033893fcb52348abb61`; no archive file was edited.
The inherited local Git store lacks 113 unrelated historical blobs, including
two archive files. Their committed references remain unchanged; Task 22 does
not claim a fresh byte-level audit of unavailable historical files. All
required WT root inputs and Task 21 witnesses are present and hash-verified.

## Bounded search

The complete recorded plan contains 27 independent starts: all ten inherited
WT roots, all seven distinct saved Task 21 resting witnesses (including its
replay-failing witness), and ten additional WT-only high-carbon starts.
Each witness stays with its inherited immutable parameter background. Its
complete cell/lumen state, capacities and membrane potentials supply the
initial coordinates; pathological values are projected into the declared
box before restoration. The original and projected coordinates are both
saved. There is no genotype input and no penalty for distance from a seed.

The numerical coordinates reuse Task 21's two charge-neutral compartments:
log Cl, pH, log TIC, log volume and log Na/K ratio in each compartment,
fifteen natural-log capacity folds, and two open membrane potentials.
All capacity folds remain inside [0.01,100] throughout the optimization.
The direct state bounds retain the contract. Na/K-ratio bounds and the
intracellular TIC upper bound follow from its Na/K and osmotic bounds.
The TIC positivity exploration floor is 1e-12 mM; it is not an added
physiological acceptance condition. Membrane potentials are unbounded
numerical coordinates. A 1e-12 inward log offset at capacity endpoints
prevents floating-point exponentiation from exceeding the hard bounds.

Bounded trust-region least squares restores the production balance/current
transcription and scaled negative physiological/thermodynamic/anti-cancellation
constraint margins. The relative weight 100 is a numerical feasibility
weight. Dependent constraints may remain violated at failed endpoints;
they are hard rejection gates, not accepted soft penalties. Each start has
at most 1200 function evaluations and `ftol=xtol=gtol=1e-11`, with Jacobian
column scaling. Four worker processes evaluate independent WT starts.
There is no bound expansion, rate rescaling of the scientific model or
acceptance-tolerance change.

All 27 searches completed and their failed endpoints were preserved. Their
smallest maximum scaled resting RHS is approximately 0.00188291 against
the inherited 1e-7 gate. Every endpoint fails the resting balance and
positive-AE4-affinity gates; four also fail intracellular Na. None violates
the hard capacity box. These outcomes alone would establish only numerical
failure, not structural impossibility.

## Production validation and analytic exclusion

Every endpoint is evaluated through the unchanged production model with
Python-float capacities and full conserved state, and saved JSON payloads
are replayed through production in focused tests. Production values control
Cl/pH/Na/K/volume, luminal physiology, osmolarity, complete RHS, current,
charge, conservation, no-slip and anti-cancellation decisions. The solute
ledger records individual cycle and species sources before cancellation;
water is reported separately with its own units.

An independent analytic necessary-condition argument establishes a
domain-wide contradiction. The luminal pH/TIC/buffer domain forces
`TA_l - TIC_l < -0.00645 mM`. The production lumen and cell carbon-minus-
alkalinity identities then force positive bath-to-cell CO2, even after all
four permitted residual errors. Intracellular CO2 must therefore be below
bath CO2. The declared cell pH/Cl/cation ranges consequently bound pooled
AE4 affinity above by a strictly negative number, less than -0.0764.
See `structural_contradiction.md` for the derivation and
`structural_certificate.json` for all ten per-root certificates.

The proof uses monotonic analytic inequalities with 60-digit evaluation and
conservative rounded margins, not optimizer flags or a numerical grid as an
infeasibility certificate. It applies to this fixed chassis and predeclared
domain. It does not disprove the pooled AE4 mechanism in general.

The lexicographic WT objective is defined in the contract and reported for
every failed endpoint for auditability. With no admissible resting state,
the admissible ranking is empty. There are no per-solution WT trajectories,
no accepted WT payloads, and no Task 22 600 s integrations. The explicit
WT dynamic status table records that all calcium conditions are unavailable
because no admissible resting WT exists. The complete empty ensemble is
frozen and remotely verified before final reporting.

## Test scope

Focused tests cover each boundary and nonfinite value, all fifteen capacity
coordinates, zero-reference and fraction handling, production osmolarity,
signed chloride-pool and stoichiometric flux definitions, proof identities
and margins, all 27 saved production replays, input preservation, and
checkpoint ordering/tamper/empty-ensemble protections.

Broad validation uses an explicit WT-only allowlist: dimensions, dynamic
WT validation, NKCC stimulation, pooled-mechanism assays, cAMP/PKA family
API/equilibrium/positivity classes, and whole-cell conserved-state,
acid-base, homeostasis/topology and parameter-domain classes. The broad
suite passes 76 tests and 10 subtests. Its short synthetic WT integrations
test existing numerical infrastructure; they are not Task 22 admissible-WT
trajectories. The full historical suite is deliberately not invoked because
it contains genotype regressions, including a knockout evaluation inside
the whole-cell integrated-balance class.

`validation_execution.json` records exact test arguments and exit codes;
the focused and broad transcripts are stored beside it. Historical genotype
regressions and Task 22 genotype evaluations both remain at zero. Python
uses NumPy 2.3.5 and SciPy 1.17.0, matching the pinned analysis requirements.

## Reproduction

Read and verify the immutable contract checkpoint first. The WT search
entry point is `python -m modern_full_model.task22_search --workers 4` with
`PYTHONPATH=src`; it refuses to run after the WT manifest exists. Existing
attempt files cannot be overwritten by the runner. The certificate entry
point is `python -m modern_full_model.task22_certificate`, and the evidence
freezer is `python -m modern_full_model.task22_checkpoint`. Reproduction of
new searches belongs in a separate unfrozen checkout with the same verified
Stage-0 contract, preserving the published evidence. A nonempty or
certificate-inconsistent result cannot be frozen by the empty-ensemble
writer and must be investigated.
