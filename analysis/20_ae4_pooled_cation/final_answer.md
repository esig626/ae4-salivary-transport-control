# Task 20 result for review

**The prescribed pooled-cation model cannot support any of the ten inherited
WT states under P0 or P10. All 20 decisions are structurally infeasible.**
The transport law is implemented and verified, and the requested experiment
is completely accounted for. No admissible pooled WT model exists from which
to obtain a matched genotype prediction under the task's constraints.

Primary classification from the required list:
**`TASK 20 NUMERICALLY INCONCLUSIVE`**.

The list has no structural-infeasibility category. This classification is used
only to denote the unavailable phenotype comparison. The precise scientific
status is **`STRUCTURALLY_INFEASIBLE_POOLED_AFFINITY`**. The WT sign contradiction
is proved analytically and independently checked at 50-digit precision; it is
not a numerical optimizer or integration failure. “Extreme shifts” would be
misleading because no finite allowed capacity shift can satisfy the targets.

## What was established

| Required question | Result |
|---|---|
| Does the new law eliminate opposing Na/K net transport? | Yes, by construction. Na and K partition one scalar cycle on its donor side. All assay and inherited-reference evaluations satisfy no slip. There are zero eligible new whole-cell trajectories, so no claim of simulated phenotype validation is made. |
| Can all ten inherited WT states be supported under P0 and P10? | No: 0/10 under each condition. Every original state remains byte-for-byte unchanged in the decision payloads. |
| How much capacity rebalancing is needed? | No finite nonnegative capacity solution exists. Canonical log-fold objectives and accepted capacities are undefined, recorded as null. |
| Does AE4 5% lower or raise secretion in P0/P10? | Unevaluable in both conditions. No genotype simulation or refit was performed. |
| Does AE4 loss lower resting Cl? | Unevaluable; pH, Na, K and compensating-flux changes are also unavailable. |
| Does AE2 remain comparatively neutral? | Unevaluable for the pooled architecture. The reused legacy comparators show less than 0.161% change. |
| Is the result robust across roots/calcium? | The structural sign failure holds at all ten roots, including both provenance families. It prevents evaluation at calcium 0.10, 0.25 and 0.50 uM; no calcium-response conclusion is inferred. |
| Is the pooled mechanism sufficient? | It eliminates slip but is incompatible with the prescribed positive-loading/fixed-state experiment. Further mechanism or fixed-state premise work is required before its secretion effect can be assessed. This does not identify a unique true microscopic mechanism. |

## Decisive constraint

At every inherited root, the specified pooled affinity is negative:
`A=-2.638392096` to `-1.650095042`. Local detailed balance requires
`k_forward/k_reverse=exp(A)<1`, hence `J<=0` for every nonnegative active AE4
capacity. P0 and P10 both require `J>0`. The target flux would therefore make
`J*A<0`, violating the prescribed thermodynamic condition.

The production state contains only 4.4299–7.1889 mM intracellular bicarbonate,
whereas the zero-affinity diagnostic at its other fixed activities is
16.3164–16.6780 mM. Those diagnostic concentrations were not substituted into
the inherited states. Membrane voltage reclosure and changes to other
capacities cannot alter this electroneutral affinity at a fixed state.

## Reused comparator results

| Condition | AE4 5% / WT secretion | AE2 loss / WT secretion |
|---|---:|---:|
| Legacy, inherited productive loading | 1.08229–1.21663 | 1.000318–1.001601 |
| Pooled P0 | unavailable | unavailable |
| Legacy, 10% productive loading | 1.18617–1.48992 | 1.000260–1.000997 |
| Pooled P10 | unavailable | unavailable |

Each legacy range covers all ten roots and three calcium values. All 60 legacy
comparison rows were reused after hash verification; no paired pooled effect
or architecture-induced direction change can be calculated.

## Artifacts and verification

`results/20_ae4_pooled_cation/` contains the contract, 14-test assay record,
reference comparison, 20 fixed-state decision payloads, frozen manifest,
remote checkpoint receipt, three complete 60-row unavailable dynamic tables,
native-state diagnostic table, 120-row four-condition comparison, reuse
manifest, post-freeze summary, and the complete test transcript.

The WT checkpoint was committed, pushed and remotely verified at
`a2902286c148c886ee2d8188dabbd020d0a2d32f` before post-freeze comparison access.
The 187 inherited frozen inputs remain unchanged. New code consists of the
separate selectable mechanism and Task 20 analysis/reporting modules. Existing
scientific equation files, prior results and `archive/` remain unchanged.
The branch is submitted for review; no automatic merge is performed.

The focused and inherited suite passed **122 tests and 10 subtests**. The one
warning is an inherited helper named `test_parameters` returning an object;
it is unrelated to the pooled law. No tests failed or were skipped.

WT stage commands used in this execution, after adding the Task 20 code:

```bash
python -m src.modern_full_model.task20_fixed_wt assays
python -m src.modern_full_model.task20_fixed_wt prepare
python -m src.modern_full_model.task20_fixed_wt solve
PYTHONPATH=src python -m pytest -q tests/test_task20_pooled_ae4.py tests/test_task20_fixed_wt.py tests/test_modern_ae4.py tests/test_task18_fixed_wt.py
python -m src.modern_full_model.task20_fixed_wt freeze
```

The freeze also includes the checked-in evidence/method reports and test
transcript. These preparation stages refuse to overwrite an existing frozen
WT manifest. Checkpoint creation and remote verification preceded the
post-freeze report. On this completed checkout, verify the immutable freeze
with the test command below. The reporting command can reproduce the derived
tables (and refreshes their reporting timestamps):

```bash
python -m src.modern_full_model.task20_postfreeze
PYTHONPATH=src python -m pytest -q tests/test_task20_pooled_ae4.py tests/test_task20_fixed_wt.py tests/test_task20_postfreeze.py tests/test_modern_ae4.py tests/test_modern_full_model_whole_cell.py tests/test_modern_nkcc_stimulation.py tests/test_modern_camp_pka.py tests/test_modern_dimensions.py tests/test_task18_fixed_wt.py tests/test_task18_postfreeze.py tests/test_task19_pump_coupling.py
```

The reporting command reproduces tables from frozen inputs and never launches
new genotype calculations for this empty feasible set. The test transcript
records the environment's test-run details.
