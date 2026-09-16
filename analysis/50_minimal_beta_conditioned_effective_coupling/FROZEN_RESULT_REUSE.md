# Task 50C frozen result reuse

**QUALIFIED FROZEN REUSE ACCEPTED. TARGET-CALIBRATED CONSTRUCTION.**

All 28 selected frozen artefacts are byte identical to Task 41
publication `86f135e075aa91a70a752b650d7bb894353fee04`. Their Git object identities and SHA-256 values
are in `output/frozen_reuse_manifest.json`. Initial states, state names,
protocols, onset conventions, solver settings and shared parameter hashes match.
The common rest hash is `a0b96ef16e5967de4dfe2bdabb882f2a28d478e211f920bc19cd735a5a143acb`. All saved production conservation gates
passed. WT's secretion summary equals Task 40's exactly.

| Frozen case | Cumulative secretion, 0 to 600 s, pL | Cumulative deficit | Flow deficit at 600 s |
| --- | ---: | ---: | ---: |
| WT | 0.992524544 | 0.000000% | 0.000000% |
| AE4 5% | 0.762622459 | 23.163365% | 20.108928% |
| AE4 null | 0.692175197 | 30.261151% | 27.038018% |

These are inherited Task 41 results, not new simulations. The full stored
coefficient differs slightly from the complement of the mandated Task 50
lambda. The user authorised numerical treatment of that difference, and 50B
passed the predeclared fixed state comparisons. Accordingly the mutant outputs
are reused with that explicit precision qualification, not claimed as bit
identical Task 50 solutions or a rigorous global trajectory error bound.

The saved cumulative observations at every 60 s from 60 through 600 s show
deficits between 23.1634% and 28.8157% at 5% AE4,
and between 30.2612% and 37.0477% at null AE4.
The endpoint flow deficits are also substantial. Thus the inherited effect is
present across the recorded window and is not supported by only one selected
instant. This is a finite window observation, not a stationary, infinite time
or parameter robustness result. No experimental timing was fitted or used as
a gate. The existing table, not a newly integrated trajectory, supplies these
observations.

WT, REST, AE4 null CCh only and AE2 knockout remain exact parent initial value
problems under the declared matching states and protocols. This does not prove
agreement with every corresponding experiment. No directly matching CCh only
or AE2 Task 40/41 production cohort was identified in this frozen set, so the
formal identity and 50B software checks suffice; no additional trajectory was
run and no different Task 48 bath or chronic resting state was substituted.

Chronic established genotype resting chloride/pH adaptation remains unresolved.
For IPR only the calcium coordinate is zero, so the added coupling is off and
cannot repair the known positive uptake discrepancy from exact null genotype
rest. The microscopic coupling and long duration physiological validity also
remain unresolved. The frozen null endpoint has retained chloride and increased
pH/volume; passing the original finite window gates is not preserved overall
homeostasis. The phenotype selected parameter is not independently identified
and supplies no evidence of a molecular AE4 to TMEM16A interaction.

New 50C model evaluations, production trajectories, stationary solves and fits:
**zero**. `reuse_frozen_results.py` imports only the standard library, reads
existing results and verifies their provenance.
