# Task 50 final reduced model report

**COMPLETE WITH AUTHORISED PRECISION QUALIFICATION.**
**Classification: TARGET-CALIBRATED CONSTRUCTION.**

The fixed beta conditioned scalar is implemented. WT, REST, AE4 null CCh only
and AE2 knockout nest their corresponding parent model exactly. The inherited
combined stimulus secretion deficits are **23.1634% at 5% AE4** and
**30.2612% at null AE4**, through qualified reuse of frozen Task 41 results.
Ten software tests pass. No new production trajectory, stationary solve,
parameter fit or mechanism search was performed. Task 50 ends at 50D.

## Fixed law and implementation

`g_Cl,eff = g_parent * [1 - lambda * a_Ca * beta * (1 - e_AE4)]`

`lambda = 0.89488127156712`

The sole new coefficient is fixed at the mandated literal. It derives from
Task 41's phenotype calibration and is not independently measured. The existing
normalised calcium coordinate is `a_Ca`; `beta` is the existing protocol input;
and `e_AE4` is the existing expression fraction. No new signalling state,
protocol timing or adjustable parameter is introduced.

The public implementation is
[`Task50EffectiveCouplingModel`](../../src/modern_full_model/task50_effective_coupling.py).
It accepts the unchanged Task 40 `MinimalNbcModel` template and inherits the
Task 41 conductance reconstruction, RHS and integration dispatch. It replaces
only the existing apical chloride conductance before the full electrical
closure. The parent calcium Hill gate remains intact. Secretion and water
outputs are not multiplied directly. All 13 inherited state coordinates remain.
No earlier scientific source file or frozen result was edited.

## Exact identities and the frozen coefficient qualification

The Task 41 expression expands as

`(1-a) + a[b + (1-b)e] = 1-a+ab+ae-abe = 1-(1-b)a(1-e)`.

| Declared condition | Factor | Result at matching state and protocol |
| --- | --- | --- |
| WT, `e=1` | Exactly `1` | Exact parent WT RHS and initial value problem |
| REST, `a=0` | Exactly `1` | Exact parent RHS for the same genotype |
| AE4 null CCh only, `beta=0` | Exactly `1` | Exact parent CCh only RHS and initial value problem |
| AE2 knockout, AE4 retained at `e=1` | Exactly `1` | Exact parent AE2 knockout RHS and initial value problem |
| Combined stimulus, `beta=1` | `1-lambda*a*(1-e)` | Exact Task 41 algebra at displayed `b=0.10511872843288`; qualified numerical comparison to the actual frozen coefficient |

Exact vector field nesting follows because a unit factor returns the original
parent object. With the same initial state and protocol, the initial value
problem is unchanged. REST nesting does not make the acute WT initial state a
null equilibrium. AE2 nesting means identity to the parent AE2 knockout, not
identity of that knockout to WT. Neither identity proves experimental agreement.

The original 50A audit found that Task 41 actually stores
`b=0.10511872843288446`, whereas the mandated lambda has exact decimal
complement `0.10511872843288`. For those decimal literals, at beta one,

`m50 - m41_frozen = -4.46e-15*a*(1-e)`.

The exact inheritance premise therefore failed and 50A stopped before source
changes. The user subsequently authorised retaining the exact lambda and
treating this rounding difference through numerical equivalence. That
authorisation lifts the operational stop only; it does not make the two
literals equal. The audit and all frozen files are preserved. See
[`EQUIVALENCE_PROOF.md`](EQUIVALENCE_PROOF.md) and decisions D50-04 to D50-07.

## Software verification

The ten tests include 21 bit identical parent comparisons and six full active
comparisons against Task 41 at the common initial state and saved final states.
The largest factor difference is `4.5102810375396984e-15`; the largest RHS
component difference is `8.090750291955828e-15` in native units. The fixed
state current, charge and chloride source comparisons and the inherited
conservation gates pass. Unrelated fixed state water, AE4, homeostasis,
speciation, CO2, outflow and state charge evaluations are unchanged.

The numerical allowances were published before implementation: factor absolute
tolerance `8e-15`; active full model relative tolerance `2e-12`, with absolute
allowances `2e-13` for RHS and amount quantities, `2e-14 V`, `2e-25 A` and
`1e-24 S`. Full model comparisons use the combined relative and absolute
criterion. Exact parent cases have no numerical allowance. No tolerance was
adjusted after seeing results. Integration dispatch was verified with a mock.

Verification used **80 core evaluations**, including constructor checks, and
**zero integrations**. The first setup attempt stopped with zero evaluations
because the historical Task 40 manifest expected an obsolete AGENTS hash.
All 149 other inputs matched; the corrected task local fixture pins the current
binding instructions and verifies unchanged scientific inputs. A later file
reporter path error also performed no model evaluation. Both issues and their
resolutions are retained. See
[`IMPLEMENTATION_VERIFICATION.md`](IMPLEMENTATION_VERIFICATION.md),
[`verification_attempt_02.json`](output/verification_attempt_02.json) and the
two original verification logs in `output/`.

These comparisons establish the declared software equivalence at the tested
states. They do not prove exact mutant solution identity or a rigorous global
trajectory error bound for the unequal coefficients.

## Frozen phenotype and persistence

All **28** selected artefacts match original Task 41 publication
`86f135e075aa91a70a752b650d7bb894353fee04` byte for byte. Initial states,
state names, protocols, onset conventions, solver settings and shared parameter
records match. The common rest hash is
`a0b96ef16e5967de4dfe2bdabb882f2a28d478e211f920bc19cd735a5a143acb`.
The saved production conservation gates passed; WT secretion equals Task 40's
summary exactly. Provenance is recorded in
[`frozen_reuse_manifest.json`](output/frozen_reuse_manifest.json).

| Frozen case | Cumulative secretion, 0 to 600 s, pL | Cumulative deficit | Flow deficit at 600 s |
| --- | ---: | ---: | ---: |
| WT | 0.992524544 | 0.000000% | 0.000000% |
| AE4 5% | 0.762622459 | 23.163365% | 20.108928% |
| AE4 null | 0.692175197 | 30.261151% | 27.038018% |

These are inherited Task 41 results under the user's numerical precision
authorisation, not newly simulated Task 50 outputs. Across every saved 60 s
cumulative observation from 60 through 600 s, the deficits range from
23.1634% to 28.8157% at 5% AE4 and 30.2612% to 37.0477% at null AE4.
The effect is substantial throughout those observations and at the flow
endpoint. This establishes finite window persistence, not stationary,
infinite time or parameter robustness. No experimental timing was fitted.
See [`FROZEN_RESULT_REUSE.md`](FROZEN_RESULT_REUSE.md) and
[`inherited_phenotype.json`](output/inherited_phenotype.json).

No directly matching CCh only or AE2 production cohort was identified in this
Task 40/41 frozen set. Their exact parent initial value problem identities and
software checks suffice without another trajectory. No Task 48 result with a
different bath or genotype resting state was substituted.

## What is established and what remains unresolved

Task 50 provides one beta conditioned algebraic scalar for the existing Task 41
constructive direction, with exact parent nesting and the explicitly qualified
inheritance of its combined stimulus phenotype. It requires no new signalling
dynamics or mechanism search. It does not independently estimate lambda or
identify a molecular AE4 to TMEM16A interaction or microscopic cAMP/PKA pathway.

Chronic established knockout resting chloride/pH adaptation remains unresolved.
In the IPR only arm `a_Ca=0`, so the added factor is off and cannot repair the
positive null uptake discrepancy from exact genotype rest. The saved mutant
states also retain the prior pH/volume limitations; finite window conservation
does not establish overall homeostasis or long duration physiological validity.
No excluded mechanism was reopened to address these limitations.

## Permitted article interpretation

Conservation driven compensation largely masks AE4 deletion in the reconstructed
transport model. A single phenotype calibrated, beta conditioned effective
coupling between AE4 availability and stimulated apical chloride secretory
capacity retains the substantial Task 41 secretion deficit through the qualified
frozen result reuse described above, while leaving WT, REST, CCh only and
AE2 control behaviour nested in their parent models. This is a
**TARGET-CALIBRATED CONSTRUCTION**. The coupling is an effective network
parameter, not an identified molecular AE4 to TMEM16A interaction or independent
validation of that interaction.

## Publication, reproducibility and stop

Work began at `8f5fefd9fc0b563a6bd963895211b4acf76972b0` and remained solely
on `analysis/task-50-minimal-beta-conditioned-effective-coupling`.

| Boundary | Published and remotely verified commit |
| --- | --- |
| Original 50A precision audit | `e4c72143fe7d71ae0034d15887312062626b4636` |
| Authorised 50A continuation before implementation | `cc69edfb88fbcf17b6d60e0cbc70c0130c0d978a` |
| 50B before frozen result reuse | `537ce3204d11a1a6cc51b1dd39032d12cb42df5e` |
| 50C before final report | `88997ef6888498e6ea6bdd9bd4acd7e5ec108dbc` |
| 50D final report and ledgers | `7cb032e8bb9372c910d0e5400712ea79beb23e8f` |

The original audit verification closeout was
`9bb9c9513a5816b872d75fd7ab86ee621effdb73`. Each scientific boundary updated
the mandatory ledger before dependent work. Machine readable checkpoint and
publication receipts are under `output/`; the final verified 50D commit is
recorded in `output/publication_50D.json`. All nine 50D artefact hashes match
the fetched remote commit. Task 50 is stopped; the final metadata closeout
records that verification without further scientific work.

`run_verification.py` records test outcomes and tested source hashes;
`reuse_frozen_results.py` reads existing files without importing a model or
solver. Both guard against optional repetition after success. The final
integrity audit checks original checkpoint receipts and unchanged source and
frozen artefacts without repeating scientific evaluation. The two original
test logs are explicitly included at 50D because the general repository ignore
rule excluded their `.log` paths from the earlier commit.

The authoritative cumulative result is R50D in
[`MANDATORY_RESEARCH_LEDGER.md`](../../docs/MANDATORY_RESEARCH_LEDGER.md).
Task 50 is complete under the stated precision qualification. No Task 51 or
Task 52 work has begun; the two remaining scientific shots are unstarted.
