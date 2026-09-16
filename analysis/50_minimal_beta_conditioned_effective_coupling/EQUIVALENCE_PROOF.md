# Task 50A algebra, source mapping and frozen coefficient audit

Current disposition: the user has authorised continuation with the fixed lambda
and numerical equivalence for the stored coefficient difference. The original
audit below is preserved; see the continuation addendum at the end.

**Disposition: STOP at 50A. The four parent nesting identities hold formally;
exact equality with the actual frozen Task 41 coefficient does not.**

Classification: FORMAL DEDUCTION and provenance qualification concerning a
**TARGET-CALIBRATED CONSTRUCTION**. This is not a rejection of the scientific
sufficiency of the effective coupling, independent validation, or evidence of
a molecular AE4 to TMEM16A interaction.

The branch was verified at `8f5fefd9fc0b563a6bd963895211b4acf76972b0`.
No scientific source was modified and no model evaluation, trajectory,
stationary solve, parameter inference or search was performed.

## Exact algebra

Let `a` be the existing normalised calcium activation, `beta` the existing
protocol input and `e` the AE4 expression coordinate. Write

`m50(a,beta,e) = 1 - lambda*a*beta*(1-e)`.

The Task 41 factor expands exactly:

`m41(a,e;b) = (1-a) + a[b + (1-b)e]`

`= 1-a+ab+ae-abe`

`= 1-a[(1-b)-(1-b)e]`

`= 1-(1-b)a(1-e)`.

Consequently `m50(a,1,e) = m41(a,e;b)` for all `a,e` if
`lambda = 1-b`. This identity is exact for the reported, shortened decimal
`b_reported = 0.10511872843288` and the mandated
`lambda = 0.89488127156712`.

The parent identities follow without using any particular value of lambda:

| Condition | Substitution | Multiplier | Meaning |
| --- | --- | --- | --- |
| WT | `e=1` | `1` | Same parent for any activation and protocol |
| REST | `a=0` | `1` | Same parent at the same state and genotype |
| AE4 null, CCh only | `beta=0`, `e=0` | `1` | Same CCh only parent initial value problem |
| AE2 knockout | `e=1` | `1` | Same AE2 knockout parent initial value problem |

REST nesting does not assert that an acute null state copied from WT is a
null equilibrium. AE2 nesting means equality with the parent AE2 knockout,
not equality of that knockout with WT.

## Code mapping and conditional vector field proof

| Symbol or operation | Existing code and precise meaning |
| --- | --- |
| `a` | `normalized_secretory_activation` in `src/modern_full_model/nbc_minimal.py`: clip `(calcium-resting_calcium)/(fully_recruited_calcium-resting_calcium)` to `[0,1]`; frozen endpoints are `0.058` and `0.25` micromolar |
| `beta` | `model.stimulus(t).beta_input`, returned by `SecretagogueProtocol` in `src/modern_full_model/validation.py`; zero without active IPR and one for the frozen standard combined stimulus |
| `e` | `genotype.ae4_expression`; AE2 deletion leaves this coordinate equal to one |
| Parent | `MinimalNbcModel(StimulatedNkcc1Model(ModernFullModel(...)))` with Task 40 equal AE4 cation routing; constructed by `analysis/40_ae4_equal_cation_routing/validation_common.py` |
| Insertion | `Ae4DependentCaccModel.effective_model` in `src/modern_full_model/ae4_cacc_recruitment.py` replaces only `parameters.membranes.g_cl_apical_S`, before evaluating the parent |
| Electrical closure | `evaluate_membrane_closure_with_nbc` in `nbc_minimal.py` forms `g_cl_a = g_cl_apical_S * gate`, solves the same coupled apical, basolateral and paracellular current equations, and computes the chloride current and both compartment sources |
| Selected driver | `Task41SelectedModel` in `src/modern_full_model/task41_selected.py`; its parameter literal is the longer frozen value below |

The parent calcium Hill gate remains in the electrical closure. The normalised
coordinate `a` is the existing extra recruitment coordinate, not a replacement
Hill gate and not an AE4 regulatory state. Beta is read from the protocol, not
from the inherited dynamic AE4 activation coordinate. No new signalling state
or timing parameter is involved in the proposed change.

When the multiplier equals one, `effective_model(1)` returns the original
parent object. Applying the intended wrapper therefore leaves every RHS
component, current, source and water evaluation identical, apart from optional
diagnostic labels. These are formal properties of the prescribed insertion;
they have not been claimed as tests of an implemented Task 50 wrapper.

If the multiplier were identical to the frozen Task 41 multiplier, the same
effective parameters, protocol and state would give the same electrical
closure and hence the same complete RHS. With identical initial states, the
same locally unique initial value problem would justify frozen trajectory
reuse. The frozen standard stimulus uses `a=beta=1` for `0<t<=600 s`, and
`a=beta=0` at exact onset. All three candidate 08 cases cite the common WT rest
hash `a0b96ef16e5967de4dfe2bdabb882f2a28d478e211f920bc19cd735a5a143acb`.
Neither protocol nor initial state is the identified obstruction.

## The actual frozen coefficient prevents the requested exact identity

Three independent repository locations agree on the full stored literal:

| Location | Field | Literal |
| --- | --- | --- |
| `src/modern_full_model/task41_selected.py` | `SELECTED_NULL_CACC_RECRUITMENT` | `0.10511872843288446` |
| `results/41_ae4_loss_algebraic_design/candidate_08/design.json` | `null_stimulated_cacc_fraction` | `0.10511872843288446` |
| `results/41_ae4_loss_algebraic_design/candidate_08/frozen_inputs.json` | `selected_parameters.null_stimulated_fraction` | `0.10511872843288446` |

The production runner reads this design field through `candidate8_common.py`;
the public driver uses the same literal. Both were published in Task 41 commit
`86f135e075aa91a70a752b650d7bb894353fee04`. All seven execution file hashes in
the frozen candidate 08 input manifest still match at the Task 50 start.

With exact arithmetic on the stored decimal literals,

`b_frozen = 0.10511872843288446`,

`1-lambda = 0.10511872843288`,

`(1-lambda)-b_frozen = -0.00000000000000446`.

Thus, at `beta=1`,

`m50 - m41_frozen = -4.46e-15 * a * (1-e)`.

In particular `a=1`, `e=0` is an exact counterexample to the requested
multiplier identity. The complement of the stored decimal would be
`0.89488127156711554`; this is shown solely to identify the provenance conflict,
not selected, proposed as a fit, substituted or implemented.

The discrepancy is extremely small. No biological or practically meaningful
secretion change is inferred from it. Floating point arithmetic and a numerical
tolerance may obscure it, but they do not make unequal declared coefficients
algebraically identical. No machine precision comparison or error bound for
the full RHS has been run or claimed in 50A.

## Gate decision and preservation of prior results

The premise that the mandated lambda is *exactly* inherited from the actual
frozen Task 41 coefficient is false. The displayed algebra is valid; applying
it to the rounded report value silently changes which Task 41 model is meant.
This qualifies the precision claim in master ledger R41 and section 7. It does
not replace or invalidate Task 41's saved secretion results.

The Task 50 prompt requires the exact audit before scientific source changes
and says to stop and publish if a required equality fails. `AGENTS.md` and
master ledger section 8 also prohibit improvising around a ledger conflict.
Accordingly 50B, 50C and 50D are not executed. A trajectory cannot resolve this
known algebraic inequality, so the minimum additional numerical workload is
zero. No frozen mutant trajectory is relabelled as an exact Task 50 result.

Continuing would require an explicit resolution of the conflict between the
fixed decimal lambda and exact inheritance from the full frozen coefficient.
This audit makes neither change. Task 51 is not started.

## Reproducibility and novelty scope

`audit_50a.py` reads Git objects, compares coefficient literals with exact
decimal arithmetic, verifies the seven frozen execution hashes and records
the remote text search. It imports no scientific model or solver. Its output
is in `output/coefficient_audit.json`, `output/audited_inputs.json`,
`output/remote_heads.json` and `output/novelty_search.json`.

All 53 fetched remote branch heads were searched under `src/`, `prompts/`,
`analysis/` and `results/`, including the completed Task 49 branch. Search terms
cover beta conditioning, AE4 with CaCC/TMEM16A/apical/secretory synonyms,
recruitment functions, and both reported numerical literals. The four distinct
matching source blobs comprise the two Task 41 modules, Task 18 fixed balance
bookkeeping and Task 21's capacity names. Source review found no already
implemented beta conditioned form of the Task 41 factor. The staged Task 50
documents already describe the intended law. This is a scoped search finding,
not proof that arbitrary renamed implementations cannot exist.

The overlapping Task 40, 41, 42, 48 and remote Task 49 final reports were read.
The original source and reports remain untouched, including historical timing
language superseded by the current phenotype convention. No excluded axis was
reopened.

## Continuation addendum: numerical comparison explicitly authorised

Following publication of the original audit, the user authorised retaining
`lambda=0.89488127156712` and using numerical equivalence for the coefficient
rounding difference. This supersedes the operational stop, not the exact
arithmetic finding. The four parent identities remain exact. At beta one,
exact equivalence holds with the displayed Task 41 equation at
`b=0.10511872843288`; comparison with the full stored coefficient is numerical.

For `0<=a<=1` and the declared expression fractions `0<=e<=1`, the exact
coefficient induced factor difference has magnitude at most `4.46e-15`.
D50-07 fixes a factor allowance of `8e-15` and unit specific full evaluation
allowances before any implementation check. The existing current closure and
initial states are unchanged. No modified model result has yet been produced
at this amended 50A boundary, and no trajectory error bound is asserted.

If the software checks pass, Task 41's mutant phenotype summaries may be reused
as numerically equivalent inherited results under the user's authorisation.
They must not be labelled newly simulated, independently validated or exactly
identical solutions for the two distinct coefficient literals. WT and the
other off cases continue to have exact parent identity.
