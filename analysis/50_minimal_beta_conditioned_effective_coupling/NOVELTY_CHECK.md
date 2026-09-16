# Task 50 novelty check

## Proposed question

Can the already successful but target-selected Task 41 apical-secretory coupling be rewritten as the **smallest beta-conditioned effective coupling**, so that it retains the already demonstrated CCh+IPR secretion deficit while exactly nesting WT, AE4-KO CCh-only, rest, and AE2-KO behaviour?

This is not a new mechanism search and is not an attempt to identify a microscopic pathway.

## Mandatory prior work checked

The following master-ledger entries and source files overlap directly and are binding:

- `docs/MANDATORY_RESEARCH_LEDGER.md`, especially R38-R42, R48-R49, sections 5-10.
- `docs/PHENOTYPE_TARGET_CONVENTION.md`.
- `analysis/40_ae4_equal_cation_routing/final_answer.md`.
- `analysis/41_ae4_loss_algebraic_design/final_answer.md`.
- `analysis/41_ae4_loss_algebraic_design/attempt_ledger.md`.
- `analysis/42_catalan_2025_ae4_mechanism_classes/final_answer.md`.
- `analysis/48_joint_experimental_constraint_reconstruction/JOINT_EXPERIMENTAL_CONSTRAINT_REPORT.md`.
- Task 49 final report on `analysis/task-49-camp-vrac-secretory-branch-reconstruction`.
- `src/modern_full_model/ae4_cacc_recruitment.py` and the Task 41 selected driver.

## What is already answered

Task 41 already proved constructively that one additional coupling acting on effective stimulated apical chloride conductance can generate the required order of secretion loss. Its selected parameter was

`b = 0.10511872843288`,

which gave approximately `23.1634%` cumulative loss at AE4=5% and `30.2612%` at AE4=null under the standard combined stimulus. WT nested Task 40 exactly.

That parameter was selected using the phenotype target. It is therefore a **TARGET-CALIBRATED CONSTRUCTION**, not independent validation.

## What Task 41 did not answer

Task 41 coupled the conductance reduction to calcium activation and AE4 expression, without conditioning the added coupling on the beta/IPR arm. Interpreted literally, that construction would also alter AE4-null CCh-only secretion whenever the calcium activation is on. The primary experiment instead shows near-preserved AE4-KO CCh-only chloride uptake.

The unresolved modelling question is therefore not whether another mechanism can produce 30% loss. Task 41 already answered that. The unresolved question is whether the same one-dimensional constructive direction can be written as a **protocol-conditioned effective law** that is exactly off when beta/IPR stimulation is absent.

## Proposed reduced law

Use exactly one added scalar coupling and no signalling dynamics:

`g_Cl,eff = g_parent * [1 - lambda * a_Ca * beta * (1 - e_AE4)]`

with fixed

`lambda = 1 - 0.10511872843288 = 0.89488127156712`.

Here:

- `a_Ca` is the existing normalized calcium/channel activation already used by Task 41;
- `beta` is the existing effective protocol input, not a new state or fitted signalling variable;
- `e_AE4` is the existing AE4 expression fraction;
- `lambda` is inherited exactly from the already phenotype-calibrated Task 41 construction and must not be searched or refit.

Expanding Task 41's selected factor gives

`(1-a) + a[b + (1-b)e] = 1 - (1-b)a(1-e)`.

Therefore at `beta=1` the proposed law is algebraically identical to Task 41. At `beta=0` the new factor is exactly one.

## Why this is not repetition

Task 50 does **not** rerun the Task 41 mechanism search or tune a new parameter. It performs a model reduction/reinterpretation:

- preserve the already demonstrated combined-stimulus construction;
- make the added effective coupling conditional on the beta/IPR arm;
- prove exact nesting of protocols/genotypes where the coupling must be absent;
- remove the unjustified microscopic claim that AE4 directly recruits TMEM16A.

Task 49's swelling-gated VRAC construction is not this law. Task 50 adds no separate channel, swelling gate, conductance calibration, beta-NKCC arm or signalling dynamics.

## Parameter/evidence classification

- `lambda`: **TARGET-CALIBRATED CONSTRUCTION**, inherited from Task 41; not independently measured.
- `beta`: existing effective protocol input; **ASSUMPTION/PROTOCOL MAP**, not cAMP/PKA kinetics.
- `a_Ca`: existing parent activation coordinate; no new calcium parameter.
- No source parameter is newly fitted.

## Falsification / stop criteria

Stop Task 50 without searching if any of the following occurs:

1. the proposed factor does not algebraically reduce to the selected Task 41 factor at `beta=1`;
2. WT (`e_AE4=1`), rest (`a_Ca=0`), CCh-only (`beta=0`) or AE2 KO (`e_AE4=1`) fails exact parent nesting;
3. the implementation cannot be inserted before the same electrical closure as Task 41 without changing unrelated model equations;
4. the `beta=1` Task 50 vector field fails exact numerical equivalence to the frozen Task 41 selected vector field at matched states/inputs;
5. charge/current/Cl source bookkeeping differs from Task 41 for reasons other than the declared beta factor.

If any stop criterion is met, report it. Do not open another mechanism family or parameter search in this run.

## 50A audit qualification

The repository and all 53 fetched remote heads were searched using beta
conditioning, AE4/CaCC/TMEM16A/apical/secretory and recruitment synonyms plus
the reported coefficient literals. The exact heads, search expression, match
counts and matching source blob identities are recorded under `output/`.
No implemented beta conditioned version of the Task 41 factor was found in
this search. The matching Task 18 and Task 21 source files concern fixed
balance bookkeeping and capacity names, respectively.

The algebra above is correct for the displayed shortened value of `b`.
However, the frozen Task 41 design, input manifest and driver all use
`b=0.10511872843288446`. The mandated Task 50 lambda has complement
`0.10511872843288`, which is different. At beta one the factor difference is
exactly `-4.46e-15*a*(1-e)` when the decimal literals are treated exactly.

This precision qualification is appended without changing the predeclared
question or frozen files. It is a failure of the required exact inheritance
premise, not a failure of the generic algebra and not a biological falsification.
See `EQUIVALENCE_PROOF.md` and `output/coefficient_audit.json`. Task 50 stops
at 50A before implementation under the stated stop rule.

## Authorised continuation after 50A

The user explicitly authorised keeping the exact Task 50 lambda and treating
the coefficient precision difference through numerical equivalence. This lifts
only the rounding related stop above. D50-06 and D50-07 record the scope and
comparison tolerances before dependent compute. The original novelty audit is
reused; no further mechanism search is needed. All Task 40/41 frozen files
remain unchanged. Exact parent nesting remains mandatory. The construction
remains a TARGET-CALIBRATED CONSTRUCTION, with the frozen mutant numerical
results eligible only for expressly qualified reuse after the checks pass.
