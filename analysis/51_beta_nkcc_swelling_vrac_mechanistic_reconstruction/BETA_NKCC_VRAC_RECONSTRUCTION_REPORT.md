# Task 51: beta-NKCC → swelling → VRAC mechanistic reconstruction

**Controlling decision: UNRESOLVED / independently unidentified full model, with a formal local initiation advance.** The authorised upstream beta-NKCC term removes the exact beta-blind initiation premise of Task 49. It has not produced an independently fixed model capable of predicting the held-out AE4 secretion phenotype. Task 50 remains the accepted quantitative proof-of-sufficiency benchmark.

The direct current-to-conductance route was unavailable. The prescribed swelling fallback was executed, but its two frozen bracket endpoints did not bracket the independent swelling target. Therefore g_V remains **null, not zero**. This is a scoped identification failure under the declared observation map and procedure. It is not proof of structural non-identifiability, global mechanism impossibility, or the old Task 49 zero-gate result.

## Source reconstruction and its limits

The source architecture was published before model edits at 51A. Primary-source checks support beta/cAMP regulation of salivary NKCC: [Chaib et al. 1998](https://pubmed.ncbi.nlm.nih.gov/9880083/) reports approximately 2.5-fold IPR-stimulated NH4 influx in rat submandibular acini, while [Paulais and Turner 1992](https://doi.org/10.1172/JCI115695) reports approximately 3-fold functional regulation in rat parotid. The latter's first author is Paulais, qualifying the earlier handoff shorthand. Phosphorylation/recruitment studies support the regulatory pathway without establishing a simple direct purified-PKA mechanism.

Following the prompt, the central effective multiplier was fixed at 2.5 and the source sensitivity at 3.0. These are **cross-species assumptions/sensitivities**, not measured mouse-SMG constants or mouse confidence bounds. The older sixfold context was not adopted as a rescue value. The existing mouse beta/PKA→AE4 arm was retained once, without a new signalling state.

[Catalán et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4343136/) supplies independent mouse-SMG evidence for IPR swelling, a chloride conductance, blocker sensitivity and TMEM16A-independent secretion. It does not establish molecular LRRC8 identity or an entirely apical whole-cell patch current. The main paper was inspected; official supplement routes returned challenge/access errors. **The supporting information was not retrieved or inspected.** This is an access limitation, not a finding that the supplement lacks the needed information.

A defensible direct map would require consistent signed current and voltage, matched intracellular/extracellular chloride and reversal potential, patch swelling, the apical current fraction and cell scaling. These were not jointly identified. Figure 2's chloride reversal cannot be transferred to Figure 5's nystatin IPR protocol. Outward rectification further qualifies an ohmic conversion. No fitted nuisance factor was introduced.

Task 49's validated extraction was reused unchanged: IPR current magnitude 208.59375+/-35.15625pA and DCPIB-blocked magnitude 107.8125+/-17.1875pA at -60 mV. The narrative swelling 12.5+/-0.2% differs from the approximately 16.03% terminal graphical proxy near 451 s; the raster envelope is not a biological SEM. The narrative magnitude remained the predeclared fallback target. No alternate extraction or source clock was fitted. Detailed evidence, access routes and qualifications are in [SOURCE_ARCHITECTURE_FREEZE.md](SOURCE_ARCHITECTURE_FREEZE.md).

## The single implemented mechanism

Only the inherited NKCC activity regulator gains the additive beta increment:

`M_NKCC = M_Ca + beta * (M_beta - 1)`.

At beta=0, the complete parent regulator is unchanged; at basal calcium and beta=1, its multiplier is M_beta. At the frozen CCh endpoint, total multipliers are 3.25 and 3.75 for the central and sensitivity gains. The Palk/Benjamin concentration-response core, reversal and exact 1Na:1K:2Cl stoichiometry remain unchanged. Both inherited layers that recompute NKCC consume the same regulator result.

The existing Task 49 current is reused without changing its law:

`G_V = g_V * beta * max(V_i / V_rest,genotype - 1, 0)`.

It participates in electrical closure and cell/lumen chloride balances as a current. It contains no AE4-expression multiplier. The Task 50/41 effective conductance wrappers are absent from the Task 51 model dependency chain. NHE1, NBC, pump, potassium pathways, calcium, AE4 routing/stoichiometry, water, paracellular transport and geometry remain frozen. No new state, signalling ODE, cap, synergy, delay or gate family was introduced.

The conservation identity applies to **total** CaCC plus VRAC chloride export:

`J_total = 6P - H + 2*dNa_i/dt - dTA_i/dt - dCl_i/dt`,

with amount derivatives, and `J_total = 6P - H` at stationarity. Added NKCC loading does not by itself prove a sustained AE4-dependent secretion deficit.

Seventeen distinct software tests pass, covering parent nesting on both baths, reversal, stoichiometry, current/source signs, current closure, independent channel masks, genotype independence of the gate and the beta-NKCC-off regression fixture. There were 18 test invocations: one original assertion failed on extra diagnostic dictionary keys; only that assertion was corrected and its targeted cases rerun. Model code did not change. Both attempts remain published. These tests establish implementation properties, not biological validity.

## Prescribed independent identification

The fallback used WT Track2 genotype rest with M_beta=2.5. Its sole residual was mean relative swelling over 300–600 s minus0.125. That broad-window amplitude operator is an explicit assumption, not a measured source window; it extends beyond the graphical source trace. The narrative SEM does not cover observation-map or species uncertainty.

One deterministic scalar procedure was frozen at 51A: nonnegative numerical bracket [0,1e-6] S, Brent tolerance 1e-13 S, maximum 40 iterations if bracketed. The bracket is a numerical device, not a source-derived conductance interval. Both endpoints were executed after remote 51B verification.

| Conductance trial | Mean 300–600 s swelling | Endpoint swelling | Residual from 12.5% target | Numerical/physical checks |
|---|---:|---:|---:|---|
| 0S | 5.750804% | 6.133046% | -6.749196 percentage points | Pass |
| 1e-6S | 2.436621% | 2.421298% | -10.063379 percentage points | Pass |

The same-sign residuals did not bracket a root; no Brent iterations or accepted fit followed. The positive-g trial activated VRAC at all 600 positive integer-second samples. Both trials therefore demonstrate responsive WT chemistry and swelling, but neither is an accepted mouse parameterisation. Two endpoints do not establish monotonicity, absence of an interior/other root or a global no-solution theorem. No extra bracket search, alternative target, phenotype calibration or rescue gain followed.

All attempted trajectories, state vectors, separate NKCC/CaCC/VRAC/AE4/AE2 fluxes, voltages, ionic/carbon/pH/volume states, water outflow, cumulative secretion and conservation diagnostics are preserved in `output/trajectories/`. The inherited solver settings and physical ranges were not relaxed. Every accepted solver endpoint and integer-second sample was checked; these are sampled numerical checks, not continuous-time bounds.

The gain-3 full-trajectory sensitivity, TMEM16A-off control, calibrated paired VRAC-blockade control and quantitative patch-current validation remain **unexecuted/unavailable without an accepted conductance or observation map**. The g=0 bracket endpoint is not promoted to paired validation. The parameter file was frozen at 51C; see [INDEPENDENT_IDENTIFICATION.md](INDEPENDENT_IDENTIFICATION.md) and `output/parameter_freeze.json`.

## Prospective local prediction and immutable reveal

The 51D prediction boundary was independently verified at `e66763952907d8d70b0c2778e08da43611ac67c1` before any held-out comparison. It contains 17 explicitly labelled case/sensitivity rows: all beta-containing quantitative predictions are unavailable; three CCh-only cases reference exact beta-zero parent initial-value problems and hashed Task 48 caches. The two simulation tracks are kept separate. No calibration trial is substituted for an accepted prediction.

There is a parameter-independent local implication. At exact AE4-null genotype rest with inward basal NKCC N0, beta stimulation adds `deltaN=(M_beta-1)*N0` with source vector `deltaN*(1,1,2,0,0)`. Initially V'=0, but

`delta[Cl]'(0+) = 2*deltaN/V0 > 0`,

`V''(0+) = 4*(L_apical+L_basolateral)*deltaN/V0 > 0`.

Using cached coefficients, the central multiplier gives an initial chloride derivative increment 0.2765211776mM/s and volume acceleration 0.0308708243pL/s2; gain 3 gives 0.3686949035 and 0.0411610990 respectively. These are approximate cached coefficients of an exact-equilibrium theorem. A tiny cached resting residual is not used as a physical seed.

Thus any finite positive g opens the specified gate locally after swelling begins. At g=0, loading and swelling can start but VRAC current remains zero. This changes the upstream premise of Task 49. It does not establish sustained chloride export, fluid secretion, SPQ recovery or the required phenotype magnitude. See [PRE_AE4_PREDICTION.md](PRE_AE4_PREDICTION.md).

## Held-out assessment

All eight required questions are addressed in [HELD_OUT_COMPARISON.md](HELD_OUT_COMPARISON.md), with machine-readable fields in `output/heldout_comparison.json`. The classification remains unresolved after reveal; no retuning occurred.

| Outcome | Task 51 conclusion |
|---|---|
| AE4-null IPR initiation | Positive local chemical loading and swelling are formally established under the stated assumptions. Sustained assay response is unpredicted. |
| AE4-null secretion deficit and persistence | Unavailable; no distance to the source 35+/-4.7% deficit or Task 50 can be computed. |
| 5%-AE4 secretion deficit | Unavailable. |
| CCh versus CCh+IPR uptake | Exact beta blindness is removed; observed direction and magnitude remain unpredicted. |
| NKCC compensation versus beta demand | Unresolved without accepted matched WT/KO trajectories. |
| AE2 secretion neutrality | Beta-protocol prediction unavailable; CCh-only parent identity is retained. |
| Physical constraints | Both numerical calibration trials pass sampled checks. Independent swelling calibration failed to bracket; accepted-model validity is unavailable. |
| Chronic rest | Unchanged experimental chloride/pH mismatch. |

The AE4-null source SPQ values are 2.30+/-0.10 for CCh, 0.90+/-0.09 for CCh+IPR and 0.20+/-0.03 for IPR, in 1e-3s^-1. Initial model chloride loading in mM/s cannot be equated to that recovery metric. The positive local result is not claimed as quantitative reproduction.

Unchanged Task 48 rests predict WT/null chloride 57.854950/56.954138 mM against 50.10+/-1.50/36.50+/-1.60mM; WT/null pH 6.883125/7.224598 against 6.91+/-0.07/6.89+/-0.02. Passing broad physical-domain checks does not resolve these source discrepancies.

The Task 50 benchmark remains **TARGET-CALIBRATED CONSTRUCTION**: cumulative secretion 0.9925245441pL WT, 0.7626224586pL at 5% AE4 and 0.6921751966pL at null; deficits 23.163365% and 30.261151%. The null endpoint-flow deficit remains 27.038018%, with a substantial saved 60–600 s cumulative deficit. These values were compared only after 51D and were never Task 51 residuals or inherited multipliers.

## Publication, integrity and work budget

All dependency boundaries have separate status, ledger and machine-readable receipt updates. Each preceding remote ref, commit parent/tree and exact checkpoint receipt was independently read before dependent work.

| Checkpoint | Verified remote commit | Boundary |
|---|---|---|
| 51A | `9935d1b1ad5a7cfad78236b29cdc70ad5eeb5be8` | Source/architecture/fallback freeze |
| 51B | `f386421bf7c69181bc8c28f73e963432524748a4` | Implementation verification |
| 51C | `f35e30150d733c4aad6867c99b3a55aa31a0b531` | Independent identification failure and parameter freeze |
| 51D | `e66763952907d8d70b0c2778e08da43611ac67c1` | Immutable pre-AE4 predictions/availability |
| 51E | `07d7a6b31652c2de10e1fe03bda0af97c31357db` | Held-out comparison without retuning |
| 51F | Containing publication commit; final remote verification follows this receipt | Final report, controlling ledger decision and stop |

The final file/Git audit passes 128 historical receipt artifact hashes at their own commits, 92 frozen input hashes and 6 cached parent-result references. Original scientific source/results/archives and binding inputs are unchanged; only Task 51 files and appended mandatory-ledger entries differ from the requested start. The old ledger is an exact prefix. The tested model, parameter freeze and 51D prediction file remain unchanged across reveal. The archive branch was independently re-read and still points to `b16c30094b95f61a73d8f1977cd58e79c7bb50f6`.

Publication used the authenticated Git Data API because direct shell transport was unavailable; remote trees match exact local Git trees. A raw-commit timezone reconstruction issue and a detected CSV newline transport mismatch were corrected without scientific changes. The unsuccessful CSV tree attempt occurred before any remote commit/ref update. This operational qualification does not alter the verified boundaries.

Total scientific work: 17 distinct software tests over 18 invocations and 102 core evaluations; 2 calibration trajectories and 3,638 core evaluations; total 3,740 core evaluations. There were zero accepted parameter fits, zero accepted-model phenotype trajectories, zero resting solves, zero post-reveal model evaluations and no mechanism/parameter-family search. The final audit itself calls no model.

Upon independent remote verification of 51F, Task 51 is complete and stopped. This run consumes the authorised Task 51 shot; one funded scientific shot remains. **Task 52 has not been started.** No new mechanism or post-reveal rescue is proposed or executed. Task 50 remains preserved as the quantitative benchmark while the microscopic coupling and independently constrained sustained reconstruction remain unresolved.
