# Task 49: cAMP/VRAC reconstruction and resting source geometry

## Decision

**The single permitted swelling-gated cAMP chloride branch is not an accepted stimulated reconstruction.** Its conserved implementation passes nine unit tests, but independent cAMP data do not identify a usable conductance in the unchanged Task 48 system. A parameter-independent prediction frozen before the held-out comparison also fails: AE4-null IPR stimulation cannot initiate any chemical response from exact own-genotype rest, whereas the experiment reports positive chloride uptake. No conductance, secretion deficit or combined-stimulus uptake was fitted, and no post-reveal retuning occurred.

The mandatory resting analysis was completed in the same run. The complete fixed source matrix has rank 48 and coefficient nullity 36. It spans the conditional measured-state residual, but an exact AE4-null alkalinity certificate excludes that residual from the native fixed-state addition cone. The obstruction is the required relation **`J_NHE1=J_AE2=J_CO2,total`** under the inherited zero-rest NBC gate. At the declared measured-state lift, NHE and reverse AE2 both add alkalinity; inward basal NBC would worsen the discrepancy.

The calculation localises the equation failure but does **not** identify a unique shared correction once unmeasured resting coordinates and source null directions are acknowledged. No resting model edit is justified, and no joint reconstruction is run. Thus this run supplies a reproducible equation-level diagnosis, not a prospectively validated reconstruction or a uniquely implementable Task 50 correction. It does not claim the stronger requested resolution has been achieved.

## Position in the repository

Tasks 12 and 38–42 established the insufficiency of the tested AE4-only kinetics/stoichiometries and the large inherited NKCC compensation. The source-fixed Palk/Benjamin core and equal cation routing did not supply the phenotype. Task 41 demonstrated a target-selected apical-export construction, but its direct AE4-dependent TMEM16A recruitment and timing/state failures were not independent validation. Task 46 showed that full reference export remains feasible without AE4; Task 47 found no binding K-recycling/pump ceiling on the frozen trajectories. Task 48 then exposed the exact architectural defect: removing AE4 removes beta input from the inherited core and forces CCh and CCh+IPR to coincide under a common observation map.

The exclusion audit also corrected a factual omission in the Task 49 premise. Task 13B had already implemented this exact positive-swelling law as the isolated **Vbeta/VRAC diagnostic**. Its final archived record contains 389 roots, 69 with positive gate, 10 passing WT rest and **zero** passing WT rest while activating the gate. Its capacity was unlicensed and it was never inserted into production electrical closure. Those results were read, preserved and not recomputed. The correct novelty here is the separate current's integration into the active closure and the constrained analysis on the Task 48 system—not the invention or first test of this gate. The original exclusion entries remain preserved beneath an explicit correction.

Task 19's bounded numerical limitation, Task 31's nonconvergence and Task 42's source-class scope are retained as qualifications, not converted into global impossibility claims. All excluded mechanism searches remained closed. See [history audit](advisory/history/history_audit.md).

## Independent physiology and one declared law

[Catalán et al., PNAS 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4343136/) establish a TMEM16A-independent IPR secretory response, inhibitor/current evidence supporting a VRAC-like pathway, and delayed swelling/secretion. This supports the existence of a separate pathway; it does not uniquely establish a molecular identity or a linear positive-swelling gate. The predeclared minimal test was

`G_V = g_V beta max(V_i/V_rest,genotype − 1, 0)`, with `g_V≥0`,

`I_V = G_V(V_a−E_Cl)` and `J_Cl,cell→lumen = −10^15 I_V/F`.

The branch has no direct AE4-expression factor. It contributes to both chloride amount balances and both coupled current equations. TMEM16A is separately labelled and unchanged. All inherited NKCC, AE4, AE2, NHE1, NBC, pump, K, calcium, geometry, water and paracellular laws and parameters remain unchanged. In particular, no beta-NKCC arm was added. REST and CCh-only nest exactly when beta is zero.

The saved source metrology distinguishes absolute whole-cell current from current density, blocked current from remaining current, and biological SEM from raster/axis uncertainty:

| Independent quantity | Source/extraction result | Qualification |
|---|---:|---|
| IPR, DCPIB, NPPB concentrations | 5, 50, 100 µM | Distinct source protocols |
| IPR cell swelling | Narrative 12.5 ± 0.2% | Figure 5C terminal proxy is approximately 1.16; discrepancy preserved |
| CCh cell shrinkage | 10 ± 1% | Opposite direction from IPR |
| IPR swelling/secretion half-time | >60 s | No fitted arbitrary delay |
| Figure 4B IPR gland output, control | 10.02 ± 0.74 µL/10 min | Raster-derived mean/SEM; n=6 |
| Figure 4B IPR gland output, TMEM16A KO | 7.17 ± 1.37 µL/10 min | n=6, reported p=.092 is not an equivalence margin |
| Figure 5E IPR-induced current | 208.6 ± 35.2 pA | Absolute whole-cell current at −60 mV; n=4 |
| Figure 5E DCPIB-blocked component | 107.8 ± 17.2 pA | Blocked magnitude, not remaining current or necessarily wholly apical |

Original figures, pixel coordinates, axis calibration, digitisation uncertainty, overlays, version correction and hashes are preserved in `input/pnas2015/`. A first Figure 5C marker corridor was corrected before scientific use; both versions and the independent review remain available. Chronic TMEM16A deletion also changed NKCC1/AQP5 protein abundance in the source, so the unchanged-model acute mask is only an idealised channel-separation check. Unsupported low-Cl/bicarbonate-free protocols were not approximated and scored as exact experiments.

## Stage I: freeze, prospective implication and reveal

49B established correct off-state nesting, chloride/current sign and reversal, conservation, channel independence and electrical closure. The arbitrary conductance and imposed swelling used by the nine tests are explicitly synthetic verification fixtures, never production parameters or initial conditions.

49C accepted **no calibrated conductance**: the frozen value is `null`, not zero. Holding voltage is not chloride driving force. The accessible source did not provide a verified matched map from blocked apical current, instantaneous swelling and chloride reversal to the model conductance. Access to the supplement failed; this is not a claim that the supplement lacks the information. Whole-gland amplitude also contains an unidentified gland-to-model scale. Scale-free shapes/ratios could in principle help, but this run did not demonstrate a finite identified scale from them under the fixed activation/protocol restrictions.

The saved Task 48 WT IPR trajectory was inspected without reintegration. All 61 saved samples have nonpositive swelling gate; terminal volume change is −0.053089%, opposite the source swelling. This is a sampled result, not a continuous-time enclosure. It was not used to justify an unproved universal WT deadlock.

49D was pushed before the new held-out comparison. It freezes unavailable numerical predictions for all beta-containing protocols requiring the missing conductance, references unchanged REST/CCh artefacts by hash, and proves one exact implication. In AE4 KO, the inherited IPR core has no beta input; at its exact rest the new gate is also zero. The constant chemical state, with an evolving but disconnected AE4 regulatory coordinate, solves the augmented locally unique ODE for every finite nonnegative conductance. Thus **AE4-KO IPR predicts zero chloride concentration and SPQ uptake slope** for any fixed optical calibration. Numerical equilibrium error is not a physical swelling seed.

After that remote boundary, 49E compared the frozen implication with [Peña-Münzenmayer et al., JBC 2015 Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC4409235/#T1): AE4-KO IPR uptake is **0.20 ± 0.03 ×10⁻³ s⁻¹**, n=7. The model–data discrepancy is −6.67 reported SE units, a descriptive discrepancy rather than a calibrated p-value. A conductance adjustment cannot repair the exact initiation failure. This rejects this one gate as a sufficient reconstruction; it does not reject the experimentally supported cAMP pathway.

| Requested stimulated conclusion | Task 49 result |
|---|---|
| Separate TMEM16A-independent current | Implemented and conservation-tested |
| Independent calibrated cAMP reconstruction | Not identified; no production-eligible parameter |
| CCh-only preservation | Exact nesting with the parent |
| AE4-KO IPR-only positive uptake | Fails by the frozen exact-rest theorem |
| Combined-stimulus KO uptake and CCh contrast | Numerical prediction unavailable |
| Early preservation then delayed secretion loss | Unavailable |
| Ten-minute AE4 deficit | Unavailable; 35% was never fitted |
| AE2 beta-containing secretion neutrality | Unavailable |
| Direct AE4-dependent TMEM16A recruitment | Absent |
| Physiological resting Cl/pH | Unchanged by Stage I; assessed independently below |

CCh can move volume under combined stimulation, so the IPR-only theorem is not a proof that all CCh and CCh+IPR trajectories remain equal after adding this branch. Neither unavailable output nor an old parent prediction is relabelled as a successful new prediction. No alternative gate, source family or parameter search followed reveal.

## Stage II: measured states as equations

The full analysis is in [RESTING_SOURCE_GEOMETRY.md](RESTING_SOURCE_GEOMETRY.md). All four measured Cl/pH pairs were used as stated in the addendum, with genotype-specific cached references for unmeasured quantities and a common WT law for both controls. No secretion phenotype enters this calculation.

The single predeclared algebraic lift enforces measured Cl/pH, both charge constraints, cell water stationarity and `J_CO2,total=J_NHE1`. It retains intracellular Na concentration and complete lumen from each cached genotype rest. These are explicit nuisance assumptions, not observations. Each Cl/pH pair originally leaves eight unmeasured core-state degrees of freedom after charge constraints; the two added stationary identities reduce this to six retained degrees. The resulting states do not solve all resting balances.

The complete 64×84 source matrix contains all 21 predeclared pathways/bookkeeping columns in each cohort, including dormant NBC/VRAC signatures and genotype deletions. Rank is48 and coefficient nullity 36. After declared row/column scaling, the target norm is 1.29069 and the orthogonal residual norm is 1.36×10⁻¹⁵. Exact source decomposition matches the native RHS within 1.12×10⁻¹⁶. There is **no missing direction orthogonal to this complete source span**. This alone cannot identify a correction, because the native vector field already has the form `F=Sj`.

The nontrivial result is the native cone obstruction. In the AE4-null conditional state:

| Balance quantity | fmol/s |
|---|---:|
| NHE, H | +0.008763263 |
| AE2 inward-Cl cycle, E | −0.000631575 |
| TA residual, H−E | +0.009394838 |
| Required TA correction | −0.009394838 |
| Required TIC correction under C=H | −0.009394838 |

The TA row of every available native source direction is nonnegative, whereas the required TA correction is negative. That row is an exact Farkas separator; the cone failure persists under additional current-closure or genotype-sharing constraints at the same fixed states/directions. Nonnegative capacity reweighting cannot fix it while maintaining nonzero NHE/AE2 activity. The full other cell and lumen residuals remain visible; the TA/TIC entries alone are not a complete corrective flux.

For all genotypes, `TA_dot=H+2B−E−2A` and `TIC_dot=C+2B−E−2A`. In the native AE4 knockout, `A=B=0`, forcing **H=E=C**. With D the total CO2 permeability, Cbar the permeability-weighted reservoir CO2, and `Ccrit=(Cl_i HCO3_o/Cl_o)10^(pKa1−pH_i)`, that condition is exactly

`H = G2 tanh(log((Cbar−H/D)/Ccrit)/w_AE2)`.

It requires `0<H<min(G2,D(Cbar−Ccrit))`. Here the right-hand upper bounds are 0.005 and 0.001001641 fmol/s, both below H=0.008763263. At the declared state, AE2 additionally has negative affinity, so increasing its capacity cannot repair the sign. This is the precise acid/base compatibility failure, not a proposal to revisit excluded NHE, AE2 or NKCC searches.

Task 36's resting NBC gate is a nesting assumption and was audited only after this residual geometry was available. The NBC inward affinity is +2.733 at the unchanged AE4-KO voltage. Opening inward NBC would add +2TA per cycle, opposing the required negative source. It is not selected as the resting correction. Exact source identities further prevent transporter inference from rank alone:

`s_AE4=2s_AE2−½s_NKCC1`,

`s_NBC=2s_NHE1+2s_CO2,b+⅓s_pump,b+⅔s_K,b`.

Measurement SEs were propagated by local derivatives, without another basis or parameter fit. The conditional KO TA residual has working independent-error SE 0.00112827 fmol/s; the unknown Cl–pH correlation gives a linearised SE upper bound 0.00123094. These are not global confidence bounds. Unmeasured intracellular Na and lumen carbon have nonzero local sensitivities, and their uncertainty is unreported. The cone certificate therefore does not establish impossibility over the entire observation-compatible state manifold.

No uniquely identified shared law follows. The requirement carried forward is the scalar established-genotype carbon/alkalinity compatibility equation, its full conserved residual and its hidden-state dependence. Those must be resolved with independent state/source information before a correction is named. A source-norm choice, unmeasured genotype multiplier or saliva-output fit cannot supply the missing identification. This is a bounded mathematical handoff, not a new list of candidate mechanisms.

## Publication, verification and limits

All work is confined to the requested branch, starting at `2bf53b773c14f3c49f539957ea3f1964bc8d227f`. Only the orchestrator edited canonical files, evaluated the model and published. Five workers supplied bounded read-only advice. The original model sources, frozen parameters, prior analyses and manuscript remain unchanged.

| Checkpoint | Decision | Verified remote commit |
|---|---|---|
| 49A | Evidence/exclusion freeze; historical Vbeta correction | `e597ae7e800317ee52e71bc77d806e7cb9aadeb4` |
| 49B | Separate current and nine passing tests | `bcf70bbc833ea59ff894f76f018652415806360d` |
| 49C | Independent calibration blocked; conductance null | `253d4b213c7e73970ff742d68da3d0422373f2a6` |
| 49D | Pre-AE4 unavailable predictions and exact invariant | `c1e7d3afa2e3c4f037c3a5f1c3bc4a5b6eaa3e3c` |
| 49E | Held-out IPR failure; no retuning | `d6e0404a9b87279b507b998edea01ed0820e2b84` |
| 49G | Fixed-basis residual/cone analysis; no unique correction | `62fb3a1cf5abd0556eb81dec3e5d1e1e6a3d4382` |
| 49H | Inapplicable: correction not uniquely identified | No model edit |
| 49I | Inapplicable: two independently accepted pieces absent | No joint prediction |
| 49F | This final report and integrity audit | Recorded by final publication closeout |

Every executed boundary has an updated status file, machine-readable receipt, commit, push and independently fetched remote SHA before dependent work. Receipts contain hashes of the artefacts at their publication commit; later edits do not silently rewrite historical receipts. The final closeout seals the post-push verification record without further scientific work.

Execution limits are explicit: **zero production trajectories, zero inferred conductances, zero constrained least-squares fits, zero resting model edits, and no mechanism/parameter search**. Nine implementation tests passed. StageII used one distinct fixed basis and one central state set; an export dictionary-key error required an identical replay after the first SVD. The actual two SVD invocations and failed/accepted RHS counts are disclosed, rather than claimed as one numerical invocation. No published calculation was recomputed. Advisory reviews accepted the algebra and uncertainty with metadata clarifications; all 72 frozen source hashes passed.

The primary numerical outputs are [resting residuals](output/resting_residuals.csv), [fixed source matrix](output/fixed_source_matrix.json), [SVD/projection/cone certificate](output/source_svd_projection.json), [uncertainty](output/resting_uncertainty.json), [parameter freeze](output/parameter_freeze.json), [pre-AE4 record](output/pre_AE4_prediction.json) and [held-out comparison](output/heldout_comparison.json). These distinguish a failed scientific gate from a missing numerical output and preserve exactly what the next task can—and cannot—infer.
