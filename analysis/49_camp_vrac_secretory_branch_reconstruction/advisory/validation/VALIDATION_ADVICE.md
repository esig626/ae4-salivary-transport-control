# Task 49 validation design — advisory only

Read against remote base `2bf53b773c14f3c49f539957ea3f1964bc8d227f`. AGENTS.md was read first; the Task 49 entry sequence, binding addendum, Task 48 ledger/protocol/operator records, and existing Vbeta diagnostic were inspected. No scientific model was evaluated, no trajectory or inference was run, and no canonical file was edited. The orchestrator must accept or reject these recommendations.

## Decision boundary before calibration

Predeclare exactly one gate/law before any new calibration. For the already implemented diagnostic, this is

`g_eff = g * beta * max(V_i/V_rest,genotype - 1, 0)`, with `g >= 0`.

The fractional-volume gate is a modelling decision, not a measured channel law. Its resting reference must be each genotype's own equilibrium volume, fixed before stimulation. Do not initialise above rest to open the gate. The inherited 30 s AE4 regulatory time constant is not an independently measured VRAC delay. Changing to a different gate after failure would be another model choice, not calibration of this single law.

Only the new conductance scale may be identified in Stage I. In particular, beta-NKCC, NHE, NBC, Ca/TMEM16A, water, geometry and the existing source laws stay frozen. The genotype-rest solver also stays frozen; report the chronic-rest mismatch separately.

## Existing diagnostic and current-state prerequisites

The claim that this topology was never represented is contradicted by `src/modern_full_model/vbeta_diagnostic.py` and `analysis/13B_modern_full_model/vbeta_diagnostic.md`. They already implement the same one-scale rectified swelling law as a non-production diagnostic, with an explicit anti-deadlock gate and an unresolved current-to-conductance observation map. Preserve that record; describe Task 49 as an updated audit/integration decision rather than first discovery.

Task 48 already saved the current scientific parent's WT IPR trajectory. Read-only extraction from `analysis/48_joint_experimental_constraint_reconstruction/output/protocol_predictions/WT__IPR_ONLY.csv` gives:

| Saved-data fact | Value |
|---|---:|
| Samples | 61 |
| Initial and maximum saved cell volume | 1.259615564408113 pL |
| Final and minimum saved cell volume | 1.2589468450620391 pL |
| Maximum saved positive fractional swelling | 0 |
| Final fractional swelling | -0.0005308916188155877 |

File SHA256: `52491822da2b6ef08987d19ded2ad326faba99cb839fedaf9f37a6e7b8155ffa`.

This is a saved-output check, not a new simulation. These samples establish no sampled activation and net shrinkage; they do not mathematically exclude an unsaved transient crossing. If the orchestrator needs continuous-path certification, use one targeted current-parent prerequisite check with a recorded crossing/event bound, not a scale sweep. The old Task 13B numerical result must not be substituted for a current-parent result.

There is also an exact conditional obstruction: under IPR only, AE4 deletion removes every inherited beta action on the conserved core. At the exact genotype equilibrium, `V_i=V_rest`, so the proposed extra conductance vanishes for every `g`. The core can therefore remain exactly stationary while the disconnected AE4 activation coordinate changes. Local uniqueness of the ODE solution makes this a zero-response result for this law, independent of `g`; numerical rest residuals must not be promoted into an activation seed. This is an equation audit, not a fitted AE4 comparison.

More generally, if a baseline trajectory remains in the gate-closed region for its whole duration, it is also the trajectory for every conductance scale: the extra current is zero along that solution. All source observables then have zero sensitivity to `g`; no finite confidence interval or best conductance can be inferred from those data.

## Independent data and identifiability

The primary source's Fig. 4 documents WT/TMEM16A-KO IPR secretion. Fig. 5 documents chloride/bicarbonate and blocker responses, volume change, and a nystatin-perforated current at -60 mV with DCPIB. These are different experiments and must retain their own observation maps. Fig. 2 concerns CaCC currents and must not be mistaken for a VRAC conductance calibration. [Catalán et al., PNAS 2015, primary figure captions](https://pubmed.ncbi.nlm.nih.gov/25646474/).

| Independent observation | Permitted use | Identification prerequisite |
|---|---|---|
| IPR whole-gland output and time course | Independent reconstruction/validation | A declared acinus-to-gland scale, or a scale-free comparison with nonzero conductance sensitivity |
| IPR WT/TMEM16A-KO output relation | Independent branch separation | Matched source protocols and genotype rests; no invented equivalence margin |
| DCPIB/NPPB effects | Evidence that a distinct anion branch contributes | Ideal branch deletion is a model intervention; drug specificity and incomplete blockade are not known conductance multipliers |
| Relative cell swelling and timing | Direct gate and delay check | Preserve fluorescence-to-volume rule, time origin, digitisation and uncertainty; do not inject measured swelling as a new forcing into the unchanged water law |
| Perforated whole-cell current | Potential absolute conductance information | Driving force/reversal, measured gate, recording normalisation, cell-to-model scale and membrane localisation must be established |
| Low-Cl or zero-HCO3 output | Source constraints | Explicitly mark blocked/conditional if substrate domains, bath access or carbonate chemistry cannot reproduce the assay |

For current `I = a*g*beta*w*(V-E_Cl)`, the data identify the product on the right unless apical fraction/model scale `a`, gate `w` and driving force are known or independently bounded. A current at one voltage is not, by itself, the coefficient `g` of a swelling law. If reported in pA/pF, capacitance introduces an additional conversion; if reported in pA, do not invent a capacitance factor. The PNAS evidence worker independently confirms that Fig. 5 reports absolute current changes and that the main text describes only partial DCPIB blockade, leaving a DCPIB-insensitive contribution possible; therefore do not equate total IPR-induced current to VRAC. Perforated recording preserves physiological ions more than conventional whole-cell dialysis but does not numerically report the intracellular chloride or volume state merely by its name. Do not pair a separate group's peak swelling with current as if it were a simultaneous measurement.

For gland output `y_j = s*Q_j(g)`, one absolute endpoint cannot separately identify `g` and unknown gland scale `s`. A time course or set of within-source ratios can potentially do so: the sensitivity columns `[s*dQ/dg, Q]` must be linearly independent after nuisance scaling. Thus lack of direct current calibration alone is not a universal impossibility proof. In this gate-closed case, however, `dQ/dg=0`, and profiling or cancelling gland scale cannot recover conductance information. Do not fit a new gland scale from the AE4 phenotype.

Digitisation must save the original figure, axes, calibration points, extracted means/error bars, script and resolution-derived error. Report source SEM separately from digitisation uncertainty and mapping uncertainty. A representative trace is not a group mean with known SEM. No biological equivalence band is licensed by non-significance.

## Predeclared calibration outcome and freeze/reveal guard

Proceed to a numerical prediction only if the single law is independently source-compatible, activated under the declared protocol, and has an independently identified finite parameter value or defensible finite envelope. Freeze parameter/protocol/operator/source hashes and the uncertainty rule. A numerical off-state diagnostic scale has no inferred value and must not be labelled a calibrated zero-conductance reconstruction.

If these gates fail, use explicit dependency receipts:

- 49C: `BLOCKED_INDEPENDENT_CALIBRATION`, with `calibrated_conductance_S: null`, named source-map/activation blockers and no parameter envelope manufactured from AE4 data.
- 49D: `PRE_AE4_PREDICTION_UNAVAILABLE`, with `prediction_generated: false`, `heldout_evaluation_performed: false`, frozen reason and dependence on the verified 49C SHA. This is an immutable null-prediction record, not a completed numerical prediction.
- 49E: `HELDOUT_NUMERICAL_COMPARISON_BLOCKED`, with `parameter_retuning_performed: false`; distinguish any analytic zero-response implication from evaluated held-out prediction residuals. Do not claim the inherited Task 48 residuals are new VRAC results.

Each remains a separately committed, pushed and remotely verified boundary. Existing knowledge of the published AE4 numbers means the experiment is procedurally predeclared, not literally blinded; the enforceable safeguard is that calibration uses only the frozen independent cAMP inputs. After publishing 49E, continue immediately to 49G. A blocked Stage I is not a reason to end the run or skip the resting calculation.

## Held-out comparison if a valid numerical prediction exists

Freeze all WT/AE4-KO/AE2-KO protocol results before generating any new comparison table. Preserve separate TMEM16A/VRAC currents, loading fluxes, conserved states, voltage, volume and water integrals. Compare source-matched collection-bin integrals; record 0-120, 120-180 and 180-600 s summaries without converting qualitative initial comparability into an invented tolerance.

Retain the complete Task 48 ledger. Total saliva ratio is conditionally mappable under a shared gland scale; absolute SPQ uptake is not presently mappable without the source optical calibration and regression windows. In particular, SPQ is `F0/F`; uptake is the recovery-window slope of net concentration after accounting for dilution, not NKCC/AE4 flux or the immediate post-stimulus derivative. Missing maps require `prediction: null` and `standardised_residual: null`, not a substitute fit. A within-genotype stimulus contrast has stronger cancellation than a genotype ratio because the latter begins from different chloride concentrations.

Quantitative rows should report signed residuals against source means and SEM, with assumptions and measurement-map status. Use reported uncertainty; do not announce a new universal pass band. Early preservation, AE2 neutrality and NHE similarity remain qualitative unless their numerical source errors are extracted. A failure to reject a difference is not equivalence. Do not score raw means and their derived ratios independently or turn correlated summary contrasts into a calibrated global chi-square p-value.

Stage I acceptance requires source-faithful independent reconstruction and the held-out stimulus/timing pattern, not merely breaking the mathematical CCh=combined identity. A scientifically valid stimulated correction also does not certify the measured resting states. If either piece is absent, 49I is inapplicable; do not construct a combined model by filling in the missing piece.
