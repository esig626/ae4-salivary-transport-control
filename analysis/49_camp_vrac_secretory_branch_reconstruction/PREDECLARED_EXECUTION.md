# Task 49 execution contract

Orchestrator acceptance of the 49a addendum at requested head `2bf53b773c14f3c49f539957ea3f1964bc8d227f`. Repository instructions take precedence over advisory recommendations. Five workers provide evidence/history, electrical and validation advice only. Only the orchestrator modifies this branch, evaluates the scientific model, accepts conclusions and publishes.

## One stimulated law

The sole law is the already documented Task 13B diagnostic

`G_V = g_V beta max(V_i/V_rest,genotype - 1, 0)`, `g_V >= 0`.

Its current is `I_V = G_V (V_a-E_Cl)` and chloride export is `-10^15 I_V/F`. Insert it inside both membrane current equations, with equal/opposite cell/lumen chloride sources. Reuse the historical gate and retain its unlicensed-capacity status until independent evidence identifies a scale. The current is separate from TMEM16A; no AE4-expression factor exists. No alternative gate, delay, threshold, Hill exponent, beta-NKCC change or unmeasured synergy is permitted.

All inherited source files and parameters remain unchanged. A Task 49 local wrapper may add the new current and consistently recalculate voltages and the consequent currents under their unchanged constitutive laws. Synthetic swollen states and arbitrary nonzero test conductances are numerical verification fixtures only. They must never become production initial conditions or inferred parameters.

## Calibration gate

Only independently extracted Catalán PNAS 2015 cAMP observations may identify `g_V`. Preserve separate gland, volume and perforated-current protocols, biological SEM and digitisation uncertainty. A total IPR-induced current is not necessarily the DCPIB-sensitive component; a whole-cell component is not automatically a known apical conductance. A separate volume group's swelling is not a simultaneous patch-current gate measurement. Whole-gland amplitude has an unidentified scale; scale-free shapes and within-source ratios can be informative only if the model actually activates the gate.

First inspect the saved current-parent WT IPR trajectory, with hashes; do not recompute published trajectories. Saved samples establish sampled gate behaviour, not a continuous-time enclosure. Use exact algebra where possible. Do not promote numerical resting residuals to a swelling trigger. Source incompatibility or a missing observation map blocks a calibrated reconstruction. A blocked result must record `calibrated_conductance_S: null`, not an invented zero or a diagnostic test scale.

49D precedes any new AE4 comparison even if 49C blocks: publish an immutable unavailable-prediction record plus any parameter-independent analytic prediction. 49E then records the consequent comparison limitations. This honours the addendum's continuation requirement without calling an unavailable numerical prediction successful. Do not use the AE4 secretion endpoint or combined uptake for calibration. Existing instructions disclose those numbers, so the safeguard is procedural data separation, not literal blinding.

## Resting geometry gate

Only after 49E is pushed and remotely verified, evaluate the measured resting Cl/pH states. Use all four source cohorts: AE4 WT, AE4 KO, AE2 control and AE2 KO. Keep their reported means and SEs distinct and use each genotype's cached equilibrium for unmeasured coordinates. A shared WT model represents both controls; do not fit strain multipliers.

Cl and pH alone are not complete conserved states. Use the single algebraic lifting in `advisory/rest_design/STAGE_II_DESIGN_ADVICE.md`: retain genotype-specific Na concentration and complete lumen from the saved rest; derive cell CO2 from the necessary stationary identity `C=H`; derive TIC from measured pH; derive cell volume from cell water stationarity together with electroneutrality/osmotic accounting; derive K and TA from exact charge and measured pH. Reject nonpositive results rather than clipping. This uses balance equations without any optimiser. The retained Na/lumen are conditional inherited values, not measurements, and the resulting states are not claimed to solve all balances. Propagate the eight reported Cl/pH SEs by a local derivative; separately expose latent-state sensitivity instead of assigning invented error bars.

Predeclare the complete fixed basis: NKCC1, NHE1, AE2, equal-routing AE4, NBC, apical/basolateral pump and K, apical TMEM16A chloride, the separate (rest-inactive) VRAC chloride signature, basolateral/apical CO2, four paracellular ion paths, three water paths and lumen outflow. Include cell and lumen Na/K/Cl/TIC/TA/volume, and derive charge/current rows from those signatures. Include dormant NBC/VRAC physical signatures and record their zero native availability separately. Deleted transporter columns are zero in their genotype blocks. Water/outflow signatures are bookkeeping directions; their presence does not authorise changing water laws.

Assemble this complete source matrix once at the four central states. Compute rank/SVD, project the missing source `-F`, and exhibit its source-coefficient nullspace. Any cone test must respect reversible directions and distinguish gross nonnegative rates from signed corrections to existing nonzero fluxes. A cancelling source representation is not a fitted shared law. If necessary use at most one constrained continuous least-squares solve with all columns; no subset selection, refits, on/off combinations or parameter grids.

The residual automatically lies in the span of a complete source basis used to assemble that same vector field. Therefore a successful span projection alone cannot identify a missing transporter or justify NBC. Audit the shared constitutive equations and exact nullspace identities. State whether any unique shared correction follows; if not, preserve the non-identifiability and do not edit the model. In particular, at rest `TA_dot=H+2B-E-2A` and `TIC_dot=C+2B-E-2A`; stationarity requires `C=H`. In the AE4 knockout with the inherited NBC gate, `B=A=0`, hence `H=E`. A positive alkalinity residual cannot be repaired by positive inward NBC at the same state.

49H requires a uniquely identified lawful correction. 49I requires both independently accepted pieces and no retuning. Otherwise mark both inapplicable. Complete 49F with a precise equation-level diagnosis and its identifiability limits.

## Publication and execution boundaries

Every milestone gets CURRENT_STATUS.md, a machine-readable receipt, a commit, a push and an independently read remote SHA before dependent work. Receipts contain parent SHA and hashes; the following commit retains the prior post-push verification. No force-push or history rewrite. Reuse published results after interruption. No prior numbered analysis, frozen result, manuscript or other branch may change.
