# Task 49 advisory: separate VRAC current with inherited electrical laws

Advisory only. No canonical files were changed, no model was evaluated or integrated, and no scientific conclusion was accepted by this worker. The orchestrator must decide what to accept and publish.

## Repository finding

`src/modern_full_model/vbeta_diagnostic.py` already implements the exact proposed minimal gate and fixed-voltage current. Its `V1_RECTIFIED_SWELLING_GATE` has

`Gv = gV * beta_input * max(V_i / V_rest - 1, 0)`.

It already derives chloride reversal, cell/lumen sources and charge conversion correctly. It deliberately does not insert its current into the membrane closure. Both its production-ineligible flag and its unlicensed-conductance provenance must survive unless new Task 49 independent evidence resolves that limitation. The new work is a coupled electrical insertion, not a newly discovered gating law. Task 13B's old deadlock trajectory is not proof about the different Task 48 active state, but its anti-deadlock test remains applicable. Reuse `vbeta_effective_conductance_S`, `evaluate_vbeta_current` and `assess_pre_vbeta_trajectory` where applicable.

## Exact sign and source derivation

The active convention is `Va = phi_cell - phi_lumen`, `Vb = phi_cell - phi_bath`, `Vt = Vb - Va`. Conventional current is positive cell to lumen at the apical membrane; positive molar apical chloride flux is chloride cell to lumen.

For chloride,

`ECl = -(RT/F) log(Cl_l / Cl_i)`;

`Iv = Gv * (Va - ECl)`;

`Jv_cell_to_lumen = -Iv * 1e15/F` in fmol/s.

Thus the separate VRAC cell chloride source is `+Iv * 1e15/F`; its lumen chloride source is `-Iv * 1e15/F`. The cell conventional-charge source is `-Iv * 1e15/F`, the lumen charge source is the opposite, and their chloride sum is exactly zero. A negative apical chloride current means chloride exit and therefore a negative cell chloride source. Do not append these two chloride sources after an already closed electrical solve: that would create nonzero charge derivatives.

VRAC contributes no direct Na, K, TIC, alkalinity or water source. Its indirect effect on NBC through Vb is required by the unchanged electrogenic NBC law.

## Small task-local adapter: voltage increments around the baseline closure

The existing Task 48 model chain is `ProtocolModel(MinimalNbcModel(StimulatedNkcc1Model(SourceBathCore)))`. Wrap the final `ProtocolModel`, retaining its intervention handling, rather than reconstructing that chain or altering shared source files. Ensure the wrapper owns `evaluate`, `rhs`, and `solve_dynamics`: forwarding `solve_dynamics` through `__getattr__` would silently integrate the unaugmented model.

1. Evaluate the inherited model at the same `(t, state, genotype)` to obtain the baseline `ModelEvaluation` and its self-consistent membrane currents.
2. Evaluate Gv with the existing diagnostic gate, using the independently fixed gV and the protocol's saved genotype-specific pre-stimulus volume. If Gv is exactly zero, return the inherited RHS literally. Any additional diagnostic should not recalculate baseline fluxes.
3. For Gv > 0, solve only the voltage increments and update membrane current/source contributions. All other laws and coefficients remain identical.

Let baseline voltages be `Va0, Vb0`; let `x=Va-Va0`, `y=Vb-Vb0`. At this same state the unchanged calcium gate `h` gives

`GKa = g_k_total*h*apical_k_fraction + g_apical_background`;

`GKb = g_k_total*h*(1-apical_k_fraction) + g_basolateral_background`;

`GCl = g_cl_apical*h`;

`Ga=GKa+GCl`, `Gb=GKb`, `Gp=sum(all four paracellular conductances)`.

Here `g_cl_apical` always remains the inherited TMEM16A parameter. Do not temporarily add VRAC to it: that obscures the separate branch and would require division by a zero calcium gate in a valid edge case.

Let `eps_a=baseline.Iapical-baseline.Ipara` and `eps_b=baseline.Ibasolateral_total+baseline.Ipara`. Preserve their finite numerical residuals in the incremental equations. The new apical equation is

`eps_a + (Ga+Gv+Gp)*x - Gp*y + Gv*(Va0-ECl) = 0`.

Therefore

`x(y) = [Gp*y - Gv*(Va0-ECl) - eps_a] / (Ga+Gv+Gp)`.

The new basolateral equation is the scalar equation

`eps_b + Gb*y + (INBC(Vb0+y)-INBC0) + Gp*(y-x(y)) = 0`.

Evaluate `INBC` with the existing `evaluate_minimal_nbc`, the existing activation and `MinimalNbcParameters`, and the same concentrations. Do not change activation, capacity or NBC source law. The inherited NBC contribution is exactly zero at its resting calcium. Read the unchanged activation and baseline NBC current/cycle from inherited diagnostics, which exist even on the zero-activation return path.

The scalar left-hand side is strictly increasing whenever the inherited passive conductances make closure nonsingular. Its derivative is

`Gb + Gp*(1-Gp/(Ga+Gv+Gp)) + dINBC/dVb > 0`,

where

`dINBC/dVb = (F/1e15) * capacity * activation / (log_width*thermal_voltage) * sech(affinity/log_width)^2 >= 0`.

Consequently this is a single physical current-closure root, not a parameter search. Use the inherited absolute voltage bracket `[-0.5,0.5] V`, translated to y by subtracting Vb0. A bracket failure is a failure, not permission to widen it repeatedly. At zero NBC activation the equations are simply one 2-by-2 linear electrical system.

## Reconstruct sources from unchanged current laws

All current increments follow directly from existing conductances at this fixed state:

- `delta_I_Ka=GKa*x`;
- `delta_I_Cl_TMEM16A=GCl*x`;
- `delta_I_Kb=GKb*y`;
- `delta_I_para_species=g_para_species*(y-x)`;
- `Iv=Gv*(Va0+x-ECl)`;
- `delta_I_NBC=INBC(Vb0+y)-INBC0`.

Pump cycles are unchanged because the active pump is voltage-independent. New currents may differ because voltage changes; their physical laws and parameters do not change. Keep the inherited `cl_apical` current as TMEM16A and add a separate `vrac_apical` key; `apical_total` must sum both. Recompute `basolateral_old`, `basolateral_total` and `paracellular_total` consistently.

Convert every delta current to its molar direction with the existing `current_to_fmol_s`. Denote the existing channel delta molar fluxes by `dJKa,dJKb,dJCl`, the VRAC molar flux by `Jv`, the four paracellular delta fluxes by `dJpNa,dJpK,dJpCl,dJpHCO3`, and the NBC inward cycle change by `dB`. The membrane-source increments in `(Na,K,Cl,TIC,TA)` order are

`delta_cell = (dB, -dJKa-dJKb, -dJCl-Jv, 2*dB, 2*dB)`;

`delta_lumen = (-dJpNa, dJKa-dJpK, dJCl+Jv-dJpCl, -dJpHCO3, -dJpHCO3)`.

Add these increments to the inherited membrane sources and to RHS rows `0:5` and `6:11`. Leave volume rows 5 and 11 and all regulatory rows untouched. Leave the inherited homeostasis, AE4, water, CO2 and outflow objects untouched at this state. The existing NHEx/NKCC/AE2 sources are concentration-dependent and receive no direct electrical update.

Recompute both membrane current residuals and both membrane charge-vs-current residuals from the reconstructed sources. Recompute full-model cell charge rate and lumen charge rate after subtracting outflow; keep the existing charge convention `Na+K-Cl-TA`.

For carbon accounting, the inherited external-carbon expectation changes by `2*dB-dJpHCO3`. Equivalently, recalculate it exactly as `MinimalNbcModel.evaluate` already does, using its unchanged homeostasis/AE4/CO2/outflow and the new membrane cell/lumen TIC sources. It is incorrect to leave the diagnostic NBC current/cycle/carbon contribution stale after changing Vb. The same applies to all `minimal_nbc_*` voltage-dependent diagnostics.

## Finite independent verification, without trajectories

The following checks can be performed on a few predeclared synthetic/reference states. They are verification cases, not a parameter panel or calibration sweep.

1. Zero conductance, beta=0, exact resting volume, and sub-rest volume each return the inherited RHS exactly; CCh-only does not change. This detects an accidental change of unrelated laws.
2. Nernst reversal and its sign: at `Va=ECl` VRAC current is exactly zero; just below it, chloride exits and the cell chloride source is negative. Check cell/lumen cancellation and charge conversion against direct valence accounting.
3. At one artificially swollen test state with nonzero beta, solve the augmented two membrane equations independently as a two-variable root from direct original-current formulas. Compare both voltages and current balances with the scalar incremental implementation. A synthetic state verifies numerics only and never licenses physiological activation.
4. In a second case with NBC activated, confirm that the new NBC current is evaluated at the new Vb; verify full cell/lumen charge, carbon and chloride identities. This catches the main plausible implementation error: appending a separate anion current while leaving electrogenic NBC stale.
5. Use one zero-calcium/TMEM16A-zero test case to establish finite VRAC current and no division by the calcium gate. This also verifies the intended TMEM16A-independent topology. It is an acute model mask, not a source-backed claim about chronic knockout compensation.
6. At the same nonzero-Gv state verify unchanged pump cycles, calcium input, NKCC capacity multiplier, homeostasis/AE4/regulatory RHS, water fluxes and CO2 exchange. Allow K/TMEM16A/para/NBC currents to change solely through voltage.
7. At beta=0 or Gv=0 compare all canonical parameter objects or their hashes before and after evaluation; there should be no parameter mutation at any Gv.

Existing unit tests `tests/test_modern_vbeta_diagnostic.py` already cover the gate and fixed-voltage signs. Existing `tests/test_task36_minimal_nbc.py` supplies current/source closure conventions and absolute residual tolerances. Reuse those checks, adding only the coupled-closure and unchanged-law tests needed for this new adapter.

## Anti-deadlock and licensing requirements

`V_rest` must be fixed from the same genotype's declared pre-stimulus rest, not a common WT volume, a fitted threshold or the instantaneous volume. A failed physiological rest remains a failed physiological rest; normalising its volume does not validate it.

For the chosen one-sided gate, examine the current active model's pre-VRAC IPR trajectory before any positive conductance can be promoted. If beta-positive volumes never exceed their starting reference, the gate is identically zero for every gV on that trajectory. The augmented field then coincides with the old field throughout it, so changing gV cannot create its own missing activating swelling. No new integration at a collection of positive conductances is necessary to prove that conditional statement. Old Task 13B outcomes cannot automatically stand in for this check of the newer active model.

This argument is conditional on the active baseline trajectory remaining at or below Vrest; it is not a biological rejection of VRAC and does not imply that another threshold, direct beta current, water change or beta-NKCC recruitment is permitted in Stage I.

Even if the gate opens, the topology and conductance are separate questions. Do not replace the inherited `UNLICENSED_NO_WHOLE_CELL_CURRENT_MAP` status with production eligibility merely because the electrical adapter conserves charge. Independent source data must identify the one fixed conductance under a declared observation map. Any new Task 49 evidence should explicitly distinguish a measured fractional volume/current quantity from an assumed model whole-cell conversion.

## Code inspected

- `AGENTS.md` and Task 49 entry, addendum, exclusion/source/plan/prompt materials.
- `src/modern_full_model/membranes.py`, `model.py`, `nbc_minimal.py`, `nkcc_stimulation.py`, `vbeta_diagnostic.py`.
- `analysis/48_joint_experimental_constraint_reconstruction/protocol_layer.py`.
- `tests/test_modern_vbeta_diagnostic.py`, relevant NBC unit tests, and the previous `analysis/13B_modern_full_model/vbeta_diagnostic.md`.

No claim here independently verifies the PNAS primary source; that evidence audit belongs to the source worker and orchestrator.
