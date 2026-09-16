# Read-only advisory review of executed 49G

Disposition: **no blocking algebraic defect found** in `resting_geometry.py` or the requested saved outputs. Scientific acceptance remains the orchestrator's responsibility. This review read code and existing JSON only; it did not execute model/RHS code, reconstruct a matrix, perform SVD, solve an inverse problem, or generate a trajectory.

The algebraic lift correctly retains six declared nuisance coordinates and derives cell CO2/TIC, volume, K and TA from measured Cl/pH, both charge constraints, `CO2 influx = NHE1`, and cell water stationarity. The fixed buffer contributes once as an osmotic particle pool and through its deprotonated fraction in TA. This is explicitly a conditional lift, not a measured full state or solved equilibrium.

The source columns, units, equal-routing AE4 stoichiometry, membrane current signs, KO masks and advective outflow are consistent with the active source conventions. Dormant NBC/VRAC have physical signatures in the relaxed span matrix but zero native cone availability. Charge/current rows are bookkeeping dependencies. The saved rank 48/nullity 36 and near-zero span residual are therefore appropriately interpreted as a relaxed source decomposition, not transporter identification or a shared-law fit.

The KO TA separator is exact for the recorded conditional state: NHE1 is positive and AE2 runs in reverse, so their native TA contributions are both positive; the remaining active cell-TA columns vanish. The saved target TA is `-0.009394838357764957 fmol/s`. A negative target cannot lie in their nonnegative addition cone. Shared-law tying and current constraints cannot undo this pointwise separator. The stronger total-capacity reweighting statement is also correct when the state and affinity signs are retained: both TA-producing activities would have to be suppressed to zero.

The native NBC affinity audit is positive (`2.733003380453671`) at the unchanged voltage. The report correctly retains its zero gate and does not claim to have run basal NBC or shown a global NBC impossibility theorem. A positive inward NBC source cannot cancel the fixed-state positive TA residual. The scalar KO compatibility equation and the bounds involving frozen AE2 capacity and `D(Cbar-Ccrit)` are algebraically correct.

The uncertainty file distinguishes measurement-SE propagation from unquantified nuisance-state uncertainty. Its six nuisance directions are consistent with the two eliminated stationary identities. No reported confidence interval or nuisance-robust separator is claimed. Accordingly the justified conclusion is conditional native-cone failure plus **no uniquely identified shared correction**. The outputs do not justify a basal NBC edit, a new transporter/species assertion, or a global impossibility claim over all hidden states.

The repeat execution is disclosed as two actual assemblies/SVD calls using one distinct basis after the diagnostic-key error. It must remain disclosed; it should not be described elsewhere as one numerical SVD in the entire session.

Two small metadata clarifications can be made without recomputation:

1. `min_norm_coefficient_representative` is minimum norm in the **column-normalised coefficient coordinates** `z=D_column*delta`, not necessarily minimum Euclidean norm of the exported raw coefficients `delta`. Add that qualifier or rename the key to `column_scaled_min_norm_coefficient_representative`.
2. `source_hashes_unchanged=True` is asserted as a literal in the export script. It should be tied to the independently verified immutable-source hash check in the checkpoint/report, rather than presented as a check performed by this script. The first-attempt failure occurred during summary/export preparation; the transparent counts themselves are adequate.

Neither clarification changes any scientific number, basis, inference, or correction decision; no model rerun is required.
