# Task52 deterministic onset projections

These equations implement D52-04 and the remotely verified 52A freeze, following R51H, R51I and the Task52-start ledger. They impose the measured chloride and pH at stimulation onset. They do not solve or approximate a stationary transport state. No pre-stimulation relaxation is applied.

Let `Na`, `K`, `Cl` and `T` denote intracellular concentrations in mM, `V` cell volume in pL, and `Q` the unchanged fixed-cell anion equivalents in fmol. At the imposed pH, the inherited carbonate fractions give

`a = alpha_HCO3 + 2*alpha_CO3`.

For unchanged finite-buffer amount `B`, define the deprotonated amount and water alkalinity

`b = B / [1 + 10^(pKa_buffer - pH)]`,

`w = [OH] - [H]`.

The exact inherited alkalinity formula is then

`TA = a*T + b/V + w`.

Let `M = other_impermeant_osmoles + B`. The finite buffer contributes one osmotic particle per molecule and is included **exactly once**. `Q` is the charge parameter and is not added as another osmotic particle. The bath osmolarity `O` is the inherited sum `Na_bath + K_bath + Cl_bath + TIC_bath + untracked_osmolyte_bath`.

Bulk charge neutrality and bath-isotonicity require

`Na + K - Cl - TA - Q/V = 0`,

`Na + K + Cl + T + M/V = O`.

At fixed pH these are linear in the two unknown coordinates after using `1/V` for the alternate projection. Eliminating `K` gives

`O = 2*Cl + (1+a)*T + w + (Q+b+M)/V`,

and then

`K = Cl + a*T + w + (Q+b)/V - Na`.

For the central projection, retain the accepted Task37 sodium concentration and cell volume:

`T = [O - 2*Cl - w - (Q+b+M)/V] / (1+a)`.

Evaluate `K` from the equation above. This uniquely fixes K and TIC without an optimisation, root iteration, transporter evaluation or secretion input.

For the alternate projection, retain the accepted Task37 sodium and TIC **concentrations**:

`V = (Q+b+M) / [O - 2*Cl - w - (1+a)*T]`.

Again evaluate `K` from the charge equation. This is the sole predeclared projection sensitivity. No alternative projection is selected if this one fails.

The implementation finally obtains TA by calling the byte-exact Task37 `total_alkalinity_mM` function, converts concentrations to amounts through `n=C*V`, and leaves the original lumen and regulatory coordinates unchanged. It returns residuals calculated from the completed state and inherited osmolarity function. Independent pH/speciation reconstruction and the unchanged physical/conservation acceptance gates are applied by the caller. Nonfinite, singular or nonpositive solutions raise `ValueError`.

The imposed source values are WT `Cl=50.10 mM`, `pH=6.91`, and complete AE4 KO `Cl=36.50 mM`, `pH=6.89`. Unmeasured projected K/TIC/volume are consequences of this explicit initialisation idealisation, not new source measurements. Bath-isotonicity does not assert zero water-flow derivative because the inherited lumen state also contributes to the water laws. A passing projector cannot establish that the unchanged chronic transport model maintains either measured genotype state.

`project_onset.py` imports only the standard-library `math` module at module load. Its exact Task37 acid-base/water imports are lazy and use the isolated `task52_pre_palk` package. This implementation subtask performs no imports of scientific modules, numerical evaluations or tests; those belong to the root's single scientific verification process.
