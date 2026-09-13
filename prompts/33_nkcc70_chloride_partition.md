# Task 33: calibrate WT chloride transport partition to NKCC1 70%

Repository: `esig626/ae4-salivary-transport-control`
Branch: `codex/task-33-nkcc70-chloride-partition`

Read `AGENTS.md` first. This task starts from the completed Task 32 branch and keeps the Task 31 Cha NHE1 mechanism fixed.

## Purpose

Stop treating absolute WT intracellular chloride as the calibration target. The experimental transport architecture is the target in this task.

For WT R09 under the central stimulated condition (Ca = 0.25 uM), calibrate the chloride-transport system so that:

1. NKCC1 supplies 70% of **positive basolateral chloride uptake**;
2. the remaining 30% of positive basolateral chloride uptake is supplied by the bicarbonate-dependent exchanger arm, AE4 plus any positive AE2 uptake;
3. the apical Ca-activated chloride channel (CaCC) carries the matching chloride efflux required by stationary whole-cell chloride balance;
4. total chloride throughput is not allowed to collapse to a trivial near-zero solution.

This is a WT transport-partition calibration. It is not a fit to the AE4-null phenotype and it is not a demand that WT Cl_i equal 50 mM.

## Definitions

At the selected stationary WT state define:

- `J_N` = NKCC1 cycle influx, so NKCC1 chloride influx is `2 J_N`;
- `J_4` = signed AE4 chloride source to the cell;
- `J_2` = signed AE2 chloride source to the cell;
- `J_a` = apical CaCC chloride efflux from cell to lumen.

Define positive basolateral chloride loading as

`L_pos = 2*max(J_N,0) + max(J_4,0) + max(J_2,0)`.

Primary target:

`2*max(J_N,0) / L_pos = 0.70`.

Therefore the positive exchanger arm contributes 0.30 of `L_pos`:

`[max(J_4,0) + max(J_2,0)] / L_pos = 0.30`.

Do not count negative AE2 or AE4 counterflux as positive uptake. Report signed counterflux separately.

At stationary intracellular chloride balance verify directly from code that apical CaCC efflux equals the net cellular chloride loading after all signed chloride sources are included. CaCC is the apical output arm, not part of the denominator used for the 70/30 basolateral uptake split.

## Calibration strategy

Use only WT R09 to calibrate. R10 is validation only after the R09 parameter set is frozen.

First compute one reference WT R09 stationary state at Ca = 0.25 uM with the Task 31/32 model. Record its total positive chloride loading and CaCC efflux. This is the throughput reference.

Allowed uncertain magnitudes for this task:

1. NKCC1 capacity;
2. AE4 activity/capacity magnitude;
3. AE2 capacity;
4. CaCC conductance.

Do not change transporter stoichiometry, thermodynamic direction laws, NHE1 mechanism or density, pump, K channels, CO2, paracellular transport, water laws, stimulus model, or genotype rules.

Use one deterministic local constrained calibration, not a grid or combinatorial search. Prefer algebra for capacity updates at a fixed state, followed by a bounded local whole-cell closure solve. If remaining nonuniqueness exists, choose the solution minimizing squared log-deviation of the four allowed magnitudes from their inherited Task 31/32 values.

The calibrated WT state must satisfy all of the following:

- stationary numerical/conservation/current/charge/water closure gates;
- NKCC1 positive chloride-loading fraction 0.70 +/- 0.02;
- positive AE4+AE2 fraction 0.30 +/- 0.02;
- total positive chloride loading within 10% of the uncalibrated WT R09 Ca=0.25 reference, so the target cannot be met by shutting the system down;
- `Na_i < 30 mM`;
- cell volume `< 3 pL`;
- pH remains physiologically credible; use the experimental WT pH as context, but do not reject an otherwise valid state merely because it is not exactly 6.91;
- do not impose a hard target of 50.1 mM on WT chloride. Report whatever WT chloride results.

If AE2 remains a signed chloride counterflux at every admissible local solution, report that explicitly. Do not force its sign with a genotype switch or hidden term. In that case the 30% positive exchanger arm may be carried predominantly or entirely by AE4, and that is a result of the model under this calibration.

## Freeze before genotype evaluation

Once one R09 WT calibration satisfying the declared transport split is obtained:

1. write and hash the frozen parameter set;
2. do not change any parameter again;
3. set AE4 expression exactly to zero;
4. solve the exact AE4-null REST state from multiple physically motivated warm starts only if needed, with a small fixed limit;
5. do not use knockout chloride, pH, secretion, or the known 30-35% secretion deficit in any calibration step.

For WT and exact null report:

- Na, K, Cl, TIC, HCO3, pH, cell/lumen volumes and voltages;
- NKCC1 cycle and chloride flux;
- AE4 and AE2 signed chloride fluxes;
- NHE1 flux;
- CaCC chloride efflux;
- pump and K fluxes;
- paracellular chloride;
- osmolarities and all water fluxes;
- positive chloride-loading fractions and net chloride balance.

The exact-null resting observations (Cl_i about 36.5 mM and pH about 6.89) are validation context only.

## Dynamics

If an admissible matched R09 WT/null REST pair exists, run only the central CCh+IPR / Ca = 0.25 uM 0-600 s WT and exact-null trajectories with the frozen parameters.

Report:

- cumulative secretion;
- null/WT secretion ratio;
- integrated NKCC1, AE4, AE2 and CaCC chloride transport;
- stimulated positive-loading fractions;
- whether NKCC1 net flux changes after AE4 deletion despite frozen capacity;
- lumen osmolarity, paracellular chloride return and water flux changes.

Only after these results are frozen compare with the experimental approximately 30-35% secretion reduction. Do not tune toward it.

If R09 succeeds, run R10 only as an out-of-sample stress check with exactly the same calibrated magnitude scalars. Do not recalibrate R10.

## Numerical limits

This is a small targeted calibration, not a parameter search.

Absolute limits:

- at most 12 stationary whole-cell solver calls;
- at most 20,000 stationary residual evaluations;
- at most 6 local calibration iterations/evaluations;
- at most 4 stimulated integrations total;
- one worker and one BLAS thread;
- at most 30 minutes numerical execution.

No Cartesian grid, random search, evolutionary/Bayesian/global optimiser, pairwise enumeration, Shapley analysis, or automatic escalation to additional mechanisms.

Stop on the first exhausted limit and preserve partial outputs.

## Required conclusion

Answer plainly:

1. Was a closed WT R09 state obtained with NKCC1 carrying about 70% of positive basolateral chloride uptake?
2. How was the remaining approximately 30% divided between AE4 and AE2?
3. What CaCC flux/conductance was required to carry the resulting apical chloride output?
4. What WT chloride concentration resulted, without fitting it?
5. After freezing the WT calibration, did exact AE4 loss produce an admissible REST state?
6. How much did chloride, pH, sodium, NKCC1 flux and water secretion change?
7. Did secretion fall in the correct direction and by what magnitude?
8. What remains unresolved?

## Publication

Use the connected GitHub integration for all remote writes. Do not depend on shell Git authentication, PATs, SSH or `git push` credentials.

Publish source edits, focused tests and compact UTF-8 result/report files directly to this branch. Do not merge to `main`. Publication trouble must never trigger scientific recomputation.
